/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "SMPGCColoring.h"
using namespace std;
using namespace ColPack;



// ============================================================================
// Author: Xin Cheng
// D2GC omp (Distance Two Graph Coloring OpenMP parallel) Speculative Iterative approach
// -------------------------------------------------------
// Pseudo code is as follows
// while(|W| is not empty)
//     TentativeColoring(W)  //TC
//     CheckConflicts(W)     //CC
//     HandleConflicts(W)    //HC
// -------------------------------------------------------
// The TC could be one of
//      do_D2GC_omp_TentativeColoroing()
//      do_D2GC_omp_NetBinTentativeColoroing()
//      do_D2GC_omp_RndTentativeColoring()
//
// The CC could be one of
//      do_D2GC_omp_CheckConflicts()
//      do_D2GC_omp_NetBinCheckConflicts()
// 
// The HC could be one of
//      *do_nothing*
//      serial_greedy_coloring()
// --------------------------------------------------------
// ============================================================================
int SMPGCColoring::D2GC_OMP_SpeculativeIterative_MemOpt(
        int nT,                     // number of threads
        int &colors,                // number of colors
        vector<int>&vtxColors,      // mapping vertex to its color
        int const LOCAL_ORDER,      // parallel local ordering
        int const NUM_RND_TC,       // number of iterations using Random TC
        int const NUM_NET_TC,       // number of iteartions using NetBin TC
        int const NUM_NET_CC,       // number of iteartions using NetBin CC
        bool const b_SEQ_HC,        // whether using sequential coloring to color the remaining graphs.
        int const nVerbose          // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    ChronoDuration         tim_PPT(0), tim_LO(0);
    vector<ChronoDuration> vtim_RndTC,   vtim_NetTC,  vtim_TC;
    vector<ChronoDuration> vtim_NetCC,   vtim_CC;
    vector<int>            vnum_NetCC,   vnum_CC;
    ChronoDuration         tim_HC;

    const int N = num_nodes();                  //number of vertex
    int const D = max_degree();
    int const DDmDp1 = D*D-D+1;
    const int BufSize = DDmDp1<=0?N:min(DDmDp1, N); //maxDegree
    const vector<int>& vtxPtr = get_CSR_ia();
    const vector<int>& vtxVal = get_CSR_ja();
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    const int ColorLB0Base = max_degree();

    vector<mt19937> mts(nT);  //random generator for each thread
    colors=0;                       
    vtxColors.assign(N, -1);

    // allocate memory
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16,-1));  
    vector<int> Qsizes(nT,0);
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16,-1));
    int const Capacity = 32;
    vector<vector<int>> Workspace(nT, vector<int>(Capacity,-1));   // used by NetBinColor 
    
    // prepartition
    {
        auto start = steady_clock::now();
        Qsizes.assign(nT,N/nT);
        vector<int> disp(nT+1,0);
        for(int tid=0; tid<N%nT; tid++) Qsizes[tid]++;
        for(int tid=0; tid<nT; tid++)   disp[tid+1]=disp[tid] + Qsizes[tid];
        for(int tid=0; tid<nT; tid++)   QQ[tid].assign(const_ordered_vertex.begin()+disp[tid], const_ordered_vertex.begin()+disp[tid+1]);
        auto end = steady_clock::now();
        tim_PPT = end - start;
    }

    // local order
    ordering("OMP","PREPARTITION", LOCAL_ORDER, &tim_LO, &QQ);

    // start coloring iterations
    int niter=0;
    int num_uncolored_vertex = N;
    do{
        // phase - Tentative Coloring
        if(niter<NUM_RND_TC){
            vtim_RndTC.emplace_back(.0);
            do_D2GC_OMP_PhaseRandomTC(vtxColors, QQ, Qsizes, mts,  ColorLB0Base, &vtim_RndTC.back());
        }
        else if(niter<(NUM_RND_TC+NUM_NET_TC)){
            vtim_NetTC.emplace_back(.0);
            do_D2GC_OMP_PhaseNetBinTC_MemOpt(vtxColors, const_ordered_vertex, N, vtxPtr, vtxVal, Workspace, Capacity, &vtim_NetTC.back());
        }
        else{
            vtim_TC.emplace_back(.0);
            do_D2GC_OMP_PhaseTC_MemOpt(vtxColors, QQ, Qsizes, vtxPtr, vtxVal, Capacity, &vtim_TC.back());
        }

        // phase - Check Conflicts
        if(niter < NUM_NET_CC){
            vtim_NetCC.emplace_back(.0);
            do_D2GC_OMP_PhaseNetBinCC_MemOpt(QQ, Qsizes, vtxColors, vtxPtr, vtxVal, const_ordered_vertex, N, Capacity, &num_uncolored_vertex, &vtim_NetCC.back());
            vnum_NetCC.emplace_back( num_uncolored_vertex );
        }
        else{
            vtim_CC.emplace_back(.0);
            do_D2GC_OMP_PhaseCC(QQ, Qsizes, vtxColors, vtxPtr, vtxVal, &num_uncolored_vertex, &vtim_CC.back());
            vnum_CC.emplace_back( num_uncolored_vertex );
        }
        



        // phase - Handle Conflcits
        if( b_SEQ_HC ){
            auto start = steady_clock::now();
            Fs[0].assign(BufSize,-1);
            for(int i=0; i<nT; i++)  
                do_D2GC_Seq_GreedyColoring_NoInitF(vtxPtr, vtxVal, vtxColors, QQ[i], Qsizes[i], Fs[0], BufSize);
            tim_HC = (start - steady_clock::now());
            num_uncolored_vertex=0;
        }

        niter++;
    }while(num_uncolored_vertex!=0); //end of while coloring iterartion

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){
        stringstream ss;
        ss<<"@D2GCMemOpt"<<(b_SEQ_HC?"3P":"MP");
        ss<<"("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<")";
        ss<<"RNN("<<NUM_RND_TC;
        ss<<","<<NUM_NET_TC;
        ss<<","<<NUM_NET_CC;
        ss<<")";

        for(int i=ss.str().size(); i<30; i++) ss<<"_";
        ss<<"\t#nT_c_T";
        ss<<"\t"<<nT;
        ss<<"\t"<<colors;
        ChronoDuration tim(.0);
        {
            //tim=tim_PPT+tim_LO;
            for(const auto& x: vtim_RndTC) tim+=x;
            for(const auto& x: vtim_NetTC) tim+=x;
            for(const auto& x: vtim_TC)    tim+=x;
            for(const auto& x: vtim_NetCC) tim+=x;
            for(const auto& x: vtim_CC)    tim+=x;
            tim+=tim_HC;
        }
        ss<<"\t"<<tim.count();
        cout<<ss.str(); ss.str("");
        // more details
        if(nVerbose>1){
            cout<<"\t#niter\t"<<vtim_CC.size()+vtim_NetCC.size();
            ChronoDuration tim_TC(.0);
            {
                for(const auto& x: vtim_RndTC) tim_TC+=x;
                for(const auto& x: vtim_NetTC) tim_TC+=x;
                for(const auto& x: vtim_TC)    tim_TC+=x;
            }
            cout<<"\t#TCtime\t"<<tim_TC.count();
            ChronoDuration tim_CC(.0);
            {
                for(const auto& x: vtim_NetCC) tim_CC+=x;
                for(const auto& x: vtim_CC)    tim_CC+=x;
            }
            cout<<"\t#CCtime\t"<<tim_CC.count();
            if(b_SEQ_HC)
                cout<<"\t*SeqHCtime\t"<<tim_HC.count();
            
            int nConf=0;
            {
                for(const auto& x: vnum_NetCC) nConf+=x;
                for(const auto& x: vnum_CC)   nConf+=x;
            }
            cout<<"\t#nConf\t"<<nConf;

            // more details for each iteration 
            if(nVerbose>2) {
                if(vtim_RndTC.size()!=0) {
                    cout<<"\n*RndTCtimes "<<vtim_RndTC.size()<<" iterts: ";
                    for(const auto& x: vtim_RndTC)
                        cout<<" "<<x.count();
                }
                if(vtim_NetTC.size()!=0) {
                    cout<<"\n*NetbinTCtimes "<<vtim_NetTC.size()<<" iterts: ";
                    for(const auto& x: vtim_NetTC)
                        cout<<" "<<x.count();
                }
                if(vtim_TC.size()!=0) {
                    cout<<"\n*PlainTCtimes "<<vtim_TC.size()<<" iterts: ";
                    for(const auto& x: vtim_TC)
                        cout<<" "<<x.count();
                }
                 if(vtim_NetCC.size()!=0) {
                    cout<<"\n*NBCCtimes "<<vtim_NetCC.size()<<" iterts: ";
                    for(const auto& x: vtim_NetCC)
                        cout<<" "<<x.count();
                }
                if(vtim_CC.size()!=0) {
                    cout<<"\n*CCtimes "<<vtim_CC.size()<<" iterts: ";
                    for(const auto& x: vtim_CC)
                        cout<<" "<<x.count();
                }
                if(b_SEQ_HC){
                    cout<<"\n*SeqHCtimes "<<tim_HC.count();
                }
                if(vnum_NetCC.size()!=0) {
                    cout<<"\n*Conflicts_find_in_NBCC "<<vnum_NetCC.size()<<" iterts: ";
                    for(const auto& x: vnum_NetCC)
                        cout<<" "<<x;
                }
                if(vnum_CC.size()!=0) {
                    cout<<"\n*Conflicts_find_in_CC "<<vnum_CC.size()<<" iterts: ";
                    for(const auto& x: vnum_CC)
                        cout<<" "<<x;
                }
            }
        }// end of if nVerbose>1 
        cout<<endl;
    }// end of if nVerbose>0
    return true;
}


// ============================================================================
// Author: xin cheng
// impements of the D2GC omp Speculative Iteration approach's 
// TentativeColoring Phase
//
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseTC_MemOpt(
        vector<int>&vtxColors,         // mapping vertex to its color
        const vector<vector<int>>& QQ, // vertex to color for each thread
        const vector<int>& Qsizes,     // number of the vertex to color for each thread
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        int const Capacity,             // forbidden color array size
        ChronoDuration *ptime       // runtime
        ) {
    chrono::steady_clock::time_point start;
    if(ptime) start = steady_clock::now();
    #pragma omp parallel
    {
        int const BUFFWIDTH = Capacity;
        unsigned int F = ~0;
        int const tid = omp_get_thread_num();
        vector<int> const &Q = QQ[tid];
        int const Qsize = Qsizes[tid];
        
        for(int iu=0; iu<Qsize; iu++) {
            int const u = Q[iu];
            int offset_mask = 0;
            while(true){
                F = ~0;
                int const LOW = (offset_mask++)*BUFFWIDTH;
                // distance one neighbors
                for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                    int const vc_loc = vtxColors[vtxVal[iv]] - LOW;
                    if(vc_loc>=0 && vc_loc<BUFFWIDTH) 
                        F &= ~(1<<(vc_loc));
                }
                // distance two neighbors
                for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                    int const v = vtxVal[iv];
                    for(int iw=vtxPtr[v]; iw!=vtxPtr[v+1]; iw++) {
                        int const w  = vtxVal[iw];
                        int const wc = vtxColors[w]; 
                        if( w==u || wc<0) continue;
                        int const wc_loc = wc - LOW;
                        if(wc_loc>=0 && wc_loc<BUFFWIDTH) 
                            F &= ~(1<<(wc_loc));
                    }
                }
                // find the first settled bit, if there is any
                if(F!=0){
                    for(int i=0; i<BUFFWIDTH; i++){
                        if(F&(1<<i)) {
                            vtxColors[u] = LOW+i;
                            break;
                        }
                    }
                    break; // break the while(true)
                }
            }// end of while(true)
        }// end of for iu in Q of Qsize
    }// end of omp parallel
    if(ptime) *ptime = steady_clock::now()-start;
}


// ============================================================================
// Author: Xin Cheng
// Implements of the D2GC omp Speculative Iterative approaches
// TentativeColoring Phase 
//     using Net-Based Algorithm
// 
// The idea is as following:
// -------------------------------------------------
// for each v in G
//     vertex_set = {v} union neighbor(v)
//     coloring each vertex x in vertex_set if it havnt been colored yet.
//     with the condition that their colors are mutual exclusive.
// -------------------------------------------------
//
// When considering the memory optimization, such that memory complexity should be 
// less than O(N) for each thread.
// The Pseduo Code is as follows.  O(1)
// -------------------------------------------------
// @input W : Workspace array with fixed length. 
// @input WCapacity : e.g. WCapacity=32 or 64
// for each v in G
//     W <- empty
//     if v.color is empty:
//         W += {v}
//     for u in neighbor(v):
//         if u.color is empty:
//             W += {u}
//             if |W| == WCapacity:
//                 Clear_W_by_Coloring_MemOpt(W)
//     if |W| is not empty:
//         Clear_W_by_Coloring_MemOpt(W)
//  -------------------------------------------------
//
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseNetBinTC_MemOpt(
        vector<int>& vtxColors,             
        vector<int> const & const_ordered_Q,
        int const N,
        const vector<int>&vtxPtr,
        const vector<int>&vtxVal,
        vector<vector<int>>& Ws, //work load vertices
        int const WCapacity,     //
        ChronoDuration *pTime   
        ){
    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int>&W = Ws[tid];
        #pragma omp for
        for(int iv=0; iv<N; iv++){
            int Wsize = 0;
            const int v = const_ordered_Q[iv];
            const int vc = vtxColors[v];
            if(vc<0)
                W[Wsize++]=v;
            for(int iu=vtxPtr[v]; iu!=vtxPtr[v+1]; iu++){
                const int u = vtxVal[iu];
                const int uc = vtxColors[u];
                if(uc>=0) continue;
                W[Wsize++] = u;
                if(Wsize>=WCapacity)
                    do_D2GC_PopAllW_by_Coloring_MemOpt(v, W, Wsize, vtxColors, WCapacity, vtxPtr, vtxVal);
            }
            if(Wsize>0) 
                do_D2GC_PopAllW_by_Coloring_MemOpt(v, W, Wsize, vtxColors, WCapacity, vtxPtr, vtxVal);
        }// end of omp for
    }// end of omp parallel
    if(pTime) *pTime = steady_clock::now()-start;
}

// ============================================================================
// Author: Xin Cheng
// given a vector of vertex for coloring. greedy color the vertex with
// the constrains that O(1) memory is used.
//
// The pheduo code is as follows
// --------------------------------------------------------
// while(w.size!=0)
//     F<-do_D2GC_build_ForbiddenArray_F(offset++);
//     for all '1' in F:
//         u = W.pop();
//         color_u_based_on_the_offset and the one's position
// }while(V is not empty)
// --------------------------------------------------------
//
// ============================================================================
void SMPGCColoring::do_D2GC_PopAllW_by_Coloring_MemOpt(
        int const v,
        vector<int>& W, 
        int & Wsize,
        vector<int>& vtxColors,
        int const WCapacity, 
        vector<int> const& vtxPtr,
        vector<int> const& vtxVal
        ) {
    int Low = 0;
    while(Wsize>0){
        unsigned int F = do_D2GC_Build_MemOpt_ForbiddenArray(v, Low, vtxColors, WCapacity, vtxPtr, vtxVal);
        // for each avaible color, assign it to a vertex in W
        for(int i=0; i<WCapacity; i++){
            if( (F&(1<<i)) ){
                vtxColors[ W[--Wsize] ] = Low + i;
                if(Wsize<=0) 
                    break;
            }
        }
        Low += WCapacity;
    }
}


// ============================================================================
// Author: Xin Cheng
// build a forbidden color array of size WCapacity only for colors range in 
// [LOW, LOW+Capacity). Instead of using 'true'/'false' bool operation, here 
// uses '0'/'1' binary operation.
//
// for each value in F, 
// '0' means this color is been forbidden
// '1' means this color is available to use
//
// The Pseudo Code should be as follows:
// ---------------------------------------------------
// F <- set_to_all_1
// for node x in {v} union neighbor(v):
//     c_loc = v.color - LOW;
//     if c_lor in range [LOW, LOW+WCapacity):
//         F[c_loc] <- '0'
// return F
// ---------------------------------------------------
//
// ===========================================================================
unsigned int SMPGCColoring::do_D2GC_Build_MemOpt_ForbiddenArray(
        int const v,  
        int const Low, 
        vector<int> const & vtxColors,
        int const WCapacity,
        vector<int> const& vtxPtr,
        vector<int> const& vtxVal
        ) {
    unsigned int F = ~0;
    int const vc_loc = vtxColors[v]-Low;
    if(vc_loc>=0 && vc_loc<WCapacity)
        F &=  ~(1<<vc_loc);
    for(int iu=vtxPtr[v]; iu<vtxPtr[v+1]; iu++){
        int const uc_loc = vtxColors[vtxVal[iu]] - Low;
        if(uc_loc>=0 && vc_loc<WCapacity) 
            F &= ~(1<<uc_loc);
    }
    return F;
}



// ============================================================================
// Author: xin cheng
// Implements of the D2GC omp Speculative Iteration approaches'
// CheckConflicts phase
//      using Net-Based Algorithm
// 
// The idea is as follows
// ----------------------------
// for each v in G
//     make sure all vertex x in {v} union neighbor(v) has mutual different colors
//     if not, marke that vertex's color to be empty.
// ----------------------------
//
// Considering the Memory Optimization requirement, the pseudo code is as follows:
// -----------------------------
// for each v in G
//     do{
//         color_range = [ Low, Low+Capacity)
//         F = set_to_all_0
//         for u in {v} union neighbor(v):
//             if u.color in the color_range :
//                 if F[u.color-LOW] == '1':
//                     u.color <- empty
//                 else
//                    F[u.color-LOW] <- '0'
//      }while( max( "{v} union neighbor(v)"'s color ) < LOW+Capacity );
// re_calculate_QQ; 
// -----------------------------
//
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseNetBinCC_MemOpt(
        vector<vector<int>> & QQ,
        vector<int>& Qsizes,
        vector<int>& vtxColors,
        const vector<int> &vtxPtr,
        const vector<int> &vtxVal,
        vector<int> const &const_ordered_Q,
        int const N,
        int const Capacity,
        int * pNum, 
        ChronoDuration* pTime
        ){

    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    // sub-phase - check conflicts
    #pragma omp parallel
    {
        //const int tid = omp_get_thread_num();
        #pragma omp for
        for(int iv=0; iv<N; iv++){
            const int v = const_ordered_Q[iv];
            const int vc = vtxColors[v];
            bool bFirst=true;
            int max_c = vc;   //max color only calculated when bFirst flag is true
            int Low = 0;
            while(true){
                unsigned int F = 0;  //binary forbidden color array 
                const int vc_loc = vc - Low;
                if(vc_loc>=0 && vc_loc<Capacity) 
                    F |= 1<<vc_loc; 
                for(int iu=vtxPtr[v]; iu!=vtxPtr[v+1]; iu++) {
                    const int u  = vtxVal[iu];
                    const int uc = vtxColors[u];
                    if(bFirst && max_c<uc)  max_c = uc;
                    const int uc_loc = uc - Low;
                    if(uc_loc<0 || uc_loc>=Capacity)
                        continue;
                    if( F&(1<<uc_loc) )
                        vtxColors[u] = -1;
                    else
                        F |= 1<<uc_loc;
                }
                Low += Capacity;
                bFirst = false;
                if( max_c < Low+Capacity) 
                    break;
            }

        }
    }
    

    // sub-phase - collect and queue the uncolored vertex
    int num_uncolored_vtx_total=0;
    #pragma omp parallel reduction(+:num_uncolored_vtx_total)
    {
        const int tid = omp_get_thread_num();
        vector<int> &Q = QQ[tid];
        int Qsize = 0 ; 
        #pragma omp for
        for(int iu=0; iu<N; iu++){
            int const u = const_ordered_Q[iu];
            if( vtxColors[u] < 0)
                Q[Qsize++] = u;
        }
        Qsizes[tid] = Qsize;
        num_uncolored_vtx_total = Qsize ;
    }
    if(pTime) *pTime = steady_clock::now() - start;
    if(pNum) *pNum = num_uncolored_vtx_total;
}






