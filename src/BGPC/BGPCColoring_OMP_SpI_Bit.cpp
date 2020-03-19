/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "BGPCColoring.h"

#include <bitset>

using namespace std;
using namespace ColPack;
using namespace chrono;




// ============================================================================
// author: xin cheng
// Bipartite Graph Partial Coloring parallel implementation
// using OpenMP, multiple phases, without using atomic operation
// ============================================================================
int BGPCColoring::BGPC_OMP_SpeculativeIterative_Bit(
        int const side,             // L or R
        int nT,                     // number of threads
        int& colors,                // number of colors
        vector<int>& vtxColors,     // mapping vertex to color
        int const LOCAL_ORDER,      // Local Order
        int const NUM_NET_TC,       // Phase: Random Tentative Coloring
        int const NUM_NET_CC,       // Phase: Netbin Check Conflicts
        bool const b_SEQ_HC,        // using Seq-Coloring to solve the conflicts?
        int const nVerbose          // verbose level
        ) {
    
    if(b_SEQ_HC) {
        cerr<<"Sorry, this function does not support your arguments for now.\n"<<endl;
        exit(1);
    }
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    // for runtime debug
    ChronoDuration tim_PPT{0};
    ChronoDuration tim_LO{0};
    ChronoDuration tim_HC{0};
    vector<ChronoDuration>  vtim_Net_TC;
    vector<ChronoDuration>  vtim_Pln_TC;
    vector<ChronoDuration>  vtim_Net_CC;
    vector<ChronoDuration>  vtim_Pln_CC;
    vector<int> vnum_Net_CC;
    vector<int> vnum_Pln_CC;
    
    // variables
    const int          N             = (side==BGPC::L)?(GetLeftVertexCount()):(GetRightVertexCount()); 
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    const int          Ndst          = (side==BGPC::L)?(GetRightVertexCount()):(GetLeftVertexCount()); 
    const vector<int>& vtxVal        = GetEdges();

    const int          srcMaxDegree  = (side==BGPC::L)?GetMaximumLeftVertexDegree():GetMaximumRightVertexDegree();
    //const int          dstMaxDegree  = (side==BGPC::L)?GetMaximumRightVertexDegree():GetMaximumLeftVertexDegree();
    //int const          DDp1          = dstMaxDegree * srcMaxDegree + 1;
    //const int          BufSize       = DDp1<0?N:min(N,DDp1);
    
    // dynamic loaded variables    
    vector<int>const& const_queue_A        = get_ordered_queue_A_const();
    vector<int> const_queue_B;
    vector<mt19937> mts;
    vector<vector<int>> Workloads;  // workload space of the netbin coloring phase 
   
    // initialization
    colors=0;                       
    vtxColors.assign(N, -1);

    // allocate memory +16 for prevent false share on some 64 bit machines
    //vector<vector<int>> ForbiddenArrays(nT, vector<int>(BufSize+16,-1)); 
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16, -1)); 
    vector<int>         Qsizes(nT,0);
    if(NUM_NET_TC || NUM_NET_CC){
        const_queue_B.assign(Ndst, 0);
        for(int i=0; i<Ndst; i++)
            const_queue_B[i]=i;
    }
    if(NUM_NET_TC) {
        for(int i=0; i<nT; i++)
            mts.emplace_back(i);
    }
    if(NUM_NET_TC) {
        Workloads.assign(nT, vector<int>(srcMaxDegree,0));
    }


    // pre-partition
    init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &tim_PPT);
    
    // local ordering   
    Ordering_PrePartition_OMP(side, QQ, Qsizes, LOCAL_ORDER, &tim_LO);
    
    int niter=0;
    int num_uncolored_vertex = N;
    bool b_NET_TC=false;
    bool b_NET_CC=false;
    //cout<<"N"<<N<<endl;
    while(num_uncolored_vertex!=0){
        
        // phase Tentative Coloring
        if(niter<NUM_NET_TC) {
            vtim_Net_TC.emplace_back(0);
            do_BGPC_Phase_NetBinTC_Bit(vtxColors,const_queue_A, N, const_queue_B, Ndst, srcPtr, dstPtr, vtxVal, Workloads, mts, &vtim_Net_TC.back());
            //do_BGPC_Phase_RandomTC(vtxColors, QQ, Qsizes, srcPtr, mts, &vtim_Net_TC.back());
            //b_NET_TC = true;
        }
        else{
            if(b_NET_CC){
                init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, vtxColors, &tim_PPT);
            }
            vtim_Pln_TC.emplace_back(0);
            do_BGPC_Phase_PlainTC_Bit(vtxColors, QQ, Qsizes, srcPtr, dstPtr, vtxVal, &vtim_Pln_TC.back());
            b_NET_TC = false;
        }

        // phase Check Conflict
        if(niter<NUM_NET_CC) {
            vtim_Net_CC.emplace_back(0);
            do_BGPC_Phase_NetBinCC_Bit(vtxColors, const_queue_B, Ndst, srcPtr, dstPtr, vtxVal, &vtim_Net_CC.back());
            b_NET_CC=true;
        }
        else{
            if(b_NET_TC){
                init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &tim_PPT);
            }
            vtim_Pln_CC.emplace_back(0);
            do_BGPC_Phase_PlainCC(vtxColors, QQ, Qsizes, srcPtr, dstPtr, vtxVal, &num_uncolored_vertex, &vtim_Pln_CC.back());
            b_NET_CC=false;
            vnum_Pln_CC.emplace_back(num_uncolored_vertex);
        }
        //cout<<niter<<" : "<<num_uncolored_vertex<<endl;
        niter += 1;
    } //end while

    colors = calc_num_colors_from_vtx_colors(vtxColors);
    
    if(nVerbose>0){
        stringstream ss;
        ss<<"@BGPCSpIMemOpt(";
        ss<<(b_SEQ_HC?"3P":"MP");
        ss<<","<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        ss<<")NN(";
        ss<<NUM_NET_TC<<","<<NUM_NET_CC;
        ss<<")";
        for(int i=ss.str().size(); i<30; i++)
            ss<<"_";
        cout<<ss.str(); ss.str("");
        cout<<"\t#nT_c_T";
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration T{0};
        //T += tim_PPT;
        //T += tim_LO;
        for(auto const x : vtim_Net_TC) T+=x;
        for(auto const x : vtim_Pln_TC) T+=x;
        for(auto const x : vtim_Net_CC) T+=x;
        for(auto const x : vtim_Pln_CC) T+=x;
        cout<<"\t"<<T.count();
        if(nVerbose>1) {
            cout<<"\t#nIter_NTC_PTC_NCC_PCC\t"<<niter;
            { //if(NUM_NET_TC){
                ChronoDuration T{0}; 
                for(auto const x : vtim_Net_TC) T+=x; 
                cout<<"\t"<<T.count();
            }
            {
                ChronoDuration T{0}; 
                for(auto const x : vtim_Pln_TC) T+=x; 
                cout<<"\t"<<T.count();
            }
            
            { //if(NUM_NET_CC){
                ChronoDuration T{0}; 
                for(auto const x : vtim_Net_CC) T+=x; 
                cout<<"\t"<<T.count();
            }
            {
                ChronoDuration T{0}; 
                for(auto const x : vtim_Pln_CC) T+=x; 
                cout<<"\t"<<T.count();
            }
        }
        cout<<endl;
    }
    return true;   
}



// ============================================================================
// Author: Xin Cheng
// Tentative coloring
// ============================================================================
void BGPCColoring::do_BGPC_Phase_PlainTC_Bit(
        vector<int>& vtxColors,
        vector<vector<int>>const& QQ,
        vector<int>const& Qsizes,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        ChronoDuration *tim
        ){
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    #pragma omp parallel
    {
        int const tid = omp_get_thread_num();
        do_BGPC_Seq_Greedy_Bit(vtxColors, QQ[tid], Qsizes[tid], srcPtr, dstPtr, vtxVal);
        //do_BGPC_Seq_Greedy_Bit_with_set(vtxColors, QQ[tid], Qsizes[tid], srcPtr, dstPtr, vtxVal);
        //do_BGPC_Seq_Greedy_Bit_with_bitset(vtxColors, QQ[tid], Qsizes[tid], srcPtr, dstPtr, vtxVal);
    }
    if(tim) *tim = steady_clock::now() - start;
}


// ============================================================================
// Author: Xin Cheng
// Tentative coloring
// ============================================================================
void BGPCColoring::do_BGPC_Seq_Greedy_Bit(
        vector<int>& vtxColors,
        vector<int>const& Q,
        int const Qsize,
        vector<int>const &srcPtr,
        vector<int>const &dstPtr,
        vector<int>const &vtxVal,
        ChronoDuration* tim ){
    std::chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();

    #ifdef MaskWide
    #undef MaskWide
    #endif
    #define MaskWide 64 
    //unsigned int Mask=~0;
    uint64_t Mask = ~0;    // 1 : available  0 : forbidden
    for(int iu=0; iu<Qsize; iu++){
        int const u = Q[iu];
        int uc = vtxColors[u]; // -1;
        int offset = 0;
        while(uc==-1){
            Mask = ~0;
            int const LOW = (offset++)*MaskWide;
            for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
                auto const v = vtxVal[iv];
                for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                    auto const w = vtxVal[iw];
                    if(u==w) continue;
                    int wc= vtxColors[w] - LOW;
                    //if(wc<0) continue;
                    //wc-=LOW;   //wc -> wc_loc
                    if(wc>=0 && wc < MaskWide)
                        Mask &= (~(1<<wc));   // clear to 0
                }// end for w
            }// end for v

            // find first settled bit
            if(Mask!=0) 
                for(int i=0; i<MaskWide; i++) 
                    if(Mask&(1<<i)) {
                        vtxColors[u] = uc = LOW + i;
                        break;
                    }
        }// end while
    }// end for u
    if(tim) *tim = steady_clock::now() - start;
    #undef MaskWide
}

// ============================================================================
// Author: Xin Cheng
// Tentative coloring
// ============================================================================
void BGPCColoring::do_BGPC_Seq_Greedy_Bit_with_bitset(
        vector<int>& vtxColors,
        vector<int>const& Q,
        int const Qsize,
        vector<int>const &srcPtr,
        vector<int>const &dstPtr,
        vector<int>const &vtxVal,
        ChronoDuration* tim ){
    
    #ifdef MaskWidth
    #undef MaskWidth
    #endif
    #define MaskWidth 64

    std::chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    std::bitset<MaskWidth>Mask;
    //Mask.flip();
    
    //uint64_t Mask = ~0;    // 1 : available  0 : forbidden
    for(int iu=0; iu<Qsize; iu++){
        int const u = Q[iu];
        int uc = vtxColors[u]; // -1;
        while(uc==-1){
            Mask = ~0;
            int Low = 0;
            for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
                auto const v = vtxVal[iv];
                for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                    auto const w = vtxVal[iw];
                    if(u==w) continue;
                    int wc= vtxColors[w] - Low;
                    //if(wc<0) continue;
                    //wc-=LOW;   //wc -> wc_loc
                    if(wc>=0 && wc < MaskWidth)
                        Mask.set(wc); // Mask &= (~(1<<wc));   // clear to 0
                }// end for w
            }// end for v

            // find first settled bit
            if(Mask.count()!=MaskWidth) 
                for(int i=0; i<MaskWidth; i++) 
                    if(!Mask[i]){ 
                        vtxColors[u] = uc = Low + i;
                        break;
                    }
            Low += MaskWidth;
        }// end while
    }// end for u
    if(tim) *tim = steady_clock::now() - start;
    #undef MaskWidth
}

// ============================================================================
// Author: Xin Cheng
// Tentative coloring
// ============================================================================
void BGPCColoring::do_BGPC_Seq_Greedy_Bit_with_set(
        vector<int>& vtxColors,
        vector<int>const& Q,
        int const Qsize,
        vector<int>const &srcPtr,
        vector<int>const &dstPtr,
        vector<int>const &vtxVal,
        ChronoDuration* tim ){
    std::chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();

    #ifdef MaskWide
    #undef MaskWide
    #endif
    #define MaskWide 32
    unsigned int Mask=~0;
    //uint64_t Mask = ~0;    // 1 : available  0 : forbidden
    for(int iu=0; iu<Qsize; iu++){
        int const u = Q[iu];
        int uc = vtxColors[u]; // -1;
        if(uc!=-1) continue;
        int offset = 0;
        unordered_set<int> uniq_d2_nbs;
        for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++){
            int const v = vtxVal[iv];
            for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++){
                uniq_d2_nbs.insert(vtxVal[iw]);
            }
        }
        uniq_d2_nbs.erase(u);
        while(uc==-1){
            Mask = ~0;
            int const LOW = (offset++)*MaskWide;
            auto iter=uniq_d2_nbs.begin();
            while(iter!=uniq_d2_nbs.end()) {
                auto w =*iter;
                auto wc=vtxColors[w]-LOW;
                if(wc<MaskWide){
                    iter=uniq_d2_nbs.erase(iter);
                    if(wc>=0)
                        Mask &= (~(1<<wc));   // clear to 0
                }
                else
                    iter++;
            }
            // find first settled bit
            if(Mask!=0) 
                for(int i=0; i<MaskWide; i++) 
                    if(Mask&(1<<i)){ 
                        vtxColors[u] = uc = LOW + i;
                        break;
                    }
        }while(uc==-1);// end while
    }// end for u
    if(tim) *tim = steady_clock::now() - start;
    #undef MaskWide
}



// ============================================================================
// Author: Xin Cheng
// Tentative coloring
// ============================================================================
void BGPCColoring::do_BGPC_Phase_NetBinTC_Bit(
        vector<int>& vtxColors,
        vector<int>const& vAs,
        int const Asize,
        vector<int>const& vBs,
        int const Bsizes,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        vector<vector<int>>& Ws,               // size of at least Delt_B
        vector<mt19937>& mts,
        ChronoDuration *tim
        ){
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    #pragma omp parallel
    {
        int const tid = omp_get_thread_num();
        mt19937& mt = mts[tid];
        //vector<int>& W = Ws[tid];
        #pragma omp for
        for(int iv=0; iv<Bsizes; iv++){
            int const v = vBs[iv];
            int nbsizem1 = -dstPtr[v] + dstPtr[v+1] -1;
            for(int iu=dstPtr[v]; iu!=dstPtr[v+1]; iu++) {
                int const u = vtxVal[iu];
                if(vtxColors[u]==-1) {
                    //vtxColors[u]=nbsizem1--;
                    uniform_int_distribution<int> dist(0,nbsizem1);
                    vtxColors[u] = dist(mt);
                }
            }
        }
        /*       #pragma omp for
        for(int iv=0; iv<Bsizes; iv++){
            int const v = vBs[iv];  // =iv
            int nbsizem1 = -dstPtr[v] + dstPtr[v+1] -1;
            //for(int i=0; i<=nbsizem1; i++) W[i]=i;
            for(int iu=dstPtr[v]; iu!=dstPtr[v+1]; iu++) {
                int const u = vtxVal[iu];
                if(vtxColors[u]==-1){
                    uniform_int_distribution<int> dist(0,nbsizem1);
                    int r = dist(mt);
                    vtxColors[vtxVal[iu]] = r;//W[r];
                    //W[r]=W[nbsizem1--];
                }
            }
        }
        */
        
  /*      #pragma omp for
        for(int iu=0; iu<Asize; iu++){
            int const u = vAs[iu];
            int nbsizem1 = -srcPtr[u] + srcPtr[u+1] -1;
            if(vtxColors[u]==-1){
                uniform_int_distribution<int> dist(0, nbsizem1);
                vtxColors[u] = dist(mt);
            }
        }
        */

    }
    if(tim) *tim = steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng
// Phase check conflicts
// ============================================================================
void BGPCColoring::do_BGPC_Phase_NetBinCC_Bit(
        vector<int>& vtxColors,
        vector<int>const vBs,
        int const& Bsizes,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        ChronoDuration *tim
        ){
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    
    #define MaskWide 64 
    #pragma omp parallel
    {
        uint64_t F = ~0; 
        #pragma omp for
        for(int iv=0; iv<Bsizes; iv++){
            int const v = vBs[iv];  // =iv
            int Low = 0;
            int const nbsize = -dstPtr[v] + dstPtr[v+1];
            while( Low<nbsize ){
                F = ~0;
                for(int iu=dstPtr[v]; iu!=dstPtr[v+1]; iu++) {
                    int const u = vtxVal[iu];
                    int const uc= vtxColors[u]-Low;
                    if(uc>=0 && uc<MaskWide) {
                        if( F&(1<<uc) )
                            F&=(~(1<<uc));
                        else
                            vtxColors[u]=-1;
                    }
                }
                Low+=MaskWide;
            }
        }
    }
    if(tim) *tim = steady_clock::now() - start;
}


