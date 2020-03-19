/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "BGPCColoring.h"
using namespace std;
using namespace ColPack;
using namespace chrono;




// ============================================================================
// author: xin cheng
// Bipartite Graph Partial Coloring parallel implementation
// using OpenMP, multiple phases, without using atomic operation
// ============================================================================
int BGPCColoring::BGPC_OMP_SpeculativeIterative(
        int const side,             // L or R
        int nT,                     // number of threads
        int& colors,                // number of colors
        vector<int>& vtxColors,     // mapping vertex to color
        int const LOCAL_ORDER,      // Local Order
        int const NUM_RND_TC,       // Phase: Random Tentative Coloring
        int const NUM_NET_TC,       // Phase: Netbin Tentative Coloring
        int const NUM_NET_CC,       // Phase: Netbin Check Conflicts
        bool const b_SEQ_HC,        // using Seq-Coloring to solve the conflicts?
        int const nVerbose          // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    // for runtime debug
    ChronoDuration tim_PPT{0};
    ChronoDuration tim_LO{0};
    ChronoDuration tim_HC{0};
    vector<ChronoDuration>  vtim_Rnd_TC;
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
    const int          dstMaxDegree  = (side==BGPC::L)?GetMaximumRightVertexDegree():GetMaximumLeftVertexDegree();
    int const          DDp1          = dstMaxDegree * srcMaxDegree + 1;
    const int          BufSize       = DDp1<0?N:min(N,DDp1);
    
    // dynamic loaded variables    
    vector<int>const& const_queue_A        = get_ordered_queue_A_const();
    vector<int> const_queue_B;
    vector<mt19937> mts;
    vector<vector<int>> Workloads;  // workload space of the netbin coloring phase 
   
    // initialization
    colors=0;                       
    vtxColors.assign(N, -1);
    

    // allocate memory +16 for prevent false share on some 64 bit machines
    vector<vector<int>> ForbiddenArrays(nT, vector<int>(BufSize+16,-1)); 
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16, -1)); 
    vector<int>         Qsizes(nT,0);
    
    if(NUM_RND_TC)   for(int i=0; i<nT; i++) mts.emplace_back(i);
    if(NUM_NET_TC)   Workloads.assign(nT, vector<int>(dstMaxDegree+1, -1));
    if(NUM_NET_TC || NUM_NET_CC) { const_queue_B.assign(Ndst,0); for(int i=0; i<Ndst; i++) const_queue_B[i]=i; }
    

    // pre-partition
    init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &tim_PPT);
    
    // local ordering   
    Ordering_PrePartition_OMP(side, QQ, Qsizes, LOCAL_ORDER, &tim_LO);
    
    int niter=0;
    int num_uncolored_vertex = N;
    bool b_TC_just_used_Netbin=false;
    bool b_CC_just_used_Netbin=false;
    while(num_uncolored_vertex!=0){

        // phase Tentative Coloring
        if(niter<NUM_RND_TC){
            if(b_CC_just_used_Netbin) {
                ChronoDuration t{.0}; init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, vtxColors, &t); vtim_Net_CC.back()+=t;
            }
            vtim_Rnd_TC.emplace_back(0);
            do_BGPC_Phase_RandomTC(vtxColors, QQ, Qsizes, srcPtr, mts, &vtim_Rnd_TC.back());
            b_TC_just_used_Netbin = false;
        }
        else if(niter<NUM_NET_TC+NUM_RND_TC) {
            vtim_Net_TC.emplace_back(0);
            do_BGPC_Phase_NetbinTC(vtxColors, const_queue_B, Ndst, dstPtr, vtxVal, Workloads, ForbiddenArrays, BufSize, &vtim_Net_TC.back());
            b_TC_just_used_Netbin = true;

        }
        else
        {
            if(b_CC_just_used_Netbin){
                ChronoDuration t{.0}; init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, vtxColors, &t); vtim_Net_CC.back()+=t;
            }
            vtim_Pln_TC.emplace_back(0);
            do_BGPC_Phase_PlainTC(vtxColors, QQ, Qsizes, srcPtr, dstPtr, vtxVal, ForbiddenArrays, BufSize, &vtim_Pln_TC.back());
            b_TC_just_used_Netbin = false;
        }


        // phase Check Conflict
        if(niter<NUM_NET_CC){
            vtim_Net_CC.emplace_back(0);
            do_BGPC_Phase_NetbinCC(vtxColors, const_queue_B, Ndst, dstPtr, vtxVal, ForbiddenArrays, BufSize, &num_uncolored_vertex, &vtim_Net_CC.back());
            vnum_Net_CC.emplace_back(num_uncolored_vertex);
            b_CC_just_used_Netbin = true;
        }
        else{
            if(b_TC_just_used_Netbin){
                ChronoDuration t{.0}; 
                init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &t); // init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, vtxColors, &t); 
                vtim_Net_TC.back()+=t;
            }
            vtim_Pln_CC.emplace_back(0);
            do_BGPC_Phase_PlainCC(vtxColors, QQ, Qsizes, srcPtr, dstPtr, vtxVal, &num_uncolored_vertex, &vtim_Pln_CC.back());
            vnum_Pln_CC.emplace_back(num_uncolored_vertex);
            b_CC_just_used_Netbin = false;
        }

        // phase seq Handle Conflicts
        if(b_SEQ_HC){
            if(b_CC_just_used_Netbin){
                ChronoDuration t{.0}; init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, vtxColors, &t); vtim_Net_CC.back()+=t;
            }
            auto start = steady_clock::now();
            ForbiddenArrays[0].assign(BufSize,-1);
            for(int i=0; i<nT; i++)
                do_BGPC_Seq_Greedy_NoInitF(vtxColors, QQ[i], Qsizes[i], srcPtr, dstPtr, vtxVal, ForbiddenArrays[0], BufSize);
            num_uncolored_vertex = 0;
            tim_HC = steady_clock::now() - start;
        }
        niter += 1;
    } //end while

    colors = calc_num_colors_from_vtx_colors(vtxColors);
    
    if(nVerbose>0){
        stringstream ss;
        ss<<"@BGPCSpI(";
        ss<<(b_SEQ_HC?"3P":"MP");
        ss<<","<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        ss<<")RNN(";
        ss<<NUM_RND_TC<<","<<NUM_NET_TC<<","<<NUM_NET_CC;
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
        for(auto const x : vtim_Rnd_TC) T+=x;
        for(auto const x : vtim_Net_TC) T+=x;
        for(auto const x : vtim_Pln_TC) T+=x;
        for(auto const x : vtim_Net_CC) T+=x;
        for(auto const x : vtim_Pln_CC) T+=x;
        cout<<"\t"<<T.count();
        if(nVerbose>1) {
            cout<<"\t#nIter_RTC_NTC_PTC_NCC_PCC\t"<<niter;
            { //if(NUM_RND_TC){
                ChronoDuration T{0}; 
                for(auto const x : vtim_Rnd_TC) T+=x; 
                cout<<"\t"<<T.count();
            }
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
void BGPCColoring::do_BGPC_Phase_PlainTC(
        vector<int>& vtxColors,
        vector<vector<int>>const& QQ,
        vector<int>const& Qsizes,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        vector<vector<int>>& Fs,
        int const BufSize,
        ChronoDuration *tim
        ){
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    #pragma omp parallel
    {
        int const tid = omp_get_thread_num();
        do_BGPC_Seq_Greedy(vtxColors, QQ[tid], Qsizes[tid], srcPtr, dstPtr, vtxVal, Fs[tid], BufSize);
    }
    if(tim) *tim = steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng
// Check conflicts
// ============================================================================
void BGPCColoring::do_BGPC_Phase_PlainCC(
        vector<int>& vtxColors,
        vector<vector<int>>& QQ,
        vector<int>& Qsizes,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        int* num,
        ChronoDuration* tim
        ){
    chrono::steady_clock::time_point start;
    if(tim)  start = steady_clock::now();
    int num_uncolored_vertex=0;
    #pragma omp parallel reduction(+:num_uncolored_vertex) 
    {
        int const tid = omp_get_thread_num();
        vector<int>& Q = QQ[tid];
        int const Qsize = Qsizes[tid];
        int cfq_size=0;
        for(int iu=0; iu<Qsize; iu++){  //u v w
            int const u = Q[iu]; 
            int const uc= vtxColors[u];
            bool b_uis_bad = false;
            for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
                if(b_uis_bad) break;
                int const v = vtxVal[iv];
                for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                    int const w = vtxVal[iw];
                    if(w >= u) continue;
                    if(vtxColors[w] == uc){
                        vtxColors[u] = -1;
                        b_uis_bad = true;
                        Q[cfq_size++] = u;
                        break;
                    }
                }// end w
            }// end v
        }// end u
        Qsizes[tid]=cfq_size;
        num_uncolored_vertex=cfq_size;
    }// end omp parallel
    if(num) *num = num_uncolored_vertex;
    if(tim) *tim = steady_clock::now() - start;
}


// ============================================================================
// Author: Xin Cheng
// implement random phase TC
// ============================================================================
void BGPCColoring::do_BGPC_Phase_RandomTC(
        vector<int>& vtxColors,
        vector<vector<int>>const& QQ,
        vector<int>const& Qsizes,
        vector<int>const& srcPtr,
        vector<mt19937>& mts,
        ChronoDuration* tim
        ){
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    #pragma omp parallel
    {
        int const tid = omp_get_thread_num();
        mt19937& mt = mts[tid];
        for(int iu=0; iu<Qsizes[tid]; iu++){
            int const u = QQ[tid][iu];
            uniform_int_distribution<int> dist_local(0, srcPtr[u+1]-srcPtr[u]);
            vtxColors[u] = dist_local(mt);
        }
    }
    if(tim) *tim = steady_clock::now() - start;
}



// ============================================================================
// Author: Xin Cheng
// implement NetBin Phase TC
// Alg 8 of the paper 'greedy is good ....' by Mustafa Kemal Tas ̧,Kamer Kaya and Erik Saule
//
// pseudo alg:
// -------------------
// for each v in B do in parallel:
//     F <- thread private empty
//     W <- thread private empty
//     for each u in neighbor(v):
//         if c[u]!=-1 and c[u]!in F:
//             F += {c[u]}
//         else
//             W += {u}
//     col <- |vtxs(v)|-1
//     for each u in W do
//         while col in F do
//             col -= 1
//         c[u] = col
//         col <- col - 1
// ------------------- 
// ============================================================================
void BGPCColoring::do_BGPC_Phase_NetbinTC(
        vector<int>& vtxColors,
        vector<int>const& vtx_set_B,
        int const Ndst,                     // |B|
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        vector<vector<int>>& Ws,
        vector<vector<int>>& Fs,
        int const BufSize,
        ChronoDuration* tim
        ) {
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    #pragma omp parallel
    {
        int const tid = omp_get_thread_num();
        vector<int>& F = Fs[tid];
        F.assign(BufSize,-1);
        vector<int>& W = Ws[tid];
        int Wsize=0;
        #pragma omp for
        for(int iv=0; iv<Ndst; iv++){
            int const v = vtx_set_B[iv];
            Wsize=0;
            for(int iu=dstPtr[v]; iu!=dstPtr[v+1]; iu++) {
                int const u = vtxVal[iu];
                int const uc= vtxColors[u];
                if( uc>=0 && F[uc]!=v )
                    F[uc] = v;
                else
                    W[Wsize++] = u;
            }
            int col = dstPtr[v+1]-dstPtr[v] - 1;
            for(int iu=0; iu<Wsize; iu++){
                int const u = W[iu];
                while( F[col]==v ){
                    col -= 1;
                }
                vtxColors[u] = col--;
            }
        }
    }
    if(tim) *tim = steady_clock::now() - start;
}




// ============================================================================
// Author: Xin Cheng
// Implement NetBin Phase CC
// Alg 8 of the paper 'greedy is good ....' by Mustafa Kemal Tas ̧,Kamer Kaya and Erik Saule
// 
// Pseudo alg:
// --------------------------
// for each v in B in parallel do
//     F <- empty
//     for each u in neighbor(v) do
//         if c[u] != -1
//             if c[u] in F:
//                 c[u] <- -1
//             else:
//                 F += {c[u]}
// --------------------------
//
// ============================================================================
void BGPCColoring::do_BGPC_Phase_NetbinCC(
        vector<int>& vtxColors,                    // vertex Colors
        vector<int>const& vtx_set_B,
        int const Ndst,
        vector<int>const& dstPtr,                // Bipartite graph in CSR format
        vector<int>const& vtxVal,                // Bipartite graph in CSR format
        vector<vector<int>>& Fs,
        int const BufSize,
        int * num,
        ChronoDuration *tim
        ){
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    
    int num_conflicts = 0;
    #pragma omp parallel reduction(+:num_conflicts)
    {
        int const tid = omp_get_thread_num();
        vector<int> &F = Fs[tid];
        F.assign(BufSize, -1);
        int n_conf = 0;
        #pragma omp for
        for(int iv=0; iv<Ndst; iv++){
            int const v = vtx_set_B[iv];
            for(int iu=dstPtr[v]; iu!=dstPtr[v+1]; iu++) {
                int const u = vtxVal[iu];
                int const uc = vtxColors[u];
                if(uc>=0){
                    if( F[uc] == v ){
                        vtxColors[u] = -1; n_conf++;  }
                    else
                        F[uc] = v;
                }
            }
        }
        num_conflicts = n_conf;
    }
   
    if(tim) *tim = steady_clock::now() - start;
    if(num) *num = num_conflicts;
    return;
}



// ============================================================================
// author: xin cheng
// Bipartite Graph Partial Coloring parallel implementation
// using OpenMP, multiple phases, may using atomic operation
// ============================================================================
int BGPCColoring::BGPC_OMP_SpeculativeIterative_withoutPrePartition(
        int const side,             // L or R
        int nT,                     // number of threads
        int& colors,                // number of colors
        vector<int>& vtxColors,     // mapping vertex to color
        int const LOCAL_ORDER,      // Local Order
        bool const b_SEQ_HC,        // using Seq-Coloring to solve the conflicts?
        int const nVerbose          // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    // for runtime debug
    //ChronoDuration tim_PPT{0};
    ChronoDuration tim_LO{0};
    ChronoDuration tim_HC{0};
    //vector<ChronoDuration>  vtim_Rnd_TC;
    //vector<ChronoDuration>  vtim_Net_TC;
    vector<ChronoDuration>  vtim_Pln_TC;
    //vector<ChronoDuration>  vtim_Net_CC;
    vector<ChronoDuration>  vtim_Pln_CC;
    //vector<int> vnum_Net_CC;
    //vector<int> vnum_Pln_CC;
    
    // variables
    const int          N             = (side==BGPC::L)?(GetLeftVertexCount()):(GetRightVertexCount()); 
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    //const int          Ndst          = (side==BGPC::L)?(GetRightVertexCount()):(GetLeftVertexCount()); 
    const vector<int>& vtxVal        = GetEdges();

    const int          srcMaxDegree  = (side==BGPC::L)?GetMaximumLeftVertexDegree():GetMaximumRightVertexDegree();
    const int          dstMaxDegree  = (side==BGPC::L)?GetMaximumRightVertexDegree():GetMaximumLeftVertexDegree();
    int const          DDp1          = dstMaxDegree * srcMaxDegree + 1;
    const int          BufSize       = DDp1<0?N:min(N,DDp1);
    
    // dynamic loaded variables    
    //vector<int>const& const_queue_A        = get_ordered_queue_A_const();
    //vector<int> const_queue_B;
    //vector<mt19937> mts;
    //vector<vector<int>> Workloads;  // workload space of the netbin coloring phase 
   
    // initialization
    colors=0;                       
    vtxColors.assign(N, -1);
    

    // allocate memory +16 for prevent false share on some 64 bit machines
    vector<vector<int>> ForbiddenArrays(nT, vector<int>(BufSize+16,-1)); 
    //vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16, -1)); 
    //vector<int>         Qsizes(nT,0);
    
    vector<int>  Q(N+1+16,0); // Queue of uncolored vertex
    vector<int>  Q2(N+1+16, 0);
    int Qsize=N;      // Qsize

    //if(NUM_RND_TC)   for(int i=0; i<nT; i++) mts.emplace_back(i); //random seed
    //if(NUM_NET_TC)   Workloads.assign(nT, vector<int>(dstMaxDegree+1, -1));
    //if(NUM_NET_TC || NUM_NET_CC) { const_queue_B.assign(Ndst,0); for(int i=0; i<Ndst; i++) const_queue_B[i]=i; }
    

    // pre-partition
    //init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &tim_PPT);
    
    // local ordering   
    //Ordering_PrePartition_OMP(side, QQ, Qsizes, LOCAL_ORDER, &tim_LO);
    Q.assign(N, 0);
    for(int i=0; i<N; i++){ 
        Q[i]=i;
    }  // Nature order, or knonw as FF
    
            
    int niter=0;
    //bool b_TC_just_used_Netbin=false;
    //bool b_CC_just_used_Netbin=false;
    while(Qsize!=0){
        vtim_Pln_TC.emplace_back(0);
        vtim_Pln_CC.emplace_back(0);
        
        // phase Tentative Coloring
        do_BGPC_Phase_PlainTC_NoPTT(vtxColors, Q, Qsize, srcPtr, dstPtr, vtxVal, ForbiddenArrays, BufSize, &vtim_Pln_TC.back());

        // phase Check Conflict
        do_BGPC_Phase_PlainCC_NoPTT(vtxColors, Q, Qsize, srcPtr, dstPtr, vtxVal, Q2, &vtim_Pln_CC.back());

        // phase seq Handle Conflicts
        if(b_SEQ_HC){
            auto start = steady_clock::now();
            ForbiddenArrays[0].assign(BufSize,-1);
            do_BGPC_Seq_Greedy_NoInitF(vtxColors, Q, Qsize, srcPtr, dstPtr, vtxVal, ForbiddenArrays[0], BufSize);
            Qsize = 0;
            tim_HC = steady_clock::now() - start;
        }
        niter += 1;
    } //end while

    colors = calc_num_colors_from_vtx_colors(vtxColors);
    
    if(nVerbose>0){
        stringstream ss;
        ss<<"@BGPCSpINoPrePTT(";
        ss<<(b_SEQ_HC?"3P":"MP");
        ss<<","<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        //ss<<")RNN(";
        //ss<<"(None, None, None)";//ss<<NUM_RND_TC<<","<<NUM_NET_TC<<","<<NUM_NET_CC;
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
        //for(auto const x : vtim_Rnd_TC) T+=x;
        //for(auto const x : vtim_Net_TC) T+=x;
        for(auto const x : vtim_Pln_TC) T+=x;
        //for(auto const x : vtim_Net_CC) T+=x;
        for(auto const x : vtim_Pln_CC) T+=x;
        cout<<"\t"<<T.count();
        //if(nVerbose>1) {
        //    cout<<"\t#nIter_RTC_NTC_PTC_NCC_PCC\t"<<niter;
        //    { //if(NUM_RND_TC){
        //        ChronoDuration T{0}; 
        //        for(auto const x : vtim_Rnd_TC) T+=x; 
        //        cout<<"\t"<<T.count();
        //    }
        //    { //if(NUM_NET_TC){
        //        ChronoDuration T{0}; 
        //        for(auto const x : vtim_Net_TC) T+=x; 
         //       cout<<"\t"<<T.count();
        //    }
        //    {
        //        ChronoDuration T{0}; 
        //        for(auto const x : vtim_Pln_TC) T+=x; 
        //        cout<<"\t"<<T.count();
        //    }
        //    
        //    { //if(NUM_NET_CC){
        //        ChronoDuration T{0}; 
        //        for(auto const x : vtim_Net_CC) T+=x; 
        //        cout<<"\t"<<T.count();
        //    }
         //   {
        //        ChronoDuration T{0}; 
        //        for(auto const x : vtim_Pln_CC) T+=x; 
        //        cout<<"\t"<<T.count();
        //    }
        //}
        cout<<endl;
    }
    return true;   
}




// ============================================================================
// author: xin cheng
// Phase Tentative Coloring of the Bipartite Graph Partial Coloring parallel implementation
// using OpenMP, multiple phases, without using atomic operation
// No Pre Partition
// ============================================================================
void BGPCColoring::do_BGPC_Phase_PlainTC_NoPTT(
        vector<int>& vtxColors, 
        vector<int>const& Q,
        int const Qsize, 
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal, 
        vector<vector<int>>& Fs, 
        int const BufSize, 
        ChronoDuration* tim ){
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    #pragma omp parallel
    {
        int const tid = omp_get_thread_num();
        vector<int>& F = Fs[tid];
        int const N = Qsize;
        #pragma omp for
        for(int i=0; i<N; i++){
            int u = Q[i];
            for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++){
                const auto v = vtxVal[iv];
                for(int iw=srcPtr[v]; iw!=srcPtr[v+1]; iw++){
                    const auto w = vtxVal[iw];
                    if(w==u) continue;
                    const auto wc = vtxColors[w];
                    if(wc<0) continue;
                    F[wc]=u;
                }
            }
            int c=0;
            for(; c<BufSize; c++){
                if(F[c]!=u)
                    break;
            }
            vtxColors[u]=c;
        }// end of omp for
    }// end of omp parallel
    if(tim) *tim = steady_clock::now() - start;
}

// ============================================================================
// author: xin cheng
// Phase Check Conflicts of the Bipartite Graph Partial Coloring parallel 
// implementation
// using OpenMP, multiple phases, using atomic operation
// No Pre Partition
// ============================================================================
void BGPCColoring::do_BGPC_Phase_PlainCC_NoPTT(
        vector<int>& vtxColors, 
        vector<int>& Q, 
        int& Qsize, 
        vector<int>const& srcPtr, 
        vector<int>const& dstPtr, 
        vector<int>const& vtxVal, 
        vector<int>& Q2,
        ChronoDuration* tim
        ){
    chrono::steady_clock::time_point start;
    if(tim)  start = steady_clock::now();
    int QTail=0;
    #pragma omp parallel 
    {
        //int const tid = omp_get_thread_num();
        int const N = Qsize;
        #pragma omp for
        for(int iu=0; iu<N; iu++) {
            int const u = Q[iu];
            int const uc =vtxColors[u];
            bool b_uis_bad = false;
            for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++){
                if(b_uis_bad) break;
                int const v = vtxVal[iv];
                for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++){
                    int const w = vtxVal[iw];
                    if(w>=u) continue;
                    if(vtxColors[w] == uc){
                        vtxColors[u]=-1;
                        b_uis_bad = true;
                        int pos = __sync_fetch_and_add(&QTail, 1);           
                        Q2[pos]=u;
                        break;
                    }
                }
            }
        }// end omp for
    }// end omp parallel
    Q.swap(Q2);
    Qsize = QTail;
    if(tim) *tim = steady_clock::now() - start;
}





