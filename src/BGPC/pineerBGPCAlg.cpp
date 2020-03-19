/*************************************************************************
    > File Name: pineerBGPCAlg.cpp
    > Author: xc
    > Descriptions: 
    > Created Time: Thu 31 Oct 2019 04:27:52 PM EDT
 ************************************************************************/

#include "BGPCColoring.h"
using namespace std;
using namespace ColPack;
using namespace chrono;


// ============================================================================
// author: xin cheng
// Bipartite Graph Partial Coloring parallel implementation
// using OpenMP, multiple phases, without using atomic operation
//
// Try idea conflict not by id but by degree
// ============================================================================
int BGPCColoring::BGPC_OMP_SpeculativeIterative_ccbydeg(
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
    
    vector<int> Gdeg(N,0);
    //vector<int> G2deg(N,0);

    chrono::steady_clock::time_point deg_start = steady_clock::now();
    #pragma omp parallel for
    for(int u=0; u<N; u++){
        Gdeg[u] = srcPtr[u+1]-srcPtr[u];
        //unordered_set<int> mset;
        //for(int iv=srcPtr[u], ivEnd=srcPtr[u+1]; iv!=ivEnd; iv++){
        //    int v = vtxVal[iv];
        //    for(int iw=dstPtr[v], iwEnd=dstPtr[v+1]; iw!=iwEnd; iw++){
        //        int w = vtxVal[iw];
        //        if(u==w) continue;
        //        mset.insert(w);
        //    }
        //}
        //G2deg[u] = mset.size();
    }
    chrono::steady_clock::time_point deg_end = steady_clock::now();
    ChronoDuration deg_t = deg_end - deg_start;

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
            do_BGPC_Phase_PlainCC_ccbydeg(vtxColors, QQ, Qsizes, srcPtr, dstPtr, vtxVal, &num_uncolored_vertex, &vtim_Pln_CC.back(), Gdeg);
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
        ss<<"@BGPCSpIccbydeg(";
        ss<<(b_SEQ_HC?"3P":"MP");
        ss<<","<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        ss<<")RNN(";
        ss<<NUM_RND_TC<<","<<NUM_NET_TC<<","<<NUM_NET_CC;
        ss<<")";
        for(int i=ss.str().size(); i<30; i++)
            ss<<"_";
        cout<<ss.str(); ss.str("");
        cout<<"\t#nT_c_T_preT";
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
        cout<<"\t"<<deg_t.count();
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
// Check conflicts
//
// try idea cc by deg instead of id
// ============================================================================
void BGPCColoring::do_BGPC_Phase_PlainCC_ccbydeg(
        vector<int>& vtxColors,
        vector<vector<int>>& QQ,
        vector<int>& Qsizes,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        int* num,
        ChronoDuration* tim
        ,vector<int>const &Gdeg
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
            int const ud= Gdeg[u]; //srcPtr[u+1] - srcPtr[u];
            int const uc= vtxColors[u];
            bool b_uis_bad = false;
            for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
                if(b_uis_bad) break;
                int const v = vtxVal[iv];
                for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                    int const w = vtxVal[iw];
                    if(w==u) continue;
                    int const wd= Gdeg[w]; //srcPtr[w+1] - srcPtr[w];   //
                    if(ud>wd) continue;
                    if(ud==wd && u<w) continue;
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


// show difference between two subgraphs G1{vc1, QQ1, Qsizes1} G2{vc2,QQ2,Qsizes2} from origianl graph G(srcPtr, dstPtr, vtxVal)
// 1. show orig graph info
// 2. show G1 info
// 3. show G2 info
// 4. show G1,G2 common info
// 5. show G1,G2 diff info
void BGPCColoring::do_BGPC_compare_subgraph(
        vector<int>const& srcPtr, vector<int>const& dstPtr,  vector<int>const& vtxVal,
        vector<int>const&vtxColors1, vector<vector<int>>const&QQ1,  vector<int>const&Qsizes1,   
        vector<int>const&vtxColors2, vector<vector<int>>const&QQ2,  vector<int>const&Qsizes2, 
        vector<int>const&vtxColors3, vector<vector<int>>const&QQ3,  vector<int>const&Qsizes3,
        vector<int>const&vtxColors4, vector<vector<int>>const&QQ4,  vector<int>const&Qsizes4,
        vector<int>const&Gdeg, vector<int>const&G2deg, vector<int>const& GRnd
        ){

    // 1. 
    int NL = ((signed int)srcPtr.size())-1;
    int NR = ((signed int)dstPtr.size())-1;
    int E  = ((signed int)vtxVal.size())/2;



    set<int> G1v, G2v, G3v, G4v;
    int const nT = QQ1.size();
    for(int i=0; i<nT; i++){
        for(int j=0, jEnd=Qsizes1[i]; j!=jEnd; j++)
            G1v.insert(QQ1[i][j]);
        for(int j=0, jEnd=Qsizes2[i]; j!=jEnd; j++)
            G2v.insert(QQ2[i][j]);
        for(int j=0, jEnd=Qsizes3[i]; j!=jEnd; j++)
            G3v.insert(QQ3[i][j]);
        for(int j=0, jEnd=Qsizes4[i]; j!=jEnd; j++)
            G4v.insert(QQ4[i][j]);
    }
    
    // 2,3,4,5
    int N1L = G1v.size();
    int N2L = G2v.size();
    int N3L = G3v.size();
    int N4L = G4v.size();
    
    int E1  = 0;
    int E2  = 0;
    int E3  = 0;
    int E4  = 0;

    int cmm_v2 = 0;
    int cmm_e2 = 0;
    
    int cmm_v3 = 0;
    int cmm_e3 = 0;

    int cmm_v4 = 0;
    int cmm_e4 = 0;
    
    for(auto x : G1v) E1+= Gdeg[x];// (srcPtr[x+1]-srcPtr[x]);
    for(auto x : G2v) E2+= Gdeg[x];// (srcPtr[x+1]-srcPtr[x]);
    for(auto x : G3v) E3+= Gdeg[x];// (srcPtr[x+1]-srcPtr[x]);
    for(auto x : G4v) E4+= Gdeg[x];// (srcPtr[x+1]-srcPtr[x]);

    for(auto x : G1v){
        if(G2v.count(x)) {
            cmm_v2++;
            cmm_e2+= Gdeg[x];
        }
        if(G3v.count(x)) {
            cmm_v3++;
            cmm_e3+= Gdeg[x];
        }
        if(G4v.count(x)) {
            cmm_v4++;
            cmm_e4+= Gdeg[x];
        }

    }
    

    cout<<"@nT,"<<nT;
    cout<<",@OG,"<<NL<<","<<NR<<","<<E;
    cout<<",@Gid,"<<N1L<<","<<E1;
    cout<<",@Gdg,"<<N2L<<","<<E2;
    cout<<",@G2d,"<<N3L<<","<<E3;
    cout<<",@Grd,"<<N4L<<","<<E4;
    cout<<",@cm2,"<<cmm_v2<<","<<cmm_e2;
    cout<<",@cm3,"<<cmm_v3<<","<<cmm_e3;
    cout<<",@cm4,"<<cmm_v4<<","<<cmm_e4;
    cout<<",@uqDeg2,"<<N2L-cmm_v2<<","<<E2-cmm_e2;
    cout<<",@uqDeg3,"<<N3L-cmm_v3<<","<<E3-cmm_e3;
    cout<<",@uqDeg4,"<<N4L-cmm_v4<<","<<E4-cmm_e4;
    cout<<endl;
}


// ============================================================================
// author: xin cheng
// Bipartite Graph Partial Coloring parallel implementation
// using OpenMP, multiple phases, without using atomic operation
//
// Try idea conflict not by id but by degree
// ============================================================================
int BGPCColoring::BGPC_OMP_SpeculativeIterative_ccbydeg_debug_show_diff(
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

    vector<int> Gdeg(N,0), G2deg(N,0);
    vector<int> GRnd(N,0);
    //chrono::steady_clock::time_point deg_start = steady_clock::now();
    #pragma omp parallel for
    for(int u=0; u<N; u++){
        GRnd[u]=u;
        Gdeg[u] = srcPtr[u+1]-srcPtr[u];
        unordered_set<int> mset;
        for(int iv=srcPtr[u], ivEnd=srcPtr[u+1]; iv!=ivEnd; iv++){
            int v = vtxVal[iv];
            for(int iw=dstPtr[v], iwEnd=dstPtr[v+1]; iw!=iwEnd; iw++){
                int w = vtxVal[iw];
                if(u==w) continue;
                mset.insert(w);
            }
        }
        G2deg[u] = mset.size();
    }

    {
        int const Nm1 = N-1;
        // make GRnd random shuffle
        for(int i=0; i<Nm1; i++){
            uniform_int_distribution<int> dist(i, Nm1);
            swap(GRnd[i], GRnd[dist(m_mt)]);
        }
    }

    //chrono::steady_clock::time_point deg_end = steady_clock::now();
    //ChronoDuration deg_t = tmp_end -tmp_start;



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
            
            // 1 bkup vtxColors, QQ, Qsizes, num_uncolored_vertex and vtx_Pln_CC
            // 2 call do_BGPC_Phase_Plain_CC()
            // 3 call do_BGPC_Phase_Plain_ccbydeg()
            // 4 call do_BGPC_Phase_Plain_ccbyG2deg()
            // 5 ccbydeg by Random weight 
            // calc difference
            
            // 1.
            vector<int> ccGDeg_vtxColors(vtxColors);
            vector<int> ccGDeg_Qsizes(Qsizes);
            vector<vector<int>> ccGDeg_QQ;
            for(auto & Q : QQ) {ccGDeg_QQ.push_back({});  ccGDeg_QQ.back().assign(Q.begin(), Q.end()); }
            

            vector<int> ccG2Deg_vtxColors(vtxColors);
            vector<int> ccG2Deg_Qsizes(Qsizes);
            vector<vector<int>> ccG2Deg_QQ;
            for(auto & Q : QQ) {ccG2Deg_QQ.push_back({});  ccG2Deg_QQ.back().assign(Q.begin(), Q.end()); }
            
            vector<int> ccGRnd_vtxColors(vtxColors);
            vector<int> ccGRnd_Qsizes(Qsizes);
            vector<vector<int>> ccGRnd_QQ;
            for(auto & Q : QQ) {ccGRnd_QQ.push_back({});  ccGRnd_QQ.back().assign(Q.begin(), Q.end()); }

            // 2. 3.
            do_BGPC_Phase_PlainCC          (      vtxColors,       QQ,       Qsizes,       srcPtr, dstPtr, vtxVal, nullptr, nullptr);
            do_BGPC_Phase_PlainCC_ccbydeg  (ccGDeg_vtxColors,  ccGDeg_QQ,  ccGDeg_Qsizes,  srcPtr, dstPtr, vtxVal, nullptr, nullptr, Gdeg);
            do_BGPC_Phase_PlainCC_ccbyG2deg(ccG2Deg_vtxColors, ccG2Deg_QQ, ccG2Deg_Qsizes, srcPtr, dstPtr, vtxVal, nullptr, nullptr, G2deg);
            
            do_BGPC_Phase_PlainCC_ccbyG2deg(ccGRnd_vtxColors, ccGRnd_QQ, ccGRnd_Qsizes, srcPtr, dstPtr, vtxVal, nullptr, nullptr, GRnd);


            // 4
            do_BGPC_compare_subgraph(
                    srcPtr, dstPtr, vtxVal ,    
                    vtxColors, QQ, Qsizes,  ccGDeg_vtxColors, ccGDeg_QQ, ccGDeg_Qsizes,   
                    ccG2Deg_vtxColors, ccG2Deg_QQ, ccG2Deg_Qsizes ,
                    ccGRnd_vtxColors, ccGRnd_QQ, ccGRnd_Qsizes ,
                    Gdeg, G2deg, GRnd );
            return true;
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
        ss<<"@BGPCSpIccbydeg(";
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
// author: xin cheng
// Bipartite Graph Partial Coloring parallel implementation
// using OpenMP, multiple phases, without using atomic operation
//
// Try idea conflict not by id, not by degree, but by G^2 degree
// ============================================================================
int BGPCColoring::BGPC_OMP_SpeculativeIterative_ccbyG2deg(
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
   

    // pre calculate G deg and G2 deg
    //vector<int> Gdeg(N,0);
    vector<int> G2deg(N,0);
 
    
    chrono::steady_clock::time_point deg_start = steady_clock::now();
    #pragma omp parallel for
    for(int u=0; u<N; u++){
        //Gdeg[u] = srcPtr[u+1]-srcPtr[u];
        unordered_set<int> mset;
        for(int iv=srcPtr[u], ivEnd=srcPtr[u+1]; iv!=ivEnd; iv++){
            int v = vtxVal[iv];
            for(int iw=dstPtr[v], iwEnd=dstPtr[v+1]; iw!=iwEnd; iw++){
                int w = vtxVal[iw];
                if(u==w) continue;
                mset.insert(w);
            }
        }
        G2deg[u] = mset.size();
    }
    chrono::steady_clock::time_point deg_end = steady_clock::now();
    ChronoDuration deg2_t = deg_end - deg_start;


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
            do_BGPC_Phase_PlainCC_ccbyG2deg(vtxColors, QQ, Qsizes, srcPtr, dstPtr, vtxVal, &num_uncolored_vertex, &vtim_Pln_CC.back(), G2deg);
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
        ss<<"@BGPCSpIccbyG2deg(";
        ss<<(b_SEQ_HC?"3P":"MP");
        ss<<","<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        ss<<")RNN(";
        ss<<NUM_RND_TC<<","<<NUM_NET_TC<<","<<NUM_NET_CC;
        ss<<")";
        for(int i=ss.str().size(); i<30; i++)
            ss<<"_";
        cout<<ss.str(); ss.str("");
        cout<<"\t#nT_c_T_preG2degT";
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
        cout<<"\t"<<deg2_t.count();
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
// Check conflicts
//
// try idea cc by deg instead of id
// ============================================================================
void BGPCColoring::do_BGPC_Phase_PlainCC_ccbyG2deg(
        vector<int>& vtxColors,
        vector<vector<int>>& QQ,
        vector<int>& Qsizes,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        int* num,
        ChronoDuration* tim
        ,vector<int>const& G2deg
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
            int const ud= G2deg[u];//do_getG2deg(u, srcPtr, dstPtr, vtxVal); //srcPtr[u+1] - srcPtr[u];
            int const uc= vtxColors[u];
            bool b_uis_bad = false;
            for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
                if(b_uis_bad) break;
                int const v = vtxVal[iv];
                for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                    int const w = vtxVal[iw];
                    if(w==u) continue;
                    int const wd= G2deg[w]; //do_getG2deg(w, srcPtr, dstPtr, vtxVal);//srcPtr[w+1] - srcPtr[w];   //
                    if(ud>wd) continue;
                    if(ud==wd && u<w) continue;
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




