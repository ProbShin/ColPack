/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "BGPCColoring.h"
using namespace ColPack;
using namespace std;
using namespace chrono;
// ============================================================================
// author: xin cheng
// Bipartite Graph Partial Coloring parallel implementation
// ============================================================================
int BGPCColoring::BGPC_OMP_Hybrid_ISI_SPI(
        int const side, 
        int nT, 
        int& colors, 
        vector<int>& vtxColors, 
        int const LOCAL_ORDER, 
        int const WEIGHT_TYPE , 
        int const K, 
        int const NUM_RND_TC, 
        int const NUM_NET_TC, 
        int const NUM_NET_CC, 
        bool const b_SEQ_HC,
        int const nVerbose ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    ChronoDuration tim_PPT{0};
    ChronoDuration tim_LO{0};
    ChronoDuration tim_WT{0};
    vector<ChronoDuration>  vtim_ISI;
    vector<int> vnum_ISI;
    ChronoDuration tim_HC{0};
    vector<ChronoDuration>  vtim_Rnd_TC;
    vector<ChronoDuration>  vtim_Net_TC;
    vector<ChronoDuration>  vtim_Pln_TC;
    vector<ChronoDuration>  vtim_Net_CC;
    vector<ChronoDuration>  vtim_Pln_CC;
    vector<int> vnum_Net_CC;
    vector<int> vnum_Pln_CC;
    
    const int          N             = (side==BGPC::L)?(GetLeftVertexCount()):(GetRightVertexCount()); 
    const int          Ndst          = (side==BGPC::L)?(GetRightVertexCount()):(GetLeftVertexCount()); 
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    const vector<int>& vtxVal        = GetEdges();
    const int          srcMaxDegree  = (side==BGPC::L)?GetMaximumLeftVertexDegree():GetMaximumRightVertexDegree();
    const int          dstMaxDegree  = (side==BGPC::L)?GetMaximumRightVertexDegree():GetMaximumLeftVertexDegree();
    int const          DDp1          = dstMaxDegree * srcMaxDegree +1;
    int                BufSize       = (DDp1<0)?N:min(N,DDp1);
    
    vector<int>const& const_queue_A  = get_ordered_queue_A_const();
    vector<int>       const_queue_B;
    vector<vector<int>> Workloads;
    vector<uint64_t>    Weights;
    vector<vector<int>> IndSets;
    vector<int>         IndSetSizes;
    vector<mt19937>     mts;

    colors=0;                       
    vtxColors.assign(N, -1);
    
    // allocate memory,  +16 for prevent false share on some 64 bit machines
    vector<vector<int>> ForbiddenArrays(nT, vector<int>(BufSize+16,-1)); 
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16, -1)); 
    vector<int> Qsizes(nT,0);
    if(NUM_NET_TC || NUM_NET_CC) { const_queue_B.assign(Ndst,0); for(int i=0; i<Ndst; i++) const_queue_B[i] = i; }    // creat memory only if Netbin algs going to execute
    if(NUM_NET_TC ) { Workloads.assign(nT, vector<int>(BufSize+16,-1)); }
    if(NUM_RND_TC ) { for(int i=0; i<nT; i++) mts.emplace_back(i); }
    if(K>0) Weights.assign(N,0);
    if(K>0) IndSets.assign(nT, vector<int>(N/nT+1+16, -1));
    if(K>0) IndSetSizes.assign(nT,0);

    // pre-partition
    init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &tim_PPT); 

    // generate the weight
    do_Generate_Weights(Weights, N, WEIGHT_TYPE, &srcPtr[0], &dstPtr[0], &vtxVal[0], &tim_WT);

    // local ordering   
    Ordering_PrePartition_OMP(side, QQ, Qsizes, LOCAL_ORDER, &tim_LO); 
    
    int niter=0;
    int num_uncolored_vertex = N;
    while(niter<K && num_uncolored_vertex!=0){
    
        auto start = steady_clock::now();
        num_uncolored_vertex = 0;
        #pragma omp parallel reduction(+:num_uncolored_vertex)
        {
            int const tid = omp_get_thread_num();
            do_BGPC_Seq_FindIndSet(IndSets[tid], IndSetSizes[tid], QQ[tid], Qsizes[tid], vtxColors, srcPtr, dstPtr, vtxVal, Weights);
            num_uncolored_vertex = Qsizes[tid];
            #pragma omp barrier
            do_BGPC_Seq_Greedy(vtxColors, IndSets[tid], IndSetSizes[tid], srcPtr, dstPtr, vtxVal, ForbiddenArrays[tid], BufSize);
        }
        auto end = steady_clock::now();
        vtim_ISI.push_back(end - start);
        vnum_ISI.push_back(num_uncolored_vertex);
        niter      += 1;
    } //end while
    
    niter=0;
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
                ChronoDuration t{.0}; init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &t); vtim_Net_CC.back()+=t;
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
        ss<<"@BGPCHyb(";
        ss<<"K"<<K;
        ss<<","<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        if( WEIGHT_TYPE == WEIGHT_RANDOM ) ss<<",RndW";
        else if( WEIGHT_TYPE == WEIGHT_DEGREE_D1 ) ss<<",D1W";
        else if( WEIGHT_TYPE == WEIGHT_DEGREE_D2 ) ss<<",D2W";
        else ss<<",??W";
        ss<<","<<(b_SEQ_HC?"3P":"MP");
        ss<<",RNN"<<NUM_RND_TC<<","<<NUM_NET_TC<<","<<NUM_NET_CC;
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
        //T += tim_WT;
        for(auto const x : vtim_ISI) T+=x;
        for(auto const x : vtim_Rnd_TC) T+=x;
        for(auto const x : vtim_Net_TC) T+=x;
        for(auto const x : vtim_Pln_TC) T+=x;
        for(auto const x : vtim_Net_CC) T+=x;
        for(auto const x : vtim_Pln_CC) T+=x;
        T += tim_HC;

        cout<<"\t"<<T.count();
        if(nVerbose>1) {
            cout<<"\tISI\t";
            {
                ChronoDuration T{0}; 
                for(auto const x : vtim_ISI) T+=x;
                cout<<"\t"<<T.count();
            }
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
// limit ISI for only 1 or 2 iteration
// ============================================================================
int BGPCColoring::BGPC_OMP_Hybrid_ISI_SPI_Katmost2(
        int const side, 
        int nT, 
        int& colors, 
        vector<int>& vtxColors, 
        int const LOCAL_ORDER, 
        int const WEIGHT_TYPE , 
        int const K, 
        int const NUM_RND_TC, 
        int const NUM_NET_TC, 
        int const NUM_NET_CC, 
        bool const b_SEQ_HC,
        int const nVerbose ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    if(K!=1 && K!=2) { printf("Warning! K=%d, but it should be either 1 or 2. ",K); return false; }
    omp_set_num_threads(nT);

    ChronoDuration tim_PPT{0};
    ChronoDuration tim_LO{0};
    ChronoDuration tim_WT{0};
    vector<ChronoDuration>  vtim_ISI;
    vector<int> vnum_ISI;
    ChronoDuration tim_HC{0};
    vector<ChronoDuration>  vtim_Rnd_TC;
    vector<ChronoDuration>  vtim_Net_TC;
    vector<ChronoDuration>  vtim_Pln_TC;
    vector<ChronoDuration>  vtim_Net_CC;
    vector<ChronoDuration>  vtim_Pln_CC;
    vector<int> vnum_Net_CC;
    vector<int> vnum_Pln_CC;
    
    const int          N             = (side==BGPC::L)?(GetLeftVertexCount()):(GetRightVertexCount()); 
    const int          Ndst          = (side==BGPC::L)?(GetRightVertexCount()):(GetLeftVertexCount()); 
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    const vector<int>& vtxVal        = GetEdges();
    const int          srcMaxDegree  = (side==BGPC::L)?GetMaximumLeftVertexDegree():GetMaximumRightVertexDegree();
    const int          dstMaxDegree  = (side==BGPC::L)?GetMaximumRightVertexDegree():GetMaximumLeftVertexDegree();
    int const          DDp1          = dstMaxDegree * srcMaxDegree +1;
    int                BufSize       = (DDp1<0)?N:min(N,DDp1);
    
    vector<int>const& const_queue_A  = get_ordered_queue_A_const();
    vector<int>       const_queue_B;
    vector<vector<int>> Workloads;
    vector<uint64_t>    Weights;
    vector<vector<int>> IndSets;
    vector<int>         IndSetSizes;
    vector<mt19937>     mts;

    colors=0;                       
    vtxColors.assign(N, -1);
    
    // allocate memory,  +16 for prevent false share on some 64 bit machines
    vector<vector<int>> ForbiddenArrays(nT, vector<int>(BufSize+16,-1)); 
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16, -1)); 
    vector<int> Qsizes(nT,0);
    if(NUM_NET_TC || NUM_NET_CC) { const_queue_B.assign(Ndst,0); for(int i=0; i<Ndst; i++) const_queue_B[i] = i; }    // creat memory only if Netbin algs going to execute
    if(NUM_NET_TC ) { Workloads.assign(nT, vector<int>(BufSize+16,-1)); }
    if(NUM_RND_TC ) { for(int i=0; i<nT; i++) mts.emplace_back(i); }
    if(K>0) Weights.assign(N,0);
    if(K>0) IndSets.assign(nT, vector<int>(N/nT+1+16, -1));
    if(K>0) IndSetSizes.assign(nT,0);

    // pre-partition
    init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &tim_PPT); 

    // generate the weight
    do_Generate_Weights(Weights, N, WEIGHT_TYPE, &srcPtr[0], &dstPtr[0], &vtxVal[0], &tim_WT);

    // local ordering   
    Ordering_PrePartition_OMP(side, QQ, Qsizes, LOCAL_ORDER, &tim_LO); 
    
    int num_uncolored_vertex = N;
    // block of D2IS
    {
        auto start = steady_clock::now();

        num_uncolored_vertex = 0;
        #pragma omp parallel reduction(+:num_uncolored_vertex)
        {
            int const tid = omp_get_thread_num();
            int num_leftover=0;

            vector<int>& Q = QQ[tid];
            // first iteration of D2IS
            {
                int const Qsize = Qsizes[tid];
                for(int iu=0; iu<Qsize; iu++) {
                    int const u = Q[iu];
                    auto const uw = Weights[u];
                    bool b_uisdomain = true;
                    for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
                        if(!b_uisdomain) break;
                        int const v = vtxVal[iv];
                        for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                            int const w  = vtxVal[iw];
                            if(u==w) continue; //if(u==w || vtxColors[w]>=0 ) continue;
                            auto const ww = Weights[w];
                            if(uw > ww) continue;
                            if(uw == ww && u>w) continue;  //if not unique weight, consider tie
                            b_uisdomain = false;
                            break;
                        }
                    }
                    if(b_uisdomain)
                        vtxColors[u]=0;
                    else
                        Q[num_leftover++] = u;
                }
                Qsizes[tid]=num_leftover;
            }            
            // second iteration of D2IS
            if(K==2){
                num_leftover=0;
                int const Qsize = Qsizes[tid];
                #pragma omp barrier
                for(int iu=0; iu<Qsize; iu++) {
                    int const u = Q[iu];
                    auto const uw = Weights[u];
                    bool b_uisdomain = true;
                    for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
                        if(!b_uisdomain) break;
                        int const v = vtxVal[iv];
                        for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                            int const w  = vtxVal[iw];
                            if(u==w || vtxColors[w]==0 ) continue; //if(u==w || vtxColors[w]>=0 ) continue;
                            auto const ww = Weights[w];
                            if(uw > ww) continue;
                            if(uw == ww && u>w) continue;  //if not unique weight, consider tie
                            b_uisdomain = false;
                            break;
                        }
                    }
                    if(b_uisdomain)
                        vtxColors[u]=1;
                    else
                        Q[num_leftover++] = u;
                }
                Qsizes[tid]=num_leftover;
            }
            num_uncolored_vertex = num_leftover;
        } // end of omp parallel
        
        auto end = steady_clock::now();
        vtim_ISI.push_back(end - start);
        vnum_ISI.push_back(num_uncolored_vertex);
    } // end of block for D2IS
    int niter=0;
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
                ChronoDuration t{.0}; init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &t); vtim_Net_CC.back()+=t;
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
        ss<<"@BGPCHyb(";
        ss<<"K"<<K;
        ss<<","<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        if( WEIGHT_TYPE == WEIGHT_RANDOM ) ss<<",RndW";
        else if( WEIGHT_TYPE == WEIGHT_DEGREE_D1 ) ss<<",D1W";
        else if( WEIGHT_TYPE == WEIGHT_DEGREE_D2 ) ss<<",D2W";
        else ss<<",??W";
        ss<<","<<(b_SEQ_HC?"3P":"MP");
        ss<<",RNN"<<NUM_RND_TC<<","<<NUM_NET_TC<<","<<NUM_NET_CC;
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
        //T += tim_WT;
        for(auto const x : vtim_ISI) T+=x;
        for(auto const x : vtim_Rnd_TC) T+=x;
        for(auto const x : vtim_Net_TC) T+=x;
        for(auto const x : vtim_Pln_TC) T+=x;
        for(auto const x : vtim_Net_CC) T+=x;
        for(auto const x : vtim_Pln_CC) T+=x;
        T += tim_HC;

        cout<<"\t"<<T.count();
        if(nVerbose>1) {
            cout<<"\tISI\t";
            {
                ChronoDuration T{0}; 
                for(auto const x : vtim_ISI) T+=x;
                cout<<"\t"<<T.count();
            }
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


