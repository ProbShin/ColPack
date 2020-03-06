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
int SMPGCColoring::D2GC_OMP_SpeculativeIterative(
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
    ChronoDuration         tim_HC(0);

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
    vector<vector<int>> Workspace(nT, vector<int>(N/nT+1+16,-1));   // used by NetBinColor and NetBinCheck
   
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
        bool b_TC_using_NETBIN=false;
        // phase - Tentative Coloring
        if(niter<NUM_RND_TC){
            vtim_RndTC.emplace_back(.0);
            do_D2GC_OMP_PhaseRandomTC(vtxColors, QQ, Qsizes, mts,  ColorLB0Base, &vtim_RndTC.back());
        }
        else if(niter<(NUM_RND_TC+NUM_NET_TC)){
            vtim_NetTC.emplace_back(.0);
            do_D2GC_OMP_PhaseNetBinTC(vtxColors, const_ordered_vertex, N, vtxPtr, vtxVal, Fs, BufSize, Workspace, &vtim_NetTC.back());
            b_TC_using_NETBIN = true;
        }
        else{
            vtim_TC.emplace_back(.0);
            do_D2GC_OMP_PhaseTC(vtxColors, QQ, Qsizes, vtxPtr, vtxVal, Fs, BufSize, &vtim_TC.back());
        }

        // phase - Check Conflicts
        if(niter < NUM_NET_CC){
            // NetBin 'Tentative coloring' would destory already colored vertex. So, we had to reconstruct the uncolored vertex set.
            // However NetBin 'Check conflict' would also destroy the already colored vertex. 
            // So, we would recovere the uncolored vertex at the last stage of the phase 'Check Conflicts'.  
            vtim_NetCC.emplace_back(.0);
            do_D2GC_OMP_PhaseNetBinCC(QQ, Qsizes, vtxColors, vtxPtr, vtxVal, const_ordered_vertex, N, Fs, BufSize, &num_uncolored_vertex, &vtim_NetCC.back());
            vnum_NetCC.emplace_back( num_uncolored_vertex );
        }
        else{
            if(b_TC_using_NETBIN) {
                // NetBin 'Tentative coloring' would destory already colored vertex. So, we had to reconstruct the uncolored vertex set.
                auto start = steady_clock::now();
                #pragma omp parallel
                {
                    int const tid = omp_get_thread_num();
                    vector<int>& Q = QQ[tid];
                    int Qsize=0;
                    #pragma omp for
                    for(int i=0; i<N; i++){
                        if(vtxColors[const_ordered_vertex[i]]<0)
                            Q[Qsize++] = const_ordered_vertex[i];
                    }
                    Qsizes[tid]=Qsize;
                }
                auto end = steady_clock::now();
                vtim_NetTC.back() += end -start;
            }
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
            tim_HC = (steady_clock::now() - start);
            num_uncolored_vertex=0;
        }

        niter++;
    }while(num_uncolored_vertex!=0); //end of while coloring iterartion

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){
        stringstream ss;
        ss<<"@D2GC"<<(b_SEQ_HC?"3P":"MP");
        ss<<"("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<")";
        ss<<"RTC"<<NUM_RND_TC;
        ss<<"NBTC"<<NUM_NET_TC;
        ss<<"NBCC"<<NUM_NET_CC;

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
// Author: Xin Cheng
// author xin cheng
// D2GC omp (Distance Two Graph Coloring OpenMP parallel) Speculative Iterative approach
// with out pre partition. but using atomic operation
// -------------------------------------------------------
// Pseudo code is as follows
// while(|W| is not empty)
//     TentativeColoring(W)  //TC
//     CheckConflicts(W)     //CC
//     HandleConflicts(W)    //HC
// -------------------------------------------------------
// The TC could be one of
//      do_D2GC_omp_PhaseTentativeColoroing_atomic()
//      do_D2GC_omp_PhaseNetBinTentativeColoroing_atomic()
//      do_D2GC_omp_PhaseRndTentativeColoring_atomic()
//
// The CC could be one of
//      do_D2GC_omp_CheckConflicts_atomic()
//      do_D2GC_omp_NetBinCheckConflicts_atomic()
// 
// The HC could be one of
//      *do_nothing*
//      serial_greedy_coloring()
// --------------------------------------------------------
// ============================================================================
int SMPGCColoring::D2GC_OMP_SpeculativeIterative_NoPrePartition(
        int nT,                     // number of threads
        int &colors,                // number of colors
        vector<int>&vtxColors,      // mapping vertex to its color
        const int LOCAL_ORDER,      // parallel local ordering
        const int NUM_RND_TC,       // number of iterations using Random TC
        const int NUM_NET_TC,       // number of iteartions using NetBin TC
        const int NUM_NET_CC,       // number of iteartions using NetBin CC
        const bool b_SEQ_HC,        // whether using sequential coloring to color the remaining graphs.
        const int nVerbose          // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    vector<ChronoDuration> vtim_RndTC,   vtim_NetTC,  vtim_TC;
    vector<ChronoDuration> vtim_NetCC,   vtim_CC;
    vector<int>            vnum_NetCC,   vnum_CC;
    ChronoDuration   tim_HC{.0};
    const int N = num_nodes();                     //number of vertex
    int const D = max_degree();
    int const DDmDp1 = D*D-D+1;
    const int BufSize = DDmDp1<=0?N:min(DDmDp1, N); //maxDegree
    const vector<int>& vtxPtr = get_CSR_ia();
    const vector<int>& vtxVal = get_CSR_ja();
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    const int ColorLB0Base = max_degree();
    vector<mt19937> mts(nT);

    colors=0;                       
    vtxColors.assign(N, -1);

    // allocate memory
    vector<int> Q(const_ordered_vertex);  
    int         Qsize=N;
    vector<int> cfQ(N,-1);  // conflict queue size
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16, -1));
    vector<vector<int>> Workspace(nT, vector<int>(N/nT+1+16,-1));   // used by NetBinColor and NetBinCheck

    // start coloring iterations
    int niter=0;
    int num_uncolored_vertex=N;
    do{


        bool b_TC_using_NETBIN=false;
        // phase - Tentative Coloring
        if(niter<NUM_RND_TC){
            vtim_RndTC.emplace_back(.0);
            do_D2GC_OMP_PhaseRandomTC_NoPrePartition(vtxColors, Q, Qsize, mts, ColorLB0Base, &vtim_RndTC.back());
        }
        else if(niter<(NUM_RND_TC+NUM_NET_TC)){
            vtim_NetTC.emplace_back(.0);
            do_D2GC_OMP_PhaseNetBinTC_NoPrePartition(vtxColors, const_ordered_vertex, N, vtxPtr, vtxVal, Fs, BufSize, Workspace, &vtim_NetTC.back());
            b_TC_using_NETBIN=true;
        }
        else{
            vtim_TC.emplace_back(.0);
            do_D2GC_OMP_PhaseTC_NoPrePartition(vtxColors, Q, Qsize, vtxPtr, vtxVal, Fs, BufSize, &vtim_TC.back());
        }


        //if(niter>0)
        //    break;

        // phase - Check Conflicts
        if(niter < NUM_NET_CC){
            vtim_NetCC.emplace_back(.0);
            do_D2GC_OMP_PhaseNetBinCC_NoPrePartition(Q, Qsize, vtxColors, vtxPtr, vtxVal, const_ordered_vertex, N, Fs, BufSize, &num_uncolored_vertex, &vtim_NetCC.back());
            vnum_NetCC.emplace_back(num_uncolored_vertex);
        }
        else{
            if(b_TC_using_NETBIN){
                auto start = steady_clock::now();
                Qsize=0;
                #pragma omp parallel for
                for(int i=0; i<N; i++){
                    if(vtxColors[const_ordered_vertex[i]]<0){
                        int enQposition = __sync_fetch_and_add(&Qsize, 1);
                        Q[enQposition]=const_ordered_vertex[i];
                    }
                }
                auto end = steady_clock::now();
                vtim_NetTC.back() +=  end - start;
            }
    
            vtim_CC.emplace_back(.0);
            do_D2GC_OMP_PhaseCC_NoPrePartition(Q, Qsize, vtxColors, vtxPtr, vtxVal, cfQ, &num_uncolored_vertex, &vtim_CC.back());
            vnum_CC.emplace_back(num_uncolored_vertex);
        }
        

        // phase - Handle Conflcits
        if( b_SEQ_HC ){
            auto start = steady_clock::now();
            do_D2GC_Seq_GreedyColoring(vtxPtr, vtxVal, vtxColors, Q, Qsize, Fs[0], BufSize);
            tim_HC = (steady_clock::now() - start);
            num_uncolored_vertex = 0;
        }
    
        niter++;
    }while(num_uncolored_vertex!=0); //end of while coloring iterartion

    colors = calc_num_colors_from_vtx_colors(vtxColors);
    
    if(nVerbose>0){
        stringstream ss;
        ss<<"@D2GCAtomic"<<(b_SEQ_HC?"3P":"MP");
        ss<<"("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<")";
        ss<<"RTC"<<NUM_RND_TC;
        ss<<"NBTC"<<NUM_NET_TC;
        ss<<"NBCC"<<NUM_NET_CC;
        for(int i=ss.str().size(); i<30; i++) ss<<"_";
        ss<<"\t#nT_c_T";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim{.0};
        {
            //tim=tim_PPT+tim_LO;
            for(const auto& x: vtim_RndTC)  tim+=x;
            for(const auto& x: vtim_NetTC) tim+=x;
            for(const auto& x: vtim_TC)   tim+=x;
            for(const auto& x: vtim_NetCC) tim+=x;
            for(const auto& x: vtim_CC)   tim+=x;
            tim+=tim_HC;
        }
        cout<<"\t"<<tim.count();

        // more details
        if(nVerbose>1){
            cout<<"\t#niter\t"<<vtim_CC.size()+vtim_NetCC.size();
            ChronoDuration tim_TC(.0);
            {
                for(const auto& x: vtim_RndTC)  tim_TC+=x;
                for(const auto& x: vtim_NetTC) tim_TC+=x;
                for(const auto& x: vtim_TC)   tim_TC+=x;
            }
            cout<<"\t#TCtime\t"<<tim_TC.count();
            ChronoDuration tim_CC(.0);
            {
                for(const auto& x: vtim_NetCC)  tim_CC+=x;
                for(const auto& x: vtim_CC) tim_CC+=x;
            }
            cout<<"\t#CCtime\t"<<tim_CC.count();
            if(b_SEQ_HC)
                cout<<"\t#SeqHCtime\t"<<tim_HC.count();
            
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
                    cout<<"\n*NBTCtimes "<<vtim_NetTC.size()<<" iterts: ";
                    for(const auto& x: vtim_NetTC)
                        cout<<" "<<x.count();
                }
                if(vtim_TC.size()!=0) {
                    cout<<"\n*TCtimes "<<vtim_TC.size()<<" iterts: ";
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
// Author xin cheng
// 
// ============================================================================


