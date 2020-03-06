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
// D2GC omp Independent Set iterative approach
//      based on the JP like algorithm
// 
// ----------------------------------------------------------
// Pseudo Code
// while ( |W| is not empty )
//      find_IndependentSet
//      Coloring_IndependentSet
// ----------------------------------------------------------
//
// the independent set is distance two independent set
// ============================================================================
int SMPGCColoring::D2GC_OMP_Hybrid_ISI_K_SpI(
        int nT,                     // number of the threads
        int& colors,                // number of the colors
        vector<int>& vtxColors,     // mapping vertex to its color
        int const LOCAL_ORDER,      // parallel local order
        int const nWeightType,      // random weight or degree weight
        int const K,                // number of ISI iterations
        int const NUM_RND_TC,
        int const NUM_NET_TC,
        int const NUM_NET_CC,
        bool const b_SEQ_HC,
        int const nVerbose          // verbose level
        ) {

    
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);
   
    ChronoDuration tim_PPT(0), tim_Weight(0), tim_LO(0);   // run time for the debug
    
    vector<ChronoDuration>  vtim_ISI;              // run time each iteration 
    vector<int>             vnum_ISI;              // independent set size each iteration 
    vector<ChronoDuration>  vtim_SPI_RndTC;
    vector<ChronoDuration>  vtim_SPI_NetTC;
    vector<ChronoDuration>  vtim_SPI_TC;

    vector<ChronoDuration>  vtim_SPI_NetCC;
    vector<int>             vnum_SPI_NetCC;
    vector<ChronoDuration>  vtim_SPI_CC;
    vector<int>             vnum_SPI_CC;
    
    ChronoDuration          tim_HC{.0};

    const int N               = num_nodes();  // number of vertex
    int const D = max_degree();
    int const DDmDp1 = D*D-D+1;
    const int BufSize = DDmDp1<=0?N:min(DDmDp1, N); //maxDegree
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const vector<int>& vtxVal = get_CSR_ja(); // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    int const ColorLB0Base = max_degree();
    
    vector<mt19937> mts(nT); 
    colors=0;
    vtxColors.assign(N, -1);

    // generate random number
    vector<int> Weight(N);
    { 
        auto start = steady_clock::now();
        if( nWeightType == ISI_WEIGHT_RAND ) {
            srand(RAND_SEED);
            for(int i=0; i<N; i++) Weight[i]=i;
            std::random_shuffle(Weight.begin(), Weight.end());
        }
        else if( nWeightType == ISI_WEIGHT_DEGREE ) {
            for(int i=0; i<N; i++) Weight[i]=vtxPtr[i+1]-vtxPtr[i];
        }
        else{
            cout<<"Error: weight type '"<<nWeightType<<"' undefined"<<endl;
            exit(1);
        }
        auto end = steady_clock::now();
        tim_Weight = end - start;
    }

    // allocate memeory
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16,-1));
    vector<int> Qsizes(nT, 0); 
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16,-1));
    vector<vector<int>> IndSets(nT, vector<int>(N/nT+1+16,-1));
    vector<vector<int>> Workspace(nT, vector<int>(N/nT+1+16,-1));   // used by NetBinColor and NetBinCheck
    // pre-partition the graph
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

    // local order, prepartition, parallel
    ordering("OMP", "PREPARTITION", LOCAL_ORDER, &tim_LO, &QQ);

    // ISI
    int niter = 0;                  // number of iteration
    int num_uncolored_vertex  = N;  // number of uncolored vertex for this iteration
    while(num_uncolored_vertex!=0 && niter<K){
        auto start = steady_clock::now();
        num_uncolored_vertex = 0;
        #pragma omp parallel reduction(+: num_uncolored_vertex)
        {
            // find an independent set
            const int tid       = omp_get_thread_num();
            vector<int>& Q      = QQ[tid];
            int Qsize           = Qsizes[tid];
            vector<int>& IndSet = IndSets[tid];
            int IndSetSize      = 0;
            do_D2GC_Seq_findIndSet(IndSet, IndSetSize, vtxPtr, vtxVal, vtxColors, Weight, Q, Qsize);
            num_uncolored_vertex = Qsize;
            Qsizes[tid] = Qsize;
            #pragma omp barrier
            do_D2GC_Seq_GreedyColoring(vtxPtr, vtxVal, vtxColors, IndSet, IndSetSize, Fs[tid], BufSize, nullptr);
        }// end omp parallel
        auto end = steady_clock::now();
        vtim_ISI.push_back(end-start);
        vnum_ISI.push_back(num_uncolored_vertex);
        niter++;
    }// end of while
   
    // Speculative iterative approach
    niter=0;
    while(num_uncolored_vertex!=0){
    
        bool b_TC_using_NETBIN = false;
        // phase - Tentative Coloring
        if(niter<NUM_RND_TC){
            vtim_SPI_RndTC.emplace_back(.0);
            do_D2GC_OMP_PhaseRandomTC(vtxColors, QQ, Qsizes, mts,  ColorLB0Base, &vtim_SPI_RndTC.back());
        }
        else if(niter<(NUM_RND_TC+NUM_NET_TC)){
            vtim_SPI_NetTC.emplace_back(.0);
            do_D2GC_OMP_PhaseNetBinTC(vtxColors, const_ordered_vertex, N, vtxPtr, vtxVal, Fs, BufSize, Workspace, &vtim_SPI_NetTC.back());
            b_TC_using_NETBIN = true;
        }
        else{
            vtim_SPI_TC.emplace_back(.0);
            do_D2GC_OMP_PhaseTC(vtxColors, QQ, Qsizes, vtxPtr, vtxVal, Fs, BufSize, &vtim_SPI_TC.back());
        }
        
        // phase - Check Conflicts
        if(niter < NUM_NET_CC){
            vtim_SPI_NetCC.emplace_back(.0);
            do_D2GC_OMP_PhaseNetBinCC(QQ, Qsizes, vtxColors, vtxPtr, vtxVal, const_ordered_vertex, N, Fs, BufSize, &num_uncolored_vertex, &vtim_SPI_NetCC.back());
            vnum_SPI_NetCC.emplace_back( num_uncolored_vertex );
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
                vtim_SPI_NetTC.back() += end -start;
            }
            vtim_SPI_CC.emplace_back(.0);
            do_D2GC_OMP_PhaseCC(QQ, Qsizes, vtxColors, vtxPtr, vtxVal, &num_uncolored_vertex, &vtim_SPI_CC.back());
            vnum_SPI_CC.emplace_back( num_uncolored_vertex );
        }
        
        // phase - Handle Conflcits
        if( b_SEQ_HC ){
            auto start = steady_clock::now();
            Fs[0].assign(BufSize,-1);
            for(int i=0; i<nT; i++)  
                do_D2GC_Seq_GreedyColoring_NoInitF(vtxPtr, vtxVal, vtxColors, QQ[i], Qsizes[i], Fs[0], BufSize);
            num_uncolored_vertex=0;
            tim_HC =(steady_clock::now() - start);
        }
        niter++;
    }

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@D2GCHyb(";
        ss<<"wt";
        if( nWeightType == ISI_WEIGHT_RAND ) ss<<"Rnd";
        else if( nWeightType == ISI_WEIGHT_DEGREE )  ss<<"Deg";
        else ss<<"?";
        ss<<",k"<<K;
        ss<<",Lo"<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        ss<<","<<(b_SEQ_HC?"3P":"MP");
        ss<<")";
        ss<<"RNN("<<NUM_RND_TC<<","<<NUM_NET_TC<<","<<NUM_NET_CC<<")";
        for(int i=ss.str().size(); i<30; i++) ss<<"_";
        cout<<ss.str(); ss.str("");
        cout<<"\t#nT_c_T";
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim_Total{.0};
        {
            for(const ChronoDuration&x: vtim_ISI) tim_Total+=x;
            for(const ChronoDuration&x: vtim_SPI_RndTC) tim_Total+=x;
            for(const ChronoDuration&x: vtim_SPI_NetTC) tim_Total+=x;
            for(const ChronoDuration&x: vtim_SPI_TC)   tim_Total+=x;
            for(const ChronoDuration&x: vtim_SPI_NetCC) tim_Total+=x;
            for(const ChronoDuration&x: vtim_SPI_CC)    tim_Total+=x;
            if(b_SEQ_HC) tim_Total += tim_HC;
        }
        cout<<"\t"<<tim_Total.count();

        if(nVerbose>1){ // show detail run time
            {
            cout<<"\t#indset_iter\t"<<vtim_ISI.size();
            cout<<"\t#indset_time\t";
            ChronoDuration t_ISI_total{.0};
            for(const auto& x: vtim_ISI)
                t_ISI_total+=x;
            cout<<t_ISI_total.count();
            cout<<"\t#indset_cnodes\t";
            if(vnum_ISI.empty()) cout<<"0";
            else cout<<N-vnum_ISI.back();
            }

            {
            cout<<"\t#spiRndTC_iter\t"<<vtim_SPI_RndTC.size();
            cout<<"\t#spiRndTC_time\t";
            ChronoDuration t_SPI_rndTime{.0};
            for(const auto&x: vtim_SPI_RndTC)
                t_SPI_rndTime+=x;
            cout<<t_SPI_rndTime.count();
            }

            {
            cout<<"\t#spiNetTC_iter\t"<<vtim_SPI_NetTC.size();
            cout<<"\t#spiNetTC_time\t";
            ChronoDuration t_SPI_NetTime{.0};
            for(const auto&x: vtim_SPI_NetTC)
                t_SPI_NetTime+=x;
            cout<<t_SPI_NetTime.count();
            }

            {
            cout<<"\t#spiPlainTC_iter\t"<<vtim_SPI_TC.size();
            cout<<"\t#spiPlainTC_time\t";
            ChronoDuration t_SPI_Time{.0};
            for(const auto&x: vtim_SPI_TC)
                t_SPI_Time+=x;
            cout<<t_SPI_Time.count();
            }

            {
            cout<<"\t#spiNetCC_iter\t"<<vtim_SPI_NetCC.size();
            cout<<"\t#spiNetCC_time\t";
            ChronoDuration t_SPI_NetCCTime{.0};
            for(const auto&x: vtim_SPI_NetCC)
                t_SPI_NetCCTime+=x;
            cout<<t_SPI_NetCCTime.count();
            cout<<"\t#spiNetCC_nodes\t";
            int t_num_NetCCnodes=0;
            for(const auto&x: vnum_SPI_NetCC)
                t_num_NetCCnodes+=x;
            cout<<t_num_NetCCnodes;
            }

            {
            cout<<"\t#spiPlnCC_iter\t"<<vtim_SPI_CC.size();
            cout<<"\t#spiPlnCC_time\t";
            ChronoDuration t_SPI_PlnTime{.0};
            for(const auto&x: vtim_SPI_CC)
                t_SPI_PlnTime+=x;
            cout<<t_SPI_PlnTime.count();
            cout<<"\t#spiPlnCC_nodes\t";
            int t_num_PlnCCnodes=0;
            for(const auto&x: vnum_SPI_CC)
                t_num_PlnCCnodes+=x;
            cout<<t_num_PlnCCnodes;
            }

            if(b_SEQ_HC){
                cout<<"\t#seqHC_time\t"<<tim_HC.count();
                cout<<"\t#seqHC_nodes\t";
                int t_num_seq_HC=0;
                for(int i=0; i<nT; i++) 
                    t_num_seq_HC += Qsizes[i];
                cout<<t_num_seq_HC;
            }

            if(nVerbose>2){
                cout<<"\n*vtim_ISI";
                cout<<"\n*vnum_ISI";
                cout<<"\n*vtim_RndTC";
                cout<<"\n*vtim_NetTC";
                cout<<"\n*vtim_PlainTC";
                cout<<"\n*vtim_NetCC";
                cout<<"\n*vnum_NetCC";
                cout<<"\n*vtim_PlainCC";
                cout<<"\n*vnum_PlainCC";
            }
        } // end if nVerbose>1
        cout<<"\n";
    }
    return true;
}













