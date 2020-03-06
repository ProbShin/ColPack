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
int SMPGCColoring::D2GC_OMP_IndependentSetIterative(
        int nT,                     // number of the threads
        int& colors,                // number of the colors
        vector<int>&vtxColors,      // mapping vertex to its color
        const int LOCAL_ORDER,      // parallel local order
        const int nWeightType,      // random weight or degree weight
        const int nVerbose          // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);
   
    ChronoDuration tim_PTT, tim_Weight, tim_LO;   // run time for the debug
    vector<ChronoDuration>  vtim_ISI;              // run time each iteration 
    vector<int>             vnum_ISI;              // independent set size each iteration 

    const int N               = num_nodes();  // number of vertex
    int const D = max_degree();
    int const DDmDp1 = D*D-D+1;
    const int BufSize = DDmDp1<=0?N:min(DDmDp1, N); //maxDegree
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const vector<int>& vtxVal = get_CSR_ja(); // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    
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
    vector<int> Qsizes(nT, N/nT); 
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16,-1));
    vector<vector<int>> IndSets(nT, vector<int>(N/nT+1+16,-1));
    
    // pre-partition the graph
    {
        auto start = steady_clock::now();
        Qsizes.assign(nT, N/nT);
        for(int i=0; i<N%nT; i++) 
            Qsizes[i]++;
        vector<int> disps(nT+1, 0); 
        for(int i=1; i<nT+1; i++) 
            disps[i]=disps[i-1]+Qsizes[i-1];
        for(int i=0; i<nT; i++)
            QQ[i].assign(const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
        auto end = steady_clock::now();
        tim_PTT = end - start;
    }

    // local order, prepartition, parallel
    ordering("OMP", "PREPARTITION", LOCAL_ORDER, &tim_LO, (void *)&QQ);

    // ISI
    int niter = 0;                  // number of iteration
    int num_uncolored_vertex = N;
    do{
        auto start = steady_clock::now();
        num_uncolored_vertex  = 0;  // number of uncolored vertex for this iteration
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
            do_D2GC_Seq_GreedyColoring(vtxPtr, vtxVal, vtxColors, IndSet, IndSetSize, Fs[tid], BufSize);
        }// end omp parallel
        auto end = steady_clock::now();

        vtim_ISI.push_back(end-start);
        vnum_ISI.push_back(num_uncolored_vertex);
        niter++;
    }while(num_uncolored_vertex!=0);
    
    colors = calc_num_colors_from_vtx_colors(vtxColors);
    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@D2GCISI("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<",wt";
        if( nWeightType == ISI_WEIGHT_RAND ) ss<<"Rnd";
        else if( nWeightType == ISI_WEIGHT_DEGREE ) ss<<"Deg";
        else ss<<"??";
        ss<<",MP)";
        
        for(int i=ss.str().size(); i<37; i++) ss<<"_";
        cout<<ss.str(); ss.str("");
        cout<<"\t#nT_c_T";
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim{.0};
        for(const ChronoDuration&x: vtim_ISI) tim+=x;
        cout<<"\t"<<tim.count();

        if(nVerbose>1) {
            cout<<"\t#Iterations\t"<<vtim_ISI.size();
            if(nVerbose>3){ // show detail run time
                cout<<"\n*Time "<<tim.count()<<" ["<<vtim_ISI.size()<<"]"; 
                for(const auto &x : vtim_ISI)
                    cout<<" "<<x.count();
                cout<<"\n*num_uncolored_vertex for ["<<vnum_ISI.size()<<"]";
                for(auto &x : vnum_ISI)
                    cout<<" "<<x;
            }
        }
        cout<<"\n";
    }
    return true;
}


// ============================================================================
// Author: Xin Cheng
// 
// ============================================================================
void SMPGCColoring::do_D2GC_Seq_findIndSet(
        vector<int>& IndSet, 
        int &IndSetSize,
        vector<int> const & vtxPtr,
        vector<int> const & vtxVal,
        vector<int> const & vtxColors,
        vector<int> const & Weight,
        vector<int> & Q,
        int &Qsize
        ) {
    
    int num_leftvoer=0;
    for(int iu=0; iu<Qsize; iu++){
        bool b_uisdomain = true;
        auto const u  = Q[iu];
        auto const uw = Weight[u];
        // check distance one neighbors
        for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
            auto const v = vtxVal[iv];
            if( vtxColors[v]>=0 ) continue;   // already colored
            const auto vw= Weight[v];         
            if(uw>vw) continue;               // u is domain
            if(uw==vw && u<v) continue;       // u is tie, but uid is small 
            b_uisdomain = false; 
            break;
        }
        // check distance two neighbor
        if(b_uisdomain){
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                auto const v = vtxVal[iv];
                for(int iw=vtxPtr[v]; iw!=vtxPtr[v+1]; iw++) {
                    auto const w  = vtxVal[iw];
                    if( u==w || vtxColors[w]>=0 ) continue;   // w is u itselt, already colored
                    const auto ww = Weight[w];        
                    if(uw>ww) continue;                       // u is domain
                    if(uw==ww && u<w) continue;               // u is tie, but uid is small
                    b_uisdomain = false;
                    break;
                }
                if(b_uisdomain==false)
                    break;
            }
        }
        // build the independent set
        if(b_uisdomain)
            IndSet[IndSetSize++] = u;
        else
            Q[num_leftvoer++] = u;
    }// end for iu
    Qsize = num_leftvoer;
}

