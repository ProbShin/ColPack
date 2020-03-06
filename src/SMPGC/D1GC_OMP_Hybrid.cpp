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
// D1GC omp Hybrid algorithm (Independent Set iterative + Speculative iterative approach
// ============================================================================
int SMPGCColoring::D1GC_OMP_Hybrid_ISI_K_SpI(
        int nT,                   // number of thread
        int&colors,               // number of colors 
        vector<int>&vtxColors,    // mapping vertex to its color
        const int K,              // numnber of ISI before switch to SpI
        const int nVerbose        // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);
   
    ChronoDuration tim_Partition, tim_Weight;   // run time for the debug
    vector<ChronoDuration>  vtim_ISI, vtim_SpI_TC, vtim_SpI_CC;            // run time each iteration 
    vector<int>             vnum_ISI, vnum_SpI ;            // independent set size each iteration 

    const int N               = num_nodes();  // number of vertex
    const int BufSize         = max_degree()+1;// forbidden color buffer size 
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const vector<int>& vtxVal = get_CSR_ja(); // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    
    colors=0;
    vtxColors.assign(N, -1);

    // allocate the memeory
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16,-1));
    vector<int> Qsizes(nT, N/nT); for(int i=0; i<N%nT; i++) Qsizes[i]++;
    vector<int> Weight(N);
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16,-1));
    vector<vector<int>> IndSets(nT, vector<int>(N/nT+1+16,-1));

    // generate the weight
    { 
        auto start = steady_clock::now();
        #pragma omp parallel for
        for(int i=0; i<N; i++) Weight[i]=vtxPtr[i+1]-vtxPtr[i];
        auto end = steady_clock::now();
        tim_Weight = end - start;
    }

    
    // pre-partition the graph
    {
        auto start = steady_clock::now();
        vector<int> disps(nT+1, 0); 
        for(int i=1; i<nT+1; i++) disps[i]=disps[i-1]+Qsizes[i-1];
        for(int i=0; i<nT; i++)   QQ[i].assign(const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
        auto end = steady_clock::now();
        tim_Partition = end - start;
    }

    // JP's colorings
    int num_input_vertex_size = N;          // number of vertex to color for this iteartion
    int n_iter = 0;                         // number of iterartions
    while(num_input_vertex_size!=0 && n_iter<K){
        auto start = steady_clock::now();
        int num_uncolored_vertex  = 0;      // number of uncolored vertex for this iteration
        #pragma omp parallel reduction(+: num_uncolored_vertex)
        {
            // find an independent set
            const int tid = omp_get_thread_num();
            vector<int> &Q = QQ[tid];
            const int Qsize = Qsizes[tid];
            int newQsize = 0;
            vector<int>& IndSet = IndSets[tid];
            int IndSetSize = 0;
            for(int iu=0; iu<Qsize; iu++){
                bool b_uisdomain = true;
                const auto u = Q[iu];
                const auto uw = Weight[u];
                for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                    const auto v = vtxVal[iv];
                    if( vtxColors[v]>=0 ) continue;
                    const auto vw= Weight[v];
                    if(uw>vw) continue;
                    if(uw==vw && u<v) continue;
                    b_uisdomain=false;
                    break;
                }
                if(b_uisdomain)
                    IndSet[IndSetSize++] = u;
                else
                    Q[newQsize++] = u;
            }// end for iu
            Qsizes[tid]=newQsize;
            num_uncolored_vertex = newQsize;
            #pragma omp barrier

            // greedy color the independent set
            vector<int> &F = Fs[tid];
            F.assign(BufSize,-1);
            for(int iu=0; iu<IndSetSize; iu++){
                const auto u = IndSet[iu];
                for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){
                    const auto vc = vtxColors[vtxVal[iv]];
                    if(vc>=0) 
                        F[vc] = u;
                }
                int c;
                for(c=0; c<BufSize; c++)
                    if( F[c]!=u )
                        break;
                vtxColors[u] = c;
            }
        }// end omp parallel
        auto end = steady_clock::now();
        
        vtim_ISI.push_back(end-start);
        vnum_ISI.push_back(num_input_vertex_size - num_uncolored_vertex);
        num_input_vertex_size = num_uncolored_vertex;
        n_iter++;
    }// end of while
   
    // Speculative iterative approach
    n_iter=0;
    while(num_input_vertex_size!=0){
        vtim_SpI_TC.emplace_back(.0);
        vtim_SpI_CC.emplace_back(.0);
        vnum_SpI.emplace_back(0);

        do_D1GC_OMP_TentativeColoring(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, Fs, BufSize, vtim_SpI_TC[n_iter]);
        do_D1GC_OMP_CheckConflicts(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, vtim_SpI_CC[n_iter], vnum_SpI[n_iter]);
        
        num_input_vertex_size = vnum_SpI[n_iter];
        n_iter++;
    }

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@HybridISISpI(K_"<<K<<"_)_nT_c_T";
        for(int i=ss.str().size(); i<27; i++) ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim_Total{.0};
        for(const ChronoDuration&x: vtim_ISI) tim_Total+=x;
        for(const ChronoDuration&x: vtim_SpI_TC) tim_Total+=x;
        for(const ChronoDuration&x: vtim_SpI_CC) tim_Total+=x;
        //tim_Total += tim_Partition + tim_Weight; 
        cout<<"\t"<<tim_Total.count();
        if(nVerbose>1){ // show detail run time
            cout<<"\n*ISI "<<vtim_ISI.size()<<" iters";
            
            ChronoDuration t_ISI_total{.0};
            for(const auto& x: vtim_ISI)
                t_ISI_total+=x;
            cout<<" takes "<<t_ISI_total.count()<<" seconds (";
            for(const auto& x: vtim_ISI)
                cout<<" "<<x.count();
            cout<<")";
            

            int n_ISI_total=0;
            for(const auto& x: vnum_ISI)
                n_ISI_total+=x;
            cout<<" colors "<<n_ISI_total<<" vertices (";
            for(const auto& x: vnum_ISI)
                cout<<" "<<x;
            cout<<")";


            cout<<"\n*SpI "<<vnum_SpI.size()<<" iters";

            ChronoDuration t_SpI_total{.0};
            for(const auto& x: vtim_SpI_TC)
                t_SpI_total+=x;
            for(const auto& x: vtim_SpI_CC)
                t_SpI_total+=x;
            cout<<" takes "<<t_SpI_total.count()<<" seconds ([";
            for(const auto& x: vtim_SpI_TC) cout<<" "<<x.count();
            cout<<"],[";
            for(const auto& x: vtim_SpI_CC) cout<<" "<<x.count();
            cout<<"])";
           
            int n_SpI_total=0;
            for(const auto& x: vnum_SpI)
                n_SpI_total+=x;

            cout<<" Conflicts "<<n_SpI_total<<" (";
            for(const auto& x: vnum_SpI)
                cout<<" "<<x;
            cout<<")";
        } // end if nVerbose>1
        cout<<"\n";
    }
    return true;
}





// ============================================================================
// Author: Xin Cheng
// D1GC omp Hybrid algorithm (Independent Set iterative + Speculative iterative approach
// limit K to be 0, 1, or 2
// 0 : no jp
// 1 : just one iter of jp
// 2 : just two iter of jp
// others, not correct
// ============================================================================
int SMPGCColoring::D1GC_OMP_Hybrid_ISI_Katmost2_SpI(
        int nT,                   // number of thread
        int&colors,               // number of colors 
        vector<int>&vtxColors,    // mapping vertex to its color
        const int K,              // numnber of ISI before switch to SpI
        const int nVerbose        // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    if(K!=1 && K!=2) { printf("Warnining, number of JP iteration is %d, (it should be in {1,2}\n", K); }
    omp_set_num_threads(nT);
   
    ChronoDuration tim_Partition, tim_Weight;   // run time for the debug
    vector<ChronoDuration>  vtim_ISI, vtim_SpI_TC, vtim_SpI_CC;            // run time each iteration 
    vector<int>             vnum_ISI, vnum_SpI ;            // independent set size each iteration 

    const int N               = num_nodes();  // number of vertex
    const int BufSize         = max_degree()+1;// forbidden color buffer size 
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const vector<int>& vtxVal = get_CSR_ja(); // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    
    colors=0;
    vtxColors.assign(N, -1);

    // allocate the memeory
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16,-1));
    vector<int> Qsizes(nT, N/nT); for(int i=0; i<N%nT; i++) Qsizes[i]++;
    vector<int> Weight(N);
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16,-1));
    vector<vector<int>> IndSets(nT, vector<int>(N/nT+1+16,-1));

    // generate the weight
    { 
        auto start = steady_clock::now();
        #pragma omp parallel for
        for(int i=0; i<N; i++) Weight[i]=vtxPtr[i+1]-vtxPtr[i];
        auto end = steady_clock::now();
        tim_Weight = end - start;
    }

    
    // pre-partition the graph
    {
        auto start = steady_clock::now();
        vector<int> disps(nT+1, 0); 
        for(int i=1; i<nT+1; i++) disps[i]=disps[i-1]+Qsizes[i-1];
        for(int i=0; i<nT; i++)   QQ[i].assign(const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
        auto end = steady_clock::now();
        tim_Partition = end - start;
    }

    // JP's colorings
    int num_input_vertex_size = N;          // number of vertex to color for this iteartion
    int n_iter = 0;                         // number of iterartions
    while(num_input_vertex_size!=0 && n_iter<K){
        auto start = steady_clock::now();
        int num_uncolored_vertex  = 0;      // number of uncolored vertex for this iteration
        if(n_iter==0){
            #pragma omp parallel reduction(+: num_uncolored_vertex)
            {
                // find an independent set
                const int tid = omp_get_thread_num();
                vector<int> &Q = QQ[tid];
                const int Qsize = Qsizes[tid];
                int newQsize = 0;
                //vector<int>& IndSet = IndSets[tid];
                //int IndSetSize = 0;
                for(int iu=0; iu<Qsize; iu++){
                    bool b_uisdomain = true;
                    const auto u = Q[iu];
                    const auto uw = Weight[u];
                    for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                        const auto v = vtxVal[iv];
                        //if( vtxColors[v]>=0 ) continue;
                        const auto vw= Weight[v];
                        if(uw>vw) continue;
                        if(uw==vw && u<v) continue;
                        b_uisdomain=false;
                        break;
                    }
                    if(b_uisdomain)
                        vtxColors[u]=0; //IndSet[IndSetSize++] = u;
                    else
                        Q[newQsize++] = u;
                }// end for iu
                Qsizes[tid]=newQsize;
                num_uncolored_vertex = newQsize;
            }// end omp parallel
        }// end if niter==0
        else{ // niter==1
            #pragma omp parallel reduction(+: num_uncolored_vertex)
            {
                // find an independent set
                const int tid = omp_get_thread_num();
                vector<int> &Q = QQ[tid];
                const int Qsize = Qsizes[tid];
                int newQsize = 0;
                //vector<int>& IndSet = IndSets[tid];
                //int IndSetSize = 0;
                for(int iu=0; iu<Qsize; iu++){
                    bool b_uisdomain = true;
                    const auto u = Q[iu];
                    const auto uw = Weight[u];
                    for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                        const auto v = vtxVal[iv];
                        if( vtxColors[v]==0 ) continue;
                        const auto vw= Weight[v];
                        if(uw>vw) continue;
                        if(uw==vw && u<v) continue;
                        b_uisdomain=false;
                        break;
                    }
                    if(b_uisdomain)
                        vtxColors[u]=1; //IndSet[IndSetSize++] = u;
                    else
                        Q[newQsize++] = u;
                }// end for iu
                Qsizes[tid]=newQsize;
                num_uncolored_vertex = newQsize;
            }// end omp parallel

        }
        auto end = steady_clock::now();

        vtim_ISI.push_back(end-start);
        vnum_ISI.push_back(num_input_vertex_size - num_uncolored_vertex);
        num_input_vertex_size = num_uncolored_vertex;
        n_iter++;
    }// end of while

    // Speculative iterative approach
    n_iter=0;
    while(num_input_vertex_size!=0){
        vtim_SpI_TC.emplace_back(.0);
        vtim_SpI_CC.emplace_back(.0);
        vnum_SpI.emplace_back(0);

        do_D1GC_OMP_TentativeColoring(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, Fs, BufSize, vtim_SpI_TC[n_iter]);
        do_D1GC_OMP_CheckConflicts(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, vtim_SpI_CC[n_iter], vnum_SpI[n_iter]);

        num_input_vertex_size = vnum_SpI[n_iter];
        n_iter++;
    }

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@HybridISISpI(K_"<<K<<"_)_nT_c_T";
        for(int i=ss.str().size(); i<27; i++) ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim_Total{.0};
        for(const ChronoDuration&x: vtim_ISI) tim_Total+=x;
        for(const ChronoDuration&x: vtim_SpI_TC) tim_Total+=x;
        for(const ChronoDuration&x: vtim_SpI_CC) tim_Total+=x;
        //tim_Total += tim_Partition + tim_Weight; 
        cout<<"\t"<<tim_Total.count();
        if(nVerbose>1){ // show detail run time
            cout<<"\n*ISI "<<vtim_ISI.size()<<" iters";

            ChronoDuration t_ISI_total{.0};
            for(const auto& x: vtim_ISI)
                t_ISI_total+=x;
            cout<<" takes "<<t_ISI_total.count()<<" seconds (";
            for(const auto& x: vtim_ISI)
                cout<<" "<<x.count();
            cout<<")";


            int n_ISI_total=0;
            for(const auto& x: vnum_ISI)
                n_ISI_total+=x;
            cout<<" colors "<<n_ISI_total<<" vertices (";
            for(const auto& x: vnum_ISI)
                cout<<" "<<x;
            cout<<")";


            cout<<"\n*SpI "<<vnum_SpI.size()<<" iters";

            ChronoDuration t_SpI_total{.0};
            for(const auto& x: vtim_SpI_TC)
                t_SpI_total+=x;
            for(const auto& x: vtim_SpI_CC)
                t_SpI_total+=x;
            cout<<" takes "<<t_SpI_total.count()<<" seconds ([";
            for(const auto& x: vtim_SpI_TC) cout<<" "<<x.count();
            cout<<"],[";
            for(const auto& x: vtim_SpI_CC) cout<<" "<<x.count();
            cout<<"])";

            int n_SpI_total=0;
            for(const auto& x: vnum_SpI)
                n_SpI_total+=x;

            cout<<" Conflicts "<<n_SpI_total<<" (";
            for(const auto& x: vnum_SpI)
                cout<<" "<<x;
            cout<<")";
        } // end if nVerbose>1
        cout<<"\n";
    }
    return true;
}










