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
// D1GC omp Independent Set iterative approach with Luby's algorhtm
// ============================================================================
int SMPGCColoring::D1GC_OMP_IndependentSet_Iterative_LB(int nT, int&colors, vector<int>&vtxColors, const int LOCAL_ORDER, const int nVerbose) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);
   
    ChronoDuration tim_Partition, tim_Weight, tim_LO;   // run time for the debug
    vector<ChronoDuration>  vtim_ISI;            // run time each iteration 
    vector<int>             vnum_ISI;            // independent set size each iteration 

    const int N               = num_nodes();  // number of vertex
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const vector<int>& vtxVal = get_CSR_ja(); // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    
    colors=0;
    vtxColors.assign(N, -1);

    // generate the random number
    vector<int> Weight(N);
    { 
        auto start = steady_clock::now();
        srand(RAND_SEED);
        for(int i=0; i<N; i++) Weight[i]=i;
        std::random_shuffle(Weight.begin(), Weight.end());
        auto end = steady_clock::now();
        tim_Weight = end - start;
    }

    // allocate memeory
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16,-1));
    vector<int> Qsizes(nT, N/nT); for(int i=0; i<N%nT; i++) Qsizes[i]++;
    vector<vector<int>> IndSets(nT, vector<int>(N/nT+1+16,-1));

    // pre-partition the graph
    {
        auto start = steady_clock::now();
        vector<int> disps(nT+1, 0); 
        for(int i=1; i<nT+1; i++) 
            disps[i]=disps[i-1]+Qsizes[i-1];
        for(int i=0; i<nT; i++)
            //QQ[i].insert(QQ[i].end(), const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
            QQ[i].assign(const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
        auto end = steady_clock::now();
        tim_Partition = end - start;
    }

    // local order, prepartition, parallel
    ordering("OMP", "PREPARTITION", LOCAL_ORDER, &tim_LO, (void *)&QQ);

    // Luby's colorings
    int num_input_vertex_size = N;
    int num_uncolored_vertex  = 0;
    int niter = 0;
    do{
        
        auto start = steady_clock::now();
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

            // color the independent set
            for(int iu=0; iu<IndSetSize; iu++){
                vtxColors[ IndSet[iu] ] = niter;
            }
        }// end omp parallel
        auto end = steady_clock::now();
        
        vtim_ISI.push_back(end-start);
        vnum_ISI.push_back(num_input_vertex_size - num_uncolored_vertex);
        num_input_vertex_size = num_uncolored_vertex;
        num_uncolored_vertex = 0;
        niter++;
    }while(num_input_vertex_size!=0);
    
    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@ISILB("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<")_nT_c_T";
        for(int i=ss.str().size(); i<27; i++) ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim_Total{.0};
        for(const ChronoDuration&x: vtim_ISI) tim_Total+=x;
        cout<<"\t"<<tim_Total.count();

        if(nVerbose>1){ // show detail run time
            cout<<"\n*T total "<<tim_Total.count()<<" ["<<vtim_ISI.size()<<"]"; 
            for(auto &x : vtim_ISI)
                cout<<" "<<x.count();

            cout<<"ISsize ["<<vnum_ISI.size()<<"]";
            for(auto &x : vnum_ISI)
                cout<<" "<<x;
        } // end if nVerbose>1
        cout<<"\n";
    }
    return true;
}




// ============================================================================
// Author: Xin Cheng
// D1GC omp Independent Set iterative approach with JP' algorhtm
// ============================================================================
int SMPGCColoring::D1GC_OMP_IndependentSet_Iterative_JP(int nT, int&colors, vector<int>&vtxColors, const int LOCAL_ORDER, const int nVerbose) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);
   
    ChronoDuration tim_Partition, tim_Weight, tim_LO;   // run time for the debug
    vector<ChronoDuration>  vtim_IS;            // run time each iteration 
    vector<int>             vnum_IS;            // independent set size each iteration 

    const int N               = num_nodes();  // number of vertex
    const int BufSize         = max_degree()+1;// forbidden color buffer size 
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const vector<int>& vtxVal = get_CSR_ja(); // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    
    colors=0;
    vtxColors.assign(N, -1);

    // generate random number
    vector<int> Weight(N);
    { 
        auto start = steady_clock::now();
#ifdef D1GC_OMP_ISI_JP_USING_RANDOM_WEIGHT
        srand(RAND_SEED);
        for(int i=0; i<N; i++) Weight[i]=i;
        std::random_shuffle(Weight.begin(), Weight.end());
#else
        #pragma omp parallel for
        for(int i=0; i<N; i++) Weight[i]=vtxPtr[i+1]-vtxPtr[i];
#endif
        auto end = steady_clock::now();
        tim_Weight = end - start;
    }

    // allocate memeory
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16,-1));
    vector<int> Qsizes(nT, N/nT); for(int i=0; i<N%nT; i++) Qsizes[i]++;
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16,-1));
    vector<vector<int>> IndSets(nT, vector<int>(N/nT+1+16,-1));
    
    // pre-partition the graph
    {
        auto start = steady_clock::now();
        vector<int> disps(nT+1, 0); 
        for(int i=1; i<nT+1; i++) 
            disps[i]=disps[i-1]+Qsizes[i-1];
        for(int i=0; i<nT; i++)
            //QQ[i].insert(QQ[i].end(), const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
            QQ[i].assign(const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
        auto end = steady_clock::now();
        tim_Partition = end - start;
    }

    // local order, prepartition, parallel
    ordering("OMP", "PREPARTITION", LOCAL_ORDER, &tim_LO, (void *)&QQ);

    // JP's colorings
    int num_input_vertex_size = N;
    int num_uncolored_vertex  = 0;
    int niter = 0;
    do{
        
        auto start = steady_clock::now();
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
       
            // greedy color the independent set
            vector<int> &F = Fs[tid];
            //F.assign(BufSize,-1);   // F does not need to be init to -1
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
        
            Qsizes[tid]=newQsize;
            num_uncolored_vertex = newQsize;
        }// end omp parallel
        auto end = steady_clock::now();
        
        vtim_IS.push_back(end-start);
        vnum_IS.push_back(num_input_vertex_size - num_uncolored_vertex);
        num_input_vertex_size = num_uncolored_vertex;
        num_uncolored_vertex = 0;
        niter++;
    }while(num_input_vertex_size!=0);
    
    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@ISIJP("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<")_nT_c_T";
        for(int i=ss.str().size(); i<27; i++) ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim_Total{.0};
        for(const ChronoDuration&x: vtim_IS) tim_Total+=x;
        //tim_Total += tim_Partition + tim_Weight;
        cout<<"\t"<<tim_Total.count();

        if(nVerbose>1){ // show detail run time
            cout<<"\n*T total "<<tim_Total.count()<<" ["<<vtim_IS.size()<<"]"; 
            for(auto &x : vtim_IS)
                cout<<" "<<x.count();

            cout<<"ISsize ["<<vnum_IS.size()<<"]";
            for(auto &x : vnum_IS)
                cout<<" "<<x;
        } // end if nVerbose>1
        cout<<"\n";
    }
    return true;
}





// ============================================================================
// Author: Xin Cheng
// D1GC omp Independent Set iterative approach with JP' algorhtm
// Find both Largest Local Domain and Smallest Local Domain
// ============================================================================
int SMPGCColoring::D1GC_OMP_IndependentSet_Iterative_JP_adv(int nT, int&colors, vector<int>&vtxColors, const int LOCAL_ORDER, const int nVerbose) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);
   
    ChronoDuration tim_Partition, tim_Weight, tim_LO;   // run time for the debug
    vector<ChronoDuration>  vtim_IS;            // run time each iteration 
    vector<int>             vnum_IS;            // independent set size each iteration 

    const int N               = num_nodes();  // number of vertex
    const int BufSize         = max_degree()+1;// forbidden color buffer size 
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const vector<int>& vtxVal = get_CSR_ja(); // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 
    
    colors=0;
    vtxColors.assign(N, -1);

    // generate random number
    vector<int> Weight(N);
    { 
        auto start = steady_clock::now();
#ifdef D1GC_OMP_ISI_JP_USING_RANDOM_WEIGHT
        srand(RAND_SEED);
        for(int i=0; i<N; i++) Weight[i]=i;
        std::random_shuffle(Weight.begin(), Weight.end());
#else
        #pragma omp parallel for
        for(int i=0; i<N; i++) Weight[i]=vtxPtr[i+1]-vtxPtr[i];
#endif
        auto end = steady_clock::now();
        tim_Weight = end - start;
    }

    // allocate memeory
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16,-1));
    vector<int> Qsizes(nT, N/nT); for(int i=0; i<N%nT; i++) Qsizes[i]++;
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16,-1));
    vector<vector<int>> IndSets(nT, vector<int>(N/nT+1+16,-1));
    vector<vector<int>> IndSets2(nT, vector<int>(N/nT+1+16,-1));
    
    // pre-partition the graph
    {
        auto start = steady_clock::now();
        vector<int> disps(nT+1, 0); 
        for(int i=1; i<nT+1; i++) 
            disps[i]=disps[i-1]+Qsizes[i-1];
        for(int i=0; i<nT; i++)
            //QQ[i].insert(QQ[i].end(), const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
            QQ[i].assign(const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
        auto end = steady_clock::now();
        tim_Partition = end - start;
    }

    // local order, prepartition, parallel
    ordering("OMP", "PREPARTITION", LOCAL_ORDER, &tim_LO, (void *)&QQ);

    // JP's colorings
    int num_input_vertex_size = N;
    int num_uncolored_vertex  = 0;
    int niter = 0;
    do{
        
        auto start = steady_clock::now();
        #pragma omp parallel reduction(+: num_uncolored_vertex)
        {
            // find an independent set
            const int tid = omp_get_thread_num();
            vector<int> &Q = QQ[tid];
            const int Qsize = Qsizes[tid];
            int newQsize = 0;
            vector<int>& IndSet = IndSets[tid];
            vector<int>& IndSet2= IndSets2[tid];
            int IndSetSize = 0;
            int IndSetSize2= 0;
            for(int iu=0; iu<Qsize; iu++){
                bool b_uisLdomain = true;
                bool b_uisSdomain = true;

                const auto u = Q[iu];
                const auto uw = Weight[u];
                for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                    if(b_uisLdomain==false && b_uisSdomain==false) { break; }
                    
                    const auto v = vtxVal[iv];
                    if( vtxColors[v]>=0 ) continue;
                    const auto vw= Weight[v];
                    
                    if(b_uisLdomain){  //Ldomain == true, Sdomain=??
                        if (       uw<vw) { b_uisLdomain=false;  
                        } else if (uw>vw) { b_uisSdomain=false;
                        } else { //uw==vw
                            if( u<v) { b_uisLdomain=false; 
                            } else   { b_uisSdomain=false; 
                            }
                        } 
                    }
                    else{   // Ldomain ==false, Sdomain==true
                        if(uw>vw || (uw==vw && u>v)) {
                            b_uisSdomain=false;
                            break;
                        }
                    }
                }


                if(b_uisLdomain)
                    IndSet[IndSetSize++] = u;
                else if(b_uisSdomain)
                    IndSet2[IndSetSize2++] = u;
                else
                    Q[newQsize++] = u;
            }// end for iu
       
            // greedy color the independent set
            // coloring part A
            vector<int> &F = Fs[tid];
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
            
            #pragma omp barrier
            // coloring part B
            for(int iu=0; iu<IndSetSize2; iu++){
                const auto u = IndSet2[iu];
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


            Qsizes[tid]=newQsize;
            num_uncolored_vertex = newQsize;
        }// end omp parallel
        auto end = steady_clock::now();
        
        vtim_IS.push_back(end-start);
        vnum_IS.push_back(num_input_vertex_size - num_uncolored_vertex);
        num_input_vertex_size = num_uncolored_vertex;
        num_uncolored_vertex = 0;
        niter++;
    }while(num_input_vertex_size!=0);
    
    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@ISIJPadv("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<")_nT_c_T";
        for(int i=ss.str().size(); i<27; i++) ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim_Total{.0};
        for(const ChronoDuration&x: vtim_IS) tim_Total+=x;
        //tim_Total += tim_Partition + tim_Weight;
        cout<<"\t"<<tim_Total.count();

        if(nVerbose>1){ // show detail run time
            cout<<"\n*T total "<<tim_Total.count()<<" ["<<vtim_IS.size()<<"]"; 
            for(auto &x : vtim_IS)
                cout<<" "<<x.count();

            cout<<"ISsize ["<<vnum_IS.size()<<"]";
            for(auto &x : vnum_IS)
                cout<<" "<<x;
        } // end if nVerbose>1
        cout<<"\n";
    }
    return true;
}









