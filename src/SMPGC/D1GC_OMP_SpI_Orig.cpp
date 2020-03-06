/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#include "SMPGCColoring.h"
#include <chrono> //c++11 system time
#include <random> //c++11 random
using namespace std;
using namespace ColPack;


// ============================================================================
// Author: Xin Cheng
// ----------------------------------------------------------------------------
// D1GC speculative iterative coloring, original version.
//
// while |W| is not empty
//     TentativeColoring
//     CheckConflicts
//     HandleConlficts
// 
//
// 3P means hadleConflicts  will sequential color the conflicts vertices
//
// original version means
//     no-prepartition,  
//     merge the conflicts 
//     does not using local orderings
// ============================================================================
int SMPGCColoring::D1GC_OMP_SpeculativeIterative_3P_Original(
        int nT,                     // number of threads 
        int&colors,                 // number of colors
        vector<int>&vtxColors,      // mapping vertex to its color
        const int nVerbose          // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);
    
    ChronoDuration tim_TC, tim_CC, tim_HC;   // run time for debug
    int num_Cf = 0;                          // number of conflicts for debug

    const int N               = num_nodes();   // number of vertex
    const int BufSize         = max_degree()+1;// forbidden color buffer size
    const vector<int>& vtxPtr = get_CSR_ia();  // csr format edges
    const vector<int>& vtxVal = get_CSR_ja();  // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); // ordered vertex 

    colors=0;                       
    vtxColors.assign(N, -1);
    
    // allocate memory
    vector<int> Q(const_ordered_vertex);
    int Qsize=N;
    vector<int> cfQ(N);
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16, -1) ); 

    // Tentative Coloring
    do_TentativeColoring_NoPrePartition(nT, vtxColors, vtxPtr, vtxVal, Q, Qsize, Fs, BufSize, tim_TC);

    // Check Conflicts
    do_CheckConflicts_NoPrePartition(nT, vtxColors, vtxPtr, vtxVal, Q, Qsize, cfQ, num_Cf, tim_CC);

    // Handle Conflicts
    {
        auto start = steady_clock::now();
        vector<int>&F = Fs[0];
        F.assign(BufSize, -1);
        for(int iu=0; iu<num_Cf; iu++){
            const int u = Q[iu];
            for(int iv=vtxPtr[u]; iv<vtxPtr[u+1]; iv++){
                const int v = vtxVal[iv];
                const int vc = vtxColors[v];
                if(vc>=0) F[vc]=u;
            }
            int c;
            for(c=0; c<BufSize; c++){
                if(F[c]!=u)
                    break;
            }
            vtxColors[u]=c;
        }// end for iu
        auto end =steady_clock::now();
        tim_HC = end - start;
    }

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@SpI3POrig(None)_nT_c_T";
        for(int i=ss.str().size(); i<27; i++) ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        auto tim_Total = tim_TC + tim_CC + tim_HC;
        cout<<"\t"<<tim_Total.count();

        if(nVerbose>1){ // show detail run time
            cout<<"\t*tTC_tCC_tHC";
            cout<<"\t"<<tim_TC.count();
            cout<<"\t"<<tim_CC.count();
            cout<<"\t"<<tim_HC.count();
        
            if(nVerbose>2){ // show conflicts information
                cout<<"\t*cnfTotal";
                cout<<"\t"<<num_Cf;
            } // end if nVerbose>2
        } // end if nVerbose>1
        cout<<"\n";
    }
    return true;
}



// ============================================================================
// Author: Xin Cheng
// ----------------------------------------------------------------------------
// D1GC speculative iterative coloring, original version.
//
// while |W| is not empty
//     TentativeColoring
//     CheckConflicts
//     HandleConlficts
// 
//
// MP means hadleConflicts  will do nothing and leave uncolore vertices for next iterations
//
// original version means
//     no-prepartition,  
//     merge the conflicts 
//     does not using local orderings
// ============================================================================
int SMPGCColoring::D1GC_OMP_SpeculativeIterative_MP_Original(
        int nT,                  // number of threads
        int&colors,              // number of colors 
        vector<int>&vtxColors,   // mapping vertex to its color
        const int nVerbose       // verbose level
        ) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);
    
    vector<ChronoDuration> vtim_TC, vtim_CC;   // run time for debug
    vector<int>  vnum_Cf;                          // number of conflicts for debug

    const int N               = num_nodes();   // number of vertex
    const int BufSize         = max_degree()+1;// forbidden color buffer size
    const vector<int>& vtxPtr = get_CSR_ia();  // csr format edges
    const vector<int>& vtxVal = get_CSR_ja();  // csr format values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); // ordered vertex 

    colors=0;                       
    vtxColors.assign(N, -1);
    
    // allocate memory
    vector<int> Q(const_ordered_vertex);
    int Qsize=N;
    vector<int> cfQ(N);
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16, -1) ); 


    int n_Iter=0;
    do{

        // Tentative Coloring
        vtim_TC.emplace_back(.0);
        do_TentativeColoring_NoPrePartition(nT, vtxColors, vtxPtr, vtxVal, Q, Qsize, Fs, BufSize, vtim_TC[n_Iter]);

        // Check Conflicts
        vtim_CC.emplace_back(.0); vnum_Cf.emplace_back(0);
        do_CheckConflicts_NoPrePartition(nT, vtxColors, vtxPtr, vtxVal, Q, Qsize, cfQ, vnum_Cf[n_Iter], vtim_CC[n_Iter]);
    
        // Handle Conflicts
        cfQ.swap(Q);
        Qsize = vnum_Cf[n_Iter];
        n_Iter++;
    }while(Qsize!=0);

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@SpIMPOrig(None)_nT_c_T";
        for(int i=ss.str().size(); i<27; i++) ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        
        ChronoDuration tim_TC((ChronoDuration).0), tim_CC((ChronoDuration).0);
        for(auto & x: vtim_TC) tim_TC+=x;
        for(auto & x: vtim_CC) tim_CC+=x;
        auto tim_Total = tim_TC + tim_CC;

        cout<<"\t"<<tim_Total.count();

        if(nVerbose>1){ // show detail run time
            cout<<"\t*tTC_tCC_nCf";
            cout<<"\t"<<tim_TC.count();
            cout<<"\t"<<tim_CC.count();
            
            int sum_cf = 0 ;
            for(auto &x :vnum_Cf) sum_cf+=x;
            cout<<"\t"<<sum_cf;


            if(nVerbose>2){ // show the details for each iteration
                cout<<"\n*nCf Total "<<sum_cf<<" ["<<vnum_Cf.size()<<"]";
                for(auto &x : vnum_Cf) cout<<" "<<x;

                cout<<"\n*tTC Total "<<tim_TC.count()<<" ["<<vtim_TC.size()<<"]";
                for(auto &x : vtim_TC) cout<<" "<<x.count();

                cout<<"\n*tCC Total "<<tim_CC.count()<<" ["<<vtim_CC.size()<<"]";
                for(auto &x : vtim_CC) cout<<" "<<x.count();
            } // end if nVerbose>2
        } // end if nVerbose>1
        cout<<"\n";
    }
    return true;
}



// ============================================================================
// Author: Xin Cheng
// do tenetative coloring with original verion
//  original means no pre-partion, merge conflicts and no local ordering
// ============================================================================
void SMPGCColoring::do_TentativeColoring_NoPrePartition(
        const int nT,                  // number of threads
        vector<int>& vtxColors,        // mapping vertex to its color
        const vector<int>& vtxPtr,     // CSR format
        const vector<int>& vtxVal,     // CSR format
        const vector<int>& Q,          // vertex to color
        const int Qsize,               // number of vertex to color
        vector<vector<int>>& Fs,       // Forbidden color buffers
        const int BufSize,             // BufSize
        ChronoDuration &time           // run time
        ){
    auto start = steady_clock::now();
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int>&F = Fs[tid];
        F.assign(BufSize, -1);
        
        #pragma omp for
        for(int iu=0; iu<Qsize; iu++){
            const auto u = Q[iu];
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                const auto vc=vtxColors[vtxVal[iv]];
                if( vc >= 0) 
                    F[vc] = u;
            } 
            int c;
            for (c=0; c!=BufSize; c++)
                if(F[c]!=u)
                    break;
            vtxColors[u] = c;
        }// end omp for
    }// end omp parallel
    auto end = steady_clock::now();
    time = end - start;
}


// ============================================================================
// Author: Xin Cheng
// do check conflicts with original verion
//  original means no pre-partion, merge conflicts and no local ordering
// ============================================================================
void SMPGCColoring::do_CheckConflicts_NoPrePartition(
        const int nT,                      // number of threads
        vector<int>& vtxColors,            // mapping vertex to its color
        const vector<int>& vtxPtr,         // CSR format
        const vector<int>& vtxVal,         // CSR format
        const vector<int> &Q,              // vertex to color
        const int Qsize,                   // number of vertex to color
        vector<int>& cfQ,                  // conflict vertex
        int& cfQsize,                      // number of conflict vertex
        ChronoDuration &time              // run time
        ){
    auto start = steady_clock::now();
    cfQsize = 0;
    #pragma omp parallel for
    for(int iu=0; iu<Qsize; iu++) {
        const auto u  = Q[iu];
        const auto uc = vtxColors[u];
        for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){ 
            const auto v = vtxVal[iv];
            if(v>u && uc == vtxColors[v]) {
                vtxColors[u] = -1;  //Will prevent v from being in conflict in another pairing
                auto position =__sync_fetch_and_add(&cfQsize, 1); //increment the counter
                cfQ[position] = u;
                break;
            }  
        }// end for neighbor v
    }// end omp parallel for iu
    auto end = steady_clock::now();
    time = end - start;
}






