/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "SMPGCColoring.h"
using namespace std;
using namespace ColPack;



// ============================================================================
// Author: xin cheng
// Implements of the D2GC omp Speculative Iteration approaches'
// CheckConflicts phase
//      pre-partition the vertex 
//      without using atomic operation
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseCC(
        vector<vector<int>>& QQ,
        vector<int>& Qsizes,
        vector<int>& vtxColors,
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        int *pNum,
        ChronoDuration *pTime
        ){
    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();

    int num_uncolored_vtx_total=0;
    #pragma omp parallel reduction(+:num_uncolored_vtx_total)
    {
        const int tid = omp_get_thread_num();
        vector<int>&Q = QQ[tid];
        const int Qsize = Qsizes[tid];
        int conflicts_size = 0;
        for(int iu=0; iu<Qsize; iu++){
            const int u = Q[iu];
            const int uc = vtxColors[u];
            bool b_uis_false = false;
            // distance 1 neighbor
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){
                const int v = vtxVal[iv];
                if(v<u && vtxColors[v]==uc){
                    vtxColors[u]=-1;
                    b_uis_false = true;
                    Q[conflicts_size++] = u;
                    break;
                }
            }
            if(b_uis_false == false){
                // distance 2 neighbor
                for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                    const int v = vtxVal[iv];
                    for(int iw=vtxPtr[v]; iw!=vtxPtr[v+1]; iw++) {
                        const int w = vtxVal[iw];
                        if(w<u && vtxColors[w]==uc){
                            vtxColors[u]=-1;
                            b_uis_false = true;
                            Q[conflicts_size++] = u;
                            break;
                        }
                    }// end for w
                    if(b_uis_false)
                        break;
                }// end for v
            }
        }// end for u
        Qsizes[tid] = conflicts_size;
        num_uncolored_vtx_total = conflicts_size;
    }//end omp parallel
    
    if(pTime) *pTime = steady_clock::now() - start;
    if(pNum ) *pNum  = num_uncolored_vtx_total;
}


// ============================================================================
// Author: xin cheng
// Implements of the D2GC omp Speculative Iteration approaches'
// CheckConflicts phase
//      using Net-Based Algorithm
//
// Note: the Net-Based Algorithm is based on the Algorithm 10 of the Paper
// "Greedy is Good : Parallel Alg..." by M. K. Tas ̧ K. Kaya and E. Saule
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseNetBinCC(
        vector<vector<int>> & QQ,
        vector<int>& Qsizes,
        vector<int>& vtxColors,
        const vector<int> &vtxPtr,
        const vector<int> &vtxVal,
        vector<int> const &const_ordered_Q,
        int const N,
        vector<vector<int>> &ForbiddenArrays,
        const int BufSizes,
        int * pNum, 
        ChronoDuration* pTime
        ){

    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    // sub-phase - check conflicts
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int> &F = ForbiddenArrays[tid];
        F.assign(BufSizes, -1);
        #pragma omp for
        for(int iv=0; iv<N; iv++){
            const int v = const_ordered_Q[iv];
            const int vc = vtxColors[v];
            if(vc>=0)
                F[vc]=v;
            for(int iu=vtxPtr[v]; iu!=vtxPtr[v+1]; iu++) {
                const int u = vtxVal[iu];
                const int uc= vtxColors[u];
                if(uc>=0){
                    if( F[uc]==v )
                        vtxColors[u] = -1;
                    else
                        F[uc]=v;
                }
            }
        }
    }
    // sub-phase - collect and queue the uncolored vertex
    int num_uncolored_vtx_total=0;
    #pragma omp parallel reduction(+:num_uncolored_vtx_total)
    {
        const int tid = omp_get_thread_num();
        vector<int> &Q = QQ[tid];
        int tmpQsize = 0; //num_uncolored_vtx_total = 0 ;
        #pragma omp for
        for(int iu=0; iu<N; iu++){
            int const u = const_ordered_Q[iu];
            if( vtxColors[u] < 0)
                Q[tmpQsize++] = u;
        }
        Qsizes[tid] = tmpQsize; //num_uncolored_vtx_total;
        num_uncolored_vtx_total = tmpQsize;
    }
    if(pTime) *pTime = steady_clock::now() - start;
    if(pNum) *pNum = num_uncolored_vtx_total;
}






// ============================================================================
// author xin cheng
// D2GC parallel check conflicts, NoRePartition
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseCC_NoPrePartition(
        vector<int>& Q,
        int& Qsize,
        vector<int>& vtxColors,
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        vector<int>& cfQ,
        int * pNum,
        ChronoDuration* pTime
        ){
    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    int cfQsize = 0;
    #pragma omp parallel
    {
        #pragma omp for
        for(int iu=0; iu<Qsize; iu++){
            const int u = Q[iu];
            const int uc = vtxColors[u];
            bool b_uis_false = false;
            // distance 1 neighbor
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){
                const int v = vtxVal[iv];
                if(v<u && vtxColors[v]==uc){
                    vtxColors[u]=-1;
                    b_uis_false = true;
                    int enQposition = __sync_fetch_and_add(&cfQsize, 1); 
                    cfQ[enQposition] = u;
                    break;
                }
            }
            if(b_uis_false == false){
                // distance 2 neighbor
                for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                    const int v = vtxVal[iv];
                    for(int iw=vtxPtr[v]; iw!=vtxPtr[v+1]; iw++) {
                        const int w = vtxVal[iw];
                        if(w<u && vtxColors[w]==uc){
                            vtxColors[u]=-1;
                            b_uis_false = true;
                            int enQposition = __sync_fetch_and_add(&cfQsize, 1);
                            cfQ[enQposition] = u;
                            break;
                        }
                    }// end for w
                    if(b_uis_false)
                        break;
                }// end for v
            }
        }// end for u
    }//end omp parallel
    
    Q.swap(cfQ);
    Qsize = cfQsize;
    
    if(pTime) *pTime = steady_clock::now() - start;
    if(pNum)  *pNum  = cfQsize;
}



// ============================================================================
// author xin cheng
// 
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseNetBinCC_NoPrePartition(
        vector<int> & Q,
        int & Qsize,
        vector<int>& vtxColors,
        const vector<int> &vtxPtr,
        const vector<int> &vtxVal,
        vector<int> const & const_ordered_Q,
        int const N,
        vector<vector<int>> &ForbiddenArrays,
        const int BufSizes,
        int * pNum,
        ChronoDuration * pTime
        ){
    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    // phase - check conflicts
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int> &F = ForbiddenArrays[tid];
        F.assign(BufSizes, -1);
        
        #pragma omp for
        for(int iv=0; iv<N; iv++){
            const int v = const_ordered_Q[iv];
            const int vc = vtxColors[v];
            if(vc>=0)
                F[vc]=v;
            for(int iu=vtxPtr[v]; iu!=vtxPtr[v+1]; iu++) {
                const int u = vtxVal[iu];
                const int uc= vtxColors[u];
                if(uc>=0){
                    if( F[uc]==v )
                        vtxColors[u] = -1;
                    else
                        F[uc]=v;
                }
            }
        }
    }
    
    // phase - handle conflicts.  generate next iteartion workload
    Qsize = 0;
    #pragma omp parallel for
    for(int i=0; i<N; i++){
        int const u = const_ordered_Q[i];
        if(vtxColors[u]<0){
            int enQposition = __sync_fetch_and_add(&Qsize, 1);
            Q[enQposition]=u;
        }
    }
    if(pTime) *pTime = steady_clock::now() - start;
    if(pNum ) *pNum = Qsize;
}


