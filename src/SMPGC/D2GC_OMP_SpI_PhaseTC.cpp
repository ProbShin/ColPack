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
// impements of the D2GC omp Speculative Iteration approach's 
// TentativeColoring Phase
//
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseTC(
        vector<int>&vtxColors,         // mapping vertex to its color
        const vector<vector<int>>& QQ, // vertex to color for each thread
        const vector<int>& Qsizes,     // number of the vertex to color for each thread
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        vector<vector<int>>& Fs,       // forbidden color array
        const int BufSize,             // forbidden color array size
        ChronoDuration *ptime       // runtime
        ) {
    chrono::steady_clock::time_point start;
    if(ptime) start = steady_clock::now();
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int> const &Q = QQ[tid];
        const int Qsize = Qsizes[tid];
        vector<int>& F = Fs[tid];
        do_D2GC_Seq_GreedyColoring(vtxPtr, vtxVal, vtxColors, Q, Qsize, F, BufSize);
    }// end of omp parallel
    if(ptime) *ptime = steady_clock::now()-start;
}


// ============================================================================
// Author: Xin Cheng
// Implements of the D2GC omp Speculative Iterative approaches
// TentativeColoring Phase 
//     using Net-Based Algorithm
// 
// Note: the Net-Based Algorithm is based on the Algorithm 9 of the Paper
// "Greedy is Good : Parallel Alg..." by M. K. Tas ̧ K. Kaya and E. Saule
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseNetBinTC(
        vector<int>& vtxColors,             
        vector<int> const & const_ordered_Q,
        int const N,
        const vector<int>&vtxPtr,
        const vector<int>&vtxVal,
        vector<vector<int>>& Fs,
        const int BufSize,
        vector<vector<int>>& Ws, //work load vertices
        ChronoDuration *pTime   
        ){
    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int>&F = Fs[tid];
        vector<int>&W = Ws[tid];
        F.assign(BufSize, -1);

        #pragma omp for
        for(int iv=0; iv<N; iv++){
            int Wsize = 0;
            const int v = const_ordered_Q[iv];
            const int vc = vtxColors[v];
            if(vc>=0)
                F[vc]=v;
            else
                W[Wsize++]=v;
            for(int iu=vtxPtr[v]; iu!=vtxPtr[v+1]; iu++){
                const int u = vtxVal[iu];
                const int uc = vtxColors[u];
                if(uc>=0 and F[uc]!=v)
                    F[uc] = v;
                else
                    W[Wsize++] = u;
            }
            int col = vtxPtr[v+1] - vtxPtr[v] ;
            for(int iu=0; iu<Wsize; iu++){
                while( F[col]== v)
                    col--;
                vtxColors[W[iu]]=col;
                col--;
            }
        }// end of omp for
    }// end of omp parallel
    if(pTime) *pTime = steady_clock::now()-start;
}


// ============================================================================
// Author: xin cheng
// Implemnets of the D2GC omp Speculative Iterative approach 
// TentativeColoring Phase
//     using random coloring
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseRandomTC(
        vector<int> &vtxColors,          // mapping vertex to its color
        const vector<vector<int>>&QQ,    // vertex to color of each thread
        const vector<int>& Qsizes,       // number of the vertex of each thread
        vector<mt19937>& mts,            // random number generator of each thread
        const int ColorLB0Base,          // random number range
        ChronoDuration *pTime            // run time
        ){
    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int> const& Q = QQ[tid];
        const int Qsize  = Qsizes[tid];
        mt19937&  mt = mts[tid];
        uniform_int_distribution<int> dist_local(0, ColorLB0Base);
        for(int iu=0; iu<Qsize; iu++){
            vtxColors[ Q[iu] ] = dist_local(mt);
        }
    }// end of omp parallel
    if(pTime) *pTime = steady_clock::now() - start;
}




// ============================================================================
// author xin cheng
// D2GC tentative coloring, seperat queue, no repartion needed 
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseTC_NoPrePartition(
        vector<int> &vtxColors,
        vector<int> & Q,
        int & Qsize,
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        vector<vector<int>> &ForbiddenArrays,
        const int BufSize,
        ChronoDuration *pTime
        ){
    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int>& F  = ForbiddenArrays[tid];
        F.assign(BufSize,-1);
        
        #pragma omp for
        for(int iu=0; iu<Qsize; iu++){
            const int u = Q[iu];
            // d1 neighbors
            for(int iv = vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                const auto vc = vtxColors[ vtxVal[iv] ];
                if(vc>=0) 
                    F[vc] = u;
            }
            // d2 neighbors
            for(int iv = vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                const auto v = vtxVal[iv];
                for(int iw = vtxPtr[v]; iw!=vtxPtr[v+1]; iw++) {
                    const auto w = vtxVal[iw];
                    if(v==w) continue;
                    const auto wc = vtxColors[w];
                    if(wc>=0)
                        F[wc] = u;
                }
            }
            // coloring
            int c=0;
            for(; c!=BufSize; c++){
                if(F[c]!=u)
                    break;
            }
            vtxColors[u]=c;
        }//end for
    }
    if(pTime) *pTime = steady_clock::now() - start;
}



// ============================================================================
// Author: Xin Cheng
//
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseNetBinTC_NoPrePartition(
        vector<int> &vtxColors,
        vector<int> const & const_ordered_Q,
        int const N,
        vector<int> const & vtxPtr,
        vector<int> const & vtxVal,
        vector<vector<int>> &Fs,
        int const BufSize,
        vector<vector<int>> &W,
        ChronoDuration *pTime
        )
{
    do_D2GC_OMP_PhaseNetBinTC(vtxColors, const_ordered_Q, N, vtxPtr, vtxVal, Fs, BufSize, W, pTime)    ;
}




// ============================================================================
// author xin cheng
// D2GC tentative coloring with random coloring, without re-partition
// ============================================================================
void SMPGCColoring::do_D2GC_OMP_PhaseRandomTC_NoPrePartition(
        vector<int> &vtxColors,
        vector<int> const & Q,
        int const QSize, 
        vector<mt19937>& mts,
        int const ColorLB0Base,
        ChronoDuration *pTime
        ){
    chrono::steady_clock::time_point start;
    if(pTime) start = steady_clock::now();
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        mt19937& mt = mts[tid];
        uniform_int_distribution<int> dist_local(0, ColorLB0Base);
        #pragma omp for
        for(int iu=0; iu<QSize; iu++){
            vtxColors[ Q[iu] ] = dist_local(mt);
        }
    }// end of omp parallel
    if(pTime) *pTime = steady_clock::now() - start;
}






