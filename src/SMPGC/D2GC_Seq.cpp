/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "SMPGCColoring.h"
#include <unordered_set>
#include <unordered_map>
using namespace std;
using namespace ColPack;

// ============================================================================
// Author: Xin Cheng
// D2GC Sequential Greedy Coloring the graph
// ============================================================================
int SMPGCColoring::D2GC_Serial(int&colors, vector<int>& vtxColors, const int nVerbose) {
    omp_set_num_threads(1);
    ChronoDuration tim_Coloring;
    
    const int N = num_nodes();                  // number of vertex
    int const D = max_degree();
    int const DDmDp1 = D*D-D+1;
    const int BufSize = DDmDp1<=0?N:min(DDmDp1, N); //maxDegree
    const vector<int>& vtxPtr = get_CSR_ia();   // csr format
    const vector<int>& vtxVal = get_CSR_ja();   // csr format
    
    colors=0;                       
    vtxColors.assign(N, -1);
    
    // allocate memory
    vector<int> Q(get_ordered_vertex());
    vector<int> F; F.reserve(BufSize+1);

    // Seq Coloring
    do_D2GC_Seq_GreedyColoring(vtxPtr, vtxVal, vtxColors, Q, N, F, BufSize, &tim_Coloring);
    
    // calculate number of colors
    colors = calc_num_colors_from_vtx_colors(vtxColors);

    // display for the debug
    if(nVerbose>0){
        stringstream ss;
        ss<<"@D2GCSeq";
        for(int i=ss.str().size(); i<30; i++) ss<<"_";
        ss<<"\t#nT_c_T\t1";
        ss<<"\t"<<colors;
        ss<<"\t"<<tim_Coloring.count();
        ss<<"\n";
        cout<<ss.str();
    }
    return true;
}


// ============================================================================
// Author: Xin Cheng
// D2GC Sequential greedy coloring implemenation
// ============================================================================
void SMPGCColoring::do_D2GC_Seq_GreedyColoring(
        vector<int> const& vtxPtr,     // csr format
        vector<int> const& vtxVal,     // csr format
        vector<int>& vtxColors,        // mapping vertex to its color
        const vector<int>& Q,          // vector to color
        const int Qsize,               // number of the vector to color
        vector<int>& F,                // forbidden color array
        const int BufSize,             // forbidden color array size
        ChronoDuration *ptime          // run time
        ) {
    chrono::steady_clock::time_point start;
    if(ptime!=nullptr) { start = steady_clock::now(); }
    F.assign(BufSize,-1);              // init buffer
    for(int iu=0; iu<Qsize; iu++) {
        const auto u = Q[iu];
        // distance one neighbor
        for(int iv=vtxPtr[u]; iv<vtxPtr[u+1]; iv++) {
            const auto vc = vtxColors[vtxVal[iv]];
            if(vc >= 0) 
                F[vc] = u;
        }
        // distance two neighbor
        for(int iv=vtxPtr[u]; iv<vtxPtr[u+1]; iv++) {
            const auto v = vtxVal[iv];
            for(int iw=vtxPtr[v]; iw<vtxPtr[v+1]; iw++) {
                const auto wc = vtxColors[ vtxVal[iw] ];
                if(wc >= 0) 
                    F[wc] = u;
            }
        }
        // color based on the forbidden color
        int c;
        for(c=0; c<BufSize; c++)
            if(F[c]!=u)
                break;
        vtxColors[u]=c;
    }
    if(ptime!=nullptr) { *ptime = steady_clock::now()-start; }
}


// ============================================================================
// Author: Xin Cheng
// D2GC Sequential greedy coloring implemenation
// ============================================================================
void SMPGCColoring::do_D2GC_Seq_GreedyColoring_NoInitF(
        vector<int> const& vtxPtr,     // csr format
        vector<int> const& vtxVal,     // csr format
        vector<int>& vtxColors,        // mapping vertex to its color
        const vector<int>& Q,          // vector to color
        const int Qsize,               // number of the vector to color
        vector<int>& F,                // forbidden color array
        const int BufSize,             // forbidden color array size
        ChronoDuration *ptime          // run time
        ) {
    chrono::steady_clock::time_point start;
    if(ptime!=nullptr) { start = steady_clock::now(); }
    for(int iu=0; iu<Qsize; iu++) {
        const auto u = Q[iu];
        // distance one neighbor
        for(int iv=vtxPtr[u]; iv<vtxPtr[u+1]; iv++) {
            const auto vc = vtxColors[vtxVal[iv]];
            if(vc >= 0) 
                F[vc] = u;
        }
        // distance two neighbor
        for(int iv=vtxPtr[u]; iv<vtxPtr[u+1]; iv++) {
            const auto v = vtxVal[iv];
            for(int iw=vtxPtr[v]; iw<vtxPtr[v+1]; iw++) {
                const auto wc = vtxColors[ vtxVal[iw] ];
                if(wc >= 0) 
                    F[wc] = u;
            }
        }
        // color based on the forbidden color
        int c;
        for(c=0; c<BufSize; c++)
            if(F[c]!=u)
                break;
        vtxColors[u]=c;
    }
    if(ptime!=nullptr) { *ptime = steady_clock::now()-start; }
}





