/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "BGPCColoring.h"
#include <chrono> //c++11 system time
#include <random> //c++11 random
using namespace std;
using namespace ColPack;
using namespace chrono;

// ============================================================================
// author xin cheng
// Bipartite Graph Partial Coloring Serial Coloring
// ============================================================================
int BGPCColoring::BGPC_SEQ_Greedy(
        const int side, 
        int& colors, 
        vector<int>& vtxColors, 
        const int nVerbose
        ) {
    ChronoDuration tim_color{.0};    //double tim_color = 0;
    
    const int          N             = (side==BGPC::L)?(GetLeftVertexCount()):(GetRightVertexCount()); 
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    const vector<int>& vtxVal        = GetEdges();
    const int          srcMaxDegree  = (side==BGPC::L)?GetMaximumLeftVertexDegree():GetMaximumRightVertexDegree();
    const int          dstMaxDegree  = (side==BGPC::L)?GetMaximumRightVertexDegree():GetMaximumLeftVertexDegree();
    const int          DDp1          = dstMaxDegree * srcMaxDegree + 1;
    const int          BufSize       = DDp1<0?N:min(N, DDp1);
   
    vector<int>const& const_queue_A  = get_ordered_queue_A_const();
    //vector<int>const& const_queue_B;

    colors=0;                       
    vtxColors.assign(N, -1);

    // allocate memory
    vector<int> Mask(BufSize,-1);
    
    // coloring
    do_BGPC_Seq_Greedy(vtxColors, const_queue_A, N, srcPtr, dstPtr, vtxVal, Mask, BufSize, &tim_color);
    
    colors = calc_num_colors_from_vtx_colors(vtxColors);
    if(nVerbose>0){
        stringstream ss;
        ss<<"@BGPC(";
        ss<<((side==BGPC::L)?"L":"R");
        ss<<")";
        for(int i=ss.str().size(); i<30; i++)
            ss<<"_";
        cout<<ss.str();
        cout<<"\t#nT_c_T\t1\t"<<colors<<"\t"<<tim_color.count();
        cout<<endl;
    }
    return true;   
}


// ============================================================================
// Author: Xin Cheng
// 
// ============================================================================
void BGPCColoring::do_BGPC_Seq_Greedy(
        vector<int>& vtxColors,
        vector<int>const& Q,
        int const Qsize,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        vector<int>& F,
        int const BufSize,
        ChronoDuration *tim
        ) {
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    F.assign(BufSize, -1);
    for(int iu=0; iu<Qsize; iu++){    //u-v-w
        int const u = Q[iu];
        for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++){
            const auto v = vtxVal[iv];
            for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++){
                const auto w = vtxVal[iw];
                if(w==u) continue;
                const auto wc = vtxColors[w];
                if(wc<0) continue;
                F[wc] = u;
            }
        }
        int c = 0;
        for(; c<BufSize; c++)
            if(F[c]!=u) 
                break;
        vtxColors[u]=c;
    }
    if(tim) *tim = steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng
// 
// ============================================================================
void BGPCColoring::do_BGPC_Seq_Greedy_NoInitF(
        vector<int>& vtxColors,
        vector<int>const& Q,
        int const Qsize,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        vector<int> &F,
        int const BufSize,
        ChronoDuration *tim
        ) {
    chrono::steady_clock::time_point start;
    if(tim) start = steady_clock::now();
    
    for(int iu=0; iu<Qsize; iu++){    //u-v-w
        int const u = Q[iu];
        for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++){
            const auto v = vtxVal[iv];
            for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++){
                const auto w = vtxVal[iw];
                if(w==u) continue;
                const auto wc = vtxColors[w];
                if(wc<0) continue;
                F[wc] = u;
            }
        }
        int c = 0;
        for(; c<BufSize; c++)
            if(F[c]!=u) 
                break;
        vtxColors[u]=c;
    }
    if(tim) *tim = steady_clock::now() - start;
}




