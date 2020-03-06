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
// D1GC serial greedy coloring
// ============================================================================
int SMPGCColoring::D1_serial(int&colors, vector<int>&vtxColors, const int nVerbose) {
    omp_set_num_threads(1);
    
    ChronoDuration tim_LocalOrder;
    ChronoDuration tim_Coloring;
    
    const int N               = num_nodes();   //number of vertex
    const int BufSize         = max_degree()+1;
    const vector<int>& vtxPtr = get_CSR_ia();
    const vector<int>& vtxVal = get_CSR_ja();
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); 

    colors=0;                       
    vtxColors.assign(N, -1);

    // allocate memory
    vector<int> Q(const_ordered_vertex);  //copied to local memory
    vector<int> F; F.reserve(BufSize);

    // coloring
    {
        auto start = steady_clock::now();
        F.assign(BufSize, -1);
        for(const auto u : Q){
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                const auto vc=vtxColors[vtxVal[iv]];
                if( vc >= 0) 
                    F[vc] = u;
            } 
            int c=0;
            for (; c!=BufSize; c++)
                if(F[c]!=u)
                    break;
            vtxColors[u] = c;
        }
        auto end = steady_clock::now();
        tim_Coloring = end - start;
    } 

    colors = calc_num_colors_from_vtx_colors(vtxColors);
    
    if(nVerbose>0){
        stringstream ss;
        ss<<"@D1Serial_1_c_T";
        for(int i=ss.str().size(); i<27; i++) ss<<"_";
        cout<<ss.str();
        cout<<"\t1";
        cout<<"\t"<<colors;
        auto tim_Total = tim_LocalOrder + tim_Coloring;
        cout<<"\t"<<tim_Total.count();
        if(nVerbose>1){
            cout<<"\t@tLO_tC";
            cout<<"\t"<<tim_LocalOrder.count();
            cout<<"\t"<<tim_Coloring.count();
        }
        cout<<endl;
    }
    return true;   
}

