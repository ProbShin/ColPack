
/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#include "SMPGCColoring.h"
#include <unordered_set> // calc_num_colors_from_vtx_colors
using namespace std;
using namespace ColPack;






// ============================================================================
// author: xin cheng
// get number of distinc colors used in vtxColor.
// in case incomplete coloring, i.e. some color is -1, ignore color with -1.
// ============================================================================
int SMPGCColoring::calc_num_colors_from_vtx_colors(const vector<int>&vtxColor){
    unordered_set<int> S(vtxColor.begin(), vtxColor.end());
    return S.count(-1)?(S.size()-1):(S.size());
}



// ============================================================================
// author: xin cheng
// do local ordering 
// for parallel coloring
// ============================================================================
/*void SMPGCColoring::do_local_ordering(const int LOCAL_ORDER, vector<vector<int>>&QQ, ChronoDuration& time) {
    switch(LOCAL_ORDER) {
        case ORDER_NONE: break;
        case ORDER_LARGEST_FIRST: omp_prepartition_ldf_ordering(QQ, &time); break;
        case ORDER_SMALLEST_LAST: omp_prepartition_sdl_ordering(QQ, &time); break;
        case ORDER_NATURAL:       omp_prepartition_ntr_ordering(QQ, &time); break;
        case ORDER_RANDOM:        omp_prepartition_rnd_ordering(QQ, &time); break;
        default:
            printf("Error! unknown local order \"%d\".\n", LOCAL_ORDER);
            exit(1);           
    }
}
*/


// ============================================================================
// check if the graph is correct colored
// ============================================================================
int SMPGCColoring::cnt_d1conflict(const vector<int>& vtxColorConst, bool bVerbose){
    vector<int>         vtxColor(vtxColorConst);
    const int N         = num_nodes();
    const vector<int>& vtxPtr = get_CSR_ia();
    const vector<int>& vtxVal = get_CSR_ja();
   
    if((signed)vtxColorConst.size()!=N){
        cout<<"Error! vColors.size()!= |E| \n";
        exit(1);
    }

    int n_uncolored=0;
    int n_conflicts=0;
    #pragma omp parallel reduction(+:n_conflicts), reduction(+:n_uncolored)
    {
        #pragma omp for
        for(int v=0; v<N; v++) {
            const int vc=vtxColor[v];
            if(vc<0) {
                n_uncolored++;
                continue;
            }
            for(int iw=vtxPtr[v]; iw!=vtxPtr[v+1]; iw++) {
                const int w = vtxVal[iw];
                if(v>=w) continue; // only check one side
                if(vc == vtxColor[w] ) {
                    vtxColor[v]=-1; // prevent further conflicts, however, since no synchronize used. May count more conflicts than actual.
                    n_conflicts++;
                    break;  
                }
            }
        }
    }
    if(bVerbose && n_uncolored) printf("There are %d vertex uncolored\nThere are %d vertex has conflicts with other nodes.\n",n_uncolored, n_conflicts);
    return n_uncolored+n_conflicts;
}

// ============================================================================
// check the graph validation
// ----------------------------------------------------------------------------
// uncolored vertex will not conflict with any other vertex
// ============================================================================
int SMPGCColoring::cnt_d2conflict(const vector<int>&vtxColorConst, bool bVerbose) {
    //do it serial
    if(0)
    {
        vector<int> vtxColor(vtxColorConst);
        const int N = num_nodes();
        const vector<int>& vtxPtr = get_CSR_ia();
        const vector<int>& vtxVal = get_CSR_ja();
        
        vector<int>                            uncolored_nodes;
        unordered_map<int, unordered_set<int>> conflicts_nodes;

        for(int v=0; v<N; v++){
            const auto vc = vtxColor[v];
            if(vc<0) { uncolored_nodes.push_back(v); continue; }
            for(int iw=vtxPtr[v]; iw!=vtxPtr[v+1]; iw++){  // check d1 neighbors
                const auto w = vtxVal[iw];
                if( vc==vtxColor[w] ) 
                    conflicts_nodes[ min(v,w) ].insert(max(v,w));
            }
            for(int iw=vtxPtr[v]; iw!=vtxPtr[v+1]; iw++) { 
                const auto w = vtxVal[iw];
                for(int iu=vtxPtr[w]; iu!=vtxPtr[w+1]; iu++){  // check d2 neighbors
                    const auto u = vtxVal[iu];
                    if(v==u) continue;
                    if( vc == vtxColor[u])  
                        conflicts_nodes[min(v,u)].insert(max(v,u));
                }
            }
        }

        if(bVerbose) {
            printf("There is %d vertex uncolored\nThere is %d vertex conflicts with other nodes.\n", 
                    (int)uncolored_nodes.size(), (int)conflicts_nodes.size());
        }
        if(!uncolored_nodes.empty()){
            printf("uncolored_nodes[%d]: ", (int)uncolored_nodes.size());
            for(int i=0; i<min((int)uncolored_nodes.size(), 10); i++) 
                printf("\t%d", uncolored_nodes[i]);
            printf("\n");
        }
        if(!conflicts_nodes.empty()){
            printf("conflicts_nodes[%d]:\n", (int)conflicts_nodes.size());
            int cnt_rows=0;
            for(const auto &x : conflicts_nodes){
                if(cnt_rows++>10) { 
                    printf("...");
                    break;
                }
                printf("[%d(%d)]:", x.first, vtxColor[x.first]);
                int cnt_cols=0;
                for(const auto &y : x.second) {
                    if(cnt_cols++>10) {
                        printf("...");
                        break;
                    }
                    printf("\t%d(%d)",y, vtxColor[y]);
                }
                printf("\n");
            }
            printf("\n");
        }
        return uncolored_nodes.size()+conflicts_nodes.size();
    }
    // do it in parallel
    vector<int>  vtxColor( vtxColorConst );
    const int N = num_nodes();
    const vector<int>& vtxPtr = get_CSR_ia();
    const vector<int>& vtxVal = get_CSR_ja();
    int   n_conflicts = 0;
    int   n_uncolored = 0;
    
    #pragma omp parallel reduction(+:n_conflicts), reduction(+: n_uncolored)
    {
        #pragma omp for
        for(int v=0; v<N; v++){
            const auto vc = vtxColor[v];
            if(vc<0) { n_uncolored++; continue; }
            bool b_visbad = false;
            for(int iw=vtxPtr[v]; iw!=vtxPtr[v+1]; iw++){  // check d1 neighbors
                const auto w = vtxVal[iw];
                if( v>=w ) continue;   // only check one side
                if( vc==vtxColor[w] ) { 
                    n_conflicts ++;
                    vtxColor[v]=-1;
                    b_visbad = true;
                    break;
                }
            }
            for(int iw=vtxPtr[v]; b_visbad==false && iw!=vtxPtr[v+1]; iw++) { 
                const auto w = vtxVal[iw];
                for(int iu=vtxPtr[w]; iu!=vtxPtr[w+1]; iu++){  // check d2 neighbors
                    const auto u = vtxVal[iu];
                    if(v >= u) continue;
                    if( vc == vtxColor[u]) {
                        n_conflicts ++;
                        vtxColor[v] =-1;
                        b_visbad=true;
                    }
                }
            }
        }
    }
    if(bVerbose) {
        printf("There is %d uncolored vertices.\nThere is %d vertices conflict with other nodes.\n", (int)n_uncolored, (int)n_conflicts);
    }
    return n_uncolored + n_conflicts;
}













