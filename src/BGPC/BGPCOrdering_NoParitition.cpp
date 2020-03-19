/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "BGPCOrdering.h"
#include <ctime>  //clock
using namespace std;
using namespace ColPack;


// ============================================================================
// Author: Xin Cheng 
// ============================================================================
void BGPCOrdering::do_Natural_Ordering_NoPartition(vector<int>&Q, const int Qsize, ChronoDuration* tim){
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    #pragma omp parallel for
    for(int i=0; i<Qsize; i++)
        Q[i] = i;
    if(tim) *tim = chrono::steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng 
// Random is shuffle to natural
// ============================================================================
void BGPCOrdering::do_Random_Ordering_NoPartition(vector<int>&Q, const int Qsize, ChronoDuration* tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    if(Qsize<=1) return;
    int const Qsizem1=Qsize-1;
    for(int i=0; i<Qsizem1; i++){
        uniform_int_distribution<int> dist(i, Qsizem1);
        swap(Q[i], Q[dist(m_mt)]);
    }
    if(tim) *tim = chrono::steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng 
// Largest Degree First
// ============================================================================
void BGPCOrdering::do_D1_LargestDegreeFirst_Ordering_NoPartition(
        vector<int>&Q, 
        const int Qsize,
        vector<int>const&srcPtr,
        ChronoDuration* tim){
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int current_max_degree = -1;
    vector<vector<int>> GroupedVertexDegree;
    for(int iu=0; iu<Qsize; iu++) {
        int const u = Q[iu];
        int const ud= -srcPtr[u] +  srcPtr[u+1];
        if(current_max_degree<ud){
            current_max_degree = ud;
            GroupedVertexDegree.resize(ud+1);
        }
        GroupedVertexDegree[ud].push_back(u);
    }
    Q.clear();
    for(int d=current_max_degree, it=current_max_degree+1; it!=0; it--, d--){
        Q.insert(Q.end(), GroupedVertexDegree[d].begin(), GroupedVertexDegree[d].end());
    }
    GroupedVertexDegree.clear();
    if(tim) *tim = chrono::steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng 
// Largest Degree First
//
// the number of unique Distance 2 neighbors
// ============================================================================
void BGPCOrdering::do_D2_LargestDegreeFirst_Ordering_NoPartition(
        vector<int>&Q, 
        const int Qsize,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        ChronoDuration* tim){
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int current_max_degree = -1;
    vector<vector<int>> GroupedVertexDegree;
    for(int iu=0; iu<Qsize; iu++) {
        int const u = Q[iu];
        unordered_set<int> uniq_d2_nbs;
        for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
            int const v = vtxVal[iv];
            for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                int const w = vtxVal[iw];
                uniq_d2_nbs.insert(w);
            }
        }
        // d2 unique degree
        int const ud = uniq_d2_nbs.size();
        if(current_max_degree<ud){
            current_max_degree = ud;
            GroupedVertexDegree.resize(ud+1);
        }
        GroupedVertexDegree[ud].push_back(u);
    }
    Q.clear();
    for(int d=current_max_degree, it=current_max_degree+1; it!=0; it--, d--){
        Q.insert(Q.end(), GroupedVertexDegree[d].begin(), GroupedVertexDegree[d].end());
    }
    GroupedVertexDegree.clear();
    if(tim) *tim = chrono::steady_clock::now() - start;
}


// ============================================================================
// Author: Xin Cheng 
// Largest Degree First
//
// Define of the Degree:
// The sum degree of neighbors in B side.
// ============================================================================
void BGPCOrdering::do_D2Rough_LargestDegreeFirst_Ordering_NoPartition(
        vector<int>&Q, 
        const int Qsize,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        ChronoDuration* tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int current_max_degree = -1;

    vector<vector<int>> GroupedVertexDegree;
    for(int iu=0; iu<Qsize; iu++) {
        int const u = Q[iu];
        int ud = 0;
        for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
            int const v = vtxVal[iv];
            int const vd = -dstPtr[v] + dstPtr[v+1];
            ud += vd;
        }
        if(current_max_degree<ud){
            current_max_degree = ud;
            GroupedVertexDegree.resize(ud+1);
        }
        GroupedVertexDegree[ud].push_back(u);
    }
    Q.clear();
    for(int d=current_max_degree, it=current_max_degree+1; it!=0; it--, d--){
        Q.insert(Q.end(), GroupedVertexDegree[d].begin(), GroupedVertexDegree[d].end());
    }
    GroupedVertexDegree.clear();
    if(tim) *tim = chrono::steady_clock::now() - start;
}


// ============================================================================
// Author: Xin Cheng 
// Largest Degree First
//
// Define of the Degree:
// max degree of the neighbors in B side, and if tie, using its own degree.
// ============================================================================
void BGPCOrdering::do_BA_LargestDegreeFirst_Ordering_NoPartition(
        vector<int>&Q, 
        const int Qsize,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        ChronoDuration* tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int current_max_degree = -1;
    vector<vector<pair<int,int>>> GroupedVertexDegree;
    for(int iu=0; iu<Qsize; iu++) {
        int const u = Q[iu];
        int ud = -srcPtr[u] + srcPtr[u+1];
        int max_vd=0;
        for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
            int const v = vtxVal[iv];
            int const vd = -dstPtr[v] + dstPtr[v+1];
            if(max_vd <vd) 
                max_vd = vd;
        }
        if(current_max_degree<max_vd){
            current_max_degree = max_vd;
            GroupedVertexDegree.resize(max_vd+1);
        }
        GroupedVertexDegree[max_vd].emplace_back(ud,u);
    }
    Q.clear();
    for(int d=current_max_degree, it=current_max_degree+1; it!=0; it--, d--){
        sort(GroupedVertexDegree[d].begin(), GroupedVertexDegree[d].end());
        Q.reserve(Q.size()+GroupedVertexDegree[d].size());
        for(int i=0, iEnd=GroupedVertexDegree[d].size(); i<iEnd; i++)
            Q.push_back( GroupedVertexDegree[d][i].second );
        //Q.insert(Q.end(), GroupedVertexDegree[d].begin(), GroupedVertexDegree[d].end());
    }
    GroupedVertexDegree.clear();
    if(tim) *tim = chrono::steady_clock::now() - start;
}






