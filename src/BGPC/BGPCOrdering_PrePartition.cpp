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
void BGPCOrdering::do_Nature_Ordering_PrePartition(
        vector<vector<int>>& QQ, 
        vector<int>const& Qsizes, 
        ChronoDuration* tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int const nT = QQ.size();
    #pragma omp parallel for
    for(int i=0; i<nT; i++)
        sort(QQ[i].begin(), QQ[i].end());
    if(tim) *tim = chrono::steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng 
// ============================================================================
void BGPCOrdering::do_Random_Ordering_PrePartition(
        vector<vector<int>>& QQ, 
        vector<int>const& Qsizes, 
        ChronoDuration* tim){
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int const nT = QQ.size();
    #pragma omp parallel for
    for(int i=0; i<nT; i++){
        mt19937 mt(i);
        int const Qsize = Qsizes[i];
        int const Qsizem1 = Qsize-1;
        for(int j=0; j<Qsizem1; j++){
            uniform_int_distribution<int> dist(j, Qsizem1);
            swap(QQ[i][j], QQ[i][dist(mt)]);
        }
    }
    if(tim) *tim = chrono::steady_clock::now() - start;
}


// ============================================================================
// Author: Xin Cheng 
// ============================================================================
void BGPCOrdering::do_D1_LargestDegreeFirst_Ordering_PrePartition(
        vector<vector<int>>& QQ, 
        vector<int>const& Qsizes, 
        vector<int>const& srcPtr, 
        ChronoDuration* tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int const nT = QQ.size();
    #pragma omp parallel for
    for(int tid=0; tid<nT; tid++){
        int const Qsize = Qsizes[tid];
        vector<int>& Q = QQ[tid];
        do_D1_LargestDegreeFirst_Ordering_NoPartition(Q, Qsize, srcPtr);
    }
    if(tim) *tim = chrono::steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng 
// ============================================================================
void BGPCOrdering::do_D2_LargestDegreeFirst_Ordering_PrePartition(
        vector<vector<int>>& QQ, 
        vector<int>const& Qsizes, 
        vector<int>const& srcPtr, 
        vector<int>const& dstPtr, 
        vector<int>const& vtxVal, 
        ChronoDuration* tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int const nT = QQ.size();
    #pragma omp parallel for
    for(int tid=0; tid<nT; tid++) {
        int const Qsize = Qsizes[tid];
        vector<int>& Q = QQ[tid];
        do_D2_LargestDegreeFirst_Ordering_NoPartition(Q, Qsize, srcPtr, dstPtr, vtxVal);
    }
    if(tim) *tim = chrono::steady_clock::now() - start;
}


// ============================================================================
// Author: Xin Cheng 
// ============================================================================
void BGPCOrdering::do_D2Rough_LargestDegreeFirst_Ordering_PrePartition(
        vector<vector<int>>& QQ, 
        vector<int>const& Qsizes, 
        vector<int>const& srcPtr, 
        vector<int>const& dstPtr, 
        vector<int>const& vtxVal, 
        ChronoDuration* tim ) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int const nT = QQ.size();
    #pragma omp parallel for
    for(int tid =0; tid<nT; tid++){
        int const Qsize = Qsizes[tid];
        vector<int>& Q = QQ[tid];
        do_D2Rough_LargestDegreeFirst_Ordering_NoPartition(Q, Qsize, srcPtr, dstPtr, vtxVal);
    }
    if(tim) *tim = chrono::steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng 
// ============================================================================
void BGPCOrdering::do_BA_LargestDegreeFirst_Ordering_PrePartition(
        vector<vector<int>>& QQ, 
        vector<int>const& Qsizes, 
        vector<int>const& srcPtr, 
        vector<int>const& dstPtr, 
        vector<int>const& vtxVal, 
        ChronoDuration* tim ) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    int const nT = QQ.size();
    #pragma omp parallel for
    for(int tid=0; tid<nT; tid++) {
        int const Qsize = Qsizes[tid];
        vector<int>& Q = QQ[tid];
        do_BA_LargestDegreeFirst_Ordering_NoPartition(Q, Qsize, srcPtr, dstPtr, vtxVal);
    }
    if(tim) *tim = chrono::steady_clock::now() - start;
}










