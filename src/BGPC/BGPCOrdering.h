/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#ifndef BGPCORDERING_H
#define BGPCORDERING_H
#include <vector>
#include <iostream>
#include "BGPC.h"
#include <random>
#include <algorithm>
#include <unordered_set> // only used in do_D2_LargestDegreeFirst
using namespace std;

namespace ColPack {

// ============================================================================
// Partial Distance Two Shared Memory Parallel Greedy/Graph Coloring Ordering wrap
// ============================================================================
class BGPCOrdering : public BGPC {
public: // construction
    BGPCOrdering(){};
    BGPCOrdering(const string& file_name, const string& fmt, ChronoDuration* iotime, const int side, const int order, ChronoDuration* ordtime);
    virtual ~BGPCOrdering();
    
public: // deplete construction
    BGPCOrdering(BGPCOrdering&&)=delete;
    BGPCOrdering(const BGPCOrdering&)=delete;
    BGPCOrdering& operator=(BGPCOrdering&&)=delete;
    BGPCOrdering& operator=(const BGPCOrdering&)=delete;

public: // user api 
    void Ordering_PrePartition_OMP                        (const int side, vector<vector<int>>& QQ, vector<int>const& Qsizes, const int LOCAL_ORDER=ORDER_NONE, ChronoDuration* tim=nullptr);
    void Ordering_NoPartition                             (const int side, vector<int>& Q,          const int Qsize,     const int LOCAL_ORDER=ORDER_NONE, ChronoDuration* tim=nullptr);

public: // developer api
    void do_Natural_Ordering_NoPartition                      (vector<int>& Q, const int Qsize, ChronoDuration* tim=nullptr);
    void do_Random_Ordering_NoPartition                      (vector<int>& Q, const int Qsize, ChronoDuration* tim=nullptr);
    void do_D1_LargestDegreeFirst_Ordering_NoPartition       (vector<int>& Q, const int Qsize, vector<int>const& srcPtr, ChronoDuration* tim=nullptr);
    void do_D2_LargestDegreeFirst_Ordering_NoPartition       (vector<int>& Q, const int Qsize, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, ChronoDuration* tim=nullptr);
    void do_D2Rough_LargestDegreeFirst_Ordering_NoPartition  (vector<int>& Q, const int Qsize, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, ChronoDuration* tim=nullptr);
    void do_BA_LargestDegreeFirst_Ordering_NoPartition       (vector<int>& Q, const int Qsize, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, ChronoDuration* tim=nullptr);

    void do_Nature_Ordering_PrePartition                     (vector<vector<int>>&QQ, vector<int>const& Qsizes, ChronoDuration* tim=nullptr);
    void do_Random_Ordering_PrePartition                     (vector<vector<int>>&QQ, vector<int>const& Qsizes, ChronoDuration* tim=nullptr);
    void do_D1_LargestDegreeFirst_Ordering_PrePartition      (vector<vector<int>>&QQ, vector<int>const& Qsizes, vector<int>const& srcPtr, ChronoDuration* tim=nullptr);
    void do_D2_LargestDegreeFirst_Ordering_PrePartition      (vector<vector<int>>&QQ, vector<int>const& Qsizes, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, ChronoDuration* tim=nullptr);
    void do_D2Rough_LargestDegreeFirst_Ordering_PrePartition (vector<vector<int>>&QQ, vector<int>const& Qsizes, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, ChronoDuration* tim=nullptr);
    void do_BA_LargestDegreeFirst_Ordering_PrePartition      (vector<vector<int>>&QQ, vector<int>const& Qsizes, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, ChronoDuration* tim=nullptr);
   
    
//     void global_natural_ordering(const int side);
//     void global_random_ordering(const int side);
//     void global_largest_degree_first_ordering(const int side);
// 
//     void local_natural_ordering(vector<int>&vtxs);
//     void local_random_ordering(vector<int>&vtxs);
//     void local_largest_degree_first_ordering(const int side, vector<int>& vtxs); 
//     void local_largest_degree_first_ordering(const int side, vector<int>& vtxs, const int beg, const int end);  
// 
//    
      vector<int>&      get_ordered_queue_A()           { return m_ordered_queue_A; }
      vector<int>const& get_ordered_queue_A_const()const{ return m_ordered_queue_A; }
      int               get_ordered_queue_A_side() const{ return m_ordered_queue_A_side; }

protected: // members
    vector<int> m_ordered_queue_A;      //ordered vertex of chosen side
    //string      m_global_ordered_method;      //method
    int         m_ordered_queue_A_side;        //chosen side
};




}// endof namespace ColPack
#endif

