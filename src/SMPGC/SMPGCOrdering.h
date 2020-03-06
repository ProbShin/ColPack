/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#ifndef SMPGCORDERING_H
#define SMPGCORDERING_H
#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <omp.h>

#include "ColPackHeaders.h" //#include "GraphOrdering.h"
#include "SMPGCGraph.h"
using namespace std;

namespace ColPack {

// ============================================================================
// Author: Xin Cheng
// Shared Memory Parallel Greedy/Graph Coloring Ordering wrap
// 
// includes 
//      [sequential, parallel], 
//      [whole vertex, subset vertex, prepartitioned vertex],
//      [natural, random, largest degree first, smallest degree last]
//
// ============================================================================
class SMPGCOrdering : public SMPGCGraph {
public: // construction
    SMPGCOrdering(const string& file_name, const string& fmt, ChronoDuration*iotime, const int order, ChronoDuration* ordtime=nullptr);
    virtual ~SMPGCOrdering();

public: // deplete construction
    SMPGCOrdering(SMPGCOrdering&&)=delete;
    SMPGCOrdering(const SMPGCOrdering&)=delete;
    SMPGCOrdering& operator=(SMPGCOrdering&&)=delete;
    SMPGCOrdering& operator=(const SMPGCOrdering&)=delete;

public: // API: interface
    void ordering(
            const string&s1="SEQ", 
            const string& s2="WHOLE", 
            const int order=ORDER_NATURAL, 
            ChronoDuration* t=nullptr, 
            void* parg=nullptr
            );
    const vector<int>& get_ordered_vertex() const { return m_ordered_vertex; }
    const int          get_ordered_method() const { return m_ordered_method; }

    void   Translate_OrderId_To_OrderTag(const int LOCAL_ORDER, string&lotag);
    string Translate_OrderId_To_OrderTag(const int LOCAL_ORDER);
protected:
    // sequential ordering on the whole graphs
    void seq_whole_ntr_ordering(ChronoDuration* ordtime=nullptr);
    void seq_whole_rnd_ordering(ChronoDuration* ordtime=nullptr);
    void seq_whole_ldf_ordering(ChronoDuration* ordtime=nullptr);
    void seq_whole_sdl_ordering(ChronoDuration* ordtime=nullptr);

    // sequential ordering on the subset of the graphs 
    void seq_subset_ntr_ordering(vector<int>& Q,ChronoDuration* ordtime=nullptr);
    void seq_subset_rnd_ordering(vector<int>& Q,ChronoDuration* ordtime=nullptr);
    void seq_subset_ldf_ordering(vector<int>& Q,ChronoDuration* ordtime=nullptr);
    void seq_subset_sdl_ordering(vector<int>& Q,ChronoDuration* ordtime=nullptr);

    // parallel ordering on the whole graphs
    void omp_whole_ntr_ordering(ChronoDuration* ordtime=nullptr);
    void omp_whole_rnd_ordering(ChronoDuration* ordtime=nullptr);
    void omp_whole_ldf_ordering(ChronoDuration* ordtime=nullptr);
    void omp_whole_sdl_ordering(ChronoDuration* ordtime=nullptr);
    
    // parallel ordering on the prepartitioned graphs
    void omp_prepartition_ntr_ordering(vector<vector<int>>&QQ,ChronoDuration* ordtime=nullptr);
    void omp_prepartition_rnd_ordering(vector<vector<int>>&QQ,ChronoDuration* ordtime=nullptr);
    void omp_prepartition_ldf_ordering(vector<vector<int>>&QQ,ChronoDuration* ordtime=nullptr);
    void omp_prepartition_sdl_ordering(vector<vector<int>>&QQ,ChronoDuration* ordtime=nullptr);

    //void DynamicLargestDegreeFirstOrdering(vector<INT>& vtxs, INT N);
    //void IncidenceDegreeOrdering(vector<INT>& vtxs, INT N);
    //void LogOrdering(vector<INT>& vtxs, INT N);

protected: // members
    vector<int> m_ordered_vertex;  // used to store the whole ordered vertex   
    int      m_ordered_method;  // tags of the whole ordered vertex
};


}// endof namespace ColPack
#endif

