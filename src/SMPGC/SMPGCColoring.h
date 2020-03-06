/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#ifndef SMPGCColoring_H
#define SMPGCColoring_H
#include <vector>
#include <omp.h>
#include "ColPackHeaders.h" //#include "GraphOrdering.h"
#include "SMPGCOrdering.h"
#include <chrono>
using namespace std;
using namespace chrono;
namespace ColPack {

//=============================================================================
// Author: Xin Cheng
// Shared Memeory Parallel (Greedy)/Graph Coloring -> SMPGC
// ============================================================================
class SMPGCColoring : public SMPGCOrdering {
public: // Constructions
    SMPGCColoring(const string& graph_name, const string& fmt, ChronoDuration*iotime=nullptr, const int ord=ORDER_NATURAL, ChronoDuration*ordtime=nullptr);
    virtual ~SMPGCColoring(){}

        // Deplete constructions
        SMPGCColoring(SMPGCColoring&&)=delete;
        SMPGCColoring(const SMPGCColoring&)=delete;
        SMPGCColoring& operator=(SMPGCColoring&&)=delete;
        SMPGCColoring& operator=(const SMPGCColoring&)=delete;

public: // API
    int Coloring(int nT, int &color, vector<int>& vtxColors, const string& method, const int nVerbose);
    
    // serial algorithm
    int D1_serial(int&colors, vector<int>&vtxColors, const int nVerbose);

    // D1GC OMP speculative iterative approaches
    int D1GC_OMP_SpeculativeIterative_3P(
            int nT,                               // number of threads
            int& colors,                          // number of colors
            vector<int>&vtxColors,                // mapping vertx to its color
            const int LOCAL_ORDER=ORDER_NATURAL,  // parallel local order
            const bool bMemOpt=false,             // using Memory Optimized Buffer? Which trade runtime for memory
            const int nVerbose=1                  // verbose level
            );

    int D1GC_OMP_SpeculativeIterative_MP(
            int nT,                               // number of threads
            int& colors,                          // number of colors
            vector<int>&vtxColors,                // mapping vertex to its color
            const int LOCAL_ORDER=ORDER_NATURAL,  // parall local order
            const bool bMemOpt=false,             // using Memory Optimized Buffer? Which trade runtime for memory
            const int nVerbose=1                  // verbose level
            );
    
    
    void do_D1GC_OMP_TentativeColoring(
        const int nT,                   // number of threads
        vector<int>& vtxColors,         // mapping vertex to its color
        const vector<int>& vtxPtr,      // CSR format edges
        const vector<int>& vtxVal,      // CSR format values
        const vector<vector<int>> &QQ,  // vertex to color
        const vector<int>& Qsizes,      // local vertex size
        vector<vector<int>>&Fs,         // Forbidden color buffers
        const int BufSize,              // Buffer size
        ChronoDuration &time            // run time
        );

    void do_D1GC_OMP_CheckConflicts(
        const int nT,                // number of threads
        vector<int> &vtxColors,      // mapping vertex to its color
        const vector<int>& vtxPtr,   // CSR format edges
        const vector<int>& vtxVal,   // CSR format values
        vector<vector<int>>& QQ,     // vertex to color
        vector<int>& Qsizes,         // local vertex size
        ChronoDuration & time,       // run time
        int& num_Conflicts           // number of conlflicts
        );

    void do_D1GC_OMP_TentativeColoring_MemOpt(
        const int nT,                   // number of threads
        vector<int>& vtxColors,         // mapping vertex to its color
        const vector<int>& vtxPtr,      // CSR format edges
        const vector<int>& vtxVal,      // CSR format values
        const vector<vector<int>> &QQ,  // vertex to color
        const vector<int>& Qsizes,      // local vertex size
        vector<vector<int>>&Fs,         // Forbidden color buffers
        const int BufSize,              // Buffer size
        ChronoDuration &time            // run time
        );           
   
    // D1GC original (without local order, prepartition and using atomic operation)
    int D1GC_OMP_SpeculativeIterative_3P_Original(
        int nT,                         // number of threads
        int&colors,                     // number of colors
        vector<int>&vtxColors,          // mapping vertex to its color
        const int nVerbose              // verbose level
        );

    int D1GC_OMP_SpeculativeIterative_MP_Original(
        int nT,                         // number of threads
        int&colors,                     // number of colors
        vector<int>&vtxColors,          // mapping vertex to its color
        const int nVerbose              // verbose level
        );

    void do_TentativeColoring_NoPrePartition(
        const int nT,                   // number of threads
        vector<int>& vtxColors,         // mapping vertex to its color
        const vector<int>& vtxPtr,      // CSR format
        const vector<int>& vtxVal,      // CSR format
        const vector<int>& Q,           // vertex to color
        const int Qsize,                // number of vertex to color
        vector<vector<int>>& Fs,        // Forbidden color buffers
        const int BufSize,              // Buffer size
        ChronoDuration& time            // run time
        );
    
    void do_CheckConflicts_NoPrePartition(
        const int nT,                   // number of threads
        vector<int>& vtxColors,         // mapping vertex to its color
        const vector<int>& vtxPtr,      // CSR format
        const vector<int>& vtxVal,      // CSR format
        const vector<int> &Q,           // vertex to color
        const int Qsize,                // number of vertex to color
        vector<int>& cfQ,               // conflict vertex queue
        int& cfQsize,                   // number of conflict vertex
        ChronoDuration &time           // run time
        );



    // D1GC OMP independent-set iterative approach
    int D1GC_OMP_IndependentSet_Iterative_LB(
        int nT,                         // number of threads
        int&colors,                     // number of colors
        vector<int>&vtxColors,          // mapping vertex to its color
        const int LOCAL_ORDER,          // parallel local ordering
        const int nVerbose              // verbose level
        );
    
    int D1GC_OMP_IndependentSet_Iterative_JP(
        int nT,                         // number of threads
        int&colors,                     // number of colors
        vector<int>&vtxColors,          // mapping vertex to its color
        const int LOCAL_ORDER,          // parallel local ordering
        const int nVerbose              // verbose level
        );


    int D1GC_OMP_IndependentSet_Iterative_JP_adv(
        int nT,                         // number of threads
        int&colors,                     // number of colors
        vector<int>&vtxColors,          // mapping vertex to its color
        const int LOCAL_ORDER,          // parallel local ordering
        const int nVerbose              // verbose level
        );



    // D1GC hybrid algorithm
    int D1GC_OMP_Hybrid_ISI_K_SpI(
            int nT,                  // number of thread
            int&colors,              // number of colors
            vector<int>&vtxColors,   // mapping vertex to its color
            const int K,             // number of ISI iterations before switch to SpI
            const int nVerbose       // verbose level
            );

// D1GC hybrid algorithm
    int D1GC_OMP_Hybrid_ISI_Katmost2_SpI(
            int nT,                  // number of thread
            int&colors,              // number of colors
            vector<int>&vtxColors,   // mapping vertex to its color
            const int K,             // number of ISI iterations before switch to SpI
            const int nVerbose       // verbose level
            );





    // D2GC serial 
    int D2GC_Serial(int &color, vector<int>&vtxColors, const int nVerbose=0);

    // D2GC OMP speculative iterative approaches
    int D2GC_OMP_SpeculativeIterative(
            int nT,                               // number of threads
            int& color,                           // number of colors
            vector<int>& vtxColors,               // mapping vertex to its color
            const int LOCAL_ORDER=ORDER_NONE,     // parall local order
            const int NUM_RND_TC=1,       // number of iterations using Random TC
            const int NUM_NET_TC=0,       // number of iteartions using NetBin TC
            const int NUM_NET_CHECK=2,    // number of iteartions using NetBin CC
            const bool b_SEQ_HC=false,
            const int nVerbose=0          // verbose level
            ); 
    
 
    // D2GC omp speculative iterative approach without prepartition and local ordering but using atomic
    int D2GC_OMP_SpeculativeIterative_NoPrePartition(
            int nT,                               // number of threads
            int& color,                           // number of colors
            vector<int>& vtxColors,               // mapping vertex to its color
            int const LOCAL_ORDER=ORDER_NONE,     // parall local order
            int const NUM_RND_TC=1,               // number of random TC phase iterations
            int const NUM_NET_TC=2,               // number of NetBin CC phase iterations
            int const NUM_NET_CC=2,
            bool const b_SEQ_HC=false,
            const int nVerbose=0                  // verbose level
            ); 

    // D2GC omp speculative iterative approach without prepartition and local ordering but using atomic
    int D2GC_OMP_SpeculativeIterative_Memopt(
            int nT,                               // number of threads
            int& color,                           // number of colors
            vector<int>& vtxColors,               // mapping vertex to its color
            int const LOCAL_ORDER=ORDER_NONE,     // parall local order
            int const NUM_RND_TC=1,               // number of random TC phase iterations
            int const NUM_NET_TC=2,               // number of NetBin CC phase iterations
            int const NUM_NET_CC=2,
            bool const b_SEQ_HC=false,
            const int nVerbose=0                  // verbose level
            ); 



    int D2GC_OMP_IndependentSetIterative(
            int nT,                               // number of threads
            int& color,                           // number of colors
            vector<int>& vtxColors,               // mapping vertex to its color
            const int LOCAL_ORDER=ORDER_NONE,     // parall local order
            const int nWeightType=ISI_WEIGHT_DEGREE, 
            const int nVerbose=0                  // verbose level

            );

    int D2GC_OMP_Hybrid_ISI_K_SpI(
            int nT,                     // number of the threads
            int& colors,                // number of the colors
            vector<int>& vtxColors,     // mapping vertex to its color
            int const LOCAL_ORDER,      // parallel local order
            int const nWeightType,      // random weight or degree weight
            int const K,                // number of ISI iterations
            int const NUM_RND_TC,
            int const NUM_NET_TC,
            int const NUM_NET_CC,
            bool const b_SEQ_HC,
            int const nVerbose          // verbose level
            );


    // D2GC utility
    void do_D2GC_Seq_GreedyColoring(
            vector<int> const & vtxPtr,          // csr format
            vector<int> const & vtxVal,          // csr format
            vector<int>&vtxColors,              // mapping vertex to its color
            const vector<int>&Q,                // vertex to color
            const int Qsize,                    // number of the vertex to color
            vector<int>&F,                      // forbidden color array
            const int BufSize,                  // forbidden color array size
            ChronoDuration* pTime=nullptr       // run time
            );

    void do_D2GC_Seq_GreedyColoring_NoInitF(
            const vector<int>& vtxPtr,          // csr format
            const vector<int>& vtxVal,          // csr format
            vector<int>&vtxColors,              // mapping vertex to its color
            const vector<int>&Q,                // vertex to color
            const int Qsize,                    // number of the vertex to color
            vector<int>&F,                      // forbidden color array
            const int BufSize,                  // forbidden color array size
            ChronoDuration* pTime=nullptr       // run time
            );


    void do_D2GC_OMP_PhaseTC(
        vector<int>&vtxColors,         // mapping vertex to its color
        const vector<vector<int>>& QQ, // vertex to color for each thread
        const vector<int>& Qsizes,     // number of the vertex to color for each thread
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        vector<vector<int>>& Fs,       // forbidden color array
        const int BufSize,             // forbidden color array size
        ChronoDuration *ptime       // runtime
        ) ;

    void do_D2GC_OMP_PhaseNetBinTC(
        vector<int>& vtxColors,             
        vector<int> const & Q,
        int const  Qsize,
        const vector<int>&vtxPtr,
        const vector<int>&vtxVal,
        vector<vector<int>>& Fs,
        const int BufSize,
        vector<vector<int>>& Ws, //work load vertices
        ChronoDuration *pTime   
        );

    void do_D2GC_OMP_PhaseRandomTC(
        vector<int> &vtxColors,          // mapping vertex to its color
        const vector<vector<int>>&QQ,    // vertex to color of each thread
        const vector<int>& Qsizes,       // number of the vertex of each thread
        vector<mt19937>& mts,            // random number generator of each thread
        const int ColorLB0Base,          // random number range
        ChronoDuration *pTime            // run time
        );

    void do_D2GC_OMP_PhaseTC_NoPrePartition(
        vector<int> &vtxColors,
        vector<int> & Q,
        int & Qsize,
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        vector<vector<int>> &ForbiddenArrays,
        const int BufSize,
        ChronoDuration *pTime
        );
    
    void do_D2GC_OMP_PhaseNetBinTC_NoPrePartition(
        vector<int> &vtxColors,
        vector<int> const & const_ordered_Q,
        int const N,
        vector<int> const & vtxPtr,
        vector<int> const & vtxVal,
        vector<vector<int>> &Fs,
        int const BufSize,
        vector<vector<int>> &W,
        ChronoDuration *pTime
        );

    void do_D2GC_OMP_PhaseRandomTC_NoPrePartition(
        vector<int> &vtxColors,
        vector<int> const & Q,
        int const QSize, 
        vector<mt19937>& mts,
        int const ColorLB0Base,
        ChronoDuration *pTime
        );

    void do_D2GC_OMP_PhaseCC(
            vector<vector<int>>& QQ,
            vector<int>& Qsizes,
            vector<int>& vtxColors,
            const vector<int>& vtxPtr,
            const vector<int>& vtxVal,
            int *pNum,
            ChronoDuration *pTime
            ); 

    void do_D2GC_OMP_PhaseNetBinCC(
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
            );

    void do_D2GC_OMP_PhaseCC_NoPrePartition(
        vector<int>& Q,
        int& Qsize,
        vector<int>& vtxColors,
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        vector<int>& cfQ,
        int * pNum,
        ChronoDuration* pTime
        );

    void do_D2GC_OMP_PhaseNetBinCC_NoPrePartition(
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
        );

    void do_D2GC_Seq_findIndSet(
        vector<int>& IndSet, 
        int &IndSetSize,
        vector<int> const & vtxPtr,
        vector<int> const & vtxVal,
        vector<int> const & vtxColors,
        vector<int> const & Weight,
        vector<int> & Q,
        int &Qsize
        );


    // D2GC MemOpt
int D2GC_OMP_SpeculativeIterative_MemOpt(
        int nT,                     // number of threads
        int &colors,                // number of colors
        vector<int>&vtxColors,      // mapping vertex to its color
        int const LOCAL_ORDER,      // parallel local ordering
        int const NUM_RND_TC,       // number of iterations using Random TC
        int const NUM_NET_TC,       // number of iteartions using NetBin TC
        int const NUM_NET_CC,       // number of iteartions using NetBin CC
        bool const b_SEQ_HC,        // whether using sequential coloring to color the remaining graphs.
        int const nVerbose          // verbose level
        );

void do_D2GC_OMP_PhaseTC_MemOpt(
        vector<int>&vtxColors,         // mapping vertex to its color
        const vector<vector<int>>& QQ, // vertex to color for each thread
        const vector<int>& Qsizes,     // number of the vertex to color for each thread
        const vector<int>& vtxPtr,
        const vector<int>& vtxVal,
        int const Capacity,             // forbidden color array size
        ChronoDuration *ptime       // runtime
        ) ;

void do_D2GC_OMP_PhaseNetBinTC_MemOpt(
        vector<int>& vtxColors,             
        vector<int> const & const_ordered_Q,
        int const N,
        const vector<int>&vtxPtr,
        const vector<int>&vtxVal,
        vector<vector<int>>& Ws, //work load vertices
        int const WCapacity,     //
        ChronoDuration *pTime   
        );

void do_D2GC_PopAllW_by_Coloring_MemOpt(
        int const v,
        vector<int>& W, 
        int & Wsize,
        vector<int>& vtxColors,
        int const WCapacity, 
        vector<int> const& vtxPtr,
        vector<int> const& vtxVal
        );

unsigned int do_D2GC_Build_MemOpt_ForbiddenArray(
        int const v,  
        int const Low, 
        vector<int> const& vtxColors,
        int const WCapacity,
        vector<int> const& vtxPtr,
        vector<int> const& vtxVal
        ) ;

void do_D2GC_OMP_PhaseNetBinCC_MemOpt(
        vector<vector<int>> & QQ,
        vector<int>& Qsizes,
        vector<int>& vtxColors,
        const vector<int> &vtxPtr,
        const vector<int> &vtxVal,
        vector<int> const &const_ordered_Q,
        int const N,
        int const Capacity,
        int * pNum, 
        ChronoDuration* pTime
        );








public: // Utilites
    int cnt_d1conflict(const vector<int>& vc, bool bVerbose=false);
    int cnt_d2conflict(const vector<int>& vc, bool bVerbose=false);
    int calc_num_colors_from_vtx_colors(const vector<int> &vtxColors);
private:

unsigned int mhash(unsigned int a, unsigned int seed){
    a ^= seed;
    a = (a + 0x7ed55d16) + (a << 12);
    a = (a ^ 0xc761c23c) + (a >> 19);
    a = (a + 0x165667b1) + (a << 5);
    a = (a ^ 0xd3a2646c) + (a << 9);
    a = (a + 0xfd7046c5) + (a << 3);
    a = (a ^ 0xb55a4f09) + (a >> 16);
    return a;
}




protected:
    //int         m_total_num_colors;
    //vector<int> m_vertex_color;
    //string      m_method;

}; // end of class SMPGCColoring


}// endof namespace ColPack
#endif

