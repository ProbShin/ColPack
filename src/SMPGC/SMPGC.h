/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
******************************************************************************/


#ifndef SMPGCDEFINE_H
#define SMPGCDEFINE_H
#include <string>
#include <chrono>
#include <random>
// ============================================================================
// Author: Xin Cheng
// SMPGC: Shared Memory Parallel Graph Coloring
// ----------------------------------------------------------------------------
// **OVERVIEW**
//
// SMPGCCore:                        Graph data, IO
//     |->SMPGCOrdering:             Ordering
//          |-> SMPGCColoring:       D1 Coloring
//          |-> D2SMPGCColoring:     D2 Coloring
//
// ----------------------------------------------------------------------------
// ...
// ...
// ============================================================================

class SMPGC{
public:    typedef std::chrono::duration<double> ChronoDuration;         //display time using double number in the second unit   
//public: // in comman Computer is LP64 Model. change the following in case not. 
    //typedef unsigned long long int uint64;
    //typedef          long long int  int64;
    //typedef unsigned int           uint32;
    //typedef          int            int32;
    //#ifndef INT32
    //    typedef uint64 UINT;
    //    typedef int64  INT ;
    //#else
    //    typedef uint32 UINT;
    //    typedef int32  INT;
    //#endif

public:
    static const int RAND_SEED           = 5489u;

    static const int HASH_SEED           = 5489u;
    static const int HASH_SHIFT          = 0XC2A50F;
    static const int HASH_NUM_HASH       = 4;

    static const std::string FORMAT_MM   ;
    static const std::string FORMAT_BINARY;

    static const int ORDER_NONE          = 0;
    static const int ORDER_NATURAL       = 1;
    static const int ORDER_RANDOM        = 2;
    static const int ORDER_LARGEST_FIRST = 3;
    static const int ORDER_SMALLEST_LAST = 4;

    static const int ISI_WEIGHT_RAND     = 0;
    static const int ISI_WEIGHT_DEGREE   = 1;

    //static const int HYBRID_GM3P         = 1;
    //static const int HYBRID_GMMP         = 2;
    //static const int HYBRID_SERIAL       = 3;
    //static const int HYBRID_STREAM       = 4;


public:
    SMPGC(): m_mt(RAND_SEED) {};
    ~SMPGC(){};
public:
    SMPGC(SMPGC&&)=delete;
    SMPGC(const SMPGC&)=delete;
    SMPGC& operator=(SMPGC&&)=delete;
    SMPGC& operator=(const SMPGC&)=delete;

protected:
    std::mt19937     m_mt;  //random machine moved from Ordering
public:
    void set_rseed(const int x){ m_mt.seed(x); }
};



#endif


