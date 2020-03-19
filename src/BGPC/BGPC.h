/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#ifndef BGPC_H
#define BGPC_H
#include <vector>
#include <iostream>
#include <omp.h>
#include <map>  // used in ColPack::BipartiteGraphCore
#include <list> // used in ColPack::BipartiteGraphInputOutput
#include <set>
#include "BipartiteGraphCore.h"
#include "BipartiteGraphInputOutput.h" //#include "ColPackHeaders.h" //#include "GraphOrdering.h"
#include <random>
#include <chrono>


using namespace std;

namespace ColPack {

//=============================================================================
// Bipartite Graph Partial (Distance Two) Coloring 
// ----------------------------------------------------------------------------
// first author: xin cheng
// ----------------------------------------------------------------------------
// Planned service for Source data structure and file io. Since have been 
// implemented in the ColPack::BipartiteGrpah* class. So, just inherit from them
// And then adds some marco defines.
//=============================================================================

class BGPC : public BipartiteGraphInputOutput {//public BipartiteGraphInputOut{
public:
    typedef chrono::duration<double> ChronoDuration;         //display time using unit second    
    //typedef chrono::duration<double,micro> ChronoDuration; //display time using unit microsecond
public:
    static const int L = 1;                        // Bipartite Graph left side,  row
    static const int R = 2;                        // Bipartite Graph right side, column
    static const int ROW = 1;                      // Bipartite Graph left side,  row
    static const int COL = 2;                      // Bipartite Graph right side, column
    static const int LR  = 3;                      // both row and column

public:
    static const int  ORDER_NONE    = 0;           // None Ordering specific
    static const int  ORDER_NATURAL = 1;           // Nature Ordering , known as FF
    static const int  ORDER_RANDOM  = 2;           // Random Ordering 
    static const int  ORDER_LARGEST_FIRST = 3;     // Largest First Ordering
    static const int  ORDER_D1_LARGEST_FIRST = 3;     // Largest First Ordering
    static const int  ORDER_SMALLEST_LAST = 4;     // Smallest Last Ordering

    static const int  ORDER_D2_LARGEST_FIRST = 5;     // Largest First Ordering
    static const int  ORDER_D2ROUGH_LARGEST_FIRST = 6;     // Largest First Ordering
    static const int  ORDER_BA_LARGEST_FIRST = 7;     // Largest First Ordering
    
public:
    static const string FORMAT_MM;                 // ="MM";
    static const string FORMAT_POTHEN;             // ="POTHEN";   sqrt(G)
    static const string ORDER_STR_NATURAL;         // = "NATURAL";
    static const string ORDER_STR_RANDOM ;         // = "RANDOM";
    static const string ORDER_STR_LARGEST_FIRST;   // = "LARGEST_FIRST";
    static const string ORDER_STR_SMALLEST_LAST;   // = "SMALLEST_LAST";

public: 
    static const int WEIGHT_RANDOM = 0;
    static const int WEIGHT_DEGREE_D1 = 1;
    static const int WEIGHT_DEGREE_D2 = 2;


public:
    static const int BIT_FORBIDDENARRAY_WIDE = 32;   // limited memory implementation, 32 bit for 32bit machine
    typedef unsigned int BIT_FORBIDDENARRAY_TYPE;

    //static const int BIT_FORBIDDENARRAY_WIDE = 64; // limited memory implementation, 64 bit for 64bit machine
    //typedef unsigned long long int BIT_FORBIDDENARRAY_TYPE;

public:
    BGPC():m_mt(5489u) {}
    virtual ~BGPC(){};
    
    // Deplete constructions
    BGPC(BGPC&&)=delete;
    BGPC(const BGPC&)=delete;
    BGPC& operator=(BGPC&&)=delete;
    BGPC& operator=(const BGPC&)=delete;

public:
    void set_rseed(const int x){ m_mt.seed(x); }

protected:
    mt19937 m_mt;  // random seed
};


}// endof namespace ColPack
#endif

