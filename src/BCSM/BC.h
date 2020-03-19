/*******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#ifndef BC_HPP
#define BC_HPP

#include<vector>
#include<iostream>
#include<omp.h>
#include<map>
#include<list>
#include<set>
#include<random>
#include<chrono>


using namespace std;
namespace ColPack{

//=========================================================================
// class of BC
// ----------------------------------------------------------------
// 
//=========================================================================
class BC {
public: chrono::duration<double> ChronoDuration;
public:
    static const int L = 1;
    static const int R = 2;
    static const int ROW = 1;
    static const int COL = 2;
    static const int LR  = 3;
    static const int S = 1;
    static const int T = 2;

public:
    static const 
} // end of class BC


} // end of namespace ColPack
#endif
