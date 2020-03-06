/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#include "SMPGCColoring.h"
#include <chrono> //c++11 system time
#include <random> //c++11 random
using namespace std;
using namespace ColPack;




// ============================================================================
// Author: Xin Cheng
// Tentative Coloring modular with Memory Optimized
// using a unsigned integer as the forbidden color buffer.
// the buffer will using bit set/clear operation.
//
// Normally the forbidden color buffer should be same width with the machine
// i.e. 32 bit for 32bit machine, 64 bit for 64 bit machine
//
// NOTICE:
// since there is no standard way (cross plantform) to tell the machine is 32/64bit
// also no standard way to tell the compiler is configuered 32/64bit model.
//
//     -  One way is to using PreDefined Macro Defination. And using #if checks.
//     -  The other way is to hard code 32 or 64.
//
// To minimize the dependent on the compiler, (automake, cmake, etc), I am using the
// method, hard code for this function. The related variables are 
//  "BUFFWIDTH" and "F"
// ============================================================================
void SMPGCColoring::do_D1GC_OMP_TentativeColoring_MemOpt(
        const int nT,                   // number of thread
        vector<int>& vtxColors,         // mapping from vertex to its color
        const vector<int>& vtxPtr,      // CSR format edges
        const vector<int>& vtxVal,      // CSR format values
        const vector<vector<int>> &QQ,  // vertex to color
        const vector<int>& Qsizes,      // local vertex size
        vector<vector<int>>&Fs,         // Forbidden color buffers
        const int BufSize,              // Buffer size
        ChronoDuration &time            // run time
        ){           
    
    auto start = steady_clock::now();
    #pragma omp parallel
    {
        //const int BUFFWIDTH = 64;
        //unsigned long long int F = ~0;
        const int BUFFWIDTH = 32;
        unsigned int F = ~0;
        
        const int tid = omp_get_thread_num();
        const vector<int>& Q = QQ[tid];
        const int Qsize = Qsizes[tid];

        for(int iu=0; iu<Qsize; iu++) {
            const int u = Q[iu];
            int offset_mask=0;
            while(true){
                F = ~0;
                const int LOW = (offset_mask++)*BUFFWIDTH;
                for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                    const int vc_local=vtxColors[vtxVal[iv]] - LOW;  //dis-regards the overflow risk.
                    if(vc_local>=0 && vc_local<BUFFWIDTH) {
                        F &= ~(1<<(vc_local));  // clear the bit
                    }
                }

                // find the first settled bit, if there is any
                if(F!=0){
                    for(int i=0; i<BUFFWIDTH; i++) {
                        if(F&(1<<i)){
                            vtxColors[u]=LOW+i;
                            break;
                        }
                    }
                    break; // break while loop
                }
            }// end while(true) 
        }// end for v
    } //end omp parallel
    auto end = steady_clock::now();
    time = end - start;
}



