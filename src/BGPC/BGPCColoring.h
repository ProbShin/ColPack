/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#ifndef BGPCCOLORING_H
#define BGPCCOLORING_H
#include <vector>
#include <omp.h>
#include <chrono>
#include <sstream>    // stringstream
#include "BGPCOrdering.h" //#include "GraphOrdering.h"

using namespace std;

namespace ColPack {

// ============================================================================
// Shared Memeory Parallel Greedy/Graph Coloring
// ============================================================================
class BGPCColoring : public BGPCOrdering {
public: // Constructions
    BGPCColoring(){};
    BGPCColoring(const string& graph_name);
    BGPCColoring(const string& graph_name, const string& fmt=FORMAT_MM, ChronoDuration*iotime=nullptr, const int side=L, const int ord=ORDER_RANDOM, ChronoDuration*ordtime=nullptr);
    virtual ~BGPCColoring(){}

        // Deplete constructions
        BGPCColoring(BGPCColoring&&)=delete;
        BGPCColoring(const BGPCColoring&)=delete;
        BGPCColoring& operator=(BGPCColoring&&)=delete;
        BGPCColoring& operator=(const BGPCColoring&)=delete;

public: // API
    int Coloring(const int side, int nT, const string& method, int& nColor, vector<int>&vColors);
    
    //int get_num_colors() const { return m_total_num_colors; } 
    //const vector<int>& get_vertex_colors() const { return m_vertex_color; }
    //void         get_vertex_colors(vector<int>& x) { x.assign(m_vertex_color.begin(), m_vertex_color.end()); }
    
    // Algorithms 
    int BGPC_SEQ_Greedy(int const side, int&colors, vector<int>&vtxColors, const int bVerbose=0);
    
    int BGPC_OMP_SpeculativeIterative   (int const side, int nT, int& colors, vector<int>& vtxColors, int const LOCAL_ORDER=ORDER_NONE, 
            int const NUM_RND_TC=0, int const NUM_NET_TC=0, int const NUM_NET_CC=0, bool const b_SEQ_HC=false, int const nVerbose=0);
    int BGPC_OMP_IndependentSetIterative(int const side, int nT, int& colors, vector<int>& vtxColors, int const LOCAL_ORDER=ORDER_NONE, 
            int const  WEIGHT_TYPE=WEIGHT_DEGREE_D2, int const nVerbose=0);
    int BGPC_OMP_Hybrid_ISI_SPI         (int const side, int nT, int& colors, vector<int>& vtxColors, int const LOCAL_ORDER=ORDER_NONE, 
            int const  WEIGHT_TYPE=WEIGHT_DEGREE_D2, 
            int const K=1,
            int const NUM_RND_TC=0, int const NUM_NET_TC=0, int const NUM_NET_CC=0, bool const b_SEQ_HC=false, int const nVerbose=0);

    int BGPC_OMP_Hybrid_ISI_SPI_Katmost2   (int const side, int nT, int& colors, vector<int>& vtxColors, int const LOCAL_ORDER=ORDER_NONE, 
            int const  WEIGHT_TYPE=WEIGHT_DEGREE_D2, 
            int const K=1,
            int const NUM_RND_TC=0, int const NUM_NET_TC=0, int const NUM_NET_CC=0, bool const b_SEQ_HC=false, int const nVerbose=0);

    int BGPC_OMP_SpeculativeIterative_Atomic();
    int BGPC_OMP_SpeculativeIterative_Bit(int const side, int nT, int& colors, vector<int>& vtxColors, int const LOCAL_ORDER=ORDER_NONE, int const NUM_NET_TC=0, int const NUM_NET_CC=0, bool const b_SEQ_HC=false, int const nVerbose=0 );
   


    int BGPC_OMP_SpeculativeIterative_withoutPrePartition(
        int const side,             // L or R
        int nT,                     // number of threads
        int& colors,                // number of colors
        vector<int>& vtxColors,     // mapping vertex to color
        int const LOCAL_ORDER,      // Local Order
        bool const b_SEQ_HC,        // using Seq-Coloring to solve the conflicts?
        int const nVerbose          // verbose level
        );


    // implementations
    void do_BGPC_Seq_Greedy(vector<int>& vtxColors, vector<int> const&Q, int const Qsize, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, vector<int>&F, int const BufSize, ChronoDuration*tim=nullptr);
    void do_BGPC_Seq_Greedy_NoInitF(vector<int>& vtxColors, vector<int>const& Q, int const Qsize, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, vector<int>&F, int const BufSize, ChronoDuration*tim=nullptr);
   
    void do_BGPC_Phase_PlainTC(vector<int>&vtxColors, vector<vector<int>>const& QQ, vector<int>const& Qsizes, vector<int>const& srcPtr, vector<int>const& dstPtr,  vector<int>const& vtxVal, vector<vector<int>>& ForbiddenArrays, int const BufSize, ChronoDuration*tim=nullptr);
    void do_BGPC_Phase_PlainCC(vector<int>&vtxColors, vector<vector<int>>& QQ, vector<int>& Qsizes, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, int* num=nullptr, ChronoDuration* tim=nullptr);
    void do_BGPC_Phase_RandomTC(vector<int>& vtxColors, vector<vector<int>>const& QQ, vector<int>const& Qsizes, vector<int>const& srcPtr, vector<mt19937>& mts, ChronoDuration* tim=nullptr);
    void do_BGPC_Phase_NetbinTC(vector<int>& vtxColors, vector<int>const& ordered_queue, int const Ndst, vector<int>const& dstPtr, vector<int>const& vtxVal, vector<vector<int>>&Ws, vector<vector<int>>&Fs, int const BufSize, ChronoDuration* tim=nullptr);
    void do_BGPC_Phase_NetbinCC(vector<int>& vtxColors, vector<int>const& ordered_queue, int const Ndst, vector<int>const& dstPtr, vector<int>const& vtxVal, vector<vector<int>>&Fs, int const BufSize, int *num, ChronoDuration* tim);

    void do_BGPC_Seq_FindIndSet(vector<int>& IndSet, int& IndSetSize, vector<int>& Q, int& Qsize, vector<int>const& vtxColors, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, vector<uint64_t>const& Weights);
    
    void do_BGPC_Phase_PlainTC_atomic();
    void do_BGPC_Phase_PlainCC_atomic();
    void do_BGPC_Phase_RandomTC_atomic();
    void do_BGPC_Phase_NetbinTC_atomic();
    void do_BGPC_Phase_NetbinCC_atomic();

    void do_BGPC_Phase_PlainTC_Bit(
        vector<int>& vtxColors,
        vector<vector<int>>const& QQ,
        vector<int>const& Qsizes,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        ChronoDuration *tim=nullptr
        );
    
    //void do_BGPC_Phase_RndTC_BIT();
    void do_BGPC_Phase_NetBinTC_Bit( 
        vector<int>& vtxColors, 
        vector<int>const& vAs, int const Asize,
        vector<int>const& vBs, int const Bsizes, vector<int>const& srcPtr,
        vector<int>const& dstPtr, vector<int>const& vtxVal, vector<vector<int>>& Ws, vector<mt19937>& mts, ChronoDuration *tim=nullptr);


    void do_BGPC_Phase_NetBinCC_Bit(
        vector<int>& vtxColors, vector<int>const vBs, int const& Bsizes, vector<int>const& srcPtr,
        vector<int>const& dstPtr, vector<int>const& vtxVal, ChronoDuration *tim=nullptr );


    void do_BGPC_Seq_Greedy_Bit(
        vector<int>& vtxColors,
        vector<int>const& Q,
        int const Qsize,
        vector<int>const &srcPtr,
        vector<int>const &dstPtr,
        vector<int>const &vtxVal,
        ChronoDuration* tim=nullptr);

    void do_BGPC_Seq_Greedy_Bit_with_set(
        vector<int>& vtxColors,
        vector<int>const& Q,
        int const Qsize,
        vector<int>const &srcPtr,
        vector<int>const &dstPtr,
        vector<int>const &vtxVal,
        ChronoDuration* tim=nullptr);

    void do_BGPC_Seq_Greedy_Bit_with_bitset(
        vector<int>& vtxColors,
        vector<int>const& Q,
        int const Qsize,
        vector<int>const &srcPtr,
        vector<int>const &dstPtr,
        vector<int>const &vtxVal,
        ChronoDuration* tim=nullptr);


    void do_BGPC_Phase_PlainTC_NoPTT(
        vector<int>& vtxColors, 
        vector<int>const& Q,
        int const Qsize, 
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal, 
        vector<vector<int>>& Fs, 
        int const BufSize, 
        ChronoDuration* tim );

    void do_BGPC_Phase_PlainCC_NoPTT(
        vector<int>& vtxColors, 
        vector<int>& Q, 
        int &Qsize, 
        vector<int>const& srcPtr, 
        vector<int>const& dstPtr, 
        vector<int>const& vtxVal, 
        vector<int>& Q2,
        ChronoDuration* tim
        );
   


    // try new idea of cc by deg not by id

    int BGPC_OMP_SpeculativeIterative_ccbydeg   (int const side, int nT, int& colors, vector<int>& vtxColors, int const LOCAL_ORDER=ORDER_NONE, 
            int const NUM_RND_TC=0, int const NUM_NET_TC=0, int const NUM_NET_CC=0, bool const b_SEQ_HC=false, int const nVerbose=0);

    void do_BGPC_Phase_PlainCC_ccbydeg(vector<int>&vtxColors, vector<vector<int>>& QQ, vector<int>& Qsizes, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, int* num=nullptr, ChronoDuration* tim=nullptr,
            vector<int>const& Gdeg={}
            );

    int BGPC_OMP_SpeculativeIterative_ccbyG2deg   (int const side, int nT, int& colors, vector<int>& vtxColors, int const LOCAL_ORDER=ORDER_NONE, 
            int const NUM_RND_TC=0, int const NUM_NET_TC=0, int const NUM_NET_CC=0, bool const b_SEQ_HC=false, int const nVerbose=0);



    void do_BGPC_Phase_PlainCC_ccbyG2deg(vector<int>&vtxColors, vector<vector<int>>& QQ, vector<int>& Qsizes, vector<int>const& srcPtr, vector<int>const& dstPtr, vector<int>const& vtxVal, int* num=nullptr, ChronoDuration* tim=nullptr, 
            vector<int>const& G2deg={});


    void do_BGPC_compare_subgraph(
        vector<int>const& srcPtr, vector<int>const& dstPtr,  vector<int>const& vtxVal,
        vector<int>const&vtxColors1, vector<vector<int>>const&QQ1,  vector<int>const&Qsizes1,   
        vector<int>const&vtxColors2, vector<vector<int>>const&QQ2,  vector<int>const&Qsizes2, 
        vector<int>const&vtxColors3, vector<vector<int>>const&QQ3,  vector<int>const&Qsizes3,
        vector<int>const&vtxColors4, vector<vector<int>>const&QQ4,  vector<int>const&Qsizes4,
        vector<int>const&Gdeg, vector<int>const& G2deg, vector<int>const&GRnd
        );

    int BGPC_OMP_SpeculativeIterative_ccbydeg_debug_show_diff(
        int const side,             // L or R
        int nT,                     // number of threads
        int& colors,                // number of colors
        vector<int>& vtxColors,     // mapping vertex to color
        int const LOCAL_ORDER,      // Local Order
        int const NUM_RND_TC,       // Phase: Random Tentative Coloring
        int const NUM_NET_TC,       // Phase: Netbin Tentative Coloring
        int const NUM_NET_CC,       // Phase: Netbin Check Conflicts
        bool const b_SEQ_HC,        // using Seq-Coloring to solve the conflicts?
        int const nVerbose          // verbose level
        );
    
    // going to remove
    //int PD2_OMP_GMMP(const int side, int nT, int&color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int nVerbose=0);
    //int PD2_OMP_GMMP_debug_localorder_randomid(const int side, int nT, int&color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int nVerbose=0);
    //int PD2_OMP_GMMP_RandColor(const int side, int nT, int&color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int nVerbose=0);
    
    //// Developing Algorithms
    //int PD2_OMP_GM3P_BIT(const int side, int nT, int &color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int nVerbose=0);
    //int PD2_OMP_GMMP_BIT(const int side, int nT, int &color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int nVerbose=0);

    //// Algorithms
    //int PD2_OMP_GM3P_atomic(const int side, int nT, int &color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int nVerbose=0);
    //int PD2_OMP_GMMP_atomic(const int side, int nT, int &color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int nVerbose=0);


    //int PD2_OMP_GMMP_atomic_Delay(const int side, int nT, int &color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int nVerbose=0);


    //void rnd_tentative_coloring(vector<int>&vtxColor, const int MaxColor, vector<int> uncolored_vtxs, const int uncolored_vtxs_size,  ChronoDuration &tim);

    //// algs from Paper "Greed is good, Parallel Algorithms for Bipartite graph partial coloring on multicore architecures"
    //int ParallelColoring_NkNk(const int side, int nT, int &color, vector<int>&vtxColors, const int local_order=ORDER_NONE, const int NUM_NET_TENTATIVE=-1, const int NUM_NET_CHECK=-1, const int nVerbose=0, const int VBSIDE_GLB_ORDER=ORDER_NONE);
    
    /*
    // algs used in Greed is good ...
    void Net_Tentative_Coloring( 
            vector<int>& vtxColor, 
            const vector<int>& dstPtr, 
            const vector<int>& vtxVal, 
            const vector<int>& vBQ, 
            vector<vector<int>>& ForbiddenArrays, 
            const int BufSize,
            vector<ChronoDuration>& vtim_tents
            );

    void Net_Tentative_Coloring_v2( 
            vector<int>& vtxColor, 
            const vector<int>& dstPtr, 
            const vector<int>& vtxVal, 
            const vector<int>& vBQ, 
            vector<vector<int>>& ForbiddenArrays, 
            const int BufSize,
            vector<ChronoDuration>& vtim_tents
            );
    
    void Net_Tentative_Coloring_v3( 
            vector<int>& vtxColor, 
            const vector<int>& dstPtr, 
            const vector<int>& vtxVal, 
            const vector<int>& vBQ, 
            vector<vector<int>>& ForbiddenArrays, 
            const int BufSize,
            vector<ChronoDuration>& vtim_tents
            );

    void Net_Tentative_Coloring_v4( 
            vector<int>& vtxColor, 
            const vector<int>& dstPtr, 
            const vector<int>& vtxVal, 
            const vector<int>& vBQ, 
            vector<vector<int>>& ForbiddenArrays, 
            const int BufSize,
            vector<ChronoDuration>& vtim_tents
            );

    void Vtx_Tentative_Coloring(
            vector<int>& vtxColor, 
            const vector<int>& srcPtr, 
            const vector<int>& dstPtr,
            const vector<int>& vtxVal, 
            vector<int>& uncolored_vAQ, 
            vector<vector<int>>& ForbiddenArrays,
            const int BufSize,
            vector<ChronoDuration>& vtim_tents
            );

    void Net_Check_Conflicts(
            vector<int>& vtxColor,
            const vector<int>& dstPtr, 
            const vector<int>& vtxVal,
            const vector<int>& vBQ,
            vector<vector<int>> &ForbiddenArrays,
            const int BufSize,
            vector<ChronoDuration>& vtim_checks,
            vector<int>& vnum_Confs
            );

    void Vtx_Check_Conflicts_Delay(
            vector<int>& vtxColor, 
            const vector<int>& srcPtr, 
            const vector<int>& dstPtr, 
            const vector<int>& vtxVal, 
            vector<int>& uncolored_vAQ, 
            vector<vector<int>>& conflicts_Qs, 
            vector<ChronoDuration>& vtim_checks, 
            vector<int>& vnums_Confs
            );


    void Vtx_Check_Conflicts(
            vector<int>& vtxColor, 
            const vector<int>& srcPtr, 
            const vector<int>& dstPtr, 
            const vector<int>& vtxVal, 
            vector<int>& uncolored_vAQ, 
            vector<int>& conflicts_Q, 
            vector<ChronoDuration>& vtim_checks, 
            vector<int>& vnum_Confs
            );

 int ParallelColoring_RandomColoring_rseq(
        const int side,                 // sides, L or R
        int nT,                         // number of threads 
        int&colors,                     // number of colors
        vector<int>&vtxColor,           // vertex colors
        const int LOCAL_ORDER,          // Local_Order: NTR, RND, LF, SL,...
        const int NUM_NET_TENTATIVE,    // number of iterations on Net_tetative_colroring before switch to Vtx_tentative_coloring. -1 means infinity
        const int NUM_NET_CHECK,        // number of iterations on Net_check before switch to Vtx_check. -1 means infinity
        const int nVerbose,              // verbose level. 0 no verbose 1 basic verbose 2 more details verbose
        const int VB_GLB_ORDER         // vB side vertex, orders
        );
 int ParallelColoring_RandomColoring_romp(
        const int side,                 // sides, L or R
        int nT,                         // number of threads 
        int&colors,                     // number of colors
        vector<int>&vtxColor,           // vertex colors
        const int LOCAL_ORDER,          // Local_Order: NTR, RND, LF, SL,...
        const int NUM_NET_TENTATIVE,    // number of iterations on Net_tetative_colroring before switch to Vtx_tentative_coloring. -1 means infinity
        const int NUM_NET_CHECK,        // number of iterations on Net_check before switch to Vtx_check. -1 means infinity
        const int nVerbose,              // verbose level. 0 no verbose 1 basic verbose 2 more details verbose
        const int VB_GLB_ORDER         // vB side vertex, orders
        );


void Random_Tentative_Coloring_seq(
        vector<int>& vtxColor,                    // vertex Colors
        vector<int>& uncolored_vAQ,
        const int ColorLB,
        vector<ChronoDuration>& vtim_tents                     // record runtime for this modular and append to the end 
        );
    
void Random_Tentative_Coloring_omp(
        vector<int>& vtxColor,                    // vertex Colors
        vector<int>& uncolored_vAQ,
        const int ColorLB,
        vector<ChronoDuration>& vtim_tents                     // record runtime for this modular and append to the end 
        );
    */

//private: 
//    bool color_in_array(const int c, const vector<int> &F){  for(const auto &x : F) { if(x==c) return true; } return false; };
//    bool color_in_array_2(const int c, const vector<int> &F, const int Fsize){  for(int i=0; i<Fsize; i++){ if(F[i]==c) return true; } return false; };



public: // Utilites
    int cnt_pd2conflict(const int side, const vector<int>& vc);
    int get_lowbound_coloring(const int side);  //return other side's max degree
    double get_std_degree(const int side);  
    int calc_num_colors_from_vtx_colors(const vector<int>& vtxColor);
    
    void Translate_OrderId_To_OrderTag(const int LOCAL_ORDER, string&lotag);
    string Translate_OrderId_To_OrderTag(const int LOCAL_ORDER);


    void init_QQ_and_Qsizes(vector<vector<int>>& QQ, vector<int>&Qsizes, int const nT, int const Nsrc,  vector<int>const& ordered_queue, chrono::duration<double>* tim=nullptr);

    void init_QQ_and_Qsizes(vector<vector<int>>& QQ, vector<int>&Qsizes, int const nT, int const Nsrc,  vector<int>const& ordered_queue, vector<int>const& vtxColors, chrono::duration<double>* tim=nullptr);

    void do_Generate_Weights(vector<uint64_t>& Weights, int const N, int const WEIGTH_TYPE, int const* srcPtr=nullptr, int const* dstPtr=nullptr, int const* vtxVal=nullptr, ChronoDuration* tim=nullptr);



    void do_BGPC_generate_random_BGPC(int const A, int const B, int const E, int const rseed=1, bool const bRepeatEdge=false, ChronoDuration* tim=nullptr);

    void do_BGPC_generate_random_seq_repeat(vector<int>& srcPtr, vector<int>& dstPtr, vector<int>& vtxVal, 
            int const rseed, int const M, int const N, int const E, ChronoDuration *tim=nullptr);
    void do_BGPC_generate_random_seq_noRepeat(vector<int>& srcPtr, vector<int>& dstPtr, vector<int>& vtxVal,
            int const rseed, int const M, int const N, int const E, ChronoDuration *tim=nullptr);

        
protected:
    //int m_total_num_colors;
    //vector<int> m_vertex_color;
    string m_method;
}; // end of class SMPGCInterface


}// endof namespace ColPack
#endif

