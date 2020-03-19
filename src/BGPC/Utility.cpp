/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "BGPCColoring.h"
#include <chrono> //c++11 system time
#include <random> //c++11 random
#include <unordered_set>
#include <unordered_map>
using namespace std;
using namespace ColPack;


// ============================================================================
// check if the graph is correct colored
// ----------------------------------------------------------------------------
// return number_of_uncolor_vertex + number_of_incorrect_colored_vertex
// ============================================================================
int BGPCColoring::cnt_pd2conflict(const int side, const vector<int>& vtxColor){
    const int          N             = (side==BGPC::L)?(GetLeftVertexCount()):(GetRightVertexCount()); 
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    const vector<int>& vtxVal        = GetEdges();
    const vector<int>& ordered_queue = get_ordered_queue_A_const(); //global_ordered_vertex(side);
    
    if(N!=(signed)vtxColor.size() || N!=(signed)ordered_queue.size()){
        printf("Error! number of vertex in %s side %d is not according to vtxColor size %d or ordered_queue size %d.(graph %s)\n",
                (side==BGPC::L?"L":"R"), N, (signed)vtxColor.size(), (signed)ordered_queue.size(), m_s_InputFile.c_str());
        exit(1);
    }

    int conflicts=0;
    for (int i=0; i<N; i++){
        const int v = ordered_queue[i];
        const int vc= vtxColor[v];
        if(vc<0) { conflicts++; continue; }
        for(int iw=srcPtr[v]; iw!=srcPtr[v+1]; iw++){
            const int w = vtxVal[iw];
            for(int iu=dstPtr[w]; iu!=dstPtr[w+1]; iu++){
                const int u = vtxVal[iu];
                if(v==u) continue;
                const int uc = vtxColor[u];
                if(v<w && vc==uc)
                    conflicts++;
            }
        }
    }
    return conflicts;
}

// ============================================================================
// get the low bound of number of coloring
// ----------------------------------------------------------------------------
// the uncolored side's max degree is the low bound of number of coloring
// ============================================================================
int BGPCColoring::get_lowbound_coloring(const int side){
    return (side==BGPC::L)?GetMaximumRightVertexDegree():GetMaximumLeftVertexDegree();
}


// ============================================================================
// author: xin cheng
// ----------------------------------------------------------------------------
// return sample stand deviation
// ============================================================================
double BGPCColoring::get_std_degree(const int side){
    if(side==BGPC::L){
        const int N = GetLeftVertexCount();
        const double mean = m_d_AverageLeftVertexDegree;
        const vector<int>& srcPtr = GetLeftVertices();
        double sum=.0;
        for(int v=0; v<N; v++){
            const double deg = (double)-(srcPtr[v]-srcPtr[v+1]);
            sum+=(deg-mean)*(deg-mean);
        }
        return sqrt(sum/(N-1));
    }
    else if(side==BGPC::R) {
        const int N = GetRightVertexCount();
        const double mean = m_d_AverageRightVertexDegree;
        const vector<int>& srcPtr = GetRightVertices();
        double sum=.0;
        for(int v=0; v<N; v++){
            const double deg = (double)-(srcPtr[v]-srcPtr[v+1]);
            sum+=(deg-mean)*(deg-mean);
        }
        return sqrt(sum/(N-1));
    }else{
        const int NL = GetLeftVertexCount();
        const int NR = GetRightVertexCount();
        const int mean = m_d_AverageVertexDegree;
        const vector<int>& srcLPtr = GetLeftVertices();
        const vector<int>& srcRPtr = GetRightVertices();
        double sum=.0;
        for(int v=0; v<NL; v++){
            const double deg = (double)-(srcLPtr[v]-srcLPtr[v+1]);
            sum += (deg-mean)*(deg-mean);
        }
        for(int v=0; v<NR; v++){
            const double deg = (double)-(srcRPtr[v]-srcRPtr[v+1]);
            sum += (deg-mean)*(deg-mean);
        }
        return sqrt(sum/(NL+NR-1));
    }
}

//
// author: xin cheng
// get number of distinc colors used in vtxColor.
//
int BGPCColoring::calc_num_colors_from_vtx_colors(const vector<int>&vtxColor){
    unordered_set<int> S(vtxColor.begin(), vtxColor.end());
    return S.count(-1)?(S.size()-1):(S.size());
}



// ============================================================================
// author: xin cheng
//
// ============================================================================
string BGPCColoring::Translate_OrderId_To_OrderTag(const int LOCAL_ORDER){
    string lotag="unknow";
    switch(LOCAL_ORDER){
        case ORDER_NONE:    lotag="None"; break;
        case ORDER_NATURAL: lotag="NT"; break;
        case ORDER_RANDOM:  lotag="RD"; break;
        case ORDER_LARGEST_FIRST: lotag="LF"; break;
        default:
            printf("Error! ColPack::BGPC tring to use local order %d. which is not supported.\n",LOCAL_ORDER); 
            exit(1);
    }
    return lotag;

}


// ============================================================================
// author: xin cheng
//
// ============================================================================
void BGPCColoring::Translate_OrderId_To_OrderTag(const int LOCAL_ORDER, string&lotag){
    switch(LOCAL_ORDER){
        case ORDER_NONE:    lotag="None"; break;
        case ORDER_NATURAL: lotag="NT"; break;
        case ORDER_RANDOM:  lotag="RD"; break;
        case ORDER_LARGEST_FIRST: lotag="LF"; break;
        default:
            printf("Error! ColPack::BGPC tring to use local order %d. which is not supported.\n",LOCAL_ORDER); 
            exit(1);
    }
}

// ============================================================================
// Author: Xin Cheng
// init QQ and Qsizes
// ============================================================================
void BGPCColoring::init_QQ_and_Qsizes(vector<vector<int>>&QQ, vector<int>& Qsizes, int const nT, int const Nsrc,  vector<int>const& ordered_queue, chrono::duration<double>* tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    
    // init Qsizes
    Qsizes.assign(nT, Nsrc/nT);
    for(int i=0; i<Nsrc%nT; i++) Qsizes[i]++;

    // calculate displacement
    vector<int> disp(nT+1, 0);
    for(int i=0; i<nT; i++) disp[i+1] = disp[i] + Qsizes[i];
    
    // init QQ
    for(int i=0; i<nT; i++)
        QQ[i].assign(ordered_queue.begin()+disp[i], ordered_queue.begin()+disp[i+1]);

    if(tim) *tim = chrono::steady_clock::now() - start;
}


// ============================================================================
// Author: Xin Cheng
// init QQ and Qsizes
// ============================================================================
void BGPCColoring::init_QQ_and_Qsizes(vector<vector<int>>& QQ, vector<int>&Qsizes, int const nT, int const Nsrc,  vector<int>const& ordered_queue, vector<int>const& vtxColors, chrono::duration<double>* tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    
    // init Qsizes
    Qsizes.assign(nT, 0);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int Qsize=0;
        #pragma omp for
        for(int iu=0; iu<Nsrc; iu++) {
            int const u = ordered_queue[iu];
            if( vtxColors[u] < 0)
                QQ[tid][Qsize++] = u;
        }
        Qsizes[tid] = Qsize;
    }
    if(tim) *tim = chrono::steady_clock::now() - start;
}




// ============================================================================
// Author: Xin Cheng
// generate weights for independet set
// ============================================================================
void BGPCColoring::do_Generate_Weights(vector<uint64_t>& Weights, int const N, int const WEIGHT_TYPE, int const* srcPtr, int const* dstPtr, int const* vtxVal, ChronoDuration* tim){
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();

    if(WEIGHT_TYPE == WEIGHT_RANDOM){
        for(int i=0; i<N; i++) Weights[i]=i;
        for(int i=0; i<N-1; i++) {
            uniform_int_distribution<int> dist(i, N-1);
            swap(Weights[i], Weights[dist(m_mt)]);
        }
    }
    else if(WEIGHT_TYPE == WEIGHT_DEGREE_D1) {
       #pragma omp parallel for
        for(int i=0; i<N; i++){ 
            Weights[i] = (((uint64_t)(srcPtr[i+1]-srcPtr[i]))<<40)|(i&0xFFFFFFFFFF) ;
        }
    }
    else if(WEIGHT_TYPE == WEIGHT_DEGREE_D2) {
        #pragma omp parallel for
        for(int u=0; u<N; u++){ 
            unordered_set<int> unbs;
            for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
                int const v = vtxVal[iv];
                for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                    unbs.insert(vtxVal[iw]);
                }
            }
            Weights[u] = unbs.size();
        }
    }
    else{
        cerr<<"Err! unknow of WEIGHT_TYPE "<<WEIGHT_TYPE<<" Only support WEIGHT_RANDOM or WEIGHT_DEGREE.\n";
        exit(1);
    }
    if(tim) *tim  = chrono::steady_clock::now() - start;
}




// ============================================================================
// Author: Xin Cheng
// generate unique weights for independet set
// ============================================================================
void BGPCColoring::do_BGPC_generate_random_BGPC(int const A, int const B, int const E, int const RSEED, bool const bRepeatEdge, ChronoDuration* tim){
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
    
    // clear the relative variable
    m_ordered_queue_A.clear();
    m_ordered_queue_A_side=BGPC::L;
   
    m_i_MaximumLeftVertexDegree=0;
    m_i_MaximumRightVertexDegree=0;
    m_i_MaximumVertexDegree=0;

    m_i_MinimumLeftVertexDegree=0;
    m_i_MinimumRightVertexDegree=0;
    m_i_MinimumVertexDegree=0;

    m_d_AverageLeftVertexDegree=.0;
    m_d_AverageRightVertexDegree=.0;
    m_d_AverageVertexDegree=.0;

    m_s_InputFile="";

    m_vi_LeftVertices.clear();
    m_vi_RightVertices.clear();
    m_vi_Edges.clear();

    m_mimi2_VertexEdgeMap.clear();

    // allocate the memory
    m_ordered_queue_A.reserve(A);
    m_vi_LeftVertices.reserve(A+1);
    m_vi_RightVertices.reserve(B+1);
    m_vi_Edges.reserve(2*E);

    if(A<=0 || B<=0) return;
    
    // begin the generate the graphs
    if(bRepeatEdge) do_BGPC_generate_random_seq_repeat  (m_vi_LeftVertices, m_vi_RightVertices, m_vi_Edges, RSEED, A, B, E, tim);
    else            do_BGPC_generate_random_seq_noRepeat(m_vi_LeftVertices, m_vi_RightVertices, m_vi_Edges, RSEED, A, B, E, tim);
    
    // fill in the degree information
    {
        int left_deg_sum  = m_vi_LeftVertices.back() - m_vi_LeftVertices.front();
        int right_deg_sum = m_vi_LeftVertices.back() - m_vi_LeftVertices.front();
        m_d_AverageLeftVertexDegree = left_deg_sum *1.0/A;
        m_d_AverageRightVertexDegree= right_deg_sum*1.0/B;
        m_d_AverageVertexDegree = (left_deg_sum+right_deg_sum)*1.0/(A+B);

        m_i_MinimumLeftVertexDegree =m_i_MaximumLeftVertexDegree =m_vi_LeftVertices [1] - m_vi_LeftVertices [0];
        m_i_MinimumRightVertexDegree=m_i_MaximumRightVertexDegree=m_vi_RightVertices[1] - m_vi_RightVertices[0];
        for(int i=1; i<A; i++){
            int deg = m_vi_LeftVertices[i+1] - m_vi_LeftVertices[i];
            if(deg<m_i_MinimumLeftVertexDegree)
                m_i_MinimumLeftVertexDegree=deg;
            else if(deg>m_i_MaximumLeftVertexDegree)
                m_i_MaximumLeftVertexDegree=deg;
        }
        for(int i=1; i<B; i++){
            int deg = m_vi_RightVertices[i+1] - m_vi_RightVertices[i];
            if(deg<m_i_MinimumRightVertexDegree)
                m_i_MinimumRightVertexDegree=deg;
            else if(deg>m_i_MaximumRightVertexDegree)
                m_i_MaximumRightVertexDegree=deg;
        }
        m_i_MinimumVertexDegree = min(m_i_MinimumLeftVertexDegree, m_i_MinimumRightVertexDegree);
        m_i_MaximumVertexDegree = max(m_i_MaximumLeftVertexDegree, m_i_MaximumRightVertexDegree);

    }

    // order the graph with the nature ordering
    m_ordered_queue_A.assign(A,0);
    #pragma omp parallel for
    for(int i=0; i<A; i++)
        m_ordered_queue_A[i]=i;

    if(tim) *tim = chrono::steady_clock::now() - start;
}


// ============================================================================
// Author: Xin Cheng
// ============================================================================
void BGPCColoring::do_BGPC_generate_random_seq_repeat(vector<int>& srcPtr, vector<int>& dstPtr, vector<int>& vtxVal, 
        int const rseed, int const M, int const N, int const E,
        ChronoDuration *tim){
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();

    srcPtr.assign(M+1,0);
    dstPtr.assign(N+1,0);
    vtxVal.clear(); vtxVal.reserve(2*E);
    mt19937 mt(rseed);
    //uniform_int_distribution<int> dist (0,M*N-1);
    uniform_int_distribution<int> distM(0,M-1);
    uniform_int_distribution<int> distN(0,N-1);
    
    // generate the edges
    {
        unordered_map<int, unordered_set<int>> GA;
        unordered_map<int, unordered_set<int>> GB;
        for(int i=0; i<E; i++) {
            //int rnd = dist(mt);
            //int row = rnd / M;
            //int col = rnd % M;
            int row = distM(mt);
            int col = distN(mt);
            GA[row].insert(col);
            GB[col].insert(row);
        }

        // GA -> CSR part 1
        for(int i=0; i<M; i++) {
            int beg = srcPtr[i] = (signed)vtxVal.size();
            auto xit = GA.find(i);
            if(xit==GA.end()) continue;
            vtxVal.insert(vtxVal.end(), (xit->second).begin(), (xit->second).end());
            sort(vtxVal.begin()+beg, vtxVal.end());
        }
        srcPtr[M] = (signed)vtxVal.size();
        GA.clear();

        // GB -> CSR part 2
        for(int i=0; i<N; i++) {
            int beg = dstPtr[i] = (signed)vtxVal.size();
            auto xit = GB.find(i);
            if(xit==GB.end()) continue;
            vtxVal.insert(vtxVal.end(), (xit->second).begin(), (xit->second).end());
            sort(vtxVal.begin()+beg, vtxVal.end());
        }
        dstPtr[N] = (signed)vtxVal.size();
        GB.clear();
    }


    if(tim) *tim =  chrono::steady_clock::now() - start;
}

// ============================================================================
// Author: Xin Cheng
// ============================================================================
void BGPCColoring::do_BGPC_generate_random_seq_noRepeat(vector<int>& srcPtr, vector<int>& dstPtr, vector<int>& vtxVal,
        int const rseed, int const M, int const N, int const E,
        ChronoDuration *tim) {
    chrono::steady_clock::time_point start;
    if(tim) start = chrono::steady_clock::now();
 
    srcPtr.assign(M+1,0);
    dstPtr.assign(N+1,0);
    vtxVal.clear(); vtxVal.reserve(2*E);
    int const Alpha = E/M;
    //double const Alpha = E*1.0/M;  // Expected average left degree
    //double const p = Alpha/N = E/(MN)
    mt19937 mt(rseed);
    //uniform_real_distribution<double> dist(0,N);
    uniform_int_distribution<int> dist(0,N-1);
    // generate the edges
    {
        unordered_map<int, vector<int>> GB;
        for(int i=0; i<M; i++){
            srcPtr[i]=(signed)vtxVal.size();
            for(int j=0; j<N; j++){
                auto rnd = dist(mt);
                if(rnd<Alpha){
                    vtxVal.push_back(j);
                    GB[j].emplace_back(i);
                }
            }
        }
        srcPtr[M]=(signed)vtxVal.size();

        for(int i=0; i<N; i++){
            dstPtr[i]=(signed)vtxVal.size();
            vtxVal.insert(vtxVal.end(), GB[i].begin(), GB[i].end());
        }
        dstPtr[N]=(signed)vtxVal.size();
        GB.clear();
    }

    if(tim) *tim =  chrono::steady_clock::now() - start;
}






