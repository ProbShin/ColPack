/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#include "BGPC.h"
#include "BGPCOrdering.h"
using namespace std;
using namespace ColPack;


// ============================================================================
// Construction
// Author: Xin Cheng
// ============================================================================
BGPCOrdering::BGPCOrdering(
        const string& graph_name, 
        const string& fmt, 
        ChronoDuration* iotim, 
        const int side, 
        const int ORDER, 
        ChronoDuration* odtim) 
{
    if( fmt==FORMAT_MM) {
        chrono::steady_clock::time_point start;
        if(iotim) start = chrono::steady_clock::now();
        ReadMMBipartiteGraphCpp11(graph_name);
        if(iotim) *iotim = chrono::steady_clock::now() - start;
    }else if( fmt == FORMAT_POTHEN){
        chrono::steady_clock::time_point start;
        if(iotim) start = chrono::steady_clock::now();
        ReadMMGeneralGraphIntoPothenBipartiteGraphCpp11(graph_name);
        if(iotim) *iotim = chrono::steady_clock::now() - start;
    }
    else{
        printf("Error!PD2SMPGCOrdering only support MatrixMarket format (MM), or pothen. \"%s\" with fmt \"%s\" is under construction... \n", graph_name.c_str(), fmt.c_str());
        exit(1);
    }
    const int          N             = (side==BGPC::L)?(GetLeftVertexCount()):(GetRightVertexCount()); 
    m_ordered_queue_A_side = side;
    m_ordered_queue_A.assign(N, 0);
    #pragma omp parallel for
    for(int i=0; i<N; i++) 
        m_ordered_queue_A[i]=i;

    if(ORDER!=BGPC::ORDER_NATURAL) 
    Ordering_NoPartition(side, m_ordered_queue_A, N, ORDER, odtim);
}

BGPCOrdering::~BGPCOrdering(){}



// ============================================================================
// Author: Xin Cheng
// Ordering
// ============================================================================
void BGPCOrdering::Ordering_NoPartition(
        const int side,
        vector<int>& Q, 
        const int Qsize,
        const int ORDER,
        ChronoDuration* tim) {
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    const vector<int>& vtxVal        = GetEdges();
    switch(ORDER){
        case ORDER_NONE:
            if(tim) *tim = ChronoDuration::zero(); //chrono::duration<double> zero(0);
            break;
        case ORDER_NATURAL: 
            {
                do_Natural_Ordering_NoPartition(Q, Qsize, tim);
            }
            break;
        case ORDER_RANDOM:  
            {
                do_Random_Ordering_NoPartition(Q, Qsize, tim);
            }
            break;
        case ORDER_LARGEST_FIRST:
            {
                do_D1_LargestDegreeFirst_Ordering_NoPartition(Q, Qsize, srcPtr, tim);
            }
            break;
        case ORDER_D2_LARGEST_FIRST:
            {
                do_D2_LargestDegreeFirst_Ordering_NoPartition(Q, Qsize, srcPtr, dstPtr, vtxVal, tim);
            }
        break;
        case ORDER_D2ROUGH_LARGEST_FIRST:
            {
                do_D2Rough_LargestDegreeFirst_Ordering_NoPartition(Q, Qsize, srcPtr, dstPtr, vtxVal, tim);
            }
            break;
        case ORDER_BA_LARGEST_FIRST:
            {   
                do_BA_LargestDegreeFirst_Ordering_NoPartition(Q, Qsize, srcPtr, dstPtr, vtxVal, tim);
            }
            break;
        default:
            cout<<"Error! BGPCOrdering::Ordering does not accept '"<<ORDER<<"'"<<endl;
            exit(1);
    }
}



// ============================================================================
// Author: xin cheng
// Order
// ============================================================================
void BGPCOrdering::Ordering_PrePartition_OMP(
        int const side,
        vector<vector<int>>& QQ, 
        vector<int>const& Qsizes,
        const int ORDER, 
        ChronoDuration* tim) {
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    const vector<int>& vtxVal        = GetEdges();
    switch(ORDER){
        case ORDER_NONE: 
            if(tim) *tim = ChronoDuration::zero();
            break;
        case ORDER_NATURAL:
            do_Nature_Ordering_PrePartition(QQ, Qsizes, tim); 
            break;
        case ORDER_RANDOM:
            do_Random_Ordering_PrePartition(QQ, Qsizes, tim); 
            break;
        case ORDER_LARGEST_FIRST:
            do_D1_LargestDegreeFirst_Ordering_PrePartition(QQ, Qsizes, srcPtr, tim);
            break;
        case ORDER_D2_LARGEST_FIRST:
            do_D2_LargestDegreeFirst_Ordering_PrePartition(QQ, Qsizes, srcPtr, dstPtr, vtxVal, tim);
            break;
        case ORDER_D2ROUGH_LARGEST_FIRST:
            do_D2Rough_LargestDegreeFirst_Ordering_PrePartition(QQ, Qsizes, srcPtr, dstPtr, vtxVal, tim);
            break;
        case ORDER_BA_LARGEST_FIRST:
            do_BA_LargestDegreeFirst_Ordering_PrePartition(QQ, Qsizes, srcPtr, dstPtr, vtxVal, tim);
            break;
        default:
            cout<<"Error! ColPack::BGPCOrdering::Ordering doest not accept '"<<ORDER<<"'"<<endl;
            exit(1);
    }
    return;
}

