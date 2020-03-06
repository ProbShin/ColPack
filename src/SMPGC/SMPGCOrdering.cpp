/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "SMPGCOrdering.h"
using namespace std;
using namespace ColPack;

//
// For graph orderings, following functions are supported
//   * seq_{whole,subset}_{ntr,rnd,ldf,sdl}
//   * omp_{whole,prepartition_{ntr,rnd,ldf,sdl}
// i.e.
//     seq_whole_ntr
//     seq_whole_rnd
//     seq_whole_ldf
//     seq_whole_sdf
//     seq_subset_ntr
//     seq_subset_rnd
//     seq_subset_ldf
//     seq_subset_sdl
//  x  omp_whole_ntr    // ambiguous, thus did not implement
//     omp_whole_rnd
//     omp_whole_ldf
//     omp_whole_sdf
//     omp_prepartition_ntr
//     omp_prepartition_rnd
//     omp_prepartition_ldf
//     omp_prepartition_sdl

// ============================================================================
// Author: Xin Cheng
// Construction
// ============================================================================
SMPGCOrdering::SMPGCOrdering(const string& graph_name, const string& fmt, ChronoDuration* iotime,  const int order=ORDER_NATURAL, ChronoDuration* ordtime) 
: SMPGCGraph(graph_name, fmt, iotime) {
    m_ordered_vertex.assign(num_nodes(),0);

    switch(order){
        case ORDER_NONE:
        case ORDER_NATURAL:       seq_whole_ntr_ordering(ordtime); break;
        case ORDER_RANDOM:        seq_whole_rnd_ordering(ordtime); break;
        case ORDER_LARGEST_FIRST: seq_whole_ldf_ordering(ordtime); break;
        case ORDER_SMALLEST_LAST: seq_whole_sdl_ordering(ordtime); break;
        default:
            cout<<"The Order of 'SEQ WHOLE "<<Translate_OrderId_To_OrderTag(order)<<"' does not supported for now."<<endl;
                    exit(1);
    }
}

// ============================================================================
// Author: Xin Cheng
// Destruction
// ============================================================================
SMPGCOrdering::~SMPGCOrdering(){}

// ============================================================================
// Author: Xin Cheng
// Ordering selected function
// ============================================================================
void SMPGCOrdering::ordering(const string& s1, const string& s2, const int order, ChronoDuration* ordtime, void *parg){
    if(s1 == "SEQ"){
        if(s2 == "WHOLE"){
            switch(order){
                case ORDER_NONE:          if(ordtime) *ordtime=ChronoDuration(.0); m_ordered_method=ORDER_NONE;  break;
                case ORDER_NATURAL:       seq_whole_ntr_ordering(ordtime); break;
                case ORDER_RANDOM:        seq_whole_rnd_ordering(ordtime); break;
                case ORDER_LARGEST_FIRST: seq_whole_ldf_ordering(ordtime); break;
                case ORDER_SMALLEST_LAST: seq_whole_sdl_ordering(ordtime); break;
                default:
                    cout<<"The Order of '"<<s1<<"', '"<<s2<<"','"<<Translate_OrderId_To_OrderTag(order)<<"' does not supported for now."<<endl;
                    exit(1);
            }
        }
        else if(s2 =="SUBSET" )
            switch(order){
                case ORDER_NONE:          if(ordtime) *ordtime=ChronoDuration(.0); break;
                case ORDER_NATURAL:       seq_subset_ntr_ordering(*(vector<int>*)parg ,ordtime); break;
                case ORDER_RANDOM:        seq_subset_rnd_ordering(*(vector<int>*)parg ,ordtime); break;
                case ORDER_LARGEST_FIRST: seq_subset_ldf_ordering(*(vector<int>*)parg ,ordtime); break;
                case ORDER_SMALLEST_LAST: seq_subset_sdl_ordering(*(vector<int>*)parg ,ordtime); break;
                default:
                    cout<<"The Order of '"<<s1<<"', '"<<s2<<"','"<<Translate_OrderId_To_OrderTag(order)<<"' does not supported for now."<<endl;
                    exit(1);
            }

        else{
            cout<<"The Order of '"<<s1<<"', '"<<s2<<"','"<<order<<"' does not supported for now."<<endl;
        }
    }
    else if(s1 == "OMP"){
        if(s2 == "WHOLE"){
            switch(order){
                case ORDER_NONE:          if(ordtime) *ordtime=ChronoDuration(.0); break;
                case ORDER_NATURAL:       omp_whole_ntr_ordering(ordtime); break;
                //case ORDER_RANDOM:        omp_whole_rnd_ordering(ordtime); break; // it is ambigous , so did not implement 
                case ORDER_LARGEST_FIRST: omp_whole_ldf_ordering(ordtime); break;
                case ORDER_SMALLEST_LAST: omp_whole_sdl_ordering(ordtime); break;
                default:
                    cout<<"The Order of '"<<s1<<"', '"<<s2<<"','"<<Translate_OrderId_To_OrderTag(order)<<"' does not supported for now."<<endl;
                    exit(1);
            }
        }
        else if(s2 =="PREPARTITION" )
            switch(order){
                case ORDER_NONE:          if(ordtime) *ordtime=ChronoDuration(.0); break;
                case ORDER_NATURAL:       omp_prepartition_ntr_ordering(*(vector<vector<int>>*)parg ,ordtime); break;
                case ORDER_RANDOM:        omp_prepartition_rnd_ordering(*(vector<vector<int>>*)parg ,ordtime); break;
                case ORDER_LARGEST_FIRST: omp_prepartition_ldf_ordering(*(vector<vector<int>>*)parg ,ordtime); break;
                case ORDER_SMALLEST_LAST: omp_prepartition_sdl_ordering(*(vector<vector<int>>*)parg ,ordtime); break;
                default:
                    cout<<"The Order of '"<<s1<<"', '"<<s2<<"','"<<Translate_OrderId_To_OrderTag(order)<<"' does not supported for now."<<endl;
                    exit(1);
            }

        else{
            cout<<"The Order of '"<<s1<<"', '"<<s2<<"','"<<Translate_OrderId_To_OrderTag(order)<<"' does not supported for now."<<endl;
        }
    }
    else{
        cout<<"The Order of '"<<s1<<"', '"<<s2<<"','"<<Translate_OrderId_To_OrderTag(order)<<"' does not supported for now."<<endl;
    }

    //if(ordtime){ *(time_t*)ordtime+=clock(); *ordtime =(double)(*(time_t*)ordtime)/CLOCKS_PER_SEC; }
}


// ============================================================================
// author: xin cheng
//
// ============================================================================
void SMPGCOrdering::Translate_OrderId_To_OrderTag(const int LOCAL_ORDER, string&lotag){
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
// author: xin cheng
//
// ============================================================================
string SMPGCOrdering::Translate_OrderId_To_OrderTag(const int LOCAL_ORDER){
    switch(LOCAL_ORDER){
        case ORDER_NONE:    return "None"; break;
        case ORDER_NATURAL: return "NT"; break;
        case ORDER_RANDOM:  return "RD"; break;
        case ORDER_LARGEST_FIRST: return "LF"; break;
        default:
            printf("Error! ColPack::BGPC tring to use local order %d. which is not supported.\n",LOCAL_ORDER); 
            exit(1);
    }
}




