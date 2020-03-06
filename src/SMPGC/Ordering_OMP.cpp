/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "SMPGCOrdering.h"
using namespace std;
using namespace ColPack;

// For graph orderings, following are supported
// seq_{whole,subset}_{ntr,rnd,ldf,sdl}
// omp_{whole,prepartition_{ntr,rnd,ldf,sdl}
// i.e.
//     seq_whole_ntr
//     seq_whole_rnd
//     seq_whole_ldf
//     seq_whole_sdf
//     seq_subset_ntr
//     seq_subset_rnd
//     seq_subset_ldf
//     seq_subset_sdl          // ambiguous, so did not implement      
//     omp_whole_ntr  
//   x omp_whole_rnd           // ambiguous, so did not implement
//     omp_whole_ldf
//     omp_whole_sdf
//     omp_prepartition_ntr
//     omp_prepartition_rnd
//     omp_prepartition_ldf
//     omp_prepartition_sdl


// ============================================================================
// Author: Xin Cheng
// omp parallel natural ordering on the whole vertex set
// store the ordered vertex into the class varibale "m_ordered_vertex"
// OrderedVertex <- [0,1,2,... N-1]
// ============================================================================
void SMPGCOrdering::omp_whole_ntr_ordering(ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const int N = num_nodes();
    m_ordered_vertex.resize(N);
    #pragma omp parallel for
    for(int i=0; i<N; i++) 
        m_ordered_vertex[i]=i;
    m_ordered_method = ORDER_NATURAL;
    if(ordtime) { end = std::chrono::steady_clock::now();  *ordtime = end - start; }
}


// ============================================================================
// Author: Xin Cheng
// omp parallel  ordering on the whole vertex set
// store the ordered vertex into the class variable "m_ordered_vertex"
// Largest Degree First
// ============================================================================
void SMPGCOrdering::omp_whole_ldf_ordering(ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const int N = num_nodes();              
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const int MaxDegreeP1 = max_degree()+1;   // max degree plus one
    int nT =0;
    #pragma omp parallel
        #pragma omp master 
        {
            nT = omp_get_num_threads();
        }
    vector<vector<vector<int>>>  ThreadGroupedVertexDegree(nT); // ListOfList per thread
    
    // clear up the memory
    m_ordered_vertex.clear();
    m_ordered_method = ORDER_LARGEST_FIRST; 
    
    // fill in the ListOfList
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<vector<int>>& GroupedVertexDegree = ThreadGroupedVertexDegree[tid];
        int maxDegreeSoFar=-1;
        #pragma omp for
        for(int u=0; u<N; u++){
            const int ud = -vtxPtr[u] + vtxPtr[u+1];
            // allocate memory only if necessary
            if(maxDegreeSoFar<ud){
                GroupedVertexDegree.resize(ud+1);
                maxDegreeSoFar=ud;
            }
            GroupedVertexDegree[ud].push_back(u);
        }
    }
   
    // fill in the orderd queue with Largest Degree First
    for(int d=MaxDegreeP1-1, i=0; i<MaxDegreeP1; d--, i++){
        for(int tid=0; tid<nT; tid++)
            if(d<(signed)((ThreadGroupedVertexDegree[tid]).size()))
                m_ordered_vertex.insert(m_ordered_vertex.end(), ThreadGroupedVertexDegree[tid][d].begin(), ThreadGroupedVertexDegree[tid][d].end());
    }
    if(ordtime) { end = std::chrono::steady_clock::now();  *ordtime = end - start; }
}

// ============================================================================
// Author: Xin Cheng
// omp parallel smallest last ordering on the whole vertex set
// store the ordered vertex into the class variable "m_ordered_vertex"
//
// the algorithm is the modifcation based on the algorithm 3 of the paper
//   "New Multithreaded Ordering and Coloring Algorithms for Multicore Architectures"
//    by M.M.A.Patwary, A.H.Gebremedhin, A.Pothen 
//    The idea is make the key varibale ListOfList distributed on each thread.
// 
// Some modifications have made for the efficience.
//
// Note:
//   the result of algorithm 3 is not idential to the Smallest Degree Last ordering
//   it is a kind of compromise on the accuracy to the parallel run time.
//   The comprimize is on the across edge vertex degree. 
// ============================================================================
void SMPGCOrdering::omp_whole_sdl_ordering(ChronoDuration* ordtime){
    //int nT=0;
    //#pragma omp parallel
    //    #pragma omp master
    //    {
    //        nT = omp_get_num_threads();
    //    }
    m_ordered_vertex.clear();
    
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    
    const int N = num_nodes();
    const vector<int>& vtxPtr = get_CSR_ia();   // csr format edges
    const vector<int>& vtxVal = get_CSR_ja();   // csr format value

    vector<int> vtxDegree(N);            // mapping from vertex to its degree,   1st index 
    vector<int> vtxLocation(N);          // mapping from vertex to its location  2nd index
    vector<int> vtxTid(N);                // mapping from vertex to its thread id.
    
    m_ordered_method=ORDER_SMALLEST_LAST;
    m_ordered_vertex.assign(N,-1);
    int outputIndex = N;     // the postion of the last written into the ordered_queue

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<vector<int>> GroupedVertexDegree;  // ListOfList  [degree][vertex]
        int num_of_local_vertex = 0;  // number of unordered vertex for this thread
        
        // init the GVDs, vtxDegree and vtxLocation in parallel 
        int maxDegreeSoFar = -1;
        #pragma omp for
        for(int u=0; u<N; u++){
            const int udeg = -vtxPtr[u] + vtxPtr[u+1];
            vtxDegree[u]   = udeg;
            // allocate memory only if necessary
            if(maxDegreeSoFar<udeg){
                GroupedVertexDegree.resize(udeg+1);
                maxDegreeSoFar=udeg;
            }
            vtxLocation[u] = GroupedVertexDegree[udeg].size();
            GroupedVertexDegree[udeg].push_back(u);
            num_of_local_vertex++;
            vtxTid[u]=tid;
        }// end omp for

        const int MaxDegreeP1 = GroupedVertexDegree.size();
        int SearchStartDegree=0; //skip the unnecessary-checks for some degrees 
        while(num_of_local_vertex!=0){ // equals to zero means output Q is full, means all vertex have been ordered.
            
            // if there are isolated vertex, move them to the output Q
            if(SearchStartDegree==0 && (GroupedVertexDegree[0].empty())==false){
                const int num_iso_vtxs = GroupedVertexDegree[0].size();
                int whereToInput = __sync_fetch_and_sub(&outputIndex, num_iso_vtxs);  // atomic operation
                for(auto u : GroupedVertexDegree[0]){ // for all isolated vertex u
                    m_ordered_vertex[--whereToInput]= u;    // move u into the output Q
                    vtxDegree[u] = -1;               // degree is -1 means the vertex have been in the ordered queue already
                }
                num_of_local_vertex -= num_iso_vtxs; 
                if(num_of_local_vertex==0) break;
                GroupedVertexDegree[0].clear();      // move u out from the ListOfList
                SearchStartDegree=1;                 // search start from the next degree
            }

            // select a vertex u of the smallest degree
            int u=-1;  // selected vertex u
            for(int d=SearchStartDegree; d<MaxDegreeP1; d++){
                if(!GroupedVertexDegree[d].empty() ){
                    u = GroupedVertexDegree[d].back();    // get one (last) of vetex u
                    GroupedVertexDegree[d].pop_back();    // move u out from the ListOfList
                    int whereToInput = __sync_fetch_and_sub(&outputIndex, 1);   // atomic operation
                    m_ordered_vertex[--whereToInput]=u;             // move u into the output Q
                    vtxDegree[u] = -1;                   // degree is -1 means the vertex have been ordered already.
                    SearchStartDegree = d-1;              // next search start from d-1
                    num_of_local_vertex--;
                }
            }

            if(num_of_local_vertex==0 || u==-1){
                break;  //terminate
            }

            // update the degree of the neighbors of u
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){
                const int v = vtxVal[iv];        // neighbor v
                if(vtxTid[v]==tid){        // we neglect the off thread vertex
                    const int vd = vtxDegree[v];
                    if(vd>0){  // if v in the ListOfList, update its vtxDegree, ListOfLis position and vtxLocation
                        const int vdm1 = vd-1;
                        vtxDegree[v]--;
                        swap( GroupedVertexDegree[vd][vtxLocation[v]], GroupedVertexDegree[vd].back());
                        GroupedVertexDegree[vd].pop_back();
                        vtxLocation[v] = GroupedVertexDegree[vdm1].size();
                        GroupedVertexDegree[vdm1].push_back(v);
                    }
                }
            }
        }// end while   
    }// end omp parallel
    if(ordtime) { end = std::chrono::steady_clock::now();  *ordtime = end - start; }
}

// ============================================================================
// Author: Xin Cheng
// omp parallel natural ordering on the pre-partitioned vertex set
// input   QQ: list of list, [thread][pre partitioned vertex]
// output  QQ: list of list, [thread][ordered pre partitioned vertex]
//
// Nature ordering
// ============================================================================
void SMPGCOrdering::omp_prepartition_ntr_ordering(vector<vector<int>>& QQ,ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        sort(QQ[tid].begin(), QQ[tid].end(), [&](int &a, int&b){return a>=b; });
    }
    if(ordtime) { end = std::chrono::steady_clock::now();  *ordtime = end - start; }
}


// ============================================================================
// Author: Xin Cheng
// omp parallel random ordering on the pre-partitioned vertex set
// input   QQ: list of list, [thread][pre partitioned vertex]
// output  QQ: list of list, [thread][ordered pre partitioned vertex]
//
// Nature ordering
// ============================================================================
void SMPGCOrdering::omp_prepartition_rnd_ordering(vector<vector<int>>& QQ,ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int>& Q = QQ[tid];
        const int QsizeM1 = ((signed)Q.size())-1;
        for(int i=0; i<QsizeM1; i++){
            uniform_int_distribution<int> dist(i,QsizeM1);
            swap(Q[i], Q[dist(m_mt)]);
        }
    }
    if(ordtime) { end = std::chrono::steady_clock::now();  *ordtime = end - start; }
}



// ============================================================================
// Author: Xin Cheng
// omp parallel largest degree first ordering on the pre-partitioned vertex set
// input   QQ: list of list, [thread][pre partitioned vertex]
// output  QQ: list of list, [thread][ordered pre partitioned vertex]
//
// Largest Degree First order
// ============================================================================
void SMPGCOrdering::omp_prepartition_ldf_ordering(vector<vector<int>>& QQ,ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    cout<<"mbd"<<endl;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const vector<int>& vtxPtr = get_CSR_ia();
    #pragma omp parallel
    { 
        const int tid = omp_get_thread_num();
        vector<int>& Q = QQ[tid];
        const int Qsize = Q.size();
        vector<vector<int>> GroupedVertexDegree; //listOflist
    
        // fill in the ListOfList
        int maxDegreeSofar = -1;
        for(int iu=0; iu<Qsize; iu++) {
            const int u = Q[iu];
            const int ud= vtxPtr[u+1]-vtxPtr[u];
            // allocate memory if necessary
            if(ud<maxDegreeSofar){ 
                GroupedVertexDegree.resize(ud+1);
                maxDegreeSofar=ud;
            }
            GroupedVertexDegree[ud].push_back(u);
        }
    
        // fill in the output queue
        Q.clear();
        for(const auto &vx : GroupedVertexDegree)
            Q.insert(Q.end(), vx.begin(), vx.end());
    }

    cout<<"wocao here"<<endl<<flush;
    if(ordtime) { end = std::chrono::steady_clock::now();  *ordtime = end - start; }
}

// ============================================================================
// Author: Xin Cheng
// omp parallel ordering on the pre-partitioned vertex set
// input   QQ: list of list, [thread][pre partitioned vertex]
// output  QQ: list of list, [thread][ordered pre partitioned vertex]
//
// Smallest Degree Last order
// ============================================================================
void SMPGCOrdering::omp_prepartition_sdl_ordering(vector<vector<int>> &QQ,ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const vector<int> & vtxPtr = get_CSR_ia();
    const vector<int> & vtxVal = get_CSR_ja();
    const int N = num_nodes();
    vector<int> vtxDegree(N);            // mapping from vertex to its degree,   1st index 
    vector<int> vtxLocation(N);          // mapping from vertex to its location  2nd index
    vector<int> vtxTid(N);                // mapping from vertex to its thread id.

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        vector<int>& Q = QQ[tid];
        int outputPosition = Q.size();
        const int Qsize = Q.size();
        vector<vector<int>> GroupedVertexDegree;  //list of list. [degree][vertex]
        // init the GVDs, vtxDegree and vtxLocation in parallel 
        int maxDegreeSoFar=-1;
        for(int iu=0; iu<Qsize; iu++){
            const int u = Q[iu];
            const int udeg = -vtxPtr[u] + vtxPtr[u+1];
            vtxDegree[u]   = udeg;
            // allocate memory only if necessary
            if(maxDegreeSoFar<udeg){
                GroupedVertexDegree.resize(udeg+1);
                maxDegreeSoFar=udeg;
            }
            vtxLocation[u] = GroupedVertexDegree[udeg].size();
            GroupedVertexDegree[udeg].push_back(u);
            vtxTid[u]=tid;
        }// end for iu
        
        #pragma omp barrier
        const int MaxDegreeP1 = GroupedVertexDegree.size();
        int SearchStartDegree=0; //skip the unnecessary-checks for some degrees 
        while(outputPosition!=0){ // equals to zero means output Q is full, means all vertex have been ordered.
             
            // if there are isolated vertex, move them to the output Q
            if(SearchStartDegree==0 && (GroupedVertexDegree[0].empty())==false){
                for(auto u : GroupedVertexDegree[0]){ // for all isolated vertex u
                    Q[--outputPosition]= u;    // move u into the output Q
                    vtxDegree[u] = -1;               // degree is -1 means the vertex have been in the ordered queue already
                }
                if(outputPosition==0) break;
                GroupedVertexDegree[0].clear();      // move u out from the ListOfList
                SearchStartDegree=1;                 // search start from the next degree
            }

            // select a vertex u of the smallest degree
            int u=-1;  // selected vertex u
            for(int d=SearchStartDegree; d<MaxDegreeP1; d++){
                if(!GroupedVertexDegree[d].empty() ){
                    u = GroupedVertexDegree[d].back();    // get one (last) of vetex u
                    GroupedVertexDegree[d].pop_back();    // move u out from the ListOfList
                    Q[--outputPosition]=u;             // move u into the output Q
                    vtxDegree[u] = -1;                   // degree is -1 means the vertex have been ordered already.
                    SearchStartDegree = d-1;              // next search start from d-1
                }
            }

            if(u==-1 || outputPosition==0){
                break;  //terminate  
            }

            // update the degree of the neighbors of u
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){
                const int v = vtxVal[iv];        // neighbor v
                if(vtxTid[v]==tid){        // we neglect the off thread vertex
                    const int vd = vtxDegree[v];
                    if(vd>0){  // if v in the ListOfList, update its vtxDegree, ListOfLis position and vtxLocation
                        const int vdm1 = vd-1;
                        vtxDegree[v]--;
                        swap( GroupedVertexDegree[vd][vtxLocation[v]], GroupedVertexDegree[vd].back());
                        GroupedVertexDegree[vd].pop_back();
                        vtxLocation[v] = GroupedVertexDegree[vdm1].size();
                        GroupedVertexDegree[vdm1].push_back(v);
                    }
                }
            }
        }// end of while outputPosition!=0
    }// end omp parallel
    if(ordtime) { end = std::chrono::steady_clock::now();  *ordtime = end - start; }
}




