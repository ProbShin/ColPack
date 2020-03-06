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
//   x seq_subset_sdl      
//     omp_whole_ntr
//     omp_whole_rnd
//     omp_whole_ldf
//     omp_whole_sdf
//     omp_prepartition_ntr
//     omp_prepartition_rnd
//     omp_prepartition_ldf
//     omp_prepartition_sdl


// ============================================================================
// Author: Xin Cheng
// sequential natrual order
// store the ordered vertex into the class varibale "m_ordered_vertex"
// OrderedVertex <- [0,1,2,... N-1]
// ============================================================================
void SMPGCOrdering::seq_whole_ntr_ordering(ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const int N = num_nodes();
    m_ordered_vertex.resize(N);
    for(int i=0; i<N; i++) 
        m_ordered_vertex[i]=i;
    m_ordered_method = ORDER_NATURAL;
    if(ordtime) { end = std::chrono::steady_clock::now(); *ordtime = end-start; }
}

// ============================================================================
// Author: Xin Cheng
// sequential random order
// store the ordered vertex into the class variable "m_ordered_vertex"
// Random is a shuffle to the Natural Ordering
// ============================================================================
void SMPGCOrdering::seq_whole_rnd_ordering(ChronoDuration* ordtime) {
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const int N = num_nodes();
    m_ordered_vertex.resize(N);
    for(int i=0; i<N; i++) m_ordered_vertex[i]=i;
    for(int i=0,iEnd=N-1; i<iEnd; i++){
        uniform_int_distribution<int> dist(i, iEnd); 
        swap(m_ordered_vertex[i], m_ordered_vertex[dist(m_mt)]);
    }
    m_ordered_method = ORDER_RANDOM;
    if(ordtime) { end = std::chrono::steady_clock::now(); *ordtime = end-start; }
}

// ============================================================================
// Author: Xin Cheng
// sequential largest degree first order
// store the ordered vertex into the class variable "m_ordered_vertex"
// Largest Degree First
// ============================================================================
void SMPGCOrdering::seq_whole_ldf_ordering(ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const int N = num_nodes();              
    const vector<int>& vtxPtr = get_CSR_ia(); // csr format edges
    const int MaxDegreeP1 = max_degree()+1;   // max degree plus one
    vector<vector<int>> GroupedVertexDegree(MaxDegreeP1); // ListOfList
    
    // clear up the memory
    m_ordered_vertex.clear();
    m_ordered_method = ORDER_LARGEST_FIRST; 
    
    // fill in the ListOfList
    for(int u=0; u<N; u++){
        GroupedVertexDegree[-vtxPtr[u]+vtxPtr[u+1]].push_back(u);
    }
   
    // Largest Degree First
    for(int d=MaxDegreeP1-1, it=MaxDegreeP1; it!=0; it--, d--){
        m_ordered_vertex.insert(m_ordered_vertex.end(), GroupedVertexDegree[d].begin(), GroupedVertexDegree[d].end());
    }
    if(ordtime) { end = std::chrono::steady_clock::now(); *ordtime = end-start; }
}

// ============================================================================
// Author: Xin Cheng
// sequential smallest degree last order
// store the ordered vertex into the class variable "m_ordered_vertex"
// Smallest Degree Last
// ============================================================================
void SMPGCOrdering::seq_whole_sdl_ordering(ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const int N = num_nodes();
    const vector<int>& vtxPtr = get_CSR_ia();   // csr format edges
    const vector<int>& vtxVal = get_CSR_ja();   // csr format value
    const int MaxDegreeP1 = max_degree()+1;     // max degree plus 1
    vector<vector<int>> GroupedVertexDegree(MaxDegreeP1); //list of list. [list of degree][list of vertex]
    vector<int> vtxDegree(N);            // mapping from vertex to its degree,   1st index 
    vector<int> vtxLocation(N);          // mapping from vertex to its location  2nd index
    
    m_ordered_vertex.assign(N, -1);
    m_ordered_method=ORDER_SMALLEST_LAST;

    int outputIndex = N;                 // the postion of the last written into the ordered_queue
    
    // init the GroupedVertexDegree, vtxDegree, vtxLocation 
    for(int u=0; u<N; u++){
        const int udeg = -vtxPtr[u] + vtxPtr[u+1];
        vtxDegree[u]   = udeg;
        vtxLocation[u] = GroupedVertexDegree[udeg].size(); 
        GroupedVertexDegree[udeg].push_back(u);
    }

    int SearchStartDegree=0; //skip the unnecessary-checks for some degrees 
    while(outputIndex!=0){ // equals to zero means output Q is full, means all vertex have been ordered.
        // if there are isolated vertex, move them to the output Q
        if(SearchStartDegree==0 && GroupedVertexDegree[0].empty()==false){
            for(auto u : GroupedVertexDegree[0]){ // for all isolated vertex u
                m_ordered_vertex[--outputIndex]= u;    // move u into the output Q
                vtxDegree[u] = -1;               // degree is -1 means the vertex have been in the ordered queue already
            }
            GroupedVertexDegree[0].clear();      // move u out from the ListOfList
            SearchStartDegree=1;                 // search start from the next degree
        }

        // select vertex u of the smallest degree
        int u=-1;  // selected vertex u
        for(int d=SearchStartDegree; d<MaxDegreeP1; d++){
            if(!GroupedVertexDegree[d].empty() ){
                u = GroupedVertexDegree[d].back();    // get one (last) of vetex u
                GroupedVertexDegree[d].pop_back();    // move u out from the ListOfList
                m_ordered_vertex[--outputIndex]=u;             // move u into the output Q
                SearchStartDegree = d-1;              // next search start from d-1
            }
        }
        
        if(u==-1 || outputIndex==0){
            break; //terminate
        }

        // update the degree of the neighbors of vetex u
        for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){
            const int v = vtxVal[iv];        // neighbor v
            if(vtxDegree[v] == -1) 
                continue;                    // vertex v have been moved into the ordered queue already
            const int vd = vtxDegree[v]--;   // update vtxDegree
            const int vdm1 = vd-1;
            // update ListOfList and vtxLocation
            swap( GroupedVertexDegree[vd][vtxLocation[v]], GroupedVertexDegree[vd].back() );
            GroupedVertexDegree[vd].pop_back();
            vtxLocation[v]=GroupedVertexDegree[vdm1].size();
            GroupedVertexDegree[vdm1].push_back(v);
        }
    }// end while ordered queue is not full
    if(ordtime) { end = std::chrono::steady_clock::now(); *ordtime = end-start; }
}


// ============================================================================
// Author: Xin Cheng
// sequential natural order on a subset of vertex Q
// input   Q: contains vertex to order
// output  Q: contains ordered vertex
//
// Nature ordering
// ============================================================================
void SMPGCOrdering::seq_subset_ntr_ordering(vector<int>& Q,ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    sort(Q.begin(), Q.end(), [&](int &a, int&b){return a>=b; });
    if(ordtime) { end = std::chrono::steady_clock::now(); *ordtime = end-start; }
}


// ============================================================================
// Author: Xin Cheng
// sequential random order on a subset of vertex Q
// input   Q: contains vertex to order
// output  Q: contains ordered vertex
//
// Nature ordering
// ============================================================================
void SMPGCOrdering::seq_subset_rnd_ordering(vector<int>& Q,ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const int QsizeM1 = ((signed)Q.size())-1;
    for(int i=0; i<QsizeM1; i++){
        uniform_int_distribution<int> dist(i,QsizeM1);
        swap(Q[i], Q[dist(m_mt)]);
    }
    if(ordtime) { end = std::chrono::steady_clock::now(); *ordtime = end-start; }
}



// ============================================================================
// Author: Xin Cheng
// sequential largest degree first order on a subset of vertex Q
// input   Q: contains vertex to order
// output  Q: contains ordered vertex
//
// Largest Degree First order
// ============================================================================
void SMPGCOrdering::seq_subset_ldf_ordering(vector<int>& Q,ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const int Qsize = Q.size();
    const vector<int>& vtxPtr = get_CSR_ia();
    vector<vector<int>> GroupedVertexDegree; //listOflist
    
    // fill in the ListOfList
    int maxDegreeSofar = -1;
    for(int iu=0; iu<Qsize; iu++) {
        const int u = Q[iu];
        const int ud= -vtxPtr[u]+vtxPtr[u+1];
        // allocate memory if necessary
        if(ud<maxDegreeSofar){ 
            GroupedVertexDegree.resize(ud+1);
            maxDegreeSofar=ud;
        }
        GroupedVertexDegree[ud].push_back(u);
    }
    
    // fill in the output queue
    Q.clear();
    for(auto &vx : GroupedVertexDegree)
        Q.insert(Q.end(), vx.begin(), vx.end());
    if(ordtime) { end = std::chrono::steady_clock::now(); *ordtime = end-start; }
}


/// ============================================================================
// Author: Xin Cheng
// sequential smallest last ordering on a subset of vertex Q
// input   Q: contains vertex to order
// output  Q: contains ordered vertex
//
// Smallest Degree Last Order
// ============================================================================
void SMPGCOrdering::seq_subset_sdl_ordering(vector<int> &Q,ChronoDuration* ordtime){
    std::chrono::time_point<std::chrono::steady_clock> start,end;
    if(ordtime) { start = std::chrono::steady_clock::now(); }
    const vector<int> & vtxPtr = get_CSR_ia();
    const vector<int> & vtxVal = get_CSR_ja();

    
    unordered_map<int, int> vtxDegree;      // compared with vector<int>, trade runtime for memory
    unordered_map<int, int> vtxLocation;    // compared with vector<int>, trade runtime for memory
    unordered_set<int> vtxTid(Q.begin(), Q.end());         // compared with vector<int>, trade runtime for memory

    //const int N = num_nodes();
    //vector<int> vtxDegree(N);        // trade memory for the runtime 
    //vector<int> vtxLocation(N);      // trade memory for the runtime
    //vector<int> vtxTid(N);           // trade memory for the runtime

    const int Qsize = Q.size();
    vector<vector<int>> GroupedVertexDegree;  //list of list. [degree][vertex]
    // init the GVDs, vtxDegree and vtxLocation in parallel 
    int maxDegreeSoFar=-1;
    for(int iu=0; iu<Qsize; iu++){
        const int u = Q[iu];
        const int udeg = vtxPtr[u+1] - vtxPtr[u];
        vtxDegree[u]   = udeg;
        // allocate memory only if necessary
        if(maxDegreeSoFar<udeg){
            GroupedVertexDegree.resize(udeg+1);
            maxDegreeSoFar=udeg;
        }
        vtxLocation[u] = GroupedVertexDegree[udeg].size();
        GroupedVertexDegree[udeg].push_back(u);
    }// end for iu

    int outputPosition = Q.size();
    const int MaxDegreeP1 = GroupedVertexDegree.size();
    int SearchStartDegree=0; //skip the unnecessary-checks for some degrees 
    while(outputPosition!=0){ // equals to zero means output Q is full, means all vertex have been ordered.
        // if there are isolated vertex, move them to the output Q
        if(SearchStartDegree==0 && GroupedVertexDegree[0].empty()==false){
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
        
        if(u==-1 || outputPosition==0) {
            break;
        }

        // update the degree of the neighbors of u
        for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){
            const int v = vtxVal[iv];        // neighbor v
            if(vtxTid.count(v)){        // we neglect the off thread vertex
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
    if(ordtime) { end = std::chrono::steady_clock::now(); *ordtime = end-start; }
}







