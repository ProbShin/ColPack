/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/
#include "SMPGCColoring.h"
using namespace std;
using namespace ColPack;

// ============================================================================
// Author: xin cheng
// D1GC parallel speculative iterative approach
// 
// while |W| is not empty
//   TentativeColoring
//   CheckConflicts
//   HandleConlficts
//
// 3P means, handle conflicts will sequential coloring the conflicts
// ============================================================================
int SMPGCColoring::D1GC_OMP_SpeculativeIterative_3P(int nT, int& colors, vector<int>&vtxColors, const int LOCAL_ORDER, const bool bMemOpt, const int nVerbose) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    ChronoDuration tim_Partition{.0}, tim_LO{.0}, tim_TC{.0}, tim_CC{.0}, tim_HC{.0}; // run time for debug
    int num_Cf=0;                                          // number of conflicts for debug
    const int N               = num_nodes();   // number of vertex
    const int BufSize         = max_degree()+1;// forbidden color buffer size 
    const vector<int>& vtxPtr = get_CSR_ia();  // csr format edges
    const vector<int>& vtxVal = get_CSR_ja();  // csr format vtx values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); // ordered vertexs 

    colors=0;                       
    vtxColors.assign(N, -1);

    // allocate memory
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16, -1) ); // for_each(QQ.begin(), QQ.end(), [&](vector<int>&x){x.reserve(N/nT+1+16); } ); //1-odd/even, 16-bus width
    vector<vector<int>> Fs(nT, vector<int>(BufSize+1+16, -1) ); // for_each(Fs.begin(), Fs.end(), [&](vector<int>&x){x.reserve(BufSize+1+16);} );
    vector<int> Qsizes(nT, -1);   // fifteen -1 fill-in to prevent false sharing

    // pre-partition the graph
    {
        auto start = steady_clock::now();
        for(int i=0; i<nT; i++) Qsizes[i]=N/nT; 
        for(int i=0; i<N%nT; i++) Qsizes[i]++;
        vector<int> disps(nT+1, 0); for(int i=1; i<nT+1; i++) disps[i]=disps[i-1]+Qsizes[(i-1)];
        for(int i=0; i<nT; i++)
            QQ[i].assign(const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
        
        auto end = steady_clock::now();
        tim_Partition = end - start;
    }

    // local order, prepartition, parallel
    ordering("OMP", "PREPARTITION", LOCAL_ORDER, &tim_LO, (void *)&QQ);

    // Tentative Coloring
    if(bMemOpt==false){
        do_D1GC_OMP_TentativeColoring(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, Fs, BufSize, tim_TC);
    }
    else{
        do_D1GC_OMP_TentativeColoring_MemOpt(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, Fs, BufSize, tim_TC);
    }

    // Check Conflicts
    do_D1GC_OMP_CheckConflicts(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, tim_CC, num_Cf);

    // Handle Conflicts
    {
        auto start = steady_clock::now();
        vector<int>&F = Fs[0];
        F.assign(BufSize, -1);
        for(int tid=0; tid<nT; tid++){
            const vector<int>& Q=QQ[tid];
            const int Qsize = Qsizes[tid];
            for(int iu=0; iu<Qsize; iu++){
                const int u = Q[iu];
                for(int iv=vtxPtr[u]; iv<vtxPtr[u+1]; iv++){
                    const int v = vtxVal[iv];
                    const int vc = vtxColors[v];
                    if(vc>=0) F[vc]=u;
                }
                int c;
                for(c=0; c<BufSize; c++){
                    if(F[c]!=u)
                        break;
                }
                vtxColors[u]=c;
            }// end for iu
        }//end for tid
        auto end =steady_clock::now();
        tim_HC = end - start;
    }

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@SpI3P";
        if(bMemOpt) ss<<"MemOpt";
        ss<<"("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<")_nT_c_T";
        for(int i=ss.str().size(); i<27; i++)
            ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        auto tim_Total = /* tim_Partition +*/ tim_LO + tim_TC + tim_CC + tim_HC;
        cout<<"\t"<<tim_Total.count();

        if(nVerbose>1){ // show detail run time
            cout<<"\t*tPtn_tLO_tTC_tCC_tHC";
            cout<<"\t"<<tim_Partition.count();
            cout<<"\t"<<tim_LO.count();
            cout<<"\t"<<tim_TC.count();
            cout<<"\t"<<tim_CC.count();
            cout<<"\t"<<tim_HC.count();
        
            if(nVerbose>2){ // show conflicts information
                cout<<"\t*cnfTotal_(cnf1,cnt2,...)";
                cout<<"\t"<<num_Cf<<"\t(";
                for(int i=0; i<nT; i++)
                    cout<<" "<<Qsizes[i];
                cout<<" )";
            } // end if nVerbose>2
        } // end if nVerbose>1
        cout<<"\n";
    }
 
    return true;
}

// ============================================================================
// Author: Xin Cheng
// Distance one graph coloring Parallel speculative iterative approach, Multiple phases
//
// while |W| is not empty
//   TentativeColoring
//   CheckConflicts
//   HandleConlficts
//
// MP means, handle conflicts will leave the uncolored vertex to the next iteration
// ============================================================================
int SMPGCColoring::D1GC_OMP_SpeculativeIterative_MP(int nT, int&colors, vector<int>&vtxColors, const int LOCAL_ORDER, const bool bMemOpt, const int nVerbose){
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    ChronoDuration tim; 
    ChronoDuration tim_Partition, tim_LO;             // run time for debug
    vector<ChronoDuration> vtim_TC, vtim_CC, vtim_HC; // run time for debug
    vector<int> vnum_Cf;                       // num conflicts for debug
    
    const int N               = num_nodes();   //number of vertex
    const int BufSize         = max_degree()+1;//forbidden color buffer size 
    const vector<int>& vtxPtr = get_CSR_ia();  //csr format edges
    const vector<int>& vtxVal = get_CSR_ja();  //csr format vtx values
    const vector<int>& const_ordered_vertex = get_ordered_vertex(); // ordered vertexs 

    colors=0;                       
    vtxColors.assign(N, -1);

    // allocate memory
    vector<vector<int>> QQ(nT); for_each(QQ.begin(), QQ.end(), [&](vector<int>&x){x.reserve(N/nT+1+16); } ); //1-odd/even, 16-bus width
    vector<vector<int>> Fs(nT); for_each(Fs.begin(), Fs.end(), [&](vector<int>&x){x.reserve(BufSize+1+16);} );
    vector<int> Qsizes(nT, 0);   // fifteen -1 fill-in to prevent false sharing

    // pre-partition the graph
    {
        auto start = steady_clock::now();
        for(int i=0; i<nT; i++) Qsizes[i]=N/nT; 
        for(int i=0; i<N%nT; i++) Qsizes[i]++;
        vector<int> disps(nT+1, 0); for(int i=1; i<nT+1; i++) disps[i]=disps[i-1]+Qsizes[(i-1)];
        for(int i=0; i<nT; i++)
            QQ[i].assign(const_ordered_vertex.begin()+disps[i], const_ordered_vertex.begin()+disps[i+1]);
        auto end = steady_clock::now();
        tim_Partition = end - start;
    }

    // local order, prepartition, parallel
    ordering("OMP", "PREPARTITION", LOCAL_ORDER, &tim_LO, (void *)&QQ);

    int n_iter=0;
    do{
        // Tentative Coloring
        vtim_TC.push_back(ChronoDuration(.0));
        if(bMemOpt==false){
            do_D1GC_OMP_TentativeColoring(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, Fs, BufSize, vtim_TC[n_iter]);
        }
        else{
            do_D1GC_OMP_TentativeColoring_MemOpt(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, Fs, BufSize, vtim_TC[n_iter]);
        }

        // Check Conflicts
        vtim_CC.push_back(ChronoDuration(.0)); vnum_Cf.push_back(0);
        do_D1GC_OMP_CheckConflicts(nT, vtxColors, vtxPtr, vtxVal, QQ, Qsizes, vtim_CC[n_iter], vnum_Cf[n_iter]);
        for(int i=0; i<nT; i++)
            vnum_Cf[n_iter] += Qsizes[i];  //for debug
        
        // Handle Conflics
        // Do nothing
    }while(vnum_Cf[n_iter++]!=0);

    colors = calc_num_colors_from_vtx_colors(vtxColors);

    if(nVerbose>0){ // show basic information
        stringstream ss;
        ss<<"@SpIMP";
        if(bMemOpt) ss<<"MemOpt";
        ss<<"("<<Translate_OrderId_To_OrderTag(LOCAL_ORDER)<<")_nT_c_T";
        for(int i=ss.str().size(); i<27; i++) 
            ss<<"_";
        cout<<ss.str();
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration tim_TC{.0}, tim_CC{.0};
        for(auto &x: vtim_TC) tim_TC+=x;
        for(auto &x: vtim_CC) tim_CC+=x;
        auto tim_Total = /* tim_Partition +*/ tim_LO + tim_TC + tim_CC;
        cout<<"\t"<<tim_Total.count();

        if(nVerbose>1){ // show detail run time
            cout<<"\t*tPtn_tLO_tTC_tCC";
            cout<<"\t"<<tim_Partition.count();
            cout<<"\t"<<tim_LO.count();
            cout<<"\t"<<tim_TC.count();
            cout<<"\t"<<tim_CC.count();
        
            if(nVerbose>2){ // show details for each iteration
                cout<<"\n*tTC "<<tim_TC.count()<<"\t=\t";
                for(auto &x: vtim_TC)
                    cout<<" "<<x.count();
                
                cout<<"\n*tCC "<<tim_CC.count()<<"\t=\t";
                for(auto &x: vtim_CC)
                    cout<<" "<<x.count();

                int totalConf=0;
                for(auto x: vnum_Cf)
                    totalConf+=x;
                cout<<"\n*cnfTotal "<<totalConf<<"\t=\t";
                for(auto x: vnum_Cf)
                    cout<<" "<<x;
            } // end if nVerbose>2
        } // end if nVerbose>1
        cout<<"\n";
    }
 
    return true;
}



// ============================================================================
// Author: Xin Cheng
// Phase Tentative Coloring  of D1GC omp
// ----------------------------------------------------------------------------
// tentative coloring the graphs without considering other threads run in parallel
// ============================================================================
void SMPGCColoring::do_D1GC_OMP_TentativeColoring(
        const int nT,                   // number of thread
        vector<int>& vtxColors,         // mapping from vertex to its color
        const vector<int>& vtxPtr,      // CSR format edges
        const vector<int>& vtxVal,      // CSR format values
        const vector<vector<int>> &QQ,  // vertex to color
        const vector<int>& Qsizes,      // local vertex size
        vector<vector<int>>&Fs,         // Forbidden color buffers
        const int BufSize,              // Buffer size
        ChronoDuration &time             // run time
        ){           
    auto start = steady_clock::now();
    #pragma omp parallel 
    {
        const int tid = omp_get_thread_num();
        const vector<int>& Q = QQ[tid];
        const int Qsize = Qsizes[tid];
        vector<int>&F = Fs[tid];
        F.assign(BufSize, -1);
        for(int iu=0; iu<Qsize; iu++){
            const int u = Q[iu];
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++) {
                const int vc = vtxColors[vtxVal[iv]];
                if(vc>=0)
                    F[vc]=u;
            }
            int c;
            for(c=0; c<BufSize; c++){
                if(F[c]!=u)
                    break;
            }
            vtxColors[u]=c;
        }// end for u
    }// end of omp parallel
    auto end = steady_clock::now();
    time = end - start;
}


// ============================================================================
// Author: Xin Cheng
// Phase CC, ( Check Conflicts )
// ----------------------------------------------------------------------------
// check conflicts, if there is, reset its color and move it to a set
// ============================================================================
void SMPGCColoring::do_D1GC_OMP_CheckConflicts(
        const int nT,                // number of thread
        vector<int> &vtxColors,      // mapping vertex to its color
        const vector<int>& vtxPtr,   // CSR format edges
        const vector<int>& vtxVal,   // CSR format values
        vector<vector<int>>& QQ,     // vertex to color
        vector<int>& Qsizes,         // local vertex size
        ChronoDuration & time,       // run time
        int& num_Conflicts                    // number of conlflicts
        ){
    int nCf = 0 ;
    auto start = steady_clock::now();
    #pragma omp parallel reduction(+:nCf)
    {
        const int tid = omp_get_thread_num();
        vector<int>& Q= QQ[tid];
        const int Qsize = Qsizes[tid];
        for(int iu=0; iu<Qsize; iu++){
            const int u  = Q[iu];
            const int uc = vtxColors[u];
            for(int iv=vtxPtr[u]; iv!=vtxPtr[u+1]; iv++){
                const int v = vtxVal[iv];
                if(v>u && vtxColors[v]==uc){
                    vtxColors[u] = -1;
                    Q[nCf++]=u;
                    break;
                }
            }
        }//end for iu
        Qsizes[tid] = nCf;
    }// end omp paralle
    auto end = steady_clock::now();
    time = end - start;
    num_Conflicts = nCf; 
}    
    


