/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "BGPCColoring.h"
using namespace std;
using namespace ColPack;
using namespace chrono;

// ============================================================================
// author: xin cheng
// Bipartite Graph Partial Coloring parallel implementation
// ============================================================================
int BGPCColoring::BGPC_OMP_IndependentSetIterative(
        const int side, 
        int nT, 
        int &colors, 
        vector<int>&vtxColors,
        int const LOCAL_ORDER, 
        int const WEIGHT_TYPE , 
        const int nVerbose) {
    if(nT<=0) { printf("Warning, number of threads changed from %d to 1\n",nT); nT=1; }
    omp_set_num_threads(nT);

    ChronoDuration tim_PPT{0};
    ChronoDuration tim_LO{0};
    ChronoDuration tim_WT{0};
    vector<ChronoDuration>  vtim_ISI;
    vector<int> vnum_ISI;
    
    const int          N             = (side==BGPC::L)?(GetLeftVertexCount()):(GetRightVertexCount()); 
    const vector<int>& srcPtr        = (side==BGPC::L)?(GetLeftVertices()   ):(GetRightVertices()   );
    const vector<int>& dstPtr        = (side==BGPC::L)?(GetRightVertices()  ):(GetLeftVertices()    );
    const vector<int>& vtxVal        = GetEdges();
    const int          srcMaxDegree  = (side==BGPC::L)?GetMaximumLeftVertexDegree():GetMaximumRightVertexDegree();
    const int          dstMaxDegree  = (side==BGPC::L)?GetMaximumRightVertexDegree():GetMaximumLeftVertexDegree();
    const int          DDp1          = dstMaxDegree * srcMaxDegree +1;
    const int          BufSize       = DDp1<0?N:min(N, DDp1);
    
    vector<int>const& const_queue_A  = get_ordered_queue_A_const();
    //vector<int>const& const_queue_B;
    vector<uint64_t> Weights;
    vector<vector<int>> IndSets;
    vector<int> IndSetSizes;
    
    colors=0;                       
    vtxColors.assign(N, -1);
    
    // allocate memory +16 for prevent false share on some 64 bit machines
    vector<vector<int>> ForbiddenArrays(nT, vector<int>(BufSize+16,-1)); 
    vector<vector<int>> QQ(nT, vector<int>(N/nT+1+16, -1)); 
    vector<int>         Qsizes(nT, 0);
    Weights.assign(N, 0);
    IndSets.assign(nT, vector<int>(N/nT+1+16));
    IndSetSizes.assign(nT,0);
    
    // pre-partition
    init_QQ_and_Qsizes(QQ, Qsizes, nT, N, const_queue_A, &tim_PPT); 

    // generate the weight
    do_Generate_Weights(Weights, N, WEIGHT_TYPE, &srcPtr[0], &dstPtr[0], &vtxVal[0], &tim_WT);

    // local ordering   
    Ordering_PrePartition_OMP(side, QQ, Qsizes, LOCAL_ORDER, &tim_LO);

    
    int niter=0;
    int num_uncolored_vertex = N;
    while(num_uncolored_vertex!=0){
       
        auto start = steady_clock::now();
        num_uncolored_vertex = 0;
        #pragma omp parallel reduction(+:num_uncolored_vertex)
        {
            int const tid = omp_get_thread_num();
            do_BGPC_Seq_FindIndSet(IndSets[tid], IndSetSizes[tid], QQ[tid], Qsizes[tid], vtxColors, srcPtr, dstPtr, vtxVal, Weights);
            num_uncolored_vertex = Qsizes[tid];
            #pragma omp barrier
            do_BGPC_Seq_Greedy(vtxColors, IndSets[tid], IndSetSizes[tid], srcPtr, dstPtr, vtxVal, ForbiddenArrays[tid], BufSize);
        }
        auto end = steady_clock::now();
        vtim_ISI.push_back(end - start);
        vnum_ISI.push_back(num_uncolored_vertex);
        //if(niter%10==0) cout<<niter<<"\t"<<vnum_ISI.back()<<endl;
        niter      += 1;
    } //end while

    colors = calc_num_colors_from_vtx_colors(vtxColors);
    if(nVerbose>0){
        stringstream ss;
        ss<<"@BGPCISI(";
        ss<<Translate_OrderId_To_OrderTag(LOCAL_ORDER);
        if( WEIGHT_TYPE == WEIGHT_RANDOM ) ss<<",RndW";
        else if( WEIGHT_TYPE == WEIGHT_DEGREE_D1 ) ss<<",D1W";
        else if( WEIGHT_TYPE == WEIGHT_DEGREE_D2 ) ss<<",D2W";
        else ss<<",??W";
        for(int i=ss.str().size(); i<30; i++)
            ss<<"_";
        cout<<ss.str(); ss.str("");
        cout<<"\t#nT_c_T";
        cout<<"\t"<<nT;
        cout<<"\t"<<colors;
        ChronoDuration T{0};
        //T += tim_PPT;
        //T += tim_LO;
        //T += tim_WT;
        for(auto const x : vtim_ISI) T+=x;
        cout<<"\t"<<T.count();
        if(nVerbose>1) {
            //TODO
        }
        cout<<endl;
    }
    return true;   
}



// ============================================================================
// Author: Xin Cheng
// find independet set
// ============================================================================
void BGPCColoring::do_BGPC_Seq_FindIndSet(
        vector<int>& IndSet,
        int& IndSetSize,
        vector<int>& Q,
        int& Qsize,
        vector<int>const& vtxColors,
        vector<int>const& srcPtr,
        vector<int>const& dstPtr,
        vector<int>const& vtxVal,
        vector<uint64_t>const& Weights
        ) {
    int num_leftover=0;
    for(int iu=0; iu<Qsize; iu++) {
        int const u = Q[iu];
        auto const uw = Weights[u];
        bool b_uisdomain = true;
        for(int iv=srcPtr[u]; iv!=srcPtr[u+1]; iv++) {
            if(!b_uisdomain) break;
            int const v = vtxVal[iv];
            for(int iw=dstPtr[v]; iw!=dstPtr[v+1]; iw++) {
                int const w  = vtxVal[iw];
                if(u==w || vtxColors[w]>=0 ) continue;
                auto const ww = Weights[w];
                if(uw > ww) continue;
                //if(uw == ww && u>w) continue;
                b_uisdomain = false;
                break;
            }
        }
        if(b_uisdomain)
            IndSet[IndSetSize++] = u;
        else
            Q[num_leftover++] = u;
    }// end iu
    Qsize = num_leftover;
}






