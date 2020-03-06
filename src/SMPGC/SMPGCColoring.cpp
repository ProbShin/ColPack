/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "SMPGCColoring.h"
using namespace std;
using namespace ColPack;


// ============================================================================
// Author: Xin Cheng
// Construction
// ============================================================================
SMPGCColoring::SMPGCColoring(const string& graph_name, const string& fmt, ChronoDuration* iotime, const int glb_order, ChronoDuration *ordtime) 
: SMPGCOrdering(graph_name, fmt, iotime, glb_order, ordtime){
}


// ============================================================================
// Interface
// ============================================================================
int SMPGCColoring::Coloring(int nT,int &total_num_colors, vector<int>& vertex_color, const string& method, const int nVerbose){
    //Method follows the following pattern:
    //
    //"  DISTANCE_ONE_OMP_                                                "
    //"                   <GM3P/GMMP/SERIAL/JP/MTJP>[_<LF/SL/NT/RD/NONE>] "
    //"                   HB[MT]JP_<GM3P/GMMP/SERIAL>[_<LF/SL/NT/RD/NONE>]"
    //"  DISTANCE_TWO_OMP_                                                "
    //"                   <GM3P/GMMP/SERIAL>[_<LF/SL/NT/RD/NONE>]         "
    //
    //For example
    //  DISTANCE_ONE_OMP_GM3P_RD
    //  DISTANCE_ONE_HBMTJP_SERIAL
    //  DISTANCE_TWO_GM3P
    //
    //
    
   return true; 

}

