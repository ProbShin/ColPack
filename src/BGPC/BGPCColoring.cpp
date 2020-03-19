/******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include "BGPCColoring.h"
#include <chrono> //c++11 system time
#include <random> //c++11 random
#include <unordered_map>
using namespace std;
using namespace ColPack;


// ============================================================================
// Interface
// ============================================================================
int BGPCColoring::Coloring(const int side, int nT, const string& method, int& nColor, vector<int>& vColors){
    m_method = method;
    unordered_map<string,int> ALGs{
        { "PD2_OMP_SERIAL", 0},
        { "PD2_SERIAL",     0},
        { "PD2_OMP_GM3P",   1 },
        { "PD2_OMP_GM3P_LF",11 },
        { "PD2_OMP_GMMP", 2},
        { "PD2_OMP_GMMP_LF", 21},
        { "PD2_OMP_GM3P_ATOMIC", 3 },
        { "PD2_OMP_GMMP_ATOMIC", 4},
        { "PD2_OMP_GM3P_BIT", 5},
        { "PD2_OMP_GMMP_BIT", 6},
        { "PD2_OMP_GMMP_NET_N1N1", 7 },
        { "PD2_OMP_GMMP_NET_N2N1", 8 }
    };

    const int alg_id = ALGs[ method ];
    switch(alg_id){
        case 0: 
            //PD2_serial(side, nColor, vColors);
            break;
        case 1: 
            //PD2_OMP_GM3P(side, nT, nColor, vColors);
            break;
        case 2: 
            //PD2_OMP_GMMP(side, nT, nColor, vColors);
            break;
        case 3: 
            
            break;
        case 4: break;
        case 5: break;
        case 6: break;
        default:
            cout<<"Unknown method "<<method<<endl;
            exit(1);
    }
    return true;
}


// ============================================================================
// Construction
// ============================================================================
BGPCColoring::BGPCColoring(const string& graph_name)
: BGPCOrdering(graph_name, FORMAT_MM, nullptr, L, ORDER_NATURAL, nullptr)/*, m_total_num_colors(0)*/ {
}

// ============================================================================
// Construction
// ============================================================================
BGPCColoring::BGPCColoring(const string& graph_name, const string& fmt, ChronoDuration* iotime,  const int side, const int glb_order, ChronoDuration *ordtime) 
: BGPCOrdering(graph_name, fmt, iotime, side, glb_order, ordtime) /*, m_total_num_colors(0)*/ {
}






