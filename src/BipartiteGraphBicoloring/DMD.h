/*******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#ifndef DMDECOMPOSITION_H
#define DMDECOMPOSITION_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
using namespace std;

class DMDecomposition{
   public: DMDecomposition(const string & DMDFileName);
public: ~DMDecomposition() {};

public:
struct {
    vector<int> E,U,O;
} L,R;

};




#endif
