/*******************************************************************************
    This file is part of ColPack, which is under its License protection.
    You should have received a copy of the License. If not, see 
    <https://github.com/CSCsw/ColPack>
*******************************************************************************/

#include<DMD.h>
using namespace std;

// ============================================================================
// construction from the file
// Author: Xin Cheng
// ============================================================================
DMDecomposition::DMDecomposition(
        string const &DMDFileName
        ) {

    // vector<int>& Lside_vtx_E = L.E;
    // vector<int>& Lside_vtx_O = L.O;
    // vector<int>& Lside_vtx_U = L.U;

    // vector<int>& Rside_vtx_E = R.E;
    // vector<int>& Rside_vtx_O = R.O;
    // vector<int>& Rside_vtx_U = R.U;

    L.E.clear(); L.O.clear(); L.U.clear();
    R.E.clear(); R.O.clear(); R.U.clear();

    if(DMDFileName.empty()){
        fprintf(stderr, "Error, the DMD file's name is empty.\n"); exit(1);
    }

    ifstream in(DMDFileName.c_str()) ;
    if(!in.is_open()){
        fprintf(stderr, "Error, the DMD file '%s' cannot be open.\n", DMDFileName.c_str()); exit(1);
    }


    {
        string line;
        while(getline(in, line)){ // first line contains LE,LO,LU and RE,RO,RU information
            if(line.empty() || line[0]=='%') continue;
            break;
        }
        stringstream ss(line);
        int LE_entries=0, LO_entries=0, LU_entries=0;
        int RE_entries=0, RO_entries=0, RU_entries=0;

        if(!(ss >> LE_entries >> LO_entries >> LU_entries >> RE_entries >>RO_entries >> RU_entries)) {
            fprintf(stderr, "Err, the DMD file's format wrong. read a line '%s'\n", line.c_str()); exit(1);
        }

        L.E.reserve(LE_entries);
        L.O.reserve(LO_entries);
        L.U.reserve(LU_entries);
        R.E.reserve(RE_entries);
        R.O.reserve(RO_entries);
        R.U.reserve(RU_entries);


        int vid=0;
        int n_read_LE_entries=0; 
        int n_read_LO_entries=0; 
        int n_read_LU_entries=0; 
        int n_read_RE_entries=0; 
        int n_read_RO_entries=0; 
        int n_read_RU_entries=0; 
        //const int n_expect_entries = LE_entries + LO_entries + LU_entries + RE_entries + RO_entries + RU_entries;


        // read LE vertices
        {
            line="";
            while(getline(in,line)){
                if(line.empty() || line[0]=='%') continue;
                break;
            }
            ss.clear(); ss.str(line);
            string tag; ss>>tag; if(tag!="LE") {fprintf(stderr,"Err, expect read line start with 'LE', but acutal read '%s'\n",tag.c_str()); exit(1);}
            while(!ss.eof() && ss>>vid){
                L.E.push_back(vid-1);
                n_read_LE_entries++;
            }
            if(n_read_LE_entries != LE_entries) {
                fprintf(stderr, "Err, expect read LE of size %d, but actual read %d entries\n", LE_entries, n_read_LE_entries); exit(1);
            }
        }

        // read LO vertices
        {
            line="";
            while(getline(in,line)){
                if(line.empty() || line[0]=='%') continue;
                break;
            }
            ss.clear(); ss.str(line);
            string tag; ss>>tag; if(tag!="LO") {fprintf(stderr,"Err, expect read line start with 'LO', but acutal read '%s'\n",tag.c_str()); exit(1); }
            while(ss>>vid){
                L.O.push_back(vid-1);
                n_read_LO_entries++;
            }
            if(n_read_LO_entries != LO_entries) {
                fprintf(stderr, "Err, expect read LO of size %d, but actual read %d entries\n", LO_entries, n_read_LO_entries); exit(1);
            }
        }

        // read LU vertices
        {
            line="";
            while(getline(in,line)){
                if(line.empty() || line[0]=='%') continue;
                break;
            }
            ss.clear(); ss.str(line);
            string tag; ss>>tag; if(tag!="LU") {fprintf(stderr,"Err, expect read line start with 'LU', but acutal read '%s'\n",tag.c_str()); exit(1); }
            while(ss>>vid){
                L.U.push_back(vid-1);
                n_read_LU_entries++;
            }
            if(n_read_LU_entries != LU_entries) {
                fprintf(stderr, "Err, expect read LU of size %d, but actual read %d entries\n", LU_entries, n_read_LU_entries); exit(1);
            }
        }

        // read RE vertices
        {
            line="";
            while(getline(in,line)){
                if(line.empty() || line[0]=='%') continue;
                break;
            }
            ss.clear(); ss.str(line);
            string tag; ss>>tag; if(tag!="RE") {fprintf(stderr,"Err, expect read line start with 'RE', but acutal read '%s'\n",tag.c_str()); exit(1); }
            while(ss>>vid){
                R.E.push_back(vid-1);
                n_read_RE_entries++;
            }
            if(n_read_RE_entries != RE_entries) {
                fprintf(stderr, "Err, expect read RE of size %d, but actual read %d entries\n", RE_entries, n_read_RE_entries); exit(1);
            }
        }

        // read RO vertices
        {
            line="";
            while(getline(in,line)){
                if(line.empty() || line[0]=='%') continue;
                break;
            }
            ss.clear(); ss.str(line);
            string tag; ss>>tag; if(tag!="RO") {fprintf(stderr,"Err, expect read line start with 'RO', but acutal read '%s'\n",tag.c_str()); exit(1); }
            while(ss>>vid){
                R.O.push_back(vid-1);
                n_read_RO_entries++;
            }
            if(n_read_RO_entries != RO_entries) {
                fprintf(stderr, "Err, expect read RO of size %d, but actual read %d entries\n", RO_entries, n_read_RO_entries); exit(1);
            }
        }

        // read RU vertices
        {
            line="";
            while(getline(in,line)) {
                if(line.empty() || line[0]=='%') continue;
                break;
            }
            ss.clear(); ss.str(line);
            string tag; ss>>tag; if(tag!="RU") {fprintf(stderr,"Err, expect read line start with 'RU', but acutal read '%s'\n",tag.c_str()); exit(1); }
            while(ss>>vid){
                R.U.push_back(vid-1);
                n_read_RU_entries++;
            }
            if(n_read_RU_entries != RU_entries) {
                fprintf(stderr, "Err, expect read RU of size %d, but actual read %d entries\n", RU_entries, n_read_RU_entries); exit(1);
            }
        }

        in.close();

    }
} // end of funtion


    //        while (n_read_entries < L_entries && getline (in,line) ) {
    //            ss.clear(); ss.str(line);
    //            ss>>vid;  // v is one-based
    //            m_vi_IncludedLeftVertices[vid-1]=0;
    //            n_read_entries++;
    //        }
    //        // read right vertices
    //        while ( getline (in,line) ) {
    //            ss.clear(); ss.str(line);
    //            ss>>vid;  // v is one-based
    //            m_vi_IncludedRightVertices[vid-1]=0;
    //            n_read_entries++;
    //        }
    //
    //        if(n_read_entries != n_expect_entries) {
    //            fprintf(stderr, "Error! The DMD file '%s' expect read %d+%d = %d entries, but acutual read %d entry. \n", DMDFileName.c_str(), L_entries, R_entries, n_expect_entries, n_read_entries );
    //           exit(1);
    //       }

    //       // update CoverLeft/RightVertices
    //       for(int i=0; i<M; i++) {
    //           if(m_vi_IncludedLeftVertices[i]==1)
    //               m_vi_CoveredLeftVertices.push_back(i);
    //       }
    //       for(int i=0; i<N; i++) {
    //           if(m_vi_IncludedRightVertices[i]==1)
    //               m_vi_CoveredRightVertices.push_back(i);
    //       }

    //   }




