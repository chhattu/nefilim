#ifndef _PDBCOORDS_
#define _PDBCOORDS_

#include <iostream>
#include <fstream>
#include "headers.h"
#include "seq.h"
#include "vector.h"
using namespace std;

class pdbcoords {

private:

  seq<vector> allCoords;
  seq<vector> allCoords1;
  seq<vector> allCoords2;
  seq<vector> allCoords1CA; // Alpha carbons in helix 1
  seq<vector> allCoords1C; // Beta carbons in helix 1
  seq<vector> allCoords1O; // Oxygen in helix 1
  seq<vector> allCoords1N; // Nitrogen in helix 1
  seq< seq<vector> > sideChains1; 
  seq<string> aminoacids1;
  seq<vector> allCoords2CA; // Alpha carbons in helix 2
  seq<vector> allCoords2C; // Beta carbons in helix 2
  seq<vector> allCoords2O; // Oxygen in helix 2
  seq<vector> allCoords2N; // Nitrogen in helix 2
  seq< seq<vector> > sideChains2; 
  seq<string> aminoacids2;
  seq< seq<int> > hindList;  // The indices of each alpha helix (start/end)
  seq<string> atomIDs;  	 // The chain IDs for each atom
  seq<string> helixIDs;  	 // The chain IDs for each helix
  seq<int> atomInds;  	 // The indices of each atom
  seq<int> aaInds1;  	 // The indices of each amino acid in helix 1
  seq<int> aaInds2;  	 // The indices of each amino acid in helix 2
 
public:
 
  pdbcoords() {
	cerr << "You gotta gimme a file there, dumbass!" << endl; 	
  }  	

  pdbcoords(string file) {			// constructor
    findCoords(file);
  }

  pdbcoords(string file,int hel1,int hel2) {			// constructor
//    if(showCAonly)
//    	findCoordsonlyCA(file,hel1,hel2);
//    else
    	findCoords(file,hel1,hel2);
  }

  ~pdbcoords() {			// destructor
  }
  
  // Function to return all the alpha carbon atoms in the protein
  void findCoords(string file);
  // Function to return the alpha carbons of only two helices
  void findCoordsonlyCA(string file,int hel1,int hel2);
  // Function to return all atoms only two helices
  void findCoords(string file,int hel1,int hel2);
  
  seq<vector> getAllCoords() {return allCoords;};
  seq<vector> getCoords1() {return allCoords1;};
  seq<vector> getCoords2() {return allCoords2;};
  seq<vector> getCoords1CA() {return allCoords1CA;};
  seq<vector> getCoords1C() {return allCoords1C;};
  seq<vector> getCoords1N() {return allCoords1N;};
  seq<vector> getCoords1O() {return allCoords1O;};
  seq<string> getAminoAcids1() {return aminoacids1;};
  seq< seq<vector> > getSideChains1() {return sideChains1;};
  seq<vector> getCoords2CA() {return allCoords2CA;};
  seq<vector> getCoords2C() {return allCoords2C;};
  seq<vector> getCoords2N() {return allCoords2N;};
  seq<vector> getCoords2O() {return allCoords2O;};
  seq<string> getAminoAcids2() {return aminoacids2;};
  seq< seq<vector> > getSideChains2() {return sideChains2;};
  seq< seq<int> > getHindList() {return hindList;};
  seq<string> getAtomIDs() {return atomIDs;};  
  seq<string> getHelixIDs() {return helixIDs;}; 
  seq<int> getAtomInds() {return atomInds;}; 
  seq<int> getAAInds1() {return aaInds1;}; 
  seq<int> getAAInds2() {return aaInds2;}; 
};

#endif //_PDBCOORDS_
