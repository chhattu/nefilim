#include "pdbcoords.h"
#include "headers.h"

// This function gets all of the alpha carbons in the protein 
void pdbcoords::findCoords(string pdbfile) {
//	string pdbfile;
	string temp;
	string line;
	string name;
	string chainID;
	string atom;
	double x;
	double y;
	double z;
	char * pEnd;
	double avgX = 0.0;
  double avgY = 0.0;
  double avgZ = 0.0;
  int minIndex1;
  int maxIndex1;
  int aacount;
  string chainid; 
	seq<int> tempSeq;
	  
	pdbfile = "./pdbfiles/" + pdbfile + ".pdb";
  ifstream iFile(pdbfile.c_str());

  if (!iFile) {
      cout << "Error opening input file" << endl;
  }

  while (getline(iFile,line)) {
    name = line.substr(0,6);
    chainID = line.substr(21,1);
    atom = line.substr(12,3);
    
    // If there are several models in the file, we would only
    // like to display the first one
    if (name=="ENDMDL") break;
    	
    if (name=="HELIX ") {
    	string minIndex = line.substr(22,3);
    	string maxIndex = line.substr(34,3);
    	minIndex1 = (int)strtod(minIndex.c_str(),&pEnd);
    	maxIndex1 = (int)strtod(maxIndex.c_str(),&pEnd);
    	tempSeq.clear();
      tempSeq.add(minIndex1);
      tempSeq.add(maxIndex1);
      hindList.add(tempSeq);
      chainid = line.substr(19,1);
      helixIDs.add(chainid);
    }
    	
    if ((name=="ATOM  ")&&(atom==" CA")) {
    	string xstr = line.substr(30,8);
    	string ystr = line.substr(38,8);
    	string zstr = line.substr(46,8);
    	x = strtod(xstr.c_str(),&pEnd);
    	y = strtod(ystr.c_str(),&pEnd);
    	z = strtod(zstr.c_str(),&pEnd);
    	avgX+=x; avgY+=y; avgZ+=z;
    	allCoords.add( vector(x,y,z));
    	chainid = line.substr(21,1);
      atomIDs.add(chainid);
      string aacountstr = line.substr(22,4);
      aacount = (int)strtod(aacountstr.c_str(),&pEnd);
      atomInds.add(aacount);
    }
  }
  // Let's translate the coordinates so the center of rotation in the 
  // window is the center of mass of the protein.
  avgX = avgX/allCoords.size();
  avgY = avgY/allCoords.size();
  avgZ = avgZ/allCoords.size();
  for (int i=0;i<allCoords.size();i++) {
  	allCoords[i].x -= avgX;
  	allCoords[i].y -= avgY;
  	allCoords[i].z -= avgZ;
  } 
}

// This function gets all of the alpha carbons in the two helices 
void pdbcoords::findCoordsonlyCA(string pdbfile, int hel1, int hel2) {
//	string pdbfile;
	string temp;
	string line;
	string name;
	string chainID;
	string atom;
	double x;
	double y;
	double z;
	char * pEnd;
	int count;
	int aacount;
    int minIndex1;
    int maxIndex1;
    int minIndex2;
    int maxIndex2;
    string chain1id;
    string chain2id;
    double totX1 = 0.0;
    double totY1 = 0.0;
    double totZ1 = 0.0;
    double totX2 = 0.0;
    double totY2 = 0.0;
    double totZ2 = 0.0;
    
	  
	pdbfile = "./pdbfiles/" + pdbfile + ".pdb";
//	pdbfile = "C:/Documents and Settings/robert/My Documents/bioinformatics/tony/Files/" + pdbfile + ".pdb";
    ifstream iFile(pdbfile.c_str());

    if (! iFile)
    {
        cout << "Error opening input file" << endl;
    }

    while (getline(iFile,line))
    {
        name = line.substr(0,6);
        chainID = line.substr(21,1);
        atom = line.substr(12,3);
        string countstr = line.substr(7,4);
        count = (int)strtod(countstr.c_str(),&pEnd);
        string aacountstr = line.substr(22,4);
        aacount = (int)strtod(aacountstr.c_str(),&pEnd);
        
        // If there are several models in the file, we would only
        // like to display the first one
        if (name=="ENDMDL")
        	break;
        
        // Find the max & min indices for each helix
        if ((name=="HELIX ")&&(count==hel1))
        {
        	chain1id = line.substr(19,1);
        	string minIndex = line.substr(22,3);
        	string maxIndex = line.substr(34,3);
        	minIndex1 = (int)strtod(minIndex.c_str(),&pEnd);
        	maxIndex1 = (int)strtod(maxIndex.c_str(),&pEnd);
        }
        if ((name=="HELIX ")&&(count==hel2))
        {
        	chain2id = line.substr(19,1);
        	string minIndex = line.substr(22,3);
        	string maxIndex = line.substr(34,3);
        	minIndex2 = (int)strtod(minIndex.c_str(),&pEnd);
        	maxIndex2 = (int)strtod(maxIndex.c_str(),&pEnd);
        }
        
        // Find the coordinates for each helix
        if ((name=="ATOM  ")&&(aacount>=minIndex1)&&(aacount<=maxIndex1)&&(chainID==chain1id)&&(atom==" CA"))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	totX1+=x;
        	totY1+=y;
        	totZ1+=z;
        	allCoords1.add( vector(x,y,z));
        	
        }
        
        if ((name=="ATOM  ")&&(aacount>=minIndex2)&&(aacount<=maxIndex2)&&(chainID==chain2id)&&(atom==" CA"))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	totX2+=x;
        	totY2+=y;
        	totZ2+=z;
        	allCoords2.add( vector(x,y,z));
        	
        }
        
     }
     
     // Let's translate the coordinates so the center of rotation in the 
     // window is the center of mass of the two helices.
     totX1 = totX1/allCoords1.size();
     totY1 = totY1/allCoords1.size();
     totZ1 = totZ1/allCoords1.size();
     totX2 = totX2/allCoords2.size();
     totY2 = totY2/allCoords2.size();
     totZ2 = totZ2/allCoords2.size();
     double avgX = (totX1+totX2)/2;
     double avgY = (totY1+totY2)/2;
     double avgZ = (totZ1+totZ2)/2;
     for (int i=0;i<allCoords1.size();i++)
     {
     	allCoords1[i].x -= avgX;
     	allCoords1[i].y -= avgY;
     	allCoords1[i].z -= avgZ;
     } 
     for (int i=0;i<allCoords2.size();i++)
     {
     	allCoords2[i].x -= avgX;
     	allCoords2[i].y -= avgY;
     	allCoords2[i].z -= avgZ;
     } 
     
}

// This function gets all of the atoms in the two helices 
void pdbcoords::findCoords(string pdbfile, int hel1, int hel2)
{
//	string pdbfile;
	string temp;
	string line;
	string name;
	string chainID;
	string atom;
	double x;
	double y;
	double z;
	char * pEnd;
	int count;
	int aacount;
    int minIndex1;
    int maxIndex1;
    int minIndex2;
    int maxIndex2;
    string chain1id;
    string chain2id;
    double totX1 = 0.0;
    double totY1 = 0.0;
    double totZ1 = 0.0;
    double totX2 = 0.0;
    double totY2 = 0.0;
    double totZ2 = 0.0;
    seq<vector> sideChainTemp;
    
	  
	pdbfile = "./pdbfiles/" + pdbfile + ".pdb";
//	pdbfile = "C:/Documents and Settings/robert/My Documents/bioinformatics/tony/Files/" + pdbfile + ".pdb";
    ifstream iFile(pdbfile.c_str());

    if (! iFile)
    {
        cout << "Error opening input file" << endl;
    }

    while (getline(iFile,line))
    {
        name = line.substr(0,6);
        chainID = line.substr(21,1);
        atom = line.substr(12,3);
        string countstr = line.substr(7,4);
        count = (int)strtod(countstr.c_str(),&pEnd);
        string aacountstr = line.substr(22,4);
        aacount = (int)strtod(aacountstr.c_str(),&pEnd);
        
        // If there are several models in the file, we would only
        // like to display the first one
        if (name=="ENDMDL")
        	break;
        
        // Find the max & min indices for each helix
        if ((name=="HELIX ")&&(count==hel1))
        {
        	chain1id = line.substr(19,1);
        	string minIndex = line.substr(22,3);
        	string maxIndex = line.substr(34,3);
        	minIndex1 = (int)strtod(minIndex.c_str(),&pEnd);
        	maxIndex1 = (int)strtod(maxIndex.c_str(),&pEnd);
        
        //	cout << minIndex << endl;
        //	cout << maxIndex << endl;
        }
        if ((name=="HELIX ")&&(count==hel2))
        {
        	chain2id = line.substr(19,1);
        	string minIndex = line.substr(22,3);
        	string maxIndex = line.substr(34,3);
        	minIndex2 = (int)strtod(minIndex.c_str(),&pEnd);
        	maxIndex2 = (int)strtod(maxIndex.c_str(),&pEnd);
        //	cout << minIndex << endl;
        //	cout << maxIndex << endl; 
        }
        
        // Find the coordinates for each helix
        //Alpha carbons:
        if ((name=="ATOM  ")&&(aacount>=minIndex1)&&(aacount<=maxIndex1)&&(chainID==chain1id)&&(atom==" CA"))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	totX1+=x;
        	totY1+=y;
        	totZ1+=z;
        	allCoords1CA.add( vector(x,y,z));
        	aminoacids1.add(line.substr(17,3));
        	aaInds1.add(aacount);
        }
        
        //Beta carbons:
        if ((name=="ATOM  ")&&(aacount>=minIndex1)&&(aacount<=maxIndex1)&&(chainID==chain1id)&&(atom==" C "))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	allCoords1C.add( vector(x,y,z));
        	
        }
        
        //Nitrogens:
        if ((name=="ATOM  ")&&(aacount>=minIndex1)&&(aacount<=maxIndex1)&&(chainID==chain1id)&&(atom==" N "))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	allCoords1N.add( vector(x,y,z));
        	
        }
        
        //Oxygens:
        if ((name=="ATOM  ")&&(aacount>=minIndex1)&&(aacount<=maxIndex1)&&(chainID==chain1id)&&(atom==" O "))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	allCoords1O.add( vector(x,y,z));
        }
         
        if ((name=="ATOM  ")&&(aacount>=minIndex2)&&(aacount<=maxIndex2)&&(chainID==chain2id)&&(atom==" CA"))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	totX2+=x;
        	totY2+=y;
        	totZ2+=z;
        	allCoords2CA.add( vector(x,y,z));
        	aminoacids2.add(line.substr(17,3));
        	aaInds2.add(aacount);
        }
        
        
        //Beta carbons:
        if ((name=="ATOM  ")&&(aacount>=minIndex2)&&(aacount<=maxIndex2)&&(chainID==chain2id)&&(atom==" C "))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	allCoords2C.add( vector(x,y,z));
        	
        }
        
        //Nitrogens:
        if ((name=="ATOM  ")&&(aacount>=minIndex2)&&(aacount<=maxIndex2)&&(chainID==chain2id)&&(atom==" N "))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	allCoords2N.add( vector(x,y,z));
        	
        }
        
        //Oxygens:
        if ((name=="ATOM  ")&&(aacount>=minIndex2)&&(aacount<=maxIndex2)&&(chainID==chain2id)&&(atom==" O "))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	allCoords2O.add( vector(x,y,z));
        }
        
        if ((name=="ATOM  ")&&(aacount>=minIndex1)&&(aacount<=maxIndex1)&&(chainID==chain1id)&&
        					  ((atom!=" N ")&&(atom!=" O ")&&(atom!=" C ")&&(atom!=" CA")))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	sideChainTemp.add(vector(x,y,z));
        	
        }
        
        if ((name=="ATOM  ")&&(aacount>=minIndex2)&&(aacount<=maxIndex2)&&(chainID==chain2id)&&
        					  ((atom!=" N ")&&(atom!=" O ")&&(atom!=" C ")&&(atom!=" CA")))
        {
        	string xstr = line.substr(30,8);
        	string ystr = line.substr(38,8);
        	string zstr = line.substr(46,8);
        	x = strtod(xstr.c_str(),&pEnd);
        	y = strtod(ystr.c_str(),&pEnd);
        	z = strtod(zstr.c_str(),&pEnd);
        	sideChainTemp.add(vector(x,y,z));
        	
        }
        //Add the side chains
        //If we've hit the Nitrogen atom on the back bone, we now have included 
        //all of the atoms on the side chain of the previous amino acid.
        if ((name=="ATOM  ")&&(aacount>minIndex1)&&(aacount<=maxIndex1+1)&&(chainID==chain1id)&&(atom==" N "))
        {
        	sideChains1.add(sideChainTemp);
        	sideChainTemp.clear();
        }
        
        if ((name=="ATOM  ")&&(aacount>minIndex2)&&(aacount<=maxIndex2+1)&&(chainID==chain2id)&&(atom==" N "))
        {
        	sideChains2.add(sideChainTemp);
        	sideChainTemp.clear();
        }
     }
     
     // Let's translate the coordinates so the center of rotation in the 
     // window is the center of mass of the two helices.
     totX1 = totX1/allCoords1CA.size();
     totY1 = totY1/allCoords1CA.size();
     totZ1 = totZ1/allCoords1CA.size();
     totX2 = totX2/allCoords2CA.size();
     totY2 = totY2/allCoords2CA.size();
     totZ2 = totZ2/allCoords2CA.size();
     double avgX = (totX1+totX2)/2;
     double avgY = (totY1+totY2)/2;
     double avgZ = (totZ1+totZ2)/2;
     for (int i=0;i<allCoords1CA.size();i++)
     {
     	allCoords1CA[i].x -= avgX;
     	allCoords1CA[i].y -= avgY;
     	allCoords1CA[i].z -= avgZ;
     } 
     for (int i=0;i<allCoords1C.size();i++)
     {
     	allCoords1C[i].x -= avgX;
     	allCoords1C[i].y -= avgY;
     	allCoords1C[i].z -= avgZ;
     } 
     for (int i=0;i<allCoords1N.size();i++)
     {
     	allCoords1N[i].x -= avgX;
     	allCoords1N[i].y -= avgY;
     	allCoords1N[i].z -= avgZ;
     } 
     for (int i=0;i<allCoords1O.size();i++)
     {
     	allCoords1O[i].x -= avgX;
     	allCoords1O[i].y -= avgY;
     	allCoords1O[i].z -= avgZ;
     } 
     for (int i=0;i<allCoords2CA.size();i++)
     {
     	allCoords2CA[i].x -= avgX;
     	allCoords2CA[i].y -= avgY;
     	allCoords2CA[i].z -= avgZ;
     } 
     for (int i=0;i<allCoords2C.size();i++)
     {
     	allCoords2C[i].x -= avgX;
     	allCoords2C[i].y -= avgY;
     	allCoords2C[i].z -= avgZ;
     } 
     for (int i=0;i<allCoords2N.size();i++)
     {
     	allCoords2N[i].x -= avgX;
     	allCoords2N[i].y -= avgY;
     	allCoords2N[i].z -= avgZ;
     } 
     for (int i=0;i<allCoords2O.size();i++)
     {
     	allCoords2O[i].x -= avgX;
     	allCoords2O[i].y -= avgY;
     	allCoords2O[i].z -= avgZ;
     } 
     for (int i=0;i<sideChains1.size();i++)
     {
     	for (int j=0;j<sideChains1[i].size();j++)
     	{
     		sideChains1[i][j].x -= avgX;
     		sideChains1[i][j].y -= avgY;
     		sideChains1[i][j].z -= avgZ;
     	}
     } 
     for (int i=0;i<sideChains2.size();i++)
     {
     	for (int j=0;j<sideChains2[i].size();j++)
     	{
     		sideChains2[i][j].x -= avgX;
     		sideChains2[i][j].y -= avgY;
     		sideChains2[i][j].z -= avgZ;
     	}
     } 
     
}
