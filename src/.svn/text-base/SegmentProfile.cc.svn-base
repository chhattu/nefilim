// Copyright (C)  2012 Zhang Initative Research Unit
// SegmentProfile.cc -description
// written by Rojan Shrestha 

// #include<SegmentProfileMain.hh>
// #include<ModelScore.hh>
#include<init_rosetta.hh>
#include<FileIO.hh>
#include<Residue.hh>
#include<FragmentGeneratorMain.hh>
#include<ProfileVectorImpl.hh>
#include<Distance.hh>


#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include<stdlib.h>
#include<iostream> 
#include<string>
#include<map>

typedef core::pose::Pose Pose;

using namespace std;

Rosetta rosetta_option_wrapper(int argc, char* argv[]) {
  if(argc < 2) return Rosetta(0);
  char** argvs = new char*[2];
  const int new_argc = 2;

  for(int i = 0; i < argc ; i++) {
    if(strcmp(argv[i],"-r") == 0) {
      argvs[0] = argv[0];
      argvs[1] = argv[i+1];
      return Rosetta(new_argc, argvs);
    }
  }
  return Rosetta(0);
}

// what could be error arisen situations?
void option_parser(int argc, char* argv[], map<string, string> &args) {
  int i = 1;
  while(i < argc) { 
    if(strcmp(argv[i],"-r") != 0) {
      args[argv[i]] = argv[i+1]; 
      i+=2;
    }
    else 
      i+=2;
  }
}

void usage() {
  cout<<"Usage:"<<endl;
  cout<<"  SegmentProfiler -s path2score -m path2models -r @flags"<<endl;
  cout<<"   -m list of path of pdb models."<<endl;
  cout<<"   -q query pdb models."<<endl;
  cout<<"   -r flags contains rosetta database and others."<<endl;
  cout<<"   -fs fragment sizes."<<endl;
  cout<<"   -bs bin sizes."<<endl;
  cout<<"   -nr number of references."<<endl;
  cout<<"   -ra radius."<<endl;
  cout<<"   --mode {DB|SE}."<<endl;
  exit(1); 
}

int main(int argc, char* argv[]) {
  map<string, string> args;  
  option_parser(argc, argv, args);
  if(args.empty()) usage();

  //start here...
  bool verbose = false; // put debug flag to make it true
  size_t frag_size = 9; //default value
  float  bin_size  = 0.5; //default value
  size_t number_of_reference = 3; //default value
  float radius = 0.3;
  string mode = "DB";
  if(args["-m"].empty()) usage();
  if(!args["-fs"].empty())  frag_size = atoi(args["-fs"].c_str());
  if(!args["-bs"].empty()) bin_size = atof(args["-bs"].c_str());
  if(!args["-nr"].empty()) number_of_reference = atoi(args["-nr"].c_str());
  if(!args["-ra"].empty()) radius = atof(args["-ra"].c_str());
  if(!args["--mode"].empty()) mode = args["--mode"]; 

  Rosetta rosetta = rosetta_option_wrapper(argc, argv);
  if(rosetta.get_init() == 0) usage();

  ProfileVectorImpl prof_vector_impl(args["-m"], frag_size, bin_size, verbose); 

  Distance distance; 
  if(mode=="DB") {
    prof_vector_impl.filter_intra_submodels(0.0);
    clock_t start = clock(); 
    const string input_path = "./coordinates.out";
    distance.get_reference_submodels(number_of_reference, input_path); 
    distance.dump_reference_submodels("./reference.out"); 
    distance.compute_and_dump_carmsd(input_path); 
    clock_t end = clock(); 
    cout<<"Total cpu time: "
        <<(end - start)/(double)CLOCKS_PER_SEC
        <<" seconds."<<endl;
  }
  else if(mode=="SB") {
    if(args["-q"].empty()) usage();
    vector<SubModel> submodels;
    prof_vector_impl.filter_intra_submodels(submodels, args["-q"] , 0.0); 
    cout<<"Query size: "<<submodels.size()<<endl;
    distance.set_reference_submodels("./reference.out"); 
    map<string, map<string, float> > found_fragments;
    clock_t start = clock(); 
    distance.compute_distance_to_references_n_search(submodels, found_fragments, radius); 
    clock_t end = clock(); 
    cout<<"Total cpu time: "
        <<(end - start)/(double)CLOCKS_PER_SEC
        <<" seconds."<<endl;

    start = clock(); 
    distance.show_found_fragments(found_fragments); 

    map<string, SubModel> selected_submodels; 
    distance.get_selected_submodels(found_fragments, "./coordinates.out", selected_submodels); 
    // cout<<"Submodel size: "<<selected_submodels.size()<<endl;
    // cout<<"Number of selected submodels "<<selected_submodels.size()<<endl;
    distance.compute_and_dump_cosine(submodels, selected_submodels, found_fragments, "./cosine.out"); 
    end = clock(); 
    cout<<"Q-score calculation cpu time: "
        <<(end - start)/(double)CLOCKS_PER_SEC
        <<" seconds."<<endl;

    //debugging...
    cout<<"Computing rmsd to query..."<<endl;
    const string input_path = "./coordinates.out";
    start = clock(); 
    distance.set_reference_submodels(args["-q"], frag_size); 
    distance.compute_and_dump_carmsd(input_path); 
    end = clock(); 
    cout<<"RMSD calculation cpu time: "
        <<(end - start)/(double)CLOCKS_PER_SEC
        <<" seconds."<<endl;
  }
}
