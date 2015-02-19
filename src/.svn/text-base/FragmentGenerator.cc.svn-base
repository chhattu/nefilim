// Copyright (C)  2012 Zhang Initative Research Unit
// SegmentProfile.cc -description
// written by Rojan Shrestha 

#include<init_rosetta.hh>
#include<FileIO.hh>
#include<Residue.hh>
// #include<FragmentGeneratorMain.hh>
#include<FragmentManager.hh>
// #include<FragmentCluster.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include<stdlib.h>
#include<iostream>
#include<string>
#include<map>

typedef core::pose::Pose Pose;

using namespace std;

Rosetta rosetta_option_wrapper(int argc, char* argv[]) { if(argc < 2) return Rosetta(0);
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
  cout<<"  Goldselector -r @flags -s path2score -m path2models -d dblist"<<endl; 
  cout<<"                  [-a path2native -n # best_fragments -w window size]"<<endl;
  cout<<"   -r flags contains rosetta dependencies."<<endl;
  cout<<"   --path_to_scores         path to apmds score files."<<endl;
  cout<<"   --path_to_models         path to pdb models."<<endl;
  cout<<"   --window_size            (default 9)."<<endl;
  cout<<"   --path_to_native         path to the native structure."<<endl;
  cout<<"   --cutoff                normalized score cut off."<<endl;
  cout<<"   --cluster_selection_mode mode of selecting models [1 for RND, 2 for SEQ]."<<endl;
  cout<<"   --proportion             [1 for equal and 2 for density ]."<<endl;
  cout<<"   --number_of_fragment"<<endl;
  cout<<"   --fragment_per_cluster"<<endl; 
  cout<<"   --number_of_template"<<endl; 
  cout<<"   --cluster_radius"<<endl; 
  exit(1); 
}

int main(int argc, char* argv[]) {
  map<string, string> args;  
  option_parser(argc, argv, args);
  if(args.empty()) usage();

  cout<<"Rosetta: (uses for parsing pdbs and others)"<<endl;
  cout<<"-------------------------------------------"<<endl;
  Rosetta rosetta = rosetta_option_wrapper(argc, argv);
  if(rosetta.get_init() == 0) usage();

  vector<ResidueProfile> residues_profiles;
  FileIO fileio;
  fileio.read_aprdss(args["--path_to_scores"], residues_profiles, args["--path_to_native"].empty());

  size_t window_size, smode, mode;
  int number_of_fragment, fragment_per_cluster, number_of_template;
  float cutoff, cluster_radius; 
  string output_filename;
  bool verbose=true;
  (args["--window_size"].empty()) ? 
  window_size=3 : window_size =atoi(args["--window_size"].c_str());
  (args["--cutoff"].empty()) ? 
  cutoff=1.00 : cutoff=atof(args["--cutoff"].c_str());
  (args["--output_filename"].empty()) ? 
  output_filename="rsd_scores.out" : output_filename=args["--output_filename"];
  (args["--cluster_selection_mode"].empty()) ? 
  mode = 1: mode = atoi(args["--cluster_selection_mode"].c_str());

  (args["--proportion"].empty()) ? 
  smode = 1: smode = atoi(args["--proportion"].c_str());
  (args["--number_of_fragment"].empty()) ? 
  number_of_fragment=25: number_of_fragment=atoi(args["--number_of_fragment"].c_str());
  (args["--fragment_per_cluster"].empty()) ? 
  fragment_per_cluster=5: fragment_per_cluster= atoi(args["--fragment_per_cluster"].c_str());
  (args["--number_of_template"].empty()) ? 
  number_of_template=25: number_of_template= atoi(args["--number_of_template"].c_str());
  (args["--cluster_radius"].empty()) ? 
  cluster_radius = 1.0 : cluster_radius = atof(args["--cluster_radius"].c_str());

  // (args["-v"].empty()) ? verbose=false : verbose=static_cast<bool)args["-v"];

  cout<<"Computing with proportion "<<smode<<" and selection "<<mode<<endl;
  cout<<"------------"<<endl;
  FragmentManager fragment_manager(window_size, residues_profiles);
  fragment_manager.set_selection_mode(mode, smode);
  fragment_manager.set_selection_number(number_of_fragment, fragment_per_cluster, number_of_template); 

  cout<<"computing sliding mean of window size "<<window_size<<endl;
  fragment_manager.compute_sliding_mean(false);

  vector<string> score_titles;
  fragment_manager.setup_score_titles(score_titles);
  if(score_titles.empty()) { 
    cout<<"fatal: no score title for filter."<<endl;
    return 0;
  }
  fragment_manager.normalize(score_titles); 

  cout<<endl<<"Analysis..."<<endl;
  cout<<"-----------"<<endl;
  // fragment_manager.get_top_models_for_each_position(sorted_profiles, groupN, score_title); 
  // fragment_manager.extend_rsd_positions(sorted_profiles, score_title, 5); 
  // vector<ResidueProfile>::iterator vit = sorted_profiles.begin();

  map<string, Pose> poses;
  cout<<"reading denovo models..."<<endl;
  fileio.read_poses_no_stout(args["--path_to_models"], poses); 

  vector<ResidueProfile> sorted_profiles;
  vector<string>::iterator sit = score_titles.begin();

  if(!args["--path_to_native"].empty()) 
    fragment_manager.set_native_pose(args["--path_to_native"]);
  
  ofstream f_output_score;
  f_output_score.open(output_filename.c_str(), ios::out);
  for(; sit != score_titles.end(); sit++) {
    f_output_score<<"scoring title: "<<*sit<<endl;
    fragment_manager.get_top_models_for_each_position(sorted_profiles, 
                                                      cutoff, 
                                                      "NORMALIZED_" + *sit); 
    if(!args["--path_to_native"].empty()) 
      fragment_manager.compute_carmsd2native(sorted_profiles, poses); 

    if(verbose)
      fragment_manager.write_scores(f_output_score, sorted_profiles, "NORMALIZED_" + *sit);
    fragment_manager.do_clustering(sorted_profiles, 
                                   cluster_radius, 
                                   poses, "N_" + *sit, 
                                   f_output_score);
    cout<<endl;
  }
  f_output_score.close();

  cout<<"done!!!."<<endl;
}
