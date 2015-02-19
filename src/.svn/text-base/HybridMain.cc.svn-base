// Copyright (C)  2012 Zhang Initative Research Unit
// SegmentProfile.cc -description
// written by Rojan Shrestha 

#include<init_rosetta.hh>
#include<RosettaFragment.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include<functional>
#include<numeric>
#include<stdlib.h>
#include<iostream>
#include<string>
#include<map>

typedef core::pose::Pose Pose;

using namespace std;

// list<RosettaFragment> parse_rosetta_fragment(const string path2fragment, int window_size) {
//   ifstream infile(path2fragment.c_str(), ifstream::in);
//   list<RosettaFragment> fragments;
//   int rsd_pos, num_of_fragments;
//   if(infile.is_open()) {
//     bool is_first_line   = false;
//     bool is_in_top25     = false;
//     int fragment_counter = 0, top25_counter = 0; 
//     RosettaFragment rosetta_fragment_;
//     string line = "";
//     // int just_for_test = 0;
//     while(getline(infile, line)) {
//       if(line.empty()) { 
//         is_first_line = true;
//         continue; 
//       }
//       else if(line.find("position:") == 1) {
//         is_in_top25 = true; 
//         top25_counter = 0;
//         stringstream streams(line);
//         string temp1, temp2;
//         streams>>temp1>>rsd_pos>>temp2>>num_of_fragments;
//       }
//       else {
//         if(is_first_line) {
//           if(++top25_counter > 25) is_in_top25 = false; 
//           RosettaFragment rosetta_fragment(rsd_pos, line, window_size, is_in_top25); 
//           rosetta_fragment_ = rosetta_fragment;
//           is_first_line = false;
//         }
//         else  
//           rosetta_fragment_.parse_line(line); 
//         ++fragment_counter;
//       }
//       if(fragment_counter == window_size) {
//         fragments.push_back(rosetta_fragment_);
//         fragment_counter = 0;
//       }
//     }
//   }
//   infile.close();
//   return fragments;
// }

list<RosettaFragment> parse_rosetta_fragment(const string path2fragment, 
                                             int window_size, 
                                             const string tag) {
  ifstream infile(path2fragment.c_str(), ifstream::in);
  list<RosettaFragment> fragments;
  int rsd_pos, num_of_fragments, fragment_id = 0;
  if(infile.is_open()) {
    bool is_first_line   = false;
    bool is_in_top25     = false;
    int fragment_counter = 0, top25_counter = 0; 
    RosettaFragment rosetta_fragment_;
    string line = "";
    // int just_for_test = 0;
    while(getline(infile, line)) {
      if(line.empty()) { 
        is_first_line = true;
        continue; 
      }
      else if(line.find("position:") == 1) {
        is_in_top25 = true; 
        top25_counter = 0;
        stringstream streams(line);
        string temp1, temp2;
        streams>>temp1>>rsd_pos>>temp2>>num_of_fragments;
      }
      else {
        if(is_first_line) {
          if(++top25_counter > 25) is_in_top25 = false; 
          RosettaFragment rosetta_fragment(rsd_pos, line, window_size, is_in_top25, tag, ++fragment_id); 
          rosetta_fragment_ = rosetta_fragment;
          is_first_line = false;
        }
        else  
          rosetta_fragment_.parse_line(line); 
        ++fragment_counter;
      }
      if(fragment_counter == window_size) {
        fragments.push_back(rosetta_fragment_);
        fragment_counter = 0;
      }
    }
  }
  infile.close();
  return fragments;
}

void dump(map<int, list<float> > carmsds, const string filename){
  ofstream f_output_score;
  f_output_score.open(filename.c_str(), ios::out);

  map<int, list<float> >::iterator mit = carmsds.begin();
  f_output_score<<"rsd_pos, min, max, avg, size"<<endl;
  for(; mit != carmsds.end(); mit++) {
    list<float> t_carmsds = mit->second;
    t_carmsds.sort();
    list<float>::iterator flit = t_carmsds.begin(); float min = *flit;
    flit = t_carmsds.end();   float max = *(--flit);
    float avg = accumulate(t_carmsds.begin(), t_carmsds.end(), 0.0) / t_carmsds.size();

    f_output_score.width(5); f_output_score<<mit->first;
    f_output_score.setf(ios::fixed, ios::floatfield); 
    f_output_score.width(7); f_output_score.precision(2); f_output_score<<" "<<min<<" ";
    f_output_score.width(7); f_output_score.precision(2); f_output_score<<" "<<max<<" ";
    f_output_score.width(7); f_output_score.precision(2); f_output_score<<" "<<avg<<" ";
    f_output_score.width(5); f_output_score<<t_carmsds.size();
    f_output_score<<endl;
  }
  f_output_score.close();
}


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
  cout<<"  FragMeasure -r @flags -f path2fragment -w window_size -d path2pdbs"<<endl; 
  cout<<"   -r flags contains rosetta dependencies."<<endl;
  cout<<"   -n path to native pdb."<<endl;
  cout<<"   -f path to sequence fragment file."<<endl;
  cout<<"   -s path to structure fragment file."<<endl;
  cout<<"   -w window size(default 9)."<<endl;
  exit(1); 
}

int main(int argc, char* argv[]) {
  map<string, string> args;  
  option_parser(argc, argv, args);
  if(args.empty()) usage();

  Rosetta rosetta = rosetta_option_wrapper(argc, argv);
  if(rosetta.get_init() == 0) usage();
  if(args["-n"].empty()) usage();
  if(args["-f"].empty()) usage();
  if(args["-w"].empty()) usage();

  int window_size = 9;
  (args["-w"].empty()) ? window_size = 9  : window_size=atoi(args["-w"].c_str());


  Residues residues;
  Pose native_pose;
  core::io::pdb::pose_from_pdb(native_pose, args["-n"]);
  size_t total_rsds = native_pose.total_residue();

  list<RosettaFragment> seq_rosetta_fragments;
  seq_rosetta_fragments = parse_rosetta_fragment(args["-f"], window_size, "SEQ"); 
  cout<<"sequence based total fragments: "<<seq_rosetta_fragments.size()<<endl;
  list<RosettaFragment>::iterator lit = seq_rosetta_fragments.begin(); 
  for(; lit != seq_rosetta_fragments.end(); lit++) {
    int rsd_pos = lit->get_rsd_pos(); 
    if((rsd_pos + window_size - 1) > total_rsds) break;
    Pose segmented_native_pose(native_pose, rsd_pos, rsd_pos + (window_size - 1));
    lit->carmsd(segmented_native_pose);

    residues.add(rsd_pos, *lit);
  }

  list<RosettaFragment> stu_rosetta_fragments;
  stu_rosetta_fragments = parse_rosetta_fragment(args["-s"], window_size, "STR"); 
  cout<<"structure based total fragments: "<<stu_rosetta_fragments.size()<<endl;
  lit = stu_rosetta_fragments.begin(); 
  for(; lit != stu_rosetta_fragments.end(); lit++) {
    int rsd_pos = lit->get_rsd_pos(); 
    if((rsd_pos + window_size - 1) > total_rsds) break;
    Pose segmented_native_pose(native_pose, rsd_pos, rsd_pos + (window_size - 1));
    lit->carmsd(segmented_native_pose);

    residues.add(rsd_pos, *lit);
  }
  residues.do_clustering(0.5);
  cout<<"Finished!!!"<<endl;
  return 0;
}


