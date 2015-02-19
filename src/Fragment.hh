// Copyright (C)  2012 Zhang Initative Research Unit
// Fragment.hh -description
// written by Rojan Shrestha 

#ifndef  FRAGMENT_HEADER_HH
#define  FRAGMENT_HEADER_HH

#include<iostream>
#include<fstream>
#include<vector>
#include<string>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

using namespace std;

class Fragment {

  typedef core::pose::Pose Pose;

  private:
    size_t rsd_pos_;
    size_t frag_size_;
    size_t max_num_frags_;
    size_t max_frag_per_model_;
    size_t frag_per_model_;

    vector<char> sec_strs_;
    vector<string> pdb_ids_;
    vector<string> chain_ids_;
    vector<string> amino_acids_;

    vector<size_t> rsd_poss_;
    vector<float> phis_, psis_, omegas_, ca_rmsds_;

    float unknown1, unknown2, unknown4; 
    size_t unknown_id;

  private:
    void inline parse_filename(string filename, string &pdb_id, string &chain_id) {
      //no thread safe
      pdb_id   = filename.substr(0,4);
      chain_id = filename.substr(5,1);
    }

  public:

    Fragment();

    Fragment(const size_t rsd_pos, const size_t frag_size);

    Fragment(const size_t rsd_pos, const size_t frag_size, 
             const size_t max_num_frags, const size_t max_frag_per_model); 

    void insert(Pose pose, const string pdb_id, const size_t rsd_pos, const float ca_rmsd, ofstream &f_output); 

    void write_rosetta_fragment_format(ofstream &f_output, Pose pose, 
                                       const string pdb_id, 
                                       const size_t rsd_pos, 
                                       const float ca_rmsd);

    void get_native_fragment(ofstream &f_output, 
                             const Pose native_pose, 
                             const string pdb_id, 
                             const size_t rsd_pos, 
                             const float ca_rmsd); 

    void show();
    
    size_t get_total_fragments(); 

    size_t get_fragment_per_model(); 

    void reset_fragment_per_model(); 

};

#endif

