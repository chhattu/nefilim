// Copyright (C)  2013 Zhang Initative Research Unit
// FragmentManager.hh -description
// written by Rojan Shrestha 

#ifndef  FRAGMENT_MANAGER_HEADER_HH
#define  FRAGMENT_MANAGER_HEADER_HH

#include<iostream>
#include<sstream>
#include<fstream>
#include<cassert>
#include<vector>
#include<map>
#include<utility>
#include<cmath>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>

#include<Residue.hh>

typedef core::pose::Pose Pose;
typedef core::conformation::Conformation Conformation;
typedef std::pair<int, float> t_pair;

using namespace std;
using namespace core::scoring;


class FragmentManager {
  private:
    vector<ResidueProfile> rsd_profiles_; 
    size_t window_size_;
    size_t num_of_rsds_;
    size_t max_template_fragments_;

    Pose native_pose_;
    size_t mode_;
    size_t smode_;

    size_t number_of_fragment_;
    size_t fragment_per_cluster_;
    size_t number_of_templates_; 

  private:
    string get_fragment_filename();

  public:
    FragmentManager();

    FragmentManager(size_t window_size, 
                    vector<ResidueProfile> rsd_profiles);

    void set_selection_mode(size_t mode, 
                            size_t smode);

    void set_selection_number(size_t number_of_fragment, 
                              size_t fragment_per_cluster,
                              size_t number_of_templates);

    void set_native_pose(const string path2native);

    void setup_score_titles(vector<string> &score_titles); 

    void append_score_titles(vector<string> &score_titles); 

    void set_max_template_fragments(const size_t max_template_fragments);

    void compute_sliding_mean(bool value); 

    void normalize(const vector<string> score_titles); 

    void write_scores(ofstream &f_output_score, 
                      vector<ResidueProfile> rsd_profiles,
                      const string title); 

    void write_selected_templates_native_RMSD(vector<ResidueProfile> selected_rsds,
                                              ofstream &outfile); 
    
    void write_cluster_info(vector<ResidueProfile> &selected_rsds, 
                            ofstream &outfile);

    void get_top_models_for_each_position(vector<ResidueProfile> &sorted_profiles, 
                                          const size_t num_of_tops, 
                                          const string score_title); 

    void get_top_models_for_each_position(vector<ResidueProfile> &sorted_profiles, 
                                          const float threshold,  
                                          string score_title); 

    void extend_rsd_positions(vector<ResidueProfile> &sorted_rsds, 
                              const string score_title, const size_t num_of_tops); 

    void compute_carmsd2native(vector<ResidueProfile> &rsd_profiles, 
                               map<string, Pose> poses); 

    void do_clustering(vector<ResidueProfile> &rsd_profiles, 
                       const float cluster_radius, 
                       const map<string, Pose> poses,
                       const string prefix, 
                       ofstream &f_output);


};

#endif

