// Copyright (C)  2013 Zhang Initative Research Unit
// FragmentManager.cc -description
// written by Rojan Shrestha 

#include<iostream>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/Energies.hh>

#include <ObjexxFCL/format.hh>
 
#include<FragmentManager.hh>
#include<FileIO.hh>


using namespace std;

FragmentManager::FragmentManager(size_t window_size, 
                                 vector<ResidueProfile> rsd_profiles) {
  rsd_profiles_ = rsd_profiles; 
  num_of_rsds_  = rsd_profiles.size();
  window_size_  = window_size;
}

void FragmentManager::set_selection_mode(size_t mode, 
                                         size_t smode) {
  mode_  = mode;
  smode_ = smode;
}

void FragmentManager::set_selection_number(size_t number_of_fragment, 
                                           size_t fragment_per_cluster,
                                           size_t number_of_templates) {
  number_of_fragment_   = number_of_fragment;
  fragment_per_cluster_ = fragment_per_cluster;
  number_of_templates_  = number_of_templates;
}

void FragmentManager::set_native_pose(const string path2native) {
  FileIO fileio;
  fileio.read_standard_pose(native_pose_, path2native);
}

void FragmentManager::set_max_template_fragments(const size_t max_template_fragments) {
  max_template_fragments_ = max_template_fragments;
}

void FragmentManager::setup_score_titles(vector<string> &score_titles) {
  score_titles.push_back("AP_RDS"); 
  // score_titles.push_back("AP_STD");
  score_titles.push_back("SLIDING_MEAN_AP_RDS"); 
  // score_titles.push_back("SLIDING_MEAN_AP_STD"); 
}

void FragmentManager::append_score_titles(vector<string> &score_titles) {
  score_titles.push_back("NORMALIZED_SLIDING_MEAN_AP_RDS"); 
  score_titles.push_back("NORMALIZED_SLIDING_MEAN_AP_STD"); 
}

void FragmentManager::compute_sliding_mean(bool value) { 
  if(rsd_profiles_.size() <= 0) return;
  vector<ResidueProfile>::iterator outer_it = rsd_profiles_.begin(), inner_it;
  for(; outer_it <= (rsd_profiles_.end() - window_size_); outer_it++) {
    inner_it = outer_it + 1; 
    while(inner_it < outer_it + window_size_) {
      outer_it->set_slide_after_normalization(value); 
      (*outer_it) += (*inner_it);
      inner_it++;
    }
  }
  outer_it = rsd_profiles_.begin();
  for(; outer_it != rsd_profiles_.end(); outer_it++) { 
    outer_it->set_slide_after_normalization(value); 
    (*outer_it) /= window_size_;
  }
}

void FragmentManager::normalize(const vector<string> score_titles) {
  try {
    if(rsd_profiles_.size() <= 0) throw " data is not available."; 
    map<string, float> max_values, min_values;

    if(score_titles.empty()) throw " score titles are missing.";

    vector<string>::const_iterator vit = score_titles.begin();
    vector<ResidueProfile>::iterator it = rsd_profiles_.begin();
    for(; it != rsd_profiles_.end(); it++) { 
      for(vit = score_titles.begin(); vit != score_titles.end(); vit++) { 
        if(!it->normalize_score(*vit)) 
          throw " normalziation of " + (*vit) +" is stopped.";
      }
    }

    it = rsd_profiles_.begin(); 
    for(; it != rsd_profiles_.end(); it++) 
      for(vit = score_titles.begin(); vit != score_titles.end(); vit++) 
        it->compare_value(*vit, min_values, max_values);

    it = rsd_profiles_.begin(); 
    for(; it != rsd_profiles_.end(); it++) { 
      for(vit = score_titles.begin(); vit != score_titles.end(); vit++) { 
        if(!(it->normalize_score(min_values[*vit], max_values[*vit], *vit))) 
          throw " global normalization of " + (*vit) + " gave an error.";
      }
    }
  }
  catch(string e) {
    std::cerr<<"warning: "<<__FILE__<<" "<<__LINE__<<" "<<e<<endl;
  }
}

void FragmentManager::write_scores(ofstream &f_output_score, 
                                   vector<ResidueProfile> rsd_profiles,
                                   const string title) {
  try {
    // dump all the scores...(it should be selected using particular scores)
    vector<ResidueProfile>::iterator it  = rsd_profiles_.begin();
    it->write_score_title(f_output_score); 
    for(;it != rsd_profiles_.end(); it++) {
      f_output_score.width(5); 
      it->write_rsd_scores(f_output_score, "");
      f_output_score<<endl;
    }
    //dum all the scores for selected models
    for(it  = rsd_profiles.begin();it != rsd_profiles.end(); it++) {
      f_output_score.width(5); 
      it->write_rsd_scores(f_output_score, title);
      f_output_score<<endl;
    }
    if(native_pose_.empty()) return; 
    // for(it  = rsd_profiles.begin(); it != rsd_profiles.end(); it++) 
    //   it->dump_average_rmsd_to_native("FRAGMENT_RMSD", f_output_score); // change in Residue  
  }
  catch(const char* error) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<error<<endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknow eror "<<endl;
    exit(1);
  }
}

void FragmentManager::get_top_models_for_each_position(vector<ResidueProfile> &sorted_profiles, 
                                                       const size_t num_of_tops,  
                                                       string score_title) {
  cout<<"# of best candidates: "<<num_of_tops<<" using "<<score_title<<endl;
  sorted_profiles.clear();
  sorted_profiles.reserve(num_of_tops);

  vector<ResidueProfile>::iterator it = rsd_profiles_.begin();
  for(; it != rsd_profiles_.end(); it++) {
    it->set_selection_number(number_of_fragment_, fragment_per_cluster_, number_of_templates_);  
    ResidueProfile rsd_profile = it->sort_by_score(score_title, num_of_tops);
    sorted_profiles.push_back(rsd_profile);
  }
}

void FragmentManager::get_top_models_for_each_position(vector<ResidueProfile> &sorted_profiles, 
                                                       const float threshold,  
                                                       string score_title) {
  cout<<"threshold for "<<score_title<<":  "<<threshold<<endl;
  sorted_profiles.clear();

  bool is_not_increased = true;
  float temp_threshold, factor = 0.05; 
  vector<ResidueProfile>::iterator it = rsd_profiles_.begin();
  while(it != rsd_profiles_.end()) { 
  // for(; it != rsd_profiles_.end() && is_not_increased; it++) {
    if(not is_not_increased) {
      temp_threshold = threshold + factor; 
      factor += 0.05;
      cerr<<"  warning: new threhold is "<<temp_threshold<<endl;
    }
    else {
      it->set_selection_number(number_of_fragment_, fragment_per_cluster_, number_of_templates_);  
      temp_threshold = threshold; 
      factor = 0.05;
    }

    ResidueProfile rsd_profile = it->sort_by_score(score_title, temp_threshold, is_not_increased);
    if(is_not_increased) {
      sorted_profiles.push_back(rsd_profile);
      it++;
    }
  }
  // cout<<"# selected fragments:             "<<sorted_profiles.size()<<endl;
}

void FragmentManager::extend_rsd_positions(vector<ResidueProfile> &sorted_rsds, 
                                           const string score_title, const size_t num_of_tops) {
  if(score_title.empty()) return;
  vector<ResidueProfile>::iterator outer_it = sorted_rsds.begin(), inner_it;
  for(; outer_it != sorted_rsds.end(); outer_it++) {
    map<string, float> means = outer_it->get_model_n_sorted_score(score_title, num_of_tops);
    map<string, float>::iterator mit = means.begin();
    for(; mit != means.end(); mit++) {
      size_t counter = 0, gap = 0; 
      const string model_id = mit->first;
      for(inner_it = outer_it + 1; inner_it != sorted_rsds.end(); inner_it++) {
        if(!inner_it->is_model_exist(model_id, score_title)) gap++; 
        counter++;
        if(counter >= window_size_) break;
      }
      //ROS:: have to look append_score functions
      // outer_it->append_score(model_id, "FRAGMENT_LENGTH", float(counter)); 
      // outer_it->append_score(model_id, "NUMBER_OF_GAPS", float(gap)); 
    }
  }
}

void FragmentManager::compute_carmsd2native(vector<ResidueProfile> &rsd_profiles, 
                                            map<string, Pose> poses) {
  if(native_pose_.empty()) return;
  cout<<"computing ca-rsmsd with crystal structures..."<<endl;
  vector<ResidueProfile>::iterator it;
  for(it = rsd_profiles.begin(); it != rsd_profiles.end(); it++) {
    size_t rsd_pos = it->get_rsd_pos();
    if((num_of_rsds_ - rsd_pos) < window_size_) break;

    vector<string> model_scores  = it->get_model_ids();
    vector<string>::iterator vit =  model_scores.begin();
    assert(num_of_rsds_ == poses[*vit].total_residue());
    
    for(; vit != model_scores.end(); vit++) {
      Pose pose = poses[*vit];
      float lcarmsd = CA_rmsd(native_pose_, pose, rsd_pos, rsd_pos + window_size_ - 1);
      it->append_score(*vit, "FRAGMENT_RMSD", lcarmsd); 
    }
  }
} 

void FragmentManager::do_clustering(vector<ResidueProfile> &rsd_profiles, 
                                    const float cluster_radius, 
                                    const map<string, Pose> poses, 
                                    const string prefix, 
                                    ofstream &f_output) {
  cout<<"clustering..."<<endl;
  try { 
    string fragment_file = prefix + "_fragment.fragK";
    string template_file = prefix + "_template.fragK";
    ofstream f_output_fragment, f_output_template;
    f_output_fragment.open(fragment_file.c_str(), ios::out);
    f_output_template.open(template_file.c_str(), ios::out);
    if(!f_output_fragment.is_open()) throw " fragment file cannot be opened.";
    if(!f_output_template.is_open()) throw " template file cannot be opened.";

    vector<ResidueProfile>::iterator vit = rsd_profiles.begin(); 
    for(; vit != rsd_profiles.end() - (window_size_ - 1); vit++) {
      vit->set_selection_mode(mode_, smode_); 
      vit->set_selection_number(number_of_fragment_, fragment_per_cluster_, number_of_templates_);  
      size_t rsd_pos = vit->get_rsd_pos();

      if(rsd_pos > num_of_rsds_ - (window_size_ - 1)) break;
      map<string, Pose> fragment_poses;
      vit->get_standard_fragment_poses(window_size_, poses, fragment_poses);
      if(fragment_poses.empty()) throw " empty fragment poses."; 
      vit->set_total_residue(rsd_profiles.size());
      vit->do_clustering(cluster_radius, fragment_poses, 
                         f_output_fragment, f_output_template, f_output);
    }
    f_output_fragment.close();
    f_output_template.close();
  }
  catch(const char* error) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<error<<endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<endl;
    exit(1);
  }
}

// private sections
string FragmentManager::get_fragment_filename() {
  string prefix = "fragments";
  string suffix = "frag";
  int window_size = window_size_; 
  stringstream ss_window_size;
  ss_window_size<<window_size;
  return prefix + "." + suffix + ss_window_size.str();
}
