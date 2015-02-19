#ifndef SPMAIN_HEADER_HH 
#define SPMAIN_HEADER_HH 

#include<iostream>
#include<fstream>
#include<cassert>
#include<vector>
#include<map>
#include<utility>
#include<cmath>

#include<Residue.hh>
#include<ModelEnergy.hh>
#include<Fragment.hh>
#include<RosettaWrapperFwd.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>

typedef core::pose::Pose Pose;
typedef std::pair<int, float> t_pair;

using namespace std;
using namespace core::scoring;

enum scoreTitle { 
  ap_rds, ap_std, ca_rmsd, rsd_energy
};

struct comp{
  bool operator()(const t_pair &pair1, const t_pair &pair2) {
    return pair1.second > pair2.second;
  }
};

// void dump_residue_score_model(ofstream &f_log_output, string model_id, float aprds) {
//   f_log_output<<model_id<<" ";
//   f_log_output.setf(ios::fixed, ios::floatfield); 
//   f_log_output.width(5); f_log_output.precision(2); f_log_output<<aprds<<" ";
// }

// string inline convert_int2string(const int input_int) {
//   ostringstream convert;   
//   convert <<input_int;      
//   return convert.str(); 
// }

bool dump_analyzed_scores(vector<ResidueProfile> rsd_profiles, string filename) {
  ofstream f_output_score;
  if(filename.empty()) filename = "analyzed_scores.out";
  f_output_score.open(filename.c_str(), ios::out);

  vector<ResidueProfile>::iterator it  = rsd_profiles.begin();
  const string score_title = "LCA_RMSD";
  for(;it != rsd_profiles.end(); it++) 
    it->dump_analyzed_scores(score_title, f_output_score); 
  f_output_score.close();
  return true;
}

bool dump_rsd_scores(vector<ResidueProfile> rsd_profiles, string filename) {
  ofstream f_output_score;
  if(filename.empty()) filename = "scores.out";
  f_output_score.open(filename.c_str(), ios::out);

  vector<ResidueProfile>::iterator it  = rsd_profiles.begin();
  it->dump_scores_title(f_output_score); 
  for(;it != rsd_profiles.end(); it++) {
    f_output_score.width(5); 
    it->dump_rsd_scores(f_output_score);
    f_output_score<<endl;
  }
  f_output_score.close();
  return true;
}

bool dump_top_rsd_scores(vector<ResidueProfile> rsd_profiles, string filename, float cutoff) {
  vector<ResidueProfile>::iterator it = rsd_profiles.begin(); 
  vector<string> score_titles = it->get_score_titles();

  map<string, vector<ResidueProfile> > final_rsd_profiles;
  vector<string>::iterator vit = score_titles.begin(); 

  ofstream f_output_score;
  if(filename.empty()) filename = "top_selected_models.out";
  f_output_score.open(filename.c_str(), ios::out);

  while(vit != score_titles.end()) {
    if((*vit) == "AP_LENGTH" || (*vit) == "CA_RMSD" ||
      (*vit) == "AP_RDS" || (*vit) == "AP_STD" || 
      (*vit) == "SLIDING_MEAN_AP_RDS" || 
      (*vit) == "SLIDING_MEAN_AP_STD" ) { vit++; continue;}

    f_output_score<<(*vit)<<"\t\tMin.\t\tMax.\t\tAvg."<<endl;
    f_output_score<<"----------------------------------"<<endl;
    it = rsd_profiles.begin(); 
    for(; it != rsd_profiles.end(); it++) { 
      ResidueProfile rsd_profile = it->sort_by_score((*vit), cutoff);
      rsd_profile.dump_analyzed_scores("LCA_RMSD", f_output_score); 
    }
    f_output_score<<endl<<endl;
    vit++;
  }
  f_output_score.close();
  return true;
}

  // void dump_residue_scores(vector<ResidueProfile> residues_profiles, string filename) {
  //   ofstream f_output, f_log_output;
  //   f_output.open(filename.c_str(), ios::out);
  //   const string score_title1 = "AP_RDS"; 
  // 
  //   const string log_filename = "model_scores.out";
  //   f_log_output.open(log_filename.c_str(), ios::out);
  // 
  //   for(vector<ResidueProfile>::iterator it  = residues_profiles.begin(); 
  //                                        it != residues_profiles.end(); 
  //                                        it++) {
  //     const size_t cur_rsd_pos = it->get_rsd_pos();
  //     f_log_output.width(4); f_log_output<<cur_rsd_pos<<" ";
  // 
  //     map<string, float> model_scores = it->get_modelnscore(score_title1);
  //     for(map<string, float>::iterator mit  = model_scores.begin(); 
  //                                      mit != model_scores.end(); 
  //                                      mit++) {
  //       f_output.setf(ios::fixed, ios::floatfield); 
  //       f_output.width(7); f_output.precision(2); f_output<<mit->second;
  //       dump_residue_score_model(f_log_output, mit->first, mit->second);
  //     }
  //     f_output<<endl;
  //     f_log_output<<endl;
  //   }
  //   f_output.close();
  //   f_log_output.close();
  // }


  vector<ResidueProfile> get_cluster_centers(const string filename, vector<ResidueProfile> residue_profiles) {
    map<size_t, vector<size_t> > m_cluster_sizes;
    map<size_t, vector<float> > m_cluster_means;
    const string score_title = "AP_RDS";

    FileIO fileio;
    cout<<"parsing started."<<endl;
    fileio.parse_r_ckmean_output(filename, m_cluster_sizes, m_cluster_means); 
    cout<<"parsing finished."<<endl;

    vector<ResidueProfile>::iterator it = residue_profiles.begin();
    map<size_t, vector<float> >::iterator mit = m_cluster_means.begin();
    vector<ResidueProfile> selected_rsd_profiles;

    while( it != residue_profiles.end() || mit != m_cluster_means.end()) {
      const size_t cur_rsd_pos = it->get_rsd_pos();
      // cout<<"Outside residue position "<<cur_rsd_pos<<" "<<mit->first<<endl;
      if(cur_rsd_pos == mit->first) {
        vector<float> cluster_means = mit->second; 
        cout<<"Residue position "<<cur_rsd_pos<<endl;
        cout<<"Number of cluster "<<cluster_means.size()<<endl;
        // for debug
        cout<<"Rsd pos: "<< mit->first<<endl;
        vector<float>::iterator tit = cluster_means.begin();
        cout<<"Cluster center: ";
        for(; tit != cluster_means.end(); tit++) {
          cout<<(*tit)<<" ";
        }
        cout<<"I am there."<<endl;
        cout<<endl;

        if(cluster_means.empty()) continue;
        vector<string> cluster_modelnames = it->get_modelnames(score_title, cluster_means);
        if(cluster_modelnames.empty()) continue;
        cout<<"I am here."<<endl;
        selected_rsd_profiles.push_back(ResidueProfile(cur_rsd_pos, score_title, cluster_modelnames, cluster_means));
      }
      mit++; it++;
    }
    vector<ResidueProfile>::iterator vit = selected_rsd_profiles.begin();
    for(; vit != selected_rsd_profiles.end(); vit++) 
      vit->show();
    return selected_rsd_profiles;
  }

  bool normalize(vector<ResidueProfile> &rsds_profiles, const int w_size) {
    if(rsds_profiles.size() <= 0) return false;

    float min_mean  = 1000000.0, max_mean  = 0.0; 
    float min_stdev = 1000000.0, max_stdev = 0.0;
    float min_sl_mean = 1000000.0, max_sl_mean = 0.0;
    float min_sl_stdev = 1000000.0, max_sl_stdev = 0.0;
    vector<ResidueProfile>::iterator it = rsds_profiles.begin();
    for(; it != rsds_profiles.end(); it++) {
      if(!it->normalize_score("AP_RDS")) { 
        cout<<" why error"<<endl;
        // return false; 
      }
      if(!it->normalize_score("AP_STD")) { 
        cout<<" why error in AP_STD "<<endl;
        // return false;
      }
      if(it > rsds_profiles.end() - w_size) continue;
      if(!it->normalize_score("SLIDING_MEAN_AP_RDS")) {
        cout<<" why error in SLIDING MEAN_AP_RDS"<<endl;
        // return false;
      }
      if(!it->normalize_score("SLIDING_MEAN_AP_STD")) {
        cout<<" why error in SLIDING MEAN_AP_STD"<<endl;
        // return false;
      }
    }
    // return false;

    cout<<"I am there "<<endl;
    it = rsds_profiles.begin();
    for(; it != rsds_profiles.end(); it++) {
      it->compare_value("AP_RDS", min_mean, max_mean);
      it->compare_value("AP_STD", min_stdev, max_stdev);
      it->compare_value("SLIDING_MEAN_AP_RDS", min_sl_mean, max_sl_mean);
      it->compare_value("SLIDING_MEAN_AP_RDS", min_sl_stdev, max_sl_stdev);
    }

    cout<<"Normalized with all data ..."<<endl;
    cout<<"APRDS (min. & max.): "<<min_mean<<" "<<max_mean<<endl;
    cout<<"APSTD (min. & max.): "<<min_stdev<<" "<<max_stdev<<endl;
    cout<<"Sliding APRDS (min. & max.): "<<min_sl_mean<<" "<<max_sl_mean<<endl;
    cout<<"Sliding APSTD (min. & max.): "<<min_sl_stdev<<" "<<max_sl_stdev<<endl;

    it = rsds_profiles.begin();
    for(; it != rsds_profiles.end(); it++) {
      if(!it->normalize_score(min_mean, max_mean, "AP_RDS")) return false;
      if(!it->normalize_score(min_stdev, max_stdev, "AP_STD")) return false; 
      if(!it->normalize_score(min_sl_mean, max_sl_mean,"SLIDING_MEAN_AP_RDS")) return false;
      if(!it->normalize_score(min_sl_stdev, max_sl_stdev, "SLIDING_MEAN_AP_STD")) return false;
    }
    return true;
}

void compute_sliding_mean(vector<ResidueProfile> &rsds_profiles, const int win_size) { 
  if(rsds_profiles.size() <= 0) return;
  cout<<"computing sliding mean: win. size "<<win_size<<endl;
  vector<ResidueProfile>::iterator outer_it = rsds_profiles.begin(), inner_it;
  for(; outer_it <= (rsds_profiles.end()-win_size); outer_it++) {
    inner_it = outer_it + 1; 
    while(inner_it < outer_it + win_size){
      (*outer_it) += (*inner_it);
      inner_it++;
    }
  }
  outer_it = rsds_profiles.begin();
  for(; outer_it != rsds_profiles.end(); outer_it++) 
    (*outer_it) /= win_size;
}


vector<ResidueProfile> getNBestFragments(vector<ResidueProfile> residues_profiles, 
                                         const size_t topN, const size_t groupN,  
                                         const size_t frameS, scoreTitle score_title) {

  string score_title1 = "", score_title2 = "", score_title3 = "";

  if(topN < 1 ) {
    cout<<"warning: no proper value for top N."<<endl;
    residues_profiles.clear();
    return residues_profiles;
  }

  if(score_title == ap_rds) {
    score_title1 = "AP_RDS"; score_title2 = "AP_LENGTH"; score_title3 = "AP_GAP";
  }
  // else if(score_title == rsd_energy) {
  //   score_title1 = "rsd_energy"; score_title2 = "rsd_length";
  // }

  cout<<"# of best candidate: "<<groupN<<endl;

  vector<ResidueProfile>::iterator it;
  vector<ResidueProfile> sorted_residues;
  sorted_residues.reserve(residues_profiles.size());
  for(it = residues_profiles.begin(); it != residues_profiles.end(); it++) {
    ResidueProfile rsd_profile = it->sort_by_score(score_title1, groupN);
    sorted_residues.push_back(rsd_profile);
  }

  vector<ResidueProfile>::iterator outer_it = sorted_residues.begin(), inner_it;
  for(; outer_it != sorted_residues.end(); outer_it++) {
    map<string, float> means = outer_it->get_model_n_sorted_score("AP_RDS", groupN);
    map<string, float>::iterator mit = means.begin();
    for(; mit != means.end(); mit++) {
      size_t counter = 0, gap = 0; 
      const string model_id = mit->first;
      for(inner_it = outer_it + 1; inner_it != sorted_residues.end(); inner_it++) {
        if(!inner_it->is_model_exist(model_id, "AP_RDS")) gap++; 
        counter++;
        if(counter>=frameS) break;
      }
      outer_it->append_score(model_id, score_title2, float(counter)); 
      outer_it->append_score(model_id, score_title3, float(gap)); 
    }
  }

  //get topN models based on less gap
  //vector<ResidueProfile> topN_residues;
  //for(outer_it = sorted_residues.begin(); outer_it != sorted_residues.end(); outer_it++) {
  //  topN_residues.push_back(outer_it->sort_by_score(score_title1, topN));
  //}

  //cout<<"dumping scores into scores.out"<<endl;
  //dump_rsd_scores(sorted_residues, "b_scores.out");
  //dump_rsd_scores(topN_residues, "n_scores.out");
  //return topN_residues;
  return sorted_residues;
}

void compute_energies(Rosetta rosetta,vector<ResidueProfile> &residues_profiles, map<string, Pose> poses) {
  ModelEnergy modelenergy(rosetta);
  vector<float> rsd_energies;
  const string score_title = "rsd_energy";
  for(map<string, Pose>::iterator mit = poses.begin(); mit != poses.end(); mit++) {
    modelenergy.compute_rosetta_energy((*mit).second, rsd_energies);  
    vector<float>::iterator it1 = rsd_energies.begin();
    vector<ResidueProfile>::iterator it2 = residues_profiles.begin();
    while(it1 != rsd_energies.end() && it2 != residues_profiles.end()) {
      it2->append_score(mit->first, score_title, (*it1));    
      it1++; it2++;
    }
    rsd_energies.clear();
  }
  vector<ResidueProfile>::iterator it2 = residues_profiles.begin();
  for(; it2 != residues_profiles.end(); it2++) {
    it2->show();
  }
}

void compute_local_carmsd(vector<ResidueProfile> &rsd_profiles, 
                          map<string, Pose> poses, Pose native_pose, 
                          scoreTitle score_title, const size_t frag_size) {

  string score_title1 = "", score_title2 = "", score_title3 = "";
  if(score_title == ap_rds) {
    score_title1 = "AP_RDS"; score_title2 = "AP_LENGTH"; score_title3 = "LCA_RMSD";
  }
  else if(score_title == rsd_energy) {
    score_title1 = "rsd_energy"; score_title2 = "rsd_length"; score_title3 = "lrsd_carmsd";
  }

  vector<ResidueProfile>::iterator it;
  for(it  = rsd_profiles.begin(); it != rsd_profiles.end(); it++) {
    size_t rsd_pos = it->get_rsd_pos();
    map<string, float> model_scores = it->get_modelnscore(score_title2);
    map<string, float>::iterator mit =  model_scores.begin();
    size_t total_rsd = poses[mit->first].total_residue();

    if((total_rsd - rsd_pos) < frag_size) break;
    if(!model_scores.empty()) {
      for( mit = model_scores.begin(); mit != model_scores.end(); mit++) {
        Pose pose = poses[mit->first];
        float lcarmsd = CA_rmsd(native_pose, pose, rsd_pos, rsd_pos + frag_size - 1);
        it->append_score(mit->first, score_title3, lcarmsd); 
      }
    }
  }
  dump_rsd_scores(rsd_profiles, "ca_rmsd.out");
  dump_analyzed_scores(rsd_profiles,""); 
  dump_top_rsd_scores(rsd_profiles,"", 0.20); 
} 


vector<float> get_fragments(vector<float> all_coords, const size_t rsd_pos, const size_t frag_size) {
  vector<float> xyz_frag;
  xyz_frag.reserve(rsd_pos*frag_size);
  if(!all_coords.empty()) {
    vector<float>::iterator fvit      = all_coords.begin() + 3*(rsd_pos - 1); 
    vector<float>::iterator frag_fvit = xyz_frag.begin();
    xyz_frag.insert(frag_fvit, fvit, fvit + (3*frag_size));

    //debugging...
    // cout<<"Before ";
    // for(int i = 0; i< xyz_frag.size(); i++) {
    //   cout<<xyz_frag[i]<<" ";
    //   if( (i + 1) % 3 == 0)
    //     cout<<endl<<"Before ";
    // }
    // cout<<endl;

    const size_t mLen = xyz_frag.size() / 3;
    float cx=0.0, cy=0.0, cz=0.0;
    int k3 = 0;

    for(size_t i = 0; i < mLen; i++) {
      cx += xyz_frag[k3];
      cy += xyz_frag[k3 + 1];
      cz += xyz_frag[k3 + 2];
      k3 += 3;
    }
    cx/=mLen; cy/=mLen; cz/=mLen; k3 = 0; 
    for(size_t i = 0; i < mLen; i++) {
      xyz_frag[k3]     -= cx;
      xyz_frag[k3 + 1] -= cy;
      xyz_frag[k3 + 2] -= cz;
      k3+=3;
    }
  }
  return xyz_frag;
}

vector<float> get_local_rsd_coord(Pose ref_pose, const size_t cur_pos, const size_t f_size) {
  using namespace RosettaWrapper;
  vector<float> ref_xyz, local_ref_xyz;
  get_ca_xyz_using_rosetta(ref_pose, core::scoring::is_protein_CA, ref_xyz);
  if(ref_xyz.empty()) return local_ref_xyz;
  local_ref_xyz = get_fragments(ref_xyz, cur_pos, f_size); 
  return local_ref_xyz;
}

void sort_by_map_value(map<int, float> scores, vector<t_pair> &sorted_map) {
  if(!sorted_map.empty()) sorted_map.clear();
  vector<t_pair> vector_pair(scores.begin(), scores.end());
  sort(vector_pair.begin(), vector_pair.end(), comp());
  sorted_map = vector_pair;
}

void compute_fragments_carmsd(vector<float> local_ref_xyz, Pose model_pose, 
                          const size_t f_size, const float threshold_rmsd, 
                          Fragment &fragment, const string pdb_id, ofstream &f_output) {

  using namespace RosettaWrapper;
  vector<float> mod_xyz;
  map<int, float> start_rsd_rmsd;
  get_ca_xyz_using_rosetta(model_pose, core::scoring::is_protein_CA, mod_xyz);

  size_t rsd_pos = 0, total_rsd = model_pose.total_residue();
  assert(total_rsd==mod_xyz.size() / 3);
  while(++rsd_pos <= (total_rsd - f_size)) {
    vector<float> mod_xyz_frag = get_fragments(mod_xyz, rsd_pos, f_size); 
    float ca_rmsd = ca_rmsd_local(local_ref_xyz, mod_xyz_frag); 
    if (ca_rmsd < threshold_rmsd) 
      start_rsd_rmsd[rsd_pos] = ca_rmsd;
    mod_xyz_frag.clear(); 
  } 

  vector<t_pair> sorted_ca_rmsd;
  sort_by_map_value(start_rsd_rmsd, sorted_ca_rmsd);  
  vector<t_pair>::iterator it; //= sorted_ca_rmsd.begin();

  //debugging 
  for(it = sorted_ca_rmsd.begin(); it != sorted_ca_rmsd.end(); it++) {
    cout<<"ROS::CA-RMSD "<<pdb_id<<", "<<it->second <<", "<<it->first<<endl;
  }
  if(!sorted_ca_rmsd.empty()) {
    it = sorted_ca_rmsd.end() - 1;
    float min_ca_rmsd = it->second;
    int   sel_rsd_pos = it->first;
    cout<<"ROS::SEL CA_RMSD "<<sel_rsd_pos<<" "<<min_ca_rmsd<<endl;
    fragment.write_rosetta_fragment_format(f_output, model_pose, pdb_id, sel_rsd_pos, min_ca_rmsd);
  }

}

void generate_fragments_with_native(vector<ResidueProfile> &residues_profiles, 
                        const string native_path, const string path2pdbs, 
                        map<string, Pose> topNPoses, const size_t frag_size, const float ca_rmsd) {
  const size_t max_frags = 200, frag_per_residue = 25;

  vector<string> path2models;
  FileIO fileio;
  fileio.read_model_path(path2pdbs, path2models);
  vector<Fragment> fragments;

  Pose native_pose; string native_model_id;
  if(!fileio.read_pose_no_stout(native_path, native_model_id, native_pose)) {
    cout<<"Native structre path is not setup."<<endl;
    return;
  }

  ofstream f_output;
  f_output.open("fragments.frag9", ios::out);
  const size_t total_residue = residues_profiles.size();
  vector<ResidueProfile>::iterator it = residues_profiles.begin(); 
  for(;it != residues_profiles.end(); it++) {
    size_t rsd_pos = it->get_rsd_pos();
    if((rsd_pos + frag_size) > total_residue) return;

    f_output<<" position: "; f_output.width(12);f_output<<rsd_pos;
    f_output<<" neighbors: "; f_output.width(12);f_output<<frag_per_residue<<endl<<endl;
    Fragment fragment(rsd_pos, frag_size, max_frags, frag_per_residue);
    // for(size_t i = 0; i < 25; i++) {
    // fragment.get_native_fragment(f_output, native_pose, native_model_id, rsd_pos, ca_rmsd); 
    // ca_rmsd does not matter here.
    //}

    vector<string> model_ids = it->get_model_ids();
    if(model_ids.empty()) continue;

    for(vector<string>::iterator vsit = model_ids.begin(); vsit != model_ids.end(); vsit++) {
      if(fragment.get_fragment_per_model() > frag_per_residue) { 
        fragment.reset_fragment_per_model();
        break;
      }

      Pose modelPose; string model_id; 
      vector<float> ref_coords = get_local_rsd_coord(topNPoses[*vsit], rsd_pos, frag_size); 

      for(vector<string>::iterator vit  = path2models.begin(); vit != path2models.end(); vit++) {
        if(fileio.read_pose_no_stout(*vit, model_id, modelPose)) {
          if(!modelPose.empty())  {
            compute_fragments_carmsd(ref_coords, modelPose, frag_size, ca_rmsd, fragment, model_id, f_output); 
          }
        }
      }
    }
  }
  f_output.close();
}

void dumpout_rsd_profile(vector<ResidueProfile> residues_profiles) {
  // ofstream f_output;
  // f_output("rsd_profile.out");
  // for(vector<ResidueProfile>::iterator it = residues_profiles.begin(); 
  //                                      it != residues_profiles.end(); 
  //                                      it++) {
  //   size_t rsd_pos = it->get_rsd_pos();
  //   vector<float> aprsd_scores=ResidueProfile::get_score("AP-RDS"); 
  //   for(vector<float>::iterator it= aprsd.begin(); it != aprsd.end(); it++) {
  //     f_output.width(9); f_output.precision(3); f_output<<*it<<" ";
  //   }
  //   f_output<<endl;
  // }
}

void generate_fragments(vector<ResidueProfile> &residues_profiles, 
                        const string path2pdbs, map<string, Pose> topNPoses, 
                        const size_t frag_size, const float ca_rmsd) {

  const size_t max_frags = 200, frag_per_residue = 25;
  vector<string> path2models;
  FileIO fileio;
  fileio.read_model_path(path2pdbs, path2models);
  vector<Fragment> fragments;

  ofstream f_output;
  f_output.open("fragments.frag9", ios::out);

  for(vector<ResidueProfile>::iterator it = residues_profiles.begin(); it != residues_profiles.end(); it++) {
    size_t rsd_pos = it->get_rsd_pos();
    f_output<<" position: "; f_output.width(12);f_output<<rsd_pos;
    f_output<<" neighbors: "; f_output.width(12);f_output<<frag_per_residue<<endl<<endl;
    Fragment fragment(rsd_pos, frag_size, max_frags, frag_per_residue);

    vector<string> model_ids = it->get_model_ids();
    if(model_ids.empty()) continue;

    for(vector<string>::iterator vsit = model_ids.begin(); vsit != model_ids.end(); vsit++) {
      if(fragment.get_fragment_per_model() > frag_per_residue) { 
        fragment.reset_fragment_per_model();
        break;
      }

      Pose modelPose; string model_id;
      vector<float> ref_coords = get_local_rsd_coord(topNPoses[*vsit], rsd_pos, frag_size); 

      for(vector<string>::iterator vit  = path2models.begin(); vit != path2models.end(); vit++) {
        if(fileio.read_pose_no_stout(*vit, model_id, modelPose)) {
          if(!modelPose.empty())  {
            compute_fragments_carmsd(ref_coords, modelPose, frag_size, ca_rmsd, fragment, model_id, f_output); 
          }
        }
      }
    }
  }
  f_output.close();
}

// void verbose_rsd_score();

  
#endif 
