#ifndef SPMAIN_HEADER_HH 
#define SPMAIN_HEADER_HH 

#include<iostream>
#include<vector>
#include<map>

#include<Residue.hh>
#include<ModelEnergy.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>


typedef core::pose::Pose Pose;

using namespace std;
using namespace core::scoring;

enum scoreTitle { 
  ap_rds,
  ap_std,
  ca_rmsd,
  rsd_energy
};

//just for this one... has to modify
vector<ResidueProfile> run(vector<ResidueProfile> residues_profiles, 
                        const size_t topN, size_t minLen, scoreTitle score_title) {
  string score_title1 = "", score_title2 = "";
  if(score_title == ap_rds) {
    score_title1 = "ap_rds";
    score_title2 = "ap_length";
  }
  else if(score_title == rsd_energy) {
    score_title1 = "rsd_energy";
    score_title2 = "rsd_length";
  }

  cout<<"Size of vector "<<residues_profiles.size()<<endl;
  vector<ResidueProfile>::iterator it;
  vector<ResidueProfile> sorted_residues;
  sorted_residues.reserve(residues_profiles.size());

  for(it = residues_profiles.begin(); it != residues_profiles.end(); it++) {
    ResidueProfile rsd_profile = it->sort_by_score(score_title1, topN); 
    sorted_residues.push_back(rsd_profile);
  }
  cout<<topN<<" models are selected for db." <<endl;
  vector<ResidueProfile>::iterator outer_it, inner_it;
  for(outer_it = sorted_residues.begin(); outer_it != sorted_residues.end(); outer_it++) {
    map<string, float> model_scores = outer_it->get_model_n_sorted_score(score_title1, topN);
    for(map<string, float>::iterator mit = model_scores.begin();
        mit != model_scores.end(); mit++) {
      size_t counter = 0; 
      const string model_id = mit->first;
      for(inner_it = outer_it + 1; inner_it != sorted_residues.end(); inner_it++) {
        if(inner_it->is_model_exist(model_id, score_title1)) counter++;
        else { 
          if(counter>=minLen) break;
          else counter++;
        }
      }
      if(counter > 2) outer_it->append_score(model_id, score_title2, float(counter)); 
    }
    cout<<endl;
  }
  cout<<endl;
  for(outer_it = sorted_residues.begin(); outer_it != sorted_residues.end(); outer_it++) {
    outer_it->show();
  }
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
  for(vector<ResidueProfile>::iterator it2 = residues_profiles.begin();
      it2 != residues_profiles.end(); it2++) {
    it2->show();
  }
}

void compute_local_carmsd(vector<ResidueProfile> &residues_profiles, 
                          map<string, Pose> poses, Pose native_pose, 
                          scoreTitle score_title) {
  string score_title1 = "", score_title2 = "", score_title3 = "";
  if(score_title == ap_rds) {
    score_title1 = "ap_rds";
    score_title2 = "ap_length";
    score_title3 = "lap_carmsd";
  }
  else if(score_title == rsd_energy) {
    score_title1 = "rsd_energy";
    score_title2 = "rsd_length";
    score_title3 = "lrsd_carmsd";
  }

  for(vector<ResidueProfile>::iterator it = residues_profiles.begin();
      it != residues_profiles.end(); it++) {
    size_t rsd_pos = it->get_rsd_pos();
    map<string, float> model_scores = it->get_modelnscore(score_title2);
    if(!model_scores.empty()) {
      for(map<string,float>::iterator mit = model_scores.begin(); 
          mit != model_scores.end(); mit++) {
        Pose pose = poses[mit->first];
        size_t length = size_t(mit->second);
        float lcarmsd = CA_rmsd(native_pose, pose, rsd_pos, rsd_pos + length);
        it->append_score(mit->first, score_title3, lcarmsd); 
      }
    }
  }

  for(vector<ResidueProfile>::iterator it2 = residues_profiles.begin();
      it2 != residues_profiles.end(); it2++) {
    it2->show();
  }
} 

//ROS:: this function is no longer useful in future it is very rough funtion no thread safe and error free.
// void compute_min_max_avg(vector<ResidueProfile> residues_profiles, const string score_title) {
//   for(vector<ResidueProfile>::iterator it = residues_profiles.begin();
//       it != residues_profiles.end(); it++) {
//     vector<float> scores = it->get_score(score_title);
//     const size_t max_index = scores.size();
//     sort(scores.begin(), scores.end());
//     float min  = scores[0];
//     float max  = scores[max_index-1];
//     float mean = 0;
//     float sum  = 0; 
//     for(size_t i = 0; i < scores.size(); i++) {
//       sum += scores[i];
//     }
//     mean = sum / max_index;
//     cout<<"Rsd pos: " <<it->get_rsd_pos()<<" min "<<min<<" max "<<max<<" mean "<<mean<<endl;
//   }
// }

#endif 
