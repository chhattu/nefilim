#include<Residue.hh>

#include<cmath>
#include<cassert>
#include<algorithm>
#include<functional>
#include<numeric>
#include<cassert>
#include<exception>

#include <core/scoring/rms_util.hh>
#include <ObjexxFCL/format.hh>

using namespace core::scoring;
using namespace protocols::moves;

using namespace core;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace protocols;

ResidueProfile::ResidueProfile() {
  slide_after_normalization_ = false; 
}

ResidueProfile::ResidueProfile(size_t rsd_pos) {
  rsd_pos_ = rsd_pos;
  slide_after_normalization_ = false; 
}

ResidueProfile::ResidueProfile(size_t rsd_pos, 
                               string model_id, 
                               string score_title, 
                               float value) {
  rsd_pos_ = rsd_pos;
  Score score(score_title, value);
  local_scores_[model_id] = score; 
  slide_after_normalization_ = false; 
}

ResidueProfile::ResidueProfile(size_t rsd_pos, 
                               string score_title,
                               vector<string> modelnames, 
                               vector<float> values) {
  rsd_pos_ = rsd_pos;
  vector<string>::iterator sit = modelnames.begin();
  vector<float>::iterator fit  = values.begin();
  while(sit != modelnames.end() && fit != values.end()) {
    string modelname = *sit;
    float value = *fit;
    Score score(score_title, value);
    local_scores_[modelname] = score; 
    sit++; fit++;
  }
  slide_after_normalization_ = false; 
}

void ResidueProfile::set_total_residue(size_t total_residue) {
  total_residue_ = total_residue;
}

void ResidueProfile::set_selection_mode(size_t mode, size_t smode) {
  mode_  = mode;
  smode_ = smode;
}

void ResidueProfile::set_selection_number(size_t number_of_fragment, 
                                          size_t fragment_per_cluster,
                                          size_t number_of_templates) {
  number_of_fragment_   = number_of_fragment;
  fragment_per_cluster_ = fragment_per_cluster;
  number_of_templates_  = number_of_templates; 
}


size_t ResidueProfile::get_rsd_pos() {
  return rsd_pos_;
}

void ResidueProfile::set_score(map<string, Score> local_scores) {
  local_scores_ = local_scores;
}

void ResidueProfile::set_slide_after_normalization(bool value) {
  slide_after_normalization_ = value; 
}

bool ResidueProfile::is_model_exist(const string model_id, const string score_title) {
  if (local_scores_.count(model_id) == 1) {
    Score scores = local_scores_[model_id]; 
    if(scores.count(score_title) == 1) 
      return true;
    else 
      return false;
  }
  return false;
}

vector<float> ResidueProfile::get_score(string score_title) {
  vector<float> scores;
  scores.clear();
  map<string, Score>::iterator mit  = local_scores_.begin(); 
  for(; mit != local_scores_.end(); mit++) {
    scores.push_back(mit->second.get_score(score_title));
  } 
  return scores;
}

void ResidueProfile::get_score(string score_title, list<float> &scores) {
  scores.clear();
  map<string, Score>::iterator mit  = local_scores_.begin(); 
  for(; mit != local_scores_.end(); mit++)
    scores.push_back(mit->second.get_score(score_title));
}

bool ResidueProfile::get_score( const string score_title, 
                                map<string, float> &models_n_scores) {
  models_n_scores.clear();
  if(score_title.empty()) return false;

  map<string, Score>::iterator mit  = local_scores_.begin();
  float sum = 0.0;
  for(; mit != local_scores_.end(); mit++) { 
    float score = mit->second.get_score(score_title);
    models_n_scores[mit->first] = score; 
    sum += score; 
  }
  if(models_n_scores.size() <= 0) return false;
  if(sum <= 0) return false;
  return true;
}

size_t ResidueProfile::get_size() {
  return local_scores_.size();
}


string inline ResidueProfile::get_modelname(const string score_title, float aprds) {
  string modelname;
  const float eplison = 0.1;
  map<string, Score>::iterator mit = local_scores_.begin(); 
  for(; mit != local_scores_.end(); mit++) {
    // cout<<"APRDS "<<mit->second.get_score(score_title)<<", "<<aprds<<endl;
    //mistake in score title
    if(fabs(mit->second.get_score(score_title)- aprds) <= eplison) {
      cout<<"APRDS is found "<<mit->first<<" "<<aprds<<endl;
      return mit->first;
    }
  } 
  return modelname; 
}

// this algorithm can be changed in more efficent way 
vector<string> ResidueProfile::get_modelnames(const string score_title, vector<float> aprdss) {
  vector<string> modelnames; 
  modelnames.reserve(aprdss.size());
  for(vector<float>::iterator it = aprdss.begin(); it != aprdss.end(); it++) {
    const string modelname = get_modelname(score_title, (*it));
    if(!modelname.empty()) 
      modelnames.push_back(modelname);
  }
  return modelnames; 
}

float ResidueProfile::get_score(const string score_title, const string modelname) {
  float score = -1.0;
  map<string, Score>::iterator mit = local_scores_.begin(); 
  for(; mit != local_scores_.end(); mit++) {
    //can be replaced this code with two more lines and efficient
    if(mit->first == modelname) {
      score = mit->second.get_score(score_title);
      break;
    }
  } 
  return score; 
}

map<string, float> ResidueProfile::get_modelnscore(string title) {
  map<string, float> scores; // <model_id, score> 
  scores.clear();
  if(title.empty()) return scores;

  map<string, Score>::iterator its;
  for(its = local_scores_.begin(); its != local_scores_.end(); its++) {
    scores[(*its).first] = its->second.get_score(title);   
  } 

  return scores;
}


//ROS:: will be changed such that searching will be done in efficient 
//way. For both are sorted at first and then searched in same manner 
// which does not need so many iteration. Right now is just brute
//force way.


map<string, float> ResidueProfile::get_model_n_sorted_score(string title, size_t topN) {
  map<string, float> scores; // <model_id, score> 
  scores.clear();
  if(title.empty()) return scores;
  map<string, Score>::iterator its;
  for(its = local_scores_.begin(); its != local_scores_.end(); its++) {
    scores[(*its).first] = its->second.get_score(title);   
  } 

  vector<t_pair> sorted_scores; 
  sort_by_value(scores, sorted_scores); 
  if(topN > sorted_scores.size() || topN < 1) {
    cout<<"warning: no proper value for top N."<<endl;
    topN = sorted_scores.size();
  }
  map<string, float> selected_scores; 
  for(vector<t_pair>::iterator it = sorted_scores.begin(); 
      it != sorted_scores.begin() + topN; it++) {
    selected_scores[it->first] = it->second;
  }
  return selected_scores;
}

void ResidueProfile::append_score(const string model_id, 
                                  string score_title, 
                                  float value) {
  if(local_scores_.count(model_id) > 0) {
    Score score = local_scores_[model_id]; 
    score.append(score_title, value); 
    local_scores_[model_id] = score; 
    // find_min_max_scores(score_title, value);
  } 
}

void ResidueProfile::update_score(const string score_title, 
                                  map<string, float> values) {
  map<string, Score>::iterator mit1 = local_scores_.begin();
  map<string, float>::iterator mit2 = values.begin();

  while(mit2 != values.end()) {
    mit1 = local_scores_.find(mit2->first);
    if(mit1 == local_scores_.end()) continue; 
    mit1->second.update(score_title, mit2->second);
    mit2++;
  }
}

void ResidueProfile::update(size_t rsd_pos, 
                            string model_id, 
                            string score_title, 
                            float value) {
  if(rsd_pos_ != rsd_pos) return; 
  if(local_scores_.count(model_id) > 0) {
    append_score(model_id, score_title, value);
  }
  else {
    Score score(score_title, value); 
    local_scores_[model_id] = score;
  }    
}

bool ResidueProfile::sort_by_value(const t_pair m1, const t_pair m2) {
  return m1.second < m2.second;
}

void ResidueProfile::sort_by_value(map<string, float> scores, vector<t_pair> & t_vector) {
  if(!t_vector.empty()) t_vector.clear();
  if(scores.empty()) {
    t_vector.clear();
    return;
  }
  vector<t_pair> s_vector(scores.begin(), scores.end());
  sort(s_vector.begin(), s_vector.end(), comp());
  t_vector = s_vector;
}

ResidueProfile ResidueProfile::sort_by_score(const string score_title, size_t topN){
  vector<t_pair> sorted_scores; 
  sort_by_value(get_modelnscore(score_title), sorted_scores); 

  //debugging...
  // cout<<"Inside sorted "<<endl;
  // for(vector<t_pair>::iterator vit=sorted_scores.begin(); 
  //                              vit != sorted_scores.end(); vit++) {
  //   cout<<vit->second<<" "<<vit->first<<endl;
  // }

  ResidueProfile rsd_profile(rsd_pos_);
  if(sorted_scores.size() < 1) return rsd_profile;
  if(topN < 1 || topN > sorted_scores.size()) {
    cout<<"warning: no proper value for top N."<<endl;
    topN = sorted_scores.size();
  }

  map<string, Score> local_scores;
  vector<t_pair>::iterator vit=sorted_scores.begin();
  for(; vit != sorted_scores.begin() + topN; vit++) 
    local_scores[vit->first] = local_scores_[vit->first];
  rsd_profile.set_score(local_scores);    
  return rsd_profile;
}

ResidueProfile ResidueProfile::sort_by_score(const string score_title, float cutoff){
  vector<t_pair> sorted_scores; 
  sort_by_value(get_modelnscore(score_title), sorted_scores); 

  //debugging...
  // cout<<"Inside sorted "<<endl;
  // for(vector<t_pair>::iterator vit=sorted_scores.begin(); 
  //                              vit != sorted_scores.end(); vit++) {
  //   cout<<vit->second<<" "<<vit->first<<endl;
  // }

  ResidueProfile rsd_profile(rsd_pos_);
  if(sorted_scores.size() < 1) return rsd_profile;
  map<string, Score> local_scores;
  vector<t_pair>::iterator vit=sorted_scores.begin();
  for(; vit != sorted_scores.end(); vit++)  {
    float score = vit->second;
    if(score > cutoff) break; 
    local_scores[vit->first] = local_scores_[vit->first];
  }
  if(local_scores.size() < number_of_fragment_) { //change this 25 values
    cerr<<"fatal: "<<__FILE__<<" "<<__LINE__
        <<" # models are not less than "<<number_of_fragment_
        <<" for residue "<<rsd_pos_<<endl;
    exit(1);
  }

  rsd_profile.set_score(local_scores);    
  return rsd_profile;
}

ResidueProfile ResidueProfile::sort_by_score(const string score_title, 
                                             float cutoff, 
                                             bool &is_not_increased){
  vector<t_pair> sorted_scores; 
  sort_by_value(get_modelnscore(score_title), sorted_scores); 
  is_not_increased = true;

  //debugging...
  // cout<<"Inside sorted "<<endl;
  // for(vector<t_pair>::iterator vit=sorted_scores.begin(); 
  //                              vit != sorted_scores.end(); vit++) {
  //   cout<<vit->second<<" "<<vit->first<<endl;
  // }

  ResidueProfile rsd_profile(rsd_pos_);
  if(sorted_scores.size() < 1) return rsd_profile;
  map<string, Score> local_scores;
  vector<t_pair>::iterator vit=sorted_scores.begin();
  for(; vit != sorted_scores.end(); vit++)  {
    float score = vit->second;
    if(score > cutoff) break; 
    local_scores[vit->first] = local_scores_[vit->first];
  }
  if(local_scores.size() < number_of_fragment_) { //change this 25 values
    cerr<<"  warning: "<<__FILE__<<" "<<__LINE__
        <<" # models are less than "<<number_of_fragment_
        <<" for residue "<<rsd_pos_<<endl;
    is_not_increased = false;
    return rsd_profile;
  }

  rsd_profile.set_score(local_scores);    
  return rsd_profile;
}


void ResidueProfile::show() {
  if(local_scores_.empty()) return;
  cout<<"rsd pos: "<<rsd_pos_<<endl;
  map<string, Score>::iterator it;
  for(it = local_scores_.begin(); it != local_scores_.end(); it++) {
    cout<<"   "<<(*it).first<<" "; (*it).second.show();
    cout<<endl;
  }
}

void ResidueProfile::write_rsd_scores(ofstream &f_output_score,
                                      const string title) { 
  if(local_scores_.empty()) return;
  f_output_score.width(5); 
  f_output_score<<"Residue number: "<<rsd_pos_<<endl;
  map<string, Score>::iterator it;
  for(it = local_scores_.begin(); it != local_scores_.end(); it++) {
    f_output_score<<"   "<<it->first<<" "; 
    it->second.write_rsd_scores(f_output_score, title);
  }
  f_output_score<<endl;
}

void ResidueProfile::write_score_title(ofstream &f_output_score) {
  if(local_scores_.empty()) return;
  map<string, Score>::iterator it = local_scores_.begin();
  for(; it != local_scores_.end(); it++) {
    it->second.write_score_title(f_output_score);
    break;
  }
}

void ResidueProfile::write_cluster_info(ofstream &outfile, 
                                        FragmentCluster clustering_) {
  map<int, vector<string> > cluster_info = clustering_.get_clusters_info();
  map<int, vector<string> >::iterator mvit = cluster_info.begin();
  size_t total_num_model = 0;
  outfile<<"clustering for residue position: "<<rsd_pos_<<std::endl;
  for(; mvit != cluster_info.end(); mvit++) {
    vector<string> members = mvit->second;
    vector<string>::iterator vit = members.begin();
    outfile<<"cluster "; outfile.width(5);  outfile<<mvit->first;
    outfile<<"  members "; outfile.width(5); outfile<<members.size()<<endl;
    total_num_model += members.size();
    for(; vit != members.end(); vit++) { 
      outfile<<"  "<<*vit<<"  ";
      outfile.setf(ios::fixed,ios::floatfield);
      outfile.width(7); outfile.precision(3); 
      outfile<<clustering_.get_distance2center(*vit);
      float rmsd2native = get_score("FRAGMENT_RMSD", *vit);
      if(rmsd2native > -1.0) 
        outfile.width(7); outfile.precision(3); outfile<<rmsd2native;
      outfile<<endl;
    }
  }
  outfile<<"cluster density rsd pos: ";
  outfile.width(5); outfile<<rsd_pos_;
  outfile.width(5); outfile<<" cluster size: "<<cluster_info.size();
  outfile.width(5); outfile<<" # models: "<<total_num_model<<std::endl;
}

void ResidueProfile::write_selected_templates_native_RMSD(ofstream &outfile,
                                                          FragmentCluster clustering_) {
  map<string, Pose> selected_templates = clustering_.get_selected_templates(); 
  vector<string> template_ids  = clustering_.get_selected_templates_id(); 
  vector<string>::iterator vit = template_ids.begin();
  outfile<<"templates for residue positon: "<<rsd_pos_
         <<" "<<template_ids.size()<<endl;
  for(; vit != template_ids.end(); vit++) {
    float rmsd2native = get_score("FRAGMENT_RMSD", *vit);
    if(rmsd2native > -1.0)  {
      outfile<<*vit; outfile.width(7); outfile.precision(3); 
      outfile<<rmsd2native; outfile<<endl;
    }
  }
}


void ResidueProfile::dump_average_rmsd_to_native(const string score_title, 
                                                 ofstream &f_output_score) {
  list<float> scores;
  get_score(score_title, scores);
  if(scores.empty()) return;
  const int num_of_data = scores.size();
  scores.sort();
  list<float>::iterator lit = scores.begin();
  float min = *lit;
  lit = scores.end();
  float max = *(--lit);
  float avg = accumulate(scores.begin(), scores.end(), 0.0) / num_of_data;
  f_output_score<<"Residue: ";
  f_output_score.width(5); f_output_score<<rsd_pos_;
  f_output_score.setf(ios::fixed, ios::floatfield); 
  f_output_score.width(7); f_output_score.precision(2); f_output_score<<" "<<min<<" ";
  f_output_score.width(7); f_output_score.precision(2); f_output_score<<" "<<max<<" ";
  f_output_score.width(7); f_output_score.precision(2); f_output_score<<" "<<avg<<" ";
  f_output_score.width(3); f_output_score.precision(0); f_output_score<<" "<<scores.size()<<" ";
  if((scores.size()) < 25) f_output_score<<"** ";
  if((scores.size()) > 25) f_output_score<<"***";
  if((scores.size()) == 25) f_output_score<<"---";
  f_output_score<<endl;
}

vector<string> ResidueProfile::get_model_ids() {
  vector<string> model_ids;
  model_ids.clear();
  if(local_scores_.empty()) return model_ids;
  model_ids.reserve(local_scores_.size());
  map<string, Score>::iterator mit  = local_scores_.begin(); 

  for(; mit != local_scores_.end(); mit++) 
    model_ids.push_back(mit->first); 
  return model_ids;
}

vector<string> ResidueProfile::get_score_titles() {
  vector<string> score_titles;
  map<string, Score>::iterator mit = local_scores_.begin(); 
  Score score = mit->second;
  map<string, float> scores = score.get_scores();
  map<string, float>::iterator  mit1 = scores.begin();
  for(; mit1 !=  scores.end(); mit1++)
    score_titles.push_back(mit1->first);
  return score_titles;
}

ResidueProfile ResidueProfile::operator+=(ResidueProfile rsd) {
  map<string, float> values1; 
  values1.clear();
  if( not slide_after_normalization_ ) { 
    add_scores("AP_RDS", values1, rsd); 
    update_score("SLIDING_MEAN_AP_RDS", values1);
    values1.clear();
    add_scores("AP_STD", values1, rsd); 
    update_score("SLIDING_MEAN_AP_STD", values1);
  }
  else {
    add_scores("NORMALIZED_AP_RDS", values1, rsd); 
    update_score("NORMALIZED_SLIDING_MEAN_AP_RDS", values1);
    values1.clear();
    add_scores("NORMALIZED_AP_STD", values1, rsd); 
    update_score("NORMALIZED_SLIDING_MEAN_AP_STD", values1);
  }

  return *this;
}

ResidueProfile ResidueProfile::operator/=(const float win_size) {
  map<string, Score>::iterator mit = local_scores_.begin(); 
  while(mit != local_scores_.end()) { 
    mit->second /= win_size;
    mit++; 
  }
  return *this;
}

bool ResidueProfile::normalize_score(const float min, const float max, string score_title) {
  try { 
    const string new_score_title = "GLOBAL_NORMALIZED_" + score_title;
    map<string, float> scores = get_modelnscore(score_title);
    const int num_of_data = scores.size(); 
    if(num_of_data <= 0) throw " no data for normalization.";
    if(max < min) throw " unequal possibilities.";

    map<string, float>::iterator mit;
    for(mit = scores.begin(); mit != scores.end(); mit++) { 
      float nor_value = (mit->second - min) / (max- min);
      append_score(mit->first, new_score_title, nor_value);
    }
    return true;
  }
  catch(const char* e) {
    std::cerr<<"warning: "<<__FILE__<<" "<<__LINE__<<" "<<e<<std::endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"warning: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<std::endl;
    exit(1);
  }
}

bool ResidueProfile::normalize_score(string score_title) {
  const string new_score_title = "NORMALIZED_" + score_title;
  map<string, float> scores = get_modelnscore(score_title);
  const int num_of_data = scores.size();
  if(num_of_data <= 0) return false;

  list<float> values; 
  map<string, float>::iterator mit;
  for(mit = scores.begin(); mit != scores.end(); mit++) 
    values.push_back(mit->second);
  
  values.sort();
  list<float>::iterator lit = values.begin(); 
  float min_value = *lit;
  lit = values.end(); 
  float max_value = *(--lit); 

  min_values[score_title] = min_value;
  max_values[score_title] = max_value;
  // cout<<score_title<<"  "<<min_value<<"  "<<max_value<<endl;
  if(max_value < min_value) return false;

  for(mit = scores.begin(); mit != scores.end(); mit++) { 
    float nor_value = (mit->second - min_value) / (max_value - min_value);
    append_score(mit->first, new_score_title, nor_value);
  }
  return true;
}

void ResidueProfile::do_clustering(const float cluster_radius, 
                                   map<string, Pose> fragment_poses, 
                                   ofstream &f_output_fragment, 
                                   ofstream &f_output_template, 
                                   ofstream &f_output) {
  try {
    cout<<"  # fragments for clustering residue "<<rsd_pos_
        <<" : "<<fragment_poses.size()<<endl;
    if(fragment_poses.size() < number_of_fragment_) 
      throw " # fragments for clustering is less than required.";
    FragmentCluster clustering_;
    clustering_.set_total_residue(total_residue_);
    clustering_.set_poses(fragment_poses);
    clustering_.do_clustering_using_durandal(cluster_radius); 
    map<string, Pose> selected_templates = clustering_.get_selected_templates(mode_, 
                                                                              smode_,
                                                                              number_of_fragment_,
                                                                              fragment_per_cluster_); 
    write_clustered_templates(f_output_fragment, clustering_, selected_templates);
    write_template_fragments_pdb(f_output_template, clustering_);
    write_cluster_info(f_output, clustering_);
    write_selected_templates_native_RMSD(f_output, clustering_);
    //delete after valgrid check 
    selected_fragment_poses_.clear();
  }
  catch(const char* s){
    cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<s<<endl;
    exit(1);
  }
  catch(...) {
    cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<endl;
    exit(1);
  }
}

void ResidueProfile::write_clustered_templates(ofstream &f_output_fragment,
                                               FragmentCluster clustering_, 
                                               map<string, Pose> selected_templates) {
  try {
    if(selected_templates.empty()) 
      throw  " templates are not found for residue. " + rsd_pos_;
    map<string, Pose> selected_extended_templates;
    selected_extended_templates = get_extended_selected_templates(selected_templates); 
    if(selected_templates.empty()) 
      throw  " templates are not found for residue. " + rsd_pos_;
    clustering_.write_fragments(f_output_fragment, selected_extended_templates, rsd_pos_); 
    //delete after valgrid check 
    // selected_fragment_poses_.clear();
  }
  catch(const char* e) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<e<<std::endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<std::endl;
    exit(1);
  }
}

void ResidueProfile::get_standard_fragment_poses(const size_t window_size, 
                                                 map<string, Pose> poses, 
                                                 map<string, Pose> &fragment_poses) {
  vector<string> models_ids = get_model_ids();
  fragment_poses.clear(); // fragment_poses.reserve(models_ids.size());
  vector<string>::iterator vit =  models_ids.begin(); 
  int counter = 0;
  const size_t start_rsd = rsd_pos_;
  const size_t end_rsd = rsd_pos_ + (window_size - 1);
  size_t m_start_rsd, m_end_rsd;
  (start_rsd == 1) ? m_start_rsd = start_rsd : m_start_rsd = start_rsd - 1;
  while(vit != models_ids.end()) {
    Pose pose = poses[*vit];

    freopen("/dev/null", "w", stdout); 
    Pose fragment_pose = get_standard_fragment_pose(pose, start_rsd, end_rsd);
    fragment_poses[*vit] = fragment_pose;

    if(end_rsd == pose.total_residue()) {
      m_end_rsd = end_rsd; 
      m_start_rsd = start_rsd - 2;
    }
    else if(start_rsd == 1) 
      m_end_rsd = end_rsd + 2;
    else
      m_end_rsd = end_rsd + 1;
    // it is just to get one more residue in each side 

    fragment_pose = get_standard_fragment_pose(pose, m_start_rsd, m_end_rsd);
    selected_fragment_poses_[*vit] = fragment_pose;
    fclose(stdout);
    freopen("/dev/tty", "w", stdout); 
    counter++;
    vit++;
  }
}

void ResidueProfile::get_centroid_fragment_poses(const size_t window_size, 
                                                 map<string, Pose> poses, 
                                                 map<string, Pose> & fragment_poses) {
  vector<string> models_ids = get_model_ids();
  vector<string>::const_iterator vit =  models_ids.begin(); 
  fragment_poses.clear();
  while(vit != models_ids.end()) {
    Pose pose = poses[*vit];
    Pose fragment_pose = get_centroid_fragment_pose(pose, rsd_pos_, rsd_pos_ + (window_size - 1));
    fragment_poses[*vit] = fragment_pose;
    vit++;
  }
}

//changed centroid to fa_standard
//have to kick out somewhere else than here...
Pose ResidueProfile::get_centroid_fragment_pose(Pose old_pose, size_t start, size_t end) {
  using namespace core::chemical;    
  using namespace core::io::pdb;    
  Pose new_pose;
  utility::vector1< core::Size > residue_indices;
  for(; start<= end; ++start){
    residue_indices.push_back(start);
  }
  ResidueTypeSetCAP residue_set(ChemicalManager::get_instance()->residue_type_set("centroid"));
  pose_from_pose( new_pose, old_pose, *residue_set,  residue_indices);
  return new_pose;
}

Pose ResidueProfile::get_standard_fragment_pose(Pose old_pose, size_t start, size_t end) {
  using namespace core::chemical;    
  using namespace core::io::pdb;    
  Pose new_pose;
  utility::vector1< core::Size > residue_indices;
  for(; start<= end; ++start)
    residue_indices.push_back(start);
  ResidueTypeSetCAP residue_set(ChemicalManager::get_instance()->residue_type_set("fa_standard"));
  pose_from_pose( new_pose, old_pose, *residue_set,  residue_indices);
  return new_pose;
}

void ResidueProfile::write_template_fragments_pdb(ofstream &out,
                                                  FragmentCluster clustering_) {
  using namespace ObjexxFCL::fmt;
  using namespace core::conformation;
  static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );
  size_t number(0), model_counter(0);
  try {
    map<string, RPose> template_poses = clustering_.get_cluster_center_templates(number_of_templates_); 
    map<string, Pose> selected_extended_templates;
    selected_extended_templates = get_extended_selected_templates(template_poses); 

    map<string, RPose>::iterator mit  = selected_extended_templates.begin();
    for(; mit != selected_extended_templates.end(); mit++) {
      RPose pose = mit->second;
      const size_t nres(pose.total_residue());
      out<<"TER "<<mit->first<<"_"<<++model_counter<<"_"<<rsd_pos_<<endl;
      for(size_t i=1; i<= nres; ++i ) {
        Residue const & rsd( pose.residue(i) );
        for ( size_t j=1; j<= rsd.natoms(); ++j ) {
          Atom const & atom( rsd.atom(j) );
          ++number;
          char const chain( chains[ rsd.chain() ] );
          out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
          rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
          F(8,3,atom.xyz()(1)) <<
          F(8,3,atom.xyz()(2)) <<
          F(8,3,atom.xyz()(3)) <<
          F(6,2,1.0) << F(6,2,1.0) << '\n';
        }
      }
    }
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<endl;
    exit(1);
  }
}


//Private functions
void ResidueProfile::add_scores(const string score_title, 
                                map<string, float> &values1, 
                                ResidueProfile rsd) {
  const string n_score_title = "SLIDING_MEAN_" + score_title;
  map<string, float> values2; 

  if(!get_score(n_score_title, values1)) 
    get_score(score_title, values1); 

  rsd.get_score(score_title, values2); 

  assert(values1.size() == values2.size());
  map<string, float>::iterator mit1 = values1.begin(); 
  map<string, float>::iterator mit2 = values2.begin(); 

  while(mit1 != values1.end()) {
    mit2 = values2.find(mit1->first);
    if(mit2 == values2.end()) continue; 
    values1[mit1->first] = mit1->second + mit2->second;  
    mit1++;
  }
}


map<string, Pose> ResidueProfile::get_extended_selected_templates(
                                  map<string, Pose> selected_templates) {
  map<string, Pose>::iterator msit = selected_templates.begin();
  map<string, Pose> selected_extended_templates;
  for(; msit != selected_templates.end(); msit++) {
    Pose pose = selected_fragment_poses_[msit->first];
    if(not pose.empty()) 
      selected_extended_templates[msit->first] = pose;
  }
  return selected_extended_templates;
}



