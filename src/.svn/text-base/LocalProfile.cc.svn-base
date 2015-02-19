#include<LocalProfile.hh>


using namespace std;

Score::Score() {
  sliding_after_normal_ = false; 
}

Score::Score(string title, float score) {
  scores_[title] = score;
  sliding_after_normal_ = false; 
}

map<string, float> Score::get_scores() {
  return scores_;
}

float Score::get_score(const string title) {
  if(scores_.empty())
    return -1.0;
  else 
    return scores_[title];
}

size_t Score::count(const string title) {
  if(title.empty()) return 0;
  return size_t(scores_.count(title));
}

void Score::set_sliding_after_normalization(bool value) {
  sliding_after_normal_ = value; 
}

void Score::append(string title, float score) {
  scores_[title] = score;
}

void Score::update(const string title, float value) {
  scores_[title] = value;
}

void Score::show() {
  if(scores_.empty()) return; 
  map<string, float>::iterator it; 
  for( it = scores_.begin(); it != scores_.end(); it++) {
    cout<<" "<<(*it).first <<": "<<(*it).second<<" ";    
  }
}

void Score::write_rsd_scores(ofstream &f_output_score,
                             const string title) {
  if(scores_.empty()) return; 
  map<string, float>::iterator it = scores_.begin(); 
  if(title == "" || title == "all" || title == "ALL") {
    for(; it != scores_.end(); it++) {
      f_output_score.setf(ios::fixed, ios::floatfield); 
      f_output_score.width(7); f_output_score.precision(2); 
      f_output_score<<(*it).second;    
    }
    f_output_score<<endl;
    return;
  }
  for(; it != scores_.end(); it++) {
    if(it->first == title || it->first == "FRAGMENT_RMSD") {
      f_output_score.setf(ios::fixed, ios::floatfield); 
      f_output_score.width(7); f_output_score.precision(2); 
      f_output_score<<(*it).second;    
    }
  }
  f_output_score<<endl;
}


void Score::write_score_title(ofstream &f_output_score) {
  if(scores_.empty()) return; 
  map<string, float>::iterator it; 
  for( it = scores_.begin(); it != scores_.end(); it++) 
    f_output_score<<it->first<<", ";    
  f_output_score<<endl;
}

Score Score::operator/=(float win_size) {
  if( not sliding_after_normal_ ) { 
    scores_["SLIDING_MEAN_AP_RDS"] /= win_size; 
    scores_["SLIDING_MEAN_AP_STD"] /= win_size; 
  }
  else if (sliding_after_normal_) {
    scores_["NORMALIZED_SLIDING_MEAN_AP_RDS"] /= win_size; 
    scores_["NORMALIZED_SLIDING_MEAN_AP_STD"] /= win_size; 
  }
  return *this;
}
