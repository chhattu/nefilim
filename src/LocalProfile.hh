#ifndef LOCAL_PROFILE_HEADER_HH 
#define LOCAL_PROFILE_HEADER_HH 

#include<iostream>
#include<fstream>
#include<map>
#include<vector>
#include<algorithm>

using namespace std;

class Score {
  private:
    typedef std::pair<string, float> t_pair;

  private:
    map<string, float> scores_; 
    bool sliding_after_normal_ ;

  public:
    Score();

    Score(string title, float score); 

    map<string, float> get_scores();

    float get_score(string title);

    size_t count(const string title); 

    void set_sliding_after_normalization(bool value);

    void append(const string title, float value); 

    void update(const string title, float value);

    void show();

    void write_rsd_scores(ofstream &f_output_score,
                          const string title); 

    void write_score_title(ofstream &f_output_score); 

    Score operator/=(float win_size); 

};

#endif
