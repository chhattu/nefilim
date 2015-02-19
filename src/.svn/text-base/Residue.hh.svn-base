#ifndef RESIDUE_VAL_HEADER_HH 
#define RESIDUE_VAL_HEADER_HH 

#include<iostream>
#include<fstream>
#include<vector>
#include<list>
#include<map>
#include<algorithm>


#include <protocols/moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/loops/Loops.hh>
#include <core/io/pose_stream/MetaPoseInputStream.hh>
#include <core/io/pose_stream/util.hh>

#include <core/scoring/rms_util.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include<LocalProfile.hh>
#include<FragmentCluster.hh>

using namespace std;


class ResidueProfile {

  private:
    typedef std::pair<string, float> t_pair;
    typedef core::pose::Pose Pose;

    float max_mean, min_mean;
    float max_stdev, min_stdev;

    map<string, float> min_values;
    map<string, float> max_values;

    size_t rsd_pos_;
    size_t total_residue_;
    static const int max_num_rsds = 25;

    bool slide_after_normalization_; 

    map<string, Score> local_scores_;  //model_id, score combination
    map<string, Pose> selected_fragment_poses_; // it has two more residues than given window size 

    size_t mode_;
    size_t smode_;
    size_t number_of_fragment_; 
    size_t fragment_per_cluster_;    
    size_t number_of_templates_;

  private:
    bool sort_by_value(const t_pair m1, const t_pair m2);

    void sort_by_value(map<string,float> scores, vector<t_pair> &t_vector); 

    struct comp{
      bool operator()(const t_pair &pair1, const t_pair &pair2) {
        return pair1.second < pair2.second;
      }
    };

    void add_scores(const string score_title, 
                    map<string, float> &val1, 
                    ResidueProfile rsd);

    void set_cluster_values(const string score_title, 
                            map<string, float> &val1, 
                            ResidueProfile rsd);

    map<string, Pose> get_extended_selected_templates(map<string, Pose> selected_templates); 

  public:
    ResidueProfile();   

    ResidueProfile(size_t rsd_pos);

    ResidueProfile(size_t rsd_pos, 
                   string model_id, 
                   string score_title, 
                   float value);

    ResidueProfile(size_t rsd_pos, 
                   string score_title, 
                   vector<string> modelnames, 
                   vector<float> values);

    size_t get_rsd_pos();

    size_t get_size(); 

    void set_score(map<string, Score> local_scores);

    void set_slide_after_normalization(bool value);  

    void set_total_residue(size_t total_residue); 

    void set_selection_mode(size_t mode, 
                            size_t smode); 

    void set_selection_number(size_t number_of_fragment, 
                              size_t fragment_per_cluster,
                              size_t number_of_templates); 

    bool is_model_exist(const string model_id, 
                        const string score_title); 

    vector<float> get_score(string score_title); 

    void get_score(string score_title, 
                   list<float> &scores); 

    float get_score(const string score_title, 
                    const string modelname);

    bool get_score(const string score_title, 
                   map<string, float> &models_n_scores); 

    string inline get_modelname(const string score_title, 
                                float aprds); 

    vector<string> get_modelnames(const string score_title, 
                                  vector<float> aprdss); 

    map<string, float> get_modelnscore(string title); 

    map<string, float> get_model_n_sorted_score(string title, size_t topN); 

    void append_score(string model_id, string score_title, float value); 

    void update(size_t rsd_pos, string model_id, string score_title, float value); 
    
    void update_score(const string score_title, map<string, float> models_n_scores); 

    ResidueProfile sort_by_score(const string score_title, 
                                 size_t topn);

    ResidueProfile sort_by_score(const string score_title, 
                                 float cutoff);

    ResidueProfile sort_by_score(const string score_title, 
                                 float cutoff, 
                                 bool &is_not_increased);

    void show();

    void write_rsd_scores(ofstream &f_output_score,
                          const string title);

    void write_score_title(ofstream &f_output_score); 


    void dump_average_rmsd_to_native(const string score_title, ofstream &f_output_score); 

    void dump_template_fragment(ofstream &out); 

    void dump_template_fragment_from_pdb(const Pose pose, ofstream &out, const int window_size); 

    vector<string> get_model_ids();

    vector<string> get_score_titles();

    bool normalize_score(string score_title); 

    bool normalize_score(const float min, const float max, string score_title); 

    ResidueProfile operator+=(ResidueProfile rsd); 

    ResidueProfile operator/=(const float win_size); 

    void do_clustering(const float cluster_radius, 
                       map<string, Pose> fragment_poses, 
                       ofstream &f_output_fragment, 
                       ofstream &f_output_template, 
                       ofstream &f_output);

    void write_clustered_templates(ofstream &f_output_fragment,
                                   FragmentCluster clustering_, 
                                   map<string, Pose> selected_templates); 

    void get_centroid_fragment_poses(const size_t window_size, 
                                     const map<string, Pose> poses, 
                                     map<string, Pose> & fragment_pose); 

    void get_standard_fragment_poses(const size_t window_size, 
                                     const map<string, Pose> poses, 
                                     map<string, Pose> & fragment_pose);

    Pose get_centroid_fragment_pose(Pose old_pose, size_t start, size_t end); 

    Pose get_standard_fragment_pose(Pose old_pose, size_t start, size_t end); 

    void write_fragment(const string pdbname, 
                        Pose subpose, 
                        ofstream &f_output); 

    void write_fragments(ofstream &outfile, 
                         map<string, Pose > selected_templates); 

    void write_template_fragments_pdb( ofstream &out, 
                                      FragmentCluster clustering); 

    void write_cluster_info(ofstream &outfile, 
                            FragmentCluster clustering); 

    void write_selected_templates_native_RMSD(ofstream &outfile,
                                              FragmentCluster clustering); 

    inline void compare_value(const string score_title, float &lowest, float &largest) {
      float min_mean = min_values[score_title];
      float max_mean = max_values[score_title];
      lowest  = (lowest < min_mean) ? lowest:min_mean;
      largest = (largest > max_mean) ? largest:max_mean;
    }

    inline void compare_value(const string score_title, 
                              map<string, float> &lowests, map<string, float> &largests) {
      float min_mean = min_values[score_title];
      float max_mean = max_values[score_title];
      float lowest  = lowests[score_title];
      float largest = largests[score_title];
      lowest  = (lowest < min_mean) ? lowest:min_mean;
      largest = (largest > max_mean) ? largest:max_mean;
      lowests[score_title]  = lowest;
      largests[score_title] = largest;
    }

    inline void compare_mean(float &lowest, float &largest) {
      lowest  = (lowest < min_mean) ? lowest:min_mean;
      largest = (largest > max_mean) ? largest:max_mean;
    }
    
    inline void compare_stdev(float &lowest, float &largest) {
      largest = (largest > max_stdev) ? largest:max_stdev;
      lowest  = (lowest < min_stdev) ? lowest:min_stdev;
    }
};

#endif
