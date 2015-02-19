#ifndef ROSETTA_FRAGMENT_HEADER_HH 
#define ROSETTA_FRAGMENT_HEADER_HH 

#include<iostream>
#include<list>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>

#include<RosettaWrapperFwd.hh>
#include<FragmentCluster.hh>


typedef core::pose::Pose Pose;

using namespace std;


class RosettaFragment {

  private:
    int rsd_pos_;
    string pdb_id_, chain_id_;
    int start_pos_, end_pos_, window_size_;
    vector<float> psis_;
    vector<float> phis_;
    vector<float> omegas_;
    string subseq_, subsec_;
    float carmsd_; 
    Pose pose_;
    bool is_in_top25_;
    string tag_;
    int rsd_count_;


  public:
    RosettaFragment();
    
    RosettaFragment(const int rsd_pos, 
                    const string line, 
                    const int window_size, 
                    const bool is_in_top25); 

    RosettaFragment(const int rsd_pos, 
                    const string line, 
                    const int window_size, 
                    const bool is_in_top25,
                    const string tag,
                    const int rsd_counter); 

    void parse_line(const string line);

    void parse_position_line(const string line, int &rsd_pos, int &num_of_fragments); 

    void print2stdoutput();

    void carmsd(Pose segmented_native_pose); 

    void collect_carmsd(map<int, list<float> > &rsd_carmsds); 

    void collect_top25_carmsd(map<int, list<float> > &rsd_carmsds); 

    int get_rsd_pos() {
      return rsd_pos_;
    }

    int get_start_pos() {
      return start_pos_;
    }

    Pose get_pose() {
      return pose_;
    }

    void set_carmsd(const float carmsd) {
      carmsd_ = carmsd;
    }

    float get_carmsd() {
      return carmsd_;  
    }

    string get_pdb_id_with_chain() {
      return pdb_id_ + "_" + chain_id_;
    }

    void dump_pdb() {
      ostringstream ostr;
      int rsd_pos = rsd_pos_;
      ostr<<rsd_pos;
      pose_.dump_pdb(pdb_id_ + "_" + ostr.str() + ".pdb");
    }

    const string get_tag() {
      return tag_;
    }

    const string get_rsd_count() {
      stringstream str;
      str << rsd_count_;
      return str.str(); 
    }
};  

class Residues {

  private:
    // Pose native_pose_;
    map<int, list<RosettaFragment> > residues_;

  public:
    Residues();

    Residues(int residue_pos, RosettaFragment rosetta_fragment); 

    void add(int residue_pos, RosettaFragment rosetta_fragment); 

    map<string, Pose> get_poses(int rsd_pos); 

    map<string, float> get_carmsd_native(int rsd_pos); 

    void do_clustering(float cluster_radius);

    void write_fragments_selected_using_clustering(ofstream &f_output_fragment,
                                                   FragmentCluster clustering, 
                                                   map<string, Pose> selected_templates); 

    void write_cluster_info(ofstream &outfile, 
                            FragmentCluster clustering,
                            map<string, float> carmsd2native, 
                            const size_t rsd_pos); 

};

#endif
