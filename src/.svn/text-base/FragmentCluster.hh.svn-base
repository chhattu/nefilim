// Copyright (C)  2013 Zhang Initative Research Unit
// Fragment.hh -description
// written by Rojan Shrestha 

#ifndef  FRAGMENT_CLUSTER_HEADER_HH
#define  FRAGMENT_CLUSTER_HEADER_HH

#include<iostream>
#include<map>
#include<list>
#include<set>
#include<vector>
#include<string>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
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

#include<DistMatrix.hh>

typedef core::pose::Pose RPose;

using namespace std;

class FragmentCluster {
  private:
    size_t total_residue_;
    size_t total_sub_residue_;

    map<string, RPose> poses_;
    map<string, RPose> selected_poses_;
    map<int, vector<string> > clusters_info_;
    vector<string> selected_templates_; //keeps indexes of map which is selected.
    map<string, float>  distances_;
    vector<string> fragment_ids_;

  public:

    FragmentCluster();

    FragmentCluster(map<string, RPose> poses);

    void set_poses(map<string, RPose> poses);

    bool is_more_than(size_t index, 
                      size_t value); 

    vector<string> get_selected_templates_id(); //keeps indexes of map which is selected.

    map<string, RPose> get_selected_templates(); 

    map<string, RPose> get_selected_templates(const size_t mode, 
                                              const size_t dmode,
                                              size_t number_of_fragment,
                                              size_t fragment_per_cluster); 

    map<string, RPose> get_cluster_center_templates(size_t number_of_template); 

    map<int, vector<string> > get_clusters_info();

    float get_distance2center(const string pose_id); 

    void set_total_residue(size_t total_residue); 

    vector<float> get_coordinates(const string pose_id);

    void do_clustering_using_durandal(const float cluster_radius); 

    void parse_clusters(DistMatrix &dm, 
                        vector<int> biggest_cluster, 
                        vector<int> &remaining_pdbs, 
                        int cluster_id); 

    void get_normalized_ca_xyz_from_pose(RPose npose, 
                                         vector<float> &fragments);

    void get_normalized_ca_xyz_from_poses(map<string, 
                                          vector<float> > &fragments);

    void divide_equally(vector<int> &selected_nums, 
                        size_t number_of_fragment, 
                        size_t fragment_per_cluster); 

    void divide_proportions(vector<int> &cluster_members,
                            size_t number_of_fragment, 
                            size_t fragment_per_cluster); 

    void make_balance(vector<int> &selected_nums, 
                      vector<size_t> &reduced_cluster, 
                      size_t &total_num_fragments,
                      size_t frags_to_be_added);

    void remove_elements_randomly(vector<int> &selected_nums, 
                                        size_t to_be_removed); 

    void select_templates_sequentially(vector<int> cluster_members); 

    void select_templates_randomly(vector<int> cluster_members); 

    void add_elements_randomly(vector<int> &selected_nums, 
                               size_t updated_total); 

    void write_fragment(const string pdbname, 
                        RPose subpose, 
                        ofstream &f_output, 
                        const size_t rsd_pos, 
                        const size_t start_rsd_pos,
                        const size_t end_rsd_pos); 

    void write_fragments(ofstream &outfile, 
                         map<string, RPose > selected_templates,
                         const size_t rsd_pos); 



  private:
    float ca_rmsd(vector<float> coor1, vector<float> coor2); 

    string getfilename(const string full_path); 
    
};

#endif

