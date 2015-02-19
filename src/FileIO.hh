#ifndef FILEIO_HEADER_HH 
#define FILEIO_HEADER_HH 

#include<fstream>
#include<sstream>
#include<string>
#include<vector>

#include<Residue.hh>

#include<SubModel.hh> // need some forward declartion.


#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

typedef core::pose::Pose Pose;

using namespace std;

class FileIO {

  private:
    typedef ProfileVector<int> IProfileVector;

    size_t get_vector_size(const string v_filepath); 

    void get_sliding_aprds(const string path2model, 
                           const size_t w_size, 
                           vector<ResidueProfile> &slided_residues); 
    
  public:

    FileIO();

    void parse_r_ckmean_output(const string r_output_filename, 
                               map<size_t, vector<size_t> > &m_cluster_sizes, 
                               map<size_t, vector<float> > &m_cluster_means); 

    void read_model_path(const string in_file, 
                         vector<string> &path2scores);

    void read_aprds_n_carmsd(const string path2model, 
                             vector<ResidueProfile> &residues); 

    void read_aprds_n_carmsd(const string path2model, 
                             const size_t w_size, 
                             vector<ResidueProfile> &residues); 

    void append_aprds_n_carmsd(const string path2model, vector<ResidueProfile> &residues); 

    void read_aprds(const string path2model, vector<ResidueProfile> &residues); 

    void read_aprdss(const string path2score, vector<ResidueProfile> &residues, bool is_carmsd);

    void read_sliding_aprdss( const string path2score, 
                              vector<ResidueProfile> &residues, 
                              size_t w_size);

    void append_aprds(const string path2model, vector<ResidueProfile> &residues); 

    void read_standard_pose(Pose &pose, const string pdb_path); 

    void read_poses_no_stout(const vector<string> path2models, map<string, Pose> &poses); 

    void read_poses_no_stout(const string path2model, map<string, Pose> &poses); 

    bool read_pose_no_stout(const string path2model, string &model_id,  Pose &pose); 

    void read_poses(const vector<string> path2models, map<string, Pose> &poses); 

    void read_poses(const string path2model, map<string, Pose> &poses); 

    void read_pose_no_stout(const string path2model, map<string, Pose> &poses); 

    void read_profile_vectors(const string path2vector, 
                              vector<ProfileVector<int> > &profile_vectors); 

    vector<ResidueProfile> read_sliding_aprds(const string path2model, const size_t w_size); 

    void read_vectors(const string path2vector, 
                      const string path2coord, 
                      vector<ProfileVector<int> > &profile_vectors); 

  private:

    void parse_file_path(const string full_path, string & filename); 

};

#endif
