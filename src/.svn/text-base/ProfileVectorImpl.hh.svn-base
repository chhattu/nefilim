#ifndef PROFILE_VECTOR_IMPL_HEADER_HH
#define PROFILE_VECTOR_IMPL_HEADER_HH

#include<iostream>
#include<vector>
#include<map>

#include<SubModel.hh>
#include<ProfileVectorIdentity.hh>
#include<EncodeDecode.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

using namespace std;

typedef core::pose::Pose Pose;

class ProfileVectorImpl {

  typedef ProfileVector<int> IProfileVector;
  typedef ProfileVector<long> LProfileVector;
  
  public:
    const char* pdb_database_path_;
    bool verbose_; 
    size_t frag_size_;
    size_t vector_size_;
    float min_dist_; 
    float max_dist_; 
    float bin_size_; 

    EncodeDecode encode_decode_;

  private:
    string path2model_list_;
    IProfileVector ref_prof_vect_; //ROS:: debugging
    LProfileVector encoded_ref_prof_vect_; //ROS:: debugging
    vector<ProfileVectorIdentity> profile_vectors_identities_; 

  public:
    ProfileVectorImpl();

    ProfileVectorImpl(const string path2model_list);

    ProfileVectorImpl(const string path2model_list, const size_t frag_size_, 
                      const float bin_size_,  const bool verbose);

    void compute_freq_vector_profiles();

    void filter_intra_submodels(float threshold); 

    void filter_intra_submodels(vector<SubModel> &submodels, const string pdb_path, float threshold); 

    void compute_freq_vector_profiles(bool debug_mode, bool overlapping); 

    // bool compute_freq_vector_profiles(map<string, Pose> poses);

    void update_prof_vect_identity(SubModel submodel); 

    void encode_prof_vect(SubModel submodel); 

    //ROS::This is just for testing case
    bool set_ref_prof_vect(const string path2pdb);

    void set_verbosity(const bool verbose);

    float ca_rmsd(vector<float> coord1, vector<float> coord2);

    void std_out_coords(ofstream &fstream_,  vector<SubModel> submodels); 

};

#endif 
