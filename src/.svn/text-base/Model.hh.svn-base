#ifndef MODEL_HEADER_HH 
#define MODEL_HEADER_HH 

#include <iostream>
#include <map>
#include <vector>

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

#include<SubModel.hh>

using namespace std;

typedef core::pose::Pose Pose;

class Model {

  private:
    string model_name_;
    string chain_id_;
    size_t num_of_rsd_;
    string sequence_;
    string sec_struct_;
    vector<float> xyz_;

    bool verbose_;

  private:
    template <class T>
    void get_ca_xyz_using_rosetta(Pose pose, T* predicate, vector<float>& xyz);

  public:
    
    Model(); 

    Model(string model_id, Pose poses); 

    Model(string model_id, Pose poses, bool verbose); 

    void set_poses(map<string, Pose> poses);

    void set_verbose(bool verbose);


    void generate_sub_models(Pose pose, vector<SubModel> &submodels, const size_t frag_size, 
                             const float bin_size, const float min_dist, const size_t vector_size); 

    void generate_overlapping_sub_models(Pose pose, vector<SubModel> &submodels, 
                                         const size_t frag_size, const float bin_size, 
                                         const float min_dist, const size_t vector_size); 
    
    void filter_similar_intra_submodels(Pose pose, vector<SubModel> &submodels, 
                                      const size_t frag_size, const float threshold); 
};

#endif
