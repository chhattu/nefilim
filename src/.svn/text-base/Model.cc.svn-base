#include<iostream>

#include<rmsd.hh>
#include<Model.hh>
// #include<RosettaWrapper.hh> //template using is bad in c++ to call source file instead of header file
// #include<RosettaWrapperFwd.hh>

using namespace std;

Model::Model() {
}

Model::Model(string model_name, Pose pose) {
  model_name_ = model_name;
  num_of_rsd_ = pose.total_residue(); 
  verbose_ = false;
}

Model::Model(string model_name, Pose pose, bool verbose) {
  model_name_ = model_name;
  num_of_rsd_ = pose.total_residue(); 
  //chain_id_ 
  // int temp = pose.chain(1); 
  sequence_ = pose.sequence();
  sec_struct_ = pose.secstruct();
  verbose_ = verbose;
} 

void Model::set_verbose(bool verbose){
  verbose_ = verbose;
}

void Model::filter_similar_intra_submodels(Pose pose, vector<SubModel> &submodels, 
                                           const size_t frag_size, const float threshold) {
  vector<float> xyz;
  if(num_of_rsd_ > 0) xyz.reserve(num_of_rsd_);
  get_ca_xyz_using_rosetta(pose, core::scoring::is_protein_CA, xyz);
  assert(3*num_of_rsd_ == xyz.size());
  submodels.clear();

  const size_t ndim = 3; // 3-dim
  vector<float> frag_xyz; frag_xyz.reserve(ndim*frag_size);
  size_t start_rsd_pos = 0, rsd_pos = 0, gap_length = 0;
  float ca_rmsd_with_front_rsd = 0.0;
  bool isFirstEntry = true;
  SubModel ref_submodel; // default copy constructor since I don't use pointer

  for(size_t i = 0; i < num_of_rsd_; i++) {
    rsd_pos++;
    if(rsd_pos > frag_size) 
      frag_xyz.erase(frag_xyz.begin(), frag_xyz.begin() + 3);
    for(size_t m = ndim*i; m < ndim*i + ndim; m++)
      frag_xyz.push_back(xyz[m]);
    if(frag_xyz.empty()) return;

    if( rsd_pos >= frag_size) { 
      start_rsd_pos++;
      if(isFirstEntry) {
        SubModel temp_ref_submodel(model_name_, chain_id_, start_rsd_pos);
        isFirstEntry = false;
        ref_submodel = temp_ref_submodel; 
        ref_submodel.set_gap_length(1);
        ref_submodel.set_ca_coord(frag_xyz);
        continue;
      }
      gap_length++;
      SubModel submodel(model_name_, chain_id_, start_rsd_pos);
      ca_rmsd_with_front_rsd = ref_submodel.get_carmsd(frag_xyz);
      if(ca_rmsd_with_front_rsd > threshold) {
        submodels.push_back(ref_submodel);
        ref_submodel = submodel;
        ref_submodel.set_ca_coord(frag_xyz);
        ref_submodel.set_gap_length(gap_length);
        gap_length = 0;
      }
    }
  }
  submodels.push_back(ref_submodel);
}


void Model::generate_overlapping_sub_models(Pose pose, vector<SubModel> &submodels, 
                                            const size_t frag_size, const float bin_size, 
                                            const float min_dist, const size_t vector_size) {
  // using namespace RosettaWrapper;
  vector<float> xyz;
  if(num_of_rsd_ > 0) xyz.reserve(num_of_rsd_);
  get_ca_xyz_using_rosetta(pose, core::scoring::is_protein_CA, xyz);
  assert(3*num_of_rsd_ == xyz.size());
  submodels.clear();
  const size_t ndim = 3; // 3-dim
  vector<float> frag_xyz; 
  frag_xyz.reserve(ndim*frag_size);
  size_t start_rsd_pos = 1, rsd_pos = 0;

  for(size_t i = 0; i < num_of_rsd_; i++) {
    rsd_pos++;
    if(rsd_pos > frag_size) 
      frag_xyz.erase(frag_xyz.begin(), frag_xyz.begin() + 3);

    for(size_t m = ndim*i; m < ndim*i + ndim; m++)
      frag_xyz.push_back(xyz[m]);
    if(frag_xyz.empty()) return;

    if( rsd_pos >= frag_size) { // && rsd_pos < num_of_rsd_ -frag_size) {
      string frag_seq = sequence_.substr(start_rsd_pos - 1, frag_size);
      string frag_ss  = sec_struct_.substr(start_rsd_pos - 1, frag_size);
      SubModel submodel(model_name_, chain_id_, frag_seq, frag_ss, start_rsd_pos);
      cout<<"Sequence "<<frag_seq<<endl;
      cout<<"SecStruc "<<frag_ss<<endl;
      cout<<"Frag size : "<<frag_xyz.size()<<endl;
      cout<<"Residue position: "<<start_rsd_pos <<", "<<rsd_pos<<endl;
      for(size_t i = 0; i < frag_xyz.size(); i++) {
        cout<<frag_xyz[i]<<" ";
      }
      cout<<endl;
      submodel.compute_frequency_profile_vector(frag_xyz, frag_size, bin_size, min_dist, vector_size);
      submodel.set_ca_coord(frag_xyz);
      submodels.push_back(submodel);
      start_rsd_pos++;
    }
  }
}

void Model::generate_sub_models(Pose pose, vector<SubModel> &submodels, const size_t frag_size, 
                                const float bin_size, const float min_dist, const size_t vector_size) {
  // using namespace RosettaWrapper;
  vector<float> xyz;
  if(num_of_rsd_ > 0) xyz.reserve(num_of_rsd_);
  get_ca_xyz_using_rosetta(pose, core::scoring::is_protein_CA, xyz);
  assert(3*num_of_rsd_ == xyz.size());
  submodels.clear();
  const size_t ndim = 3; // 3-dim
  vector<float> frag_xyz; 
  frag_xyz.reserve(ndim*frag_size);
  int start_rsd_pos = 0;
  for(size_t rsd_pos = 0; rsd_pos < num_of_rsd_; rsd_pos++) {
    for(size_t m = ndim*rsd_pos; m < ndim*rsd_pos + ndim; m++) { 
      frag_xyz.push_back(xyz[m]);
    }
    if(frag_xyz.empty()) return;
    if((rsd_pos + 1) % frag_size == 0) {
      //cout<<"Frag size "<<frag_size<<" = "<<rsd_pos+1<<" mod "<<(rsd_pos+1)%frag_size<<endl;
      string frag_seq = sequence_.substr(rsd_pos, frag_size);
      string frag_ss  = sec_struct_.substr(rsd_pos, frag_size);
      SubModel submodel(model_name_, chain_id_, frag_seq, frag_ss, start_rsd_pos + 1);
      cout<<"Residue Position: "<<start_rsd_pos + 1<<", "<<rsd_pos + 1<<endl;
      submodel.compute_frequency_profile_vector(frag_xyz, frag_size, bin_size, min_dist, vector_size);
      submodel.set_ca_coord(frag_xyz);
      submodels.push_back(submodel);
      frag_xyz.clear();
      frag_xyz.reserve(ndim*frag_size);
      start_rsd_pos = rsd_pos;
    }
  }
}

template <class T>
void Model::get_ca_xyz_using_rosetta(Pose pose, T* predicate, vector<float>& xyz){
  using namespace core::scoring;
  const core::Size nres = pose.total_residue();
  vector<core::Vector> p_coords;
  for(core::Size i=1; i<=nres; ++i) {
    for(core::Size j=1; j<=pose.residue(i).natoms();++j) {
      if((*predicate) (pose, pose, i, j)) 
        p_coords.push_back(pose.residue(i).xyz(j));
    }
  }
  const size_t natoms = p_coords.size();
  xyz.clear();
  xyz.reserve(natoms);
  for(size_t i=0; i<natoms; ++i) {
    for(int k=0; k<3; ++k) 
      xyz.push_back(float(p_coords[i][k]));
  }
  if(verbose_) {
    cout<<"Size of natoms "<<natoms<<endl;
    for(size_t i=0; i<natoms; ++i) {
      cout<<"Residue "<<i+1<<" CA coord ";
      for(int k=0; k<3; ++k) 
        cout<<p_coords[i][k]<<",";
      cout<<endl;
    }
  }
}

