// Copyright (C)  2012 Zhang Initative Research Unit
// Fragment.cc -description
// written by Rojan Shrestha 

#include<Fragment.hh>

#include<iostream>
#include<assert.h>

using namespace std;

Fragment::Fragment(){
}

Fragment::Fragment(const size_t rsd_pos, const size_t frag_size) {
  rsd_pos_   = rsd_pos;
  frag_size_ = frag_size;
}

Fragment::Fragment( const size_t rsd_pos, const size_t frag_size, 
                    const size_t max_num_frags, const size_t max_frag_per_model) {
  unknown1   = 0.506; unknown2 = 4.707; unknown4 = 22.614; unknown_id = 5; //strange to me 
  rsd_pos_   = rsd_pos;
  frag_size_ = frag_size;
  max_num_frags_ = max_num_frags;
  max_frag_per_model_ = max_frag_per_model;
  frag_per_model_ = 0;
}


void Fragment::insert(Pose pose, const string pdb_id, const size_t rsd_pos, 
                      const float ca_rmsd, ofstream &f_output) {
  if(!pose.empty()) {
    if(pdb_ids_.size() > max_num_frags_) return;
    if(frag_per_model_ > max_frag_per_model_) return;
    const string sequence = pose.sequence();
    cout<<"ROS::SEQUENCE "<<sequence<<endl;
    for(size_t i = rsd_pos; i <= rsd_pos + frag_size_; i++) {
      if(frag_per_model_ > max_frag_per_model_) return;
      f_output<<pdb_id<<"  "<<i<<" "<<sequence.substr(i, 1)<<" "<<pose.secstruct(i)<<" "<<pose.phi(i)<<" "
              <<pose.psi(i)<<" "<<pose.omega(i)<<" "<<ca_rmsd<<endl;
      rsd_poss_.push_back(i);
      pdb_ids_.push_back(pdb_id);
      cout<<"ROS::CHAIN"<<pose.chain(2)<<endl;
      // chain_ids_.push_back(pose.chain_id);
      amino_acids_.push_back(sequence.substr(i, 1));
      sec_strs_.push_back(pose.secstruct(i));
      phis_.push_back(pose.phi(i));
      psis_.push_back(pose.psi(i));
      omegas_.push_back(pose.omega(i));
      ca_rmsds_.push_back(ca_rmsd);
    }
    f_output<<endl;
    frag_per_model_++;
  } 
}

void Fragment::get_native_fragment(ofstream &f_output, 
                                   const Pose native_pose, 
                                   const string pdb_id, 
                                   const size_t rsd_pos, 
                                   const float ca_rmsd) {
  if(native_pose.empty()) return;
  if(pdb_ids_.size() > max_num_frags_) return;
  if(frag_per_model_ > max_frag_per_model_) return;
  const string sequence = native_pose.sequence();
  const size_t total_residue   = native_pose.total_residue();
  if((rsd_pos + frag_size_) > total_residue) return;
  for(size_t i = rsd_pos; i <rsd_pos + frag_size_; i++) {
    string four_pdb_id, chain_id;
    parse_filename(pdb_id, four_pdb_id, chain_id);
    f_output.width(5); f_output<<four_pdb_id;
    f_output.width(2); f_output<<chain_id;
    f_output.width(6); f_output<<i;
    f_output.width(2); f_output<<sequence.substr(i-1, 1);
    f_output.width(2); f_output<<native_pose.secstruct(i);
    f_output.setf(ios::fixed,ios::floatfield);
    f_output.width(9); f_output.precision(3); f_output<<native_pose.phi(i);
    f_output.width(9); f_output.precision(3); f_output<<native_pose.psi(i);
    f_output.width(9); f_output.precision(3); f_output<<native_pose.omega(i);
    f_output.width(9); f_output.precision(3); f_output<<ca_rmsd;
    f_output<<endl;
  }
  f_output<<endl;
  frag_per_model_++;
}

void Fragment::write_rosetta_fragment_format(ofstream &f_output, Pose pose, const string pdb_id, 
                                              const size_t rsd_pos, const float ca_rmsd) {
  if(pose.empty()) return;
  if(pdb_ids_.size() > max_num_frags_) return;
  if(frag_per_model_ > max_frag_per_model_) return;
  const string sequence = pose.sequence();
  for(size_t i = rsd_pos; i <rsd_pos + frag_size_; i++) {
    string four_pdb_id, chain_id;
    parse_filename(pdb_id, four_pdb_id, chain_id);
    f_output.width(5); f_output<<four_pdb_id;
    f_output.width(2); f_output<<chain_id;
    f_output.width(6); f_output<<i;
    f_output.width(2); f_output<<sequence.substr(i-1, 1);
    f_output.width(2); f_output<<pose.secstruct(i);
    f_output.setf(ios::fixed,ios::floatfield);
    f_output.width(9); f_output.precision(3); f_output<<pose.phi(i);
    f_output.width(9); f_output.precision(3); f_output<<pose.psi(i);
    f_output.width(9); f_output.precision(3); f_output<<pose.omega(i);
    f_output.width(9); f_output.precision(3); f_output<<ca_rmsd;
    f_output<<endl;
  }
  f_output<<endl;
  frag_per_model_++;
}

size_t Fragment::get_total_fragments() {
  return pdb_ids_.size();
}

size_t Fragment::get_fragment_per_model() {
  return frag_per_model_;
}

void Fragment::reset_fragment_per_model() {
  frag_per_model_ = 0;
} 

void Fragment::show() {
  cout<<"Position: "<<rsd_pos_<<endl;
  for(size_t i =0; i < rsd_poss_.size(); i++) {
    cout<<pdb_ids_[i]<<" "<<chain_ids_[i]<<" "<<sec_strs_[i]<<" ";
    cout<<phis_[i]<<" "<<psis_[i]<<" "<<omegas_[i]<<endl;
    if((i+1) % frag_size_ == 0)
      cout<<endl;
  }
}

