// Copyright (C)  2012 Zhang Initative Research Unit
// RosettaWrapper.cc -description
// written by Rojan Shrestha 

#include<iostream>

#include<RosettaWrapper.hh>

using namespace std;

template <class T>
void RosettaWrapper::get_ca_xyz_using_rosetta(Pose pose, T* predicate, vector<float>& xyz){
  using namespace core::scoring;
  const core::Size nres = pose.total_residue();
  vector<core::Vector> p_coords;
  for(core::Size i=1; i<=nres; ++i) {
    for(core::Size j=1; j<=pose.residue(i).natoms();++j) {
      if((*predicate) (pose, pose, i, j)) {
        p_coords.push_back(pose.residue(i).xyz(j));
      }
    }
  }

  const size_t natoms = p_coords.size();
  xyz.clear(); xyz.reserve(natoms);
  for(size_t i=0; i<natoms; ++i) {
    for(int k=0; k<3; ++k) {
      xyz.push_back(float(p_coords[i][k]));
    }
  }

  // if(verbose_) {
  //   cout<<"Size of natoms "<<natoms<<endl;
  //   for(size_t i=0; i<natoms; ++i) {
  //     cout<<"Residue "<<i+1<<" CA coord ";
  //     for(int k=0; k<3; ++k) {
  //       cout<<p_coords[i][k]<<",";
  //     }
  //     cout<<endl;
  //   }
  // }
}

