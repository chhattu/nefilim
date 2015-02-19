// Copyright (C)  2012 Zhang Initative Research Unit
// RosettaWrapper.hh -description
// written by Rojan Shrestha 

#ifndef  ROSETTAWRAPPER_HEADER_HH
#define  ROSETTAWRAPPER_HEADER_HH

#include<vector>
#include<iostream>

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


using namespace std;


namespace RosettaWrapper {

  typedef core::pose::Pose Pose;
  
  template <class T>
  void get_ca_xyz_using_rosetta(Pose pose, T* predicate, vector<float>& xyz);

}

#endif

