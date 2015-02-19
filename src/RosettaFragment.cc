
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<functional>
#include<numeric>

#include <core/scoring/rms_util.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include<RosettaFragment.hh>
#include<rmsd.hh>

using namespace std;
using namespace core::scoring;
using namespace core;

RosettaFragment::RosettaFragment() {
} 

RosettaFragment::RosettaFragment(const int rsd_pos, 
                                 const string line, 
                                 const int window_size, 
                                 const bool is_in_top25) { 
  if(!line.empty()) { 
    rsd_pos_ = rsd_pos;
    string aminoacid, secstruct, pdb_id, chain_id;
    float phi, psi, omega;
    stringstream streams(line);
    streams>>pdb_id>>chain_id>>start_pos_>>aminoacid>>secstruct>>phi>>psi>>omega;
    end_pos_ = start_pos_ + (window_size - 1);
    transform(pdb_id.begin(), pdb_id.end(), pdb_id.begin(),::toupper);
    if(chain_id=="-" || chain_id=="_") 
      chain_id_ = "A";
    else
      chain_id_ = chain_id;
    pdb_id_ = pdb_id;
    subseq_ = aminoacid; subsec_ = secstruct;
    phis_.push_back(phi); psis_.push_back(psi); omegas_.push_back(omega);
    window_size_ = window_size;
    is_in_top25_ = is_in_top25;
  }
} 

RosettaFragment::RosettaFragment(const int rsd_pos, 
                                 const string line, 
                                 const int window_size, 
                                 const bool is_in_top25,
                                 const string tag,
                                 const int rsd_count) { 
  if(!line.empty()) { 
    rsd_pos_ = rsd_pos;
    string aminoacid, secstruct, pdb_id, chain_id;
    float phi, psi, omega;
    stringstream streams(line);
    streams>>pdb_id>>chain_id>>start_pos_>>aminoacid>>secstruct>>phi>>psi>>omega;
    end_pos_ = start_pos_ + (window_size - 1);
    transform(pdb_id.begin(), pdb_id.end(), pdb_id.begin(),::toupper);
    if(chain_id=="-" || chain_id=="_") 
      chain_id_ = "A";
    else
      chain_id_ = chain_id;
    pdb_id_ = pdb_id;
    subseq_ = aminoacid; subsec_ = secstruct;
    phis_.push_back(phi); psis_.push_back(psi); omegas_.push_back(omega);
    window_size_ = window_size;
    is_in_top25_ = is_in_top25;
    tag_ = tag;
    rsd_count_ = rsd_count;
  }
} 

void RosettaFragment::parse_line(const string line) {
  if(line.empty()) return; 
  int cur_rsd_pos;
  float phi, psi, omega;
  string aminoacid, secstruct, pdb_id, chain_id;
  stringstream streams(line);
  streams>>pdb_id>>chain_id>>cur_rsd_pos>>aminoacid>>secstruct>>phi>>psi>>omega;

  subseq_ += aminoacid; subsec_ += secstruct;
  phis_.push_back(phi); psis_.push_back(psi); omegas_.push_back(omega);
  if(phis_.size() == (size_t)window_size_ && 
     psis_.size() == (size_t)window_size_ && 
     omegas_.size() == (size_t)window_size_) {
    Pose pose;
    chemical::make_pose_from_sequence(pose, subseq_,
             *(chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID )));
    assert(pose.total_residue() == (size_t)window_size_);
    for(size_t pos = 1; pos <= pose.total_residue(); pos++) {
      if (!pose.residue(pos).is_protein() ) continue;
      pose.set_phi(pos, phis_[pos-1]); 
      pose.set_psi(pos, psis_[pos-1]); 
      pose.set_omega(pos, omegas_[pos-1]); 
    }
    pose_ = pose;
  }
}

void RosettaFragment::print2stdoutput() {
  cout<<"pdb_id: "<<pdb_id_<<" "
      <<chain_id_ <<" "<<rsd_pos_<<" "
      <<start_pos_<<" "<<end_pos_<<" "
      <<window_size_<<" "<<subseq_<<" "
      <<subsec_<<" "<<carmsd_<<"  "
      <<is_in_top25_<<endl;
  // int total_rsd = pose_.total_residue();
  // for(int i = 1; i <= total_rsd; i++) {
  //   cout<<"pdb id "<<i<<" "<<phis_[i-1]<<" "<<psis_[i-1]<<" "<<omegas_[i-1]<<endl;
  //   cout<<"pdb id "<<i<<" "<<pose_.phi(i)<<" "<<pose_.psi(i)<<" "<<pose_.omega(i)<<endl;
  // } 
}

void RosettaFragment::carmsd(Pose segmented_native_pose) {
  carmsd_ = CA_rmsd(segmented_native_pose, pose_, 1, window_size_); 
}

void RosettaFragment::collect_carmsd(map<int, list<float> > &rsd_carmsds) {
  list<float> carmsds = rsd_carmsds[rsd_pos_]; 
  carmsds.push_back(carmsd_);
  rsd_carmsds[rsd_pos_] = carmsds;
}

void RosettaFragment::collect_top25_carmsd(map<int, list<float> > &rsd_carmsds) {
  if(!is_in_top25_) return;
  list<float> carmsds = rsd_carmsds[rsd_pos_]; 
  carmsds.push_back(carmsd_);
  rsd_carmsds[rsd_pos_] = carmsds;
} 

Residues::Residues() {
}

Residues::Residues(int residue_pos, RosettaFragment rosetta_fragment) {
  list<RosettaFragment> fragments;
  fragments.clear();
  fragments.push_back(rosetta_fragment);
  residues_[residue_pos] = fragments;
}

void Residues::add(int residue_pos, RosettaFragment rosetta_fragment) {
  list<RosettaFragment> fragments = residues_[residue_pos];
  if(!fragments.empty()) 
    fragments.push_back(rosetta_fragment);
  else 
    fragments.push_back(rosetta_fragment);
  residues_[residue_pos] = fragments;
} 

map<string, Pose> Residues::get_poses(int rsd_pos) {
  map<string, Pose> poses;
  list<RosettaFragment> rosetta_fragments = residues_[rsd_pos];
  list<RosettaFragment>::iterator loit = rosetta_fragments.begin();
  while(loit != rosetta_fragments.end()) {
    Pose pose = loit->get_pose();
    poses[loit->get_rsd_count()] = pose;
  }
  return poses;
}

map<string, float> Residues::get_carmsd_native(int rsd_pos) {
  map<string, float> carmsds;
  list<RosettaFragment> rosetta_fragments = residues_[rsd_pos];
  list<RosettaFragment>::iterator loit = rosetta_fragments.begin();
  while(loit != rosetta_fragments.end()) {
    carmsds[loit->get_rsd_count()] = loit->get_carmsd();
  }
  return carmsds;
}

void Residues::do_clustering(float cluster_radius) { 
  // map<int, list<RosettaFragment> > residues;
  map<int, list<RosettaFragment> >::iterator milt = residues_.begin(); 
  int mode_  = 1, smode_ = 2;
  int number_of_fragment_  = 0, fragment_per_cluster_ = 0; 
  int total_residue = residues_.size();

  try {
    ofstream f_output_fragment, f_output;
    f_output_fragment.open("hybrid_fragment.fragK", ios::out);
    if(!f_output_fragment.is_open()) throw " fragment file cannot be opened.";

    while(milt != residues_.end()) {
      int rsd_pos = milt->first;
      map<string, Pose> fragment_poses = get_poses(rsd_pos); 
      FragmentCluster clustering;
      clustering.set_total_residue(total_residue);
      clustering.set_poses(fragment_poses);
      clustering.do_clustering_using_durandal(cluster_radius); 
      map<string, Pose> selected_templates = clustering.get_selected_templates(mode_, 
                                                                                smode_,
                                                                                number_of_fragment_,
                                                                                fragment_per_cluster_); 
      clustering.write_fragments(f_output_fragment, selected_templates, rsd_pos); 

      // write_fragments_selected_using_clustering(f_output_fragment, clustering, selected_templates);

      map<string, float> carmsd2native= get_carmsd_native(rsd_pos); 
      // write_cluster_info(f_output, clustering, carmsd2native);
      //delete after valgrid check 
      // selected_fragment_poses_.clear();
      milt++;
    }  
    f_output_fragment.close();
    f_output.close();
  }
  catch(const char* error) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<error<<endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<endl;
    exit(1);
  }
}

// void Residues::write_fragments_selected_using_clustering(ofstream &f_output_fragment,
//                                                          FragmentCluster clustering_, 
//                                                          map<string, Pose> selected_templates) {
//   try {
//     // if(selected_templates.empty()) 
//     //   throw  " templates are not found for residue. " + rsd_pos_;
//     // map<string, Pose> selected_extended_templates;
//     // if(selected_templates.empty()) 
//     //   throw  " templates are not found for residue. " + rsd_pos_;
//     // clustering_.write_fragments(f_output_fragment, selected_templates, rsd_pos_); 
//   }
//   catch(const char* e) {
//     std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<e<<std::endl;
//     exit(1);
//   }
//   catch(...) {
//     std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<std::endl;
//     exit(1);
//   }
// }

void Residues::write_cluster_info(ofstream &outfile, 
                                  FragmentCluster clustering_,
                                  map<string, float> carmsd2native,
                                  const size_t rsd_pos) {
  map<int, vector<string> > cluster_info = clustering_.get_clusters_info();
  map<int, vector<string> >::iterator mvit = cluster_info.begin();
  size_t total_num_model = 0;
  outfile<<"clustering for residue position: "<<rsd_pos<<std::endl;
  for(; mvit != cluster_info.end(); mvit++) {
    vector<string> members = mvit->second;
    vector<string>::iterator vit = members.begin();
    outfile<<"cluster "; outfile.width(5);  outfile<<mvit->first;
    outfile<<"  members "; outfile.width(5); outfile<<members.size()<<endl;
    total_num_model += members.size();
    for(; vit != members.end(); vit++) { 
      outfile<<"  "<<*vit<<"  ";
      outfile.setf(ios::fixed,ios::floatfield);
      outfile.width(7); outfile.precision(3); 
      outfile<<clustering_.get_distance2center(*vit);
      if(not carmsd2native.empty()) {
        float rmsd2native = carmsd2native["1"];
        outfile.width(7); outfile.precision(3); outfile<<rmsd2native;
      }
      outfile<<endl;
    }
  }
  outfile<<"cluster density rsd pos: ";
  outfile.width(5); outfile<<rsd_pos;
  outfile.width(5); outfile<<" cluster size: "<<cluster_info.size();
  outfile.width(5); outfile<<" # models: "<<total_num_model<<std::endl;
}


