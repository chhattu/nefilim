// Copyright (C)  2013 Zhang Initative Research Unit
// Fragment.cc -description
// written by Rojan Shrestha 

#include<fstream>
#include<exception>
#include<cassert>

#include<rmsd.hh>
#include<qcprot.hh>
#include<FragmentCluster.hh>

FragmentCluster::FragmentCluster() {
  poses_.clear();
}

FragmentCluster::FragmentCluster(map<string, RPose> poses) {
  set_poses(poses);
}

void FragmentCluster::set_total_residue(size_t total_residue) {
  total_residue_ = total_residue;
}

void FragmentCluster::set_poses(map<string, RPose> poses) {
  poses_ = poses;
  map<string, RPose>::iterator msit = poses_.begin();
  total_sub_residue_ = msit->second.total_residue(); 
  for(; msit != poses_.end(); msit++) 
    fragment_ids_.push_back(msit->first);
}

bool FragmentCluster::is_more_than(size_t index, 
                                   size_t value) {
  return clusters_info_[index].size() > value; 
}

vector<string> FragmentCluster::get_selected_templates_id() {
  return selected_templates_; 
}

map<string, RPose> FragmentCluster::get_selected_templates() {
  return selected_poses_;
}

map<string, RPose> FragmentCluster::get_selected_templates(const size_t mode,
                                                           const size_t smode,
                                                           size_t number_of_fragment,
                                                           size_t fragment_per_cluster) {
  try {
    vector<int> selected_nums;
    if(smode == 1)  //based on the proportion of cluster
      divide_proportions(selected_nums, number_of_fragment, fragment_per_cluster);
    else if(smode == 2) //each cluster gives you five.
      divide_equally(selected_nums, number_of_fragment, fragment_per_cluster); 

    if(selected_nums.empty()) throw "empty selection.";
    if(mode == 1) // random selection
      select_templates_randomly(selected_nums);
    else if(mode == 2) //random selection 
      select_templates_sequentially(selected_nums);

    if(selected_templates_.size() < number_of_fragment) 
      throw " # of template is less than required.";

    vector<string>::iterator vit = selected_templates_.begin();
    for(; vit != selected_templates_.end(); vit++) {
      selected_poses_[*vit] = poses_[*vit];
    }
    return selected_poses_; 
  }
  catch(const char* e) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<e<<endl;
    exit(1);
  }
}

map<string, RPose> FragmentCluster::get_cluster_center_templates(size_t number_of_template) {
  const size_t clust_size = clusters_info_.size();
  size_t counter = 0;
  map<string, RPose> clust_center_templates;
  map<int, vector<string> >::iterator miit = clusters_info_.begin();
  for(; miit != clusters_info_.end(); miit++) {
    if(++counter > number_of_template) break;
    if(miit->second.size() < 2) break;
    vector<string> centers       = miit->second;
    vector<string>::iterator vit = centers.begin();
    clust_center_templates[*vit] = poses_[*vit];  
  }

  if(clust_size < number_of_template) {
    std::cerr<<"  warning: "<<__FILE__<<" "<<__LINE__ 
             <<" # of clusters: "<<clust_size<<endl;
  }

  return clust_center_templates;
  
  //make more smart
  // cout<<"  warning: remaining templates are choosen from first cluster."<<endl;
  // vector<string> centers = miit->second;
  // vector<string>::iterator vit = centers.begin();
  // while(clust_center_templates.size() < number_of_template) {
  //   for(; vit != centers.end(); vit++) 
  //     clust_center_templates[*vit] = poses_[*vit];  
  // } 
  // return clust_center_templates;
}

map<int, vector<string> > FragmentCluster::get_clusters_info() {
  return clusters_info_;
}

float FragmentCluster::get_distance2center(const string pose_id) {
  return distances_[pose_id];
}

vector<float> FragmentCluster::get_coordinates(const string pose_id) {
  RPose pose = poses_[pose_id];
  vector<float> fragments; 
  get_normalized_ca_xyz_from_pose(pose, fragments); 
  return fragments;
}

void FragmentCluster::get_normalized_ca_xyz_from_pose(RPose pose, vector<float> &fragments){
  const core::Size nres = pose.total_residue();
  vector<core::Vector> p_coords;
  for(core::Size i=1; i<=nres; ++i) {
    for(core::Size j=1; j<=pose.residue(i).natoms();++j) {
      if((core::scoring::is_protein_CA) (pose, pose, i, j)) 
        p_coords.push_back(pose.residue(i).xyz(j));
    }
  }
  const size_t natoms = p_coords.size();
  float cx = 0, cy = 0, cz = 0;
  for(size_t i=0; i<natoms; ++i) { 
    cx += float(p_coords[i][0]);
    cy += float(p_coords[i][1]);
    cz += float(p_coords[i][2]);
  }
  cx /= natoms; cy /= natoms; cz /= natoms;
  vector<float> xyzs; xyzs.clear(); xyzs.reserve(natoms);
  for(size_t i=0; i<natoms; ++i) { 
    xyzs.push_back(float(p_coords[i][0]) - cx);
    xyzs.push_back(float(p_coords[i][1]) - cy);
    xyzs.push_back(float(p_coords[i][2]) - cz);
  }
  fragments = xyzs;
}

void FragmentCluster::get_normalized_ca_xyz_from_poses(map<string, vector<float> > &fragments){
  try {
    if(poses_.empty()) {
      std::cerr<<"fatal: "<<__FILE__<<" line: "<<__LINE__<<" empty models."<<std::endl;
      exit(1);
    } 
    map<string, RPose>::iterator mit = poses_.begin();
    for(; mit != poses_.end(); mit++) {
      vector<float> coords;
      get_normalized_ca_xyz_from_pose(mit->second, coords);
      fragments[mit->first] = coords;
    }
  }
  catch(std::exception &error) {
    std::cerr<<"fatal: "<<__FILE__<<" line: "<<__LINE__<<" "<<error.what()<<std::endl;
    exit(1);
  }
}

void FragmentCluster::do_clustering_using_durandal(const float cluster_radius) { 
  try {
    map<string, vector<float> > fragments;
    get_normalized_ca_xyz_from_poses(fragments);
    if(fragments.empty()) 
      throw "no fragments are selected."; 
    DistMatrix dm(fragments, cluster_radius);
    dm.brute_init();
    vector<int> remaining_pdbs = dm.get_all_PDBs_index();
    for (int num_of_cluster = 0; remaining_pdbs.size() > 0; ++num_of_cluster) {
      vector<int> biggest_cluster;
      vector< vector<int> > pole_position_clusters;
      dm.get_biggest_cluster(remaining_pdbs, biggest_cluster, pole_position_clusters);
      parse_clusters(dm, biggest_cluster, remaining_pdbs, num_of_cluster); 
    }
  }
  catch(const char* e) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<e<<endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<"unknown error."<<endl;
    exit(1);
  }
}

void FragmentCluster::parse_clusters(DistMatrix &dm, 
                                     vector<int> biggest_cluster, 
                                     vector<int> &remaining_pdbs, 
                                     int cluster_id) {
  try {
    int center_idx = biggest_cluster[0];
    vector<string> cluster_members;
    for(int i = 0; i < static_cast<int>(biggest_cluster.size()); i++) {
      int current_idx   = biggest_cluster[i];
      string pose_id    = fragment_ids_[current_idx]; // class variable 
      if(pose_id.empty()) throw "empty pose name.";
      float dist2center = 0;
      DistRange &range  = dm.get(center_idx, current_idx);
      if(range._maxi == range._mini) 
        dist2center = range._maxi;
      else throw "error in distance matrix.";
      cluster_members.push_back(pose_id);
      distances_[pose_id] = dist2center;
      remove(remaining_pdbs, current_idx);
    }
    clusters_info_[cluster_id] = cluster_members;
  }
  catch(const char* e) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<e<<endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<"unknown error."<<endl;
    exit(1);
  }
}

void FragmentCluster::divide_equally(vector<int> &selected_nums,
                                     size_t number_of_fragment,
                                     size_t fragment_per_cluster) {
  try {
    map<int, vector<string> >::iterator miit = clusters_info_.begin();
    size_t counter = 0,total_num_fragments = 0;
    vector<size_t> reduced_cluster;

    for(;miit != clusters_info_.end(); miit++, counter++ ) {
      if(total_num_fragments >= number_of_fragment) break; 
      const size_t num_of_members = miit->second.size();
      if(num_of_members >= fragment_per_cluster) {
        selected_nums.push_back(fragment_per_cluster);
        total_num_fragments += fragment_per_cluster;
        if((num_of_members - fragment_per_cluster) > 0)
          reduced_cluster.push_back(num_of_members - fragment_per_cluster);
      }
      else if(num_of_members > 2 && num_of_members < fragment_per_cluster) {
        selected_nums.push_back(num_of_members);
        total_num_fragments += num_of_members;
      }
      else if(num_of_members == 2 && num_of_members < fragment_per_cluster) {
        cout<<"  warning: template is selected from singleton cluster."<<endl;
        selected_nums.push_back(num_of_members);
        total_num_fragments += num_of_members;
      }
    }
    if(counter > fragment_per_cluster) 
      cout<<"  warning: templates are selected from " <<counter<< " clusters."<<endl;

    if(total_num_fragments < number_of_fragment) {
      cout<<"  warning: fragments are choosen repeatedly."<<endl;
      while(total_num_fragments < number_of_fragment) {
        size_t frags_to_be_added = number_of_fragment - total_num_fragments;
        make_balance(selected_nums, reduced_cluster, total_num_fragments, frags_to_be_added);
        if(reduced_cluster.size()==0) 
          throw " # fragments are not enough and increase # of fragments for cluster";
      }
      size_t frags_to_be_removed = number_of_fragment - total_num_fragments;
      if(frags_to_be_removed > 0)
        remove_elements_randomly(selected_nums, frags_to_be_removed); 
    }
  }
  catch(const char* error) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<error<<std::endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<std::endl;
    exit(1);
  }
}

void FragmentCluster::make_balance(vector<int> &selected_nums, 
                                   vector<size_t> &reduced_cluster, 
                                   size_t &total_num_fragments,
                                   size_t frags_to_be_added) { 
  size_t counter = 0; 
  for(; counter < reduced_cluster.size(); counter++) {
    selected_nums[counter]++; 
    total_num_fragments++;
    reduced_cluster[counter]--;
    frags_to_be_added--;
    if(frags_to_be_added <= 0) return; 
  }
  vector<size_t> temp_cluster;
  for(counter = 0; counter < reduced_cluster.size(); counter++) {
    if(reduced_cluster[counter] > 0)
      temp_cluster.push_back(reduced_cluster[counter]); 
  }
  reduced_cluster = temp_cluster;
}

void FragmentCluster::remove_elements_randomly(vector<int> &selected_nums, 
                                                     size_t to_be_removed) {
  while(to_be_removed-- >  0) {
    int index = rand()% (selected_nums.size() - 1);
    std::cout<<"  warning: removed element from index: "<<index<<std::endl;
    selected_nums[index--]--; 
  }
}


void FragmentCluster::divide_proportions(vector<int> &selected_nums, 
                                         size_t number_of_fragment,
                                         size_t fragment_per_cluster) {
  try {
    selected_nums.clear();
    size_t counter = 0; float total_num_fragments = 0.0;
    vector<int> upper_limits; 
    vector<float> proportions; 
    float max_fragment = static_cast<float>(number_of_fragment);

    map<int, vector<string> >::iterator miit = clusters_info_.begin();
    while(miit != clusters_info_.end() ) {
      if(total_num_fragments > max_fragment && counter >= fragment_per_cluster) {
        break;
      }
      float value = static_cast<float>(miit->second.size());
      total_num_fragments += value;
      proportions.push_back(value);
      upper_limits.push_back(miit->second.size());
      miit++; ++counter;
    }

    if(counter > fragment_per_cluster) {
      cout<<"  warning: templates are selected from "<<counter<< " clusters."<<endl;
    }

    if(total_num_fragments < max_fragment)
      throw "# template fragments are less than required.";

    int updated_total = number_of_fragment;
    for(size_t i = 0; i < proportions.size(); i++) {
      float nume = static_cast<float>(number_of_fragment) / total_num_fragments;
      int current_total = static_cast<int>(round((nume * proportions[i])));
      selected_nums.push_back(current_total);
      updated_total -= current_total; 
    }
    if(updated_total < 0) {
      while(updated_total++ < 0) {
        int index = rand()% (selected_nums.size() - 1);
        std::cout<<"  warning: removed element from index: "<<index<<std::endl;
        selected_nums[index--]--; 
      }
    }
    else if(updated_total > 0) 
      add_elements_randomly(selected_nums, updated_total); 
  }
  catch(const char* error) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<error<<std::endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<std::endl;
    exit(1);
  }
}

void FragmentCluster::add_elements_randomly(vector<int> &selected_nums, 
                                            size_t updated_total) {
  vector<size_t> upper_limits;
  for(size_t i = 0; i< selected_nums.size(); i++) {
    vector<string> members = clusters_info_[i]; 
    size_t remain = (members.size() - selected_nums[i]);
    if(remain > 0) upper_limits.push_back(remain);
  }
  while(updated_total > 0) {
    int index = rand() % (selected_nums.size() - 1);
    if(upper_limits[index] > 0) { 
      selected_nums[index]++; 
      std::cout<<"  warning: added element at index: "<<index<<"; "
               <<"current value: "<<selected_nums[index]<<"; "
               <<"# max limit: "<<upper_limits[index]<<std::endl;
      upper_limits[index]--;
      updated_total--;
    }
  }
}


void FragmentCluster::select_templates_sequentially(vector<int> selected_nums) {
                                                    // size_t number_of_fragment) {
  try {
    if(selected_nums.empty()) throw "template fragment selection failed";
    vector<int>::iterator vit = selected_nums.begin();
    map<int, vector<string> >::iterator miit = clusters_info_.begin();
    while(vit != selected_nums.end() && miit != clusters_info_.end()) {
      map<float, string> sorted_distances; 
      vector<string> members = miit->second;
      vector<string>::iterator vsit = members.begin(); 
      while(vsit != members.end()) {
        sorted_distances[distances_[*vsit]] = *vsit;
        vsit++;
      }
      map<float, string>::iterator mfit = sorted_distances.begin();
      int allowed_num_members = static_cast<int>(*vit);
      while(allowed_num_members-- > 0) {
        selected_templates_.push_back(mfit->second);
        mfit++;
      }
      miit++; vit++;
    }
  }
  catch(const char* error) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<error<<endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<std::endl;
    exit(1);
  }
}

//here is the starting point
void FragmentCluster::select_templates_randomly(vector<int> selected_nums) {
  try {
    vector<int>::iterator vit = selected_nums.begin();
    map<int, vector<string> >::iterator miit = clusters_info_.begin();
    selected_templates_.clear();
    while(vit != selected_nums.end() && miit != clusters_info_.end()) {
      vector<string> members = miit->second;
      vector<string>::iterator vsit = members.begin(); 
      int min = 0, max = members.size();
      srand(time(NULL));
      set<int> unique_checkers;
      pair<set<int>::iterator, bool> ssit; 
      int counter = *vit;
      if(counter > max) throw " indexes are not correctly selected.";
      while(counter > 0) {
        int index = rand()%(max-min);
        pair<set<int>::iterator, bool> ssit = unique_checkers.insert(index);
        if(ssit.second) {
          selected_templates_.push_back(members[index]);
          counter--;
        }
      }
      miit++; vit++;
    }
  }
  catch(const char* error) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" "<<error<<endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<std::endl;
    exit(1);
  }
}

void FragmentCluster::write_fragment(const string pdbname, 
                                     RPose subpose, 
                                     ofstream &f_output, 
                                     const size_t rsd_pos, 
                                     const size_t start_rsd_pos,
                                     const size_t end_rsd_pos) {
  try {
    if(subpose.empty()) throw std::exception();
    // if(pdbname.length() == 11) throw std::exception(); //"invalid pdbname [????_?]");
    const string pdb_id   = pdbname.substr(0,4);
    const string chain_id = pdbname.substr(5,1);
    const string sequence = subpose.sequence();
    for(size_t i = start_rsd_pos; i <= end_rsd_pos; i++) {
      f_output.width(5); f_output<<pdb_id;
      f_output.width(2); f_output<<chain_id;
      if(rsd_pos==1) {
        f_output.width(6); f_output<<rsd_pos + (i - 1);
      }
      else if((rsd_pos + total_sub_residue_ - 1) == total_residue_) {
        f_output.width(6); f_output<<rsd_pos + (i - 2);
      }
      else {
        f_output.width(6); f_output<<rsd_pos + (i - 2);
      }
      f_output.width(2); f_output<<sequence.substr(i-1, 1);
      f_output.width(2); f_output<<subpose.secstruct(i);
      f_output.setf(ios::fixed,ios::floatfield);
      f_output.width(9); f_output.precision(3); f_output<<subpose.phi(i);
      f_output.width(9); f_output.precision(3); f_output<<subpose.psi(i);
      f_output.width(9); f_output.precision(3); f_output<<subpose.omega(i);
      f_output<<endl;
    }
    f_output<<endl;
  }
  catch(const char* error) {
    std::cerr<<"fatal: "<<"FILE: "<<__FILE__
             <<" LINE: "<<__LINE__<<" "<<error<<endl;
    exit(1);
  }
  catch(...) {
    std::cerr<<"fatal: "<<"FILE: "<<__FILE__
             <<" LINE: "<<__LINE__<<" unknown error."<<endl;
    exit(1);
  }
}

void FragmentCluster::write_fragments(ofstream &outfile, 
                                      map<string, RPose > selected_templates,
                                      const size_t rsd_pos) {
  map<string, RPose>::iterator msit = selected_templates.begin(); 
  size_t total_sub_residue = msit->second.total_residue(); 
  size_t start_rsd_pos = 1, end_rsd_pos = total_sub_residue;
  size_t diff = total_sub_residue - total_sub_residue_;
  if(diff > 0) { 
    start_rsd_pos = 2;  end_rsd_pos = total_sub_residue - 1; 
    if((rsd_pos-1)==0) {
      // cout<<"First residue "<<total_sub_residue<<endl;
      start_rsd_pos = 1; end_rsd_pos = total_sub_residue - 2; 
    }
    else if(((rsd_pos + total_sub_residue_ - 1) - total_residue_)==0) {
      // cout<<"last residue "<<total_sub_residue<<endl;
      start_rsd_pos = 3; end_rsd_pos = total_sub_residue; 
    }
    // cout<<"  fragment debug rsd pos "<<rsd_pos<<" diff "<<diff
    //     <<" start "<<start_rsd_pos<<" end "<<end_rsd_pos<<endl;
    // cout<<"  subsequence "<<msit->second.sequence()<<endl;
  }
  outfile<<" position: ";  outfile.width(12);outfile<<rsd_pos;
  outfile<<" neighbors: "; outfile.width(12);outfile<<selected_templates.size()<<endl<<endl;
  for(; msit != selected_templates.end(); msit++) 
    write_fragment(msit->first, msit->second, outfile, rsd_pos, start_rsd_pos, end_rsd_pos);
}

