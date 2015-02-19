#include<FileIO.hh>
#include<stdlib.h>
#include<algorithm>
#include<cstdio>

FileIO::FileIO() {
}

void FileIO::read_model_path(const string in_file, vector<string> &path2scores) {
 string line;
 ifstream infile(in_file.c_str(), ifstream::in);
 if(infile.is_open()) {
   while(getline(infile, line)) {
     path2scores.push_back(line);
   }
   infile.close();
 }
 else {
   cout <<"Fatal: unable to open file." << endl;
   exit(1);
 }
}

void FileIO::parse_r_ckmean_output(const string r_output_filename, 
                                   map<size_t, vector<size_t> > &m_cluster_sizes, 
                                   map<size_t, vector<float> > &m_cluster_centers) {
  string line;
  bool first = true;
  vector<size_t> cluster_sizes; 
  vector<float> cluster_centers;
  size_t rsd_pos = 0;
  ifstream infile(r_output_filename.c_str(), ifstream::in);
  if(infile.is_open()) {
    while(getline(infile, line)) {
      size_t found = line.find("[1]"); 
      if(found != string::npos && found == 0) {
        stringstream ss(line);
        string element;
        while(ss>>element) {
          if(element != "[1]") {
            (first) ? cluster_sizes.push_back(atoi(element.c_str())):cluster_centers.push_back(atof(element.c_str()));
          }
        }
        if (first) { 
          ++rsd_pos;
          m_cluster_sizes[rsd_pos] = cluster_sizes; 
        }
        else 
          m_cluster_centers[rsd_pos] = cluster_centers; 

        cout<<"Residue position: "<<rsd_pos<<endl;

        (first) ? cluster_sizes.clear(): cluster_centers.clear(); 
        first = (first) ? false: true;
      }
    }
  }
  infile.close();
}

void FileIO::get_sliding_aprds(const string path2model, 
                               const size_t w_size, 
                               vector<ResidueProfile> &slided_residues) {
  vector<ResidueProfile> residue_profiles; 
  residue_profiles.clear();

  if(w_size < 2) return; // but try to redirect to one dimensional 
  string modelname = "";
  parse_file_path(path2model, modelname);
  read_aprds_n_carmsd(path2model, residue_profiles); 
  if(residue_profiles.empty()) return; //redirect to one d 

  size_t counter = 0;
  const string s_aprds = "AP_RDS"; 
  const size_t slided_tot_rsd = residue_profiles.size() - (w_size - 1);
  vector<ResidueProfile>::iterator it, in_it, sl_it;
  
  for(it = residue_profiles.begin(); it != residue_profiles.end(); it++) {
    size_t cur_rsd_pos = it->get_rsd_pos(); 
    float  sum_aprds = 0.0;
    if(cur_rsd_pos > slided_tot_rsd) break;
    for(in_it = it; in_it < it + w_size; in_it++) {
       float f_aprds = in_it->get_score(s_aprds, modelname);
       sum_aprds += f_aprds;
       // cout<<f_aprds<<" ";
    }
    
    if(cur_rsd_pos == 1 && !slided_residues.empty()) {  
      sl_it = slided_residues.begin(); 
    }

    if(slided_residues.size() < slided_tot_rsd) {
      ResidueProfile rsd_profile(cur_rsd_pos, modelname, s_aprds, sum_aprds/float(w_size));
      slided_residues.push_back(rsd_profile);
      // cout<<" initial insertion: "<<sum_aprds/float(w_size)<<endl;
    }
    else { 
      sl_it->update(cur_rsd_pos, modelname, s_aprds, sum_aprds/float(w_size));
      sl_it++;
      //cout<<" appended value: "<<sum_aprds/float(w_size)<<endl;
    }
    counter++;
  }
} 


void FileIO::read_aprds_n_carmsd(const string path2model, vector<ResidueProfile> &residues) {
  const string ap_rds="AP_RDS", ap_rsd="AP_STD", s_ca_rmsd="CA_RMSD"; 
  string line;
  ifstream infile(path2model.c_str(), ifstream::in);
  if(infile.is_open()) {
    string modelname = "";
    parse_file_path(path2model, modelname);
    if (modelname.empty()) return; 
    while(getline(infile, line)) {
      size_t rsd_pos = 0;
      float mean = 0, stdev = 0, ca_rmsd = -1; 
      stringstream ss(line);
      ss >> rsd_pos>> mean >>stdev>>ca_rmsd; 
      ResidueProfile rsd_profile(rsd_pos, modelname, ap_rds, mean);
      rsd_profile.append_score(modelname, ap_rsd, stdev); 
      rsd_profile.append_score(modelname, s_ca_rmsd, ca_rmsd); 
      residues.push_back(rsd_profile);
    }
    infile.close();
  }
}

void FileIO::append_aprds_n_carmsd(const string path2model, vector<ResidueProfile> &residues) {
  const string ap_rds="AP_RDS", ap_rsd="AP_STD", s_ca_rmsd="CA_RMSD"; 
  string line;
  if (residues.empty())  return; 
  ifstream infile(path2model.c_str(), ifstream::in);
  if(infile.is_open()) {
    string modelname = "";
    parse_file_path(path2model, modelname);
    if (modelname.empty()) return; 
    vector<ResidueProfile>::iterator it = residues.begin();
    while(getline(infile, line) && it != residues.end()) {
      size_t rsd_pos = 0;
      float mean = 0, stdev = 0, ca_rmsd = -1; 
      stringstream ss(line);
      ss >> rsd_pos>> mean >>stdev>>ca_rmsd; 
      it->update(rsd_pos, modelname, ap_rds, mean);
      it->update(rsd_pos, modelname, ap_rsd, stdev);
      it->update(rsd_pos, modelname, s_ca_rmsd, ca_rmsd);
      it++;
    }
    infile.close();
  }
}

void FileIO::read_aprds(const string path2model, vector<ResidueProfile> &residues) {
  const string ap_rds = "AP_RDS", ap_std = "AP_STD"; 
  string line;
  ifstream infile(path2model.c_str(), ifstream::in);
  if(infile.is_open()) {
    string modelname = "";
    parse_file_path(path2model, modelname);
    if (modelname.empty()) return; 
    while(getline(infile, line)) {
      size_t rsd_pos = 0;
      float mean = 0, stdev = 0; //, ca_rmsd2native;
      stringstream ss(line);
      ss >> rsd_pos>> mean >>stdev; 
      ResidueProfile rsd_profile(rsd_pos, modelname, ap_rds, mean);
      rsd_profile.append_score(modelname, ap_std, stdev); 
      residues.push_back(rsd_profile);
    }
    infile.close();
  }
}

void FileIO::append_aprds(const string path2model, 
                        vector<ResidueProfile> &residues) {
  const string ap_rds = "AP_RDS", ap_rsd = "AP_STD"; 
  string line;
  if (residues.empty())  return; 
  ifstream infile(path2model.c_str(), ifstream::in);

  if(infile.is_open()) {
    string modelname = "";
    parse_file_path(path2model, modelname);
    if (modelname.empty()) return; 
    vector<ResidueProfile>::iterator it = residues.begin();

    while(getline(infile, line) && it != residues.end()) {
      size_t rsd_pos = 0;
      float mean = 0, stdev = 0; //, ca_rmsd2native;
      stringstream ss(line);
      ss >> rsd_pos>> mean >>stdev; 
      it->update(rsd_pos, modelname, ap_rds, mean);
      it->update(rsd_pos, modelname, ap_rsd, stdev);
      it++;
    }
    infile.close();
  }
}

void FileIO::read_sliding_aprdss( const string path2score, 
                                  vector<ResidueProfile> &slided_residues,
                                  size_t w_size){
  vector<string> path2models;
  path2models.clear();
  read_model_path(path2score, path2models);
  if(path2models.empty()) return; 

  for(vector<string>::iterator it=path2models.begin(); it != path2models.end(); it++) {
    // cout<<"Model "<<(*it)<<endl;
    get_sliding_aprds((*it), w_size, slided_residues); 
    // cout<<"F...."<<endl;
  }
}


void FileIO::read_aprdss(const string path2score, vector<ResidueProfile> &residues, 
                         bool with_not_carmsd) {
  cout<<endl;
  cout<<"Reading scores..."<<endl;
  cout<<"-----------------"<<endl;
  cout<<"reading score file(s) from "<<path2score<<endl;
  vector<string> path2models;
  path2models.clear();
  read_model_path(path2score, path2models);
  if(path2models.empty()) return; 
  residues.clear();
  vector<string>::iterator it = path2models.begin();
  for(;it != path2models.end(); it++) {
    if(residues.empty()) {
      if(with_not_carmsd) 
        read_aprds((*it), residues); 
      else 
        read_aprds_n_carmsd((*it), residues); 
    }
    else { 
      if(with_not_carmsd) 
        append_aprds((*it), residues); 
      else 
        append_aprds_n_carmsd((*it), residues); 
    }
  }
  cout<<"# of residues:     "<<residues.size()<<endl;
  cout<<"# of data file(s): "<<path2models.size()<<endl<<endl;
}

void FileIO::read_standard_pose(Pose &pose, const string pdb_path) {
  freopen("/dev/null", "w", stdout); //console.log for just debugging 
  core::io::pdb::pose_from_pdb(pose, pdb_path);
  fclose(stdout);
  freopen("/dev/tty", "a", stdout); //console.log for just debugging 
}

void FileIO::read_poses_no_stout(const string path2model, map<string, Pose> &poses) {
  freopen("/dev/null", "w", stdout); //console.log for just debugging 
  read_poses(path2model, poses); 
  // fclose(stdout);
  freopen("/dev/tty", "a", stdout); //console.log for just debugging 
  cout<<"# of denovo models:               "<<poses.size()<<endl;
}

void FileIO::read_poses_no_stout(const vector<string> path2models, map<string, Pose> &poses) {
  freopen("/dev/null", "w", stdout); //console.log for just debugging 
  read_poses(path2models, poses); 
  // fclose(stdout);
  freopen("/dev/tty", "a", stdout); //console.log for just debugging 
  cout<<"# of denovo models: "<<poses.size()<<endl;
}

void FileIO::read_poses(const vector<string> path2models, map<string, Pose> &poses) {
  if(path2models.empty()) return; 
  poses.clear();
  vector<string>::const_iterator it;
  for(it=path2models.begin(); it != path2models.end(); it++) {
    Pose pose;
    string modelname = "";
    parse_file_path((*it), modelname);
    cout<<"Model:: "<<*it<<endl;
    // core::io::pdb::centroid_pose_from_pdb(pose, *it);
    core::io::pdb::pose_from_pdb(pose, *it);
    poses[modelname] = pose;
  }
}

void FileIO::read_pose_no_stout(const string path2model, map<string, Pose> &poses) {
  Pose pose;
  string modelname;
  parse_file_path(path2model, modelname);
  core::io::pdb::centroid_pose_from_pdb(pose, path2model);
  poses[modelname] = pose;
}


void FileIO::read_poses(const string path2model, map<string, Pose> &poses) {
  vector<string> path2models;
  path2models.clear();
  read_model_path(path2model, path2models);
  if(path2models.empty()) return; 
  poses.clear();
  for(vector<string>::iterator it=path2models.begin(); it != path2models.end(); it++) {
    Pose pose;
    string modelname = "";
    parse_file_path((*it), modelname);
    core::io::pdb::centroid_pose_from_pdb(pose, *it);
    // core::io::pdb::pose_from_pdb(pose, *it);
    poses[modelname] = pose;
  }
}

bool FileIO::read_pose_no_stout(const string path2model, string &model_id,  Pose &pose) {
  if(path2model.empty()) return false; 
  parse_file_path(path2model, model_id); 
  std::transform(model_id.begin(), model_id.end(),model_id.begin(), ::toupper);
  if(model_id.empty()) return false;
  core::io::pdb::centroid_pose_from_pdb(pose, path2model);
  return true;
}

void FileIO::parse_file_path(const string full_path, string &filename) {
  if(full_path.empty()) { 
    filename.empty();  
    return;
  }

  size_t file_start_pos = full_path.find_last_of("/\\");
  size_t ext_pos        = full_path.find_last_of(".");
  size_t filename_size  = ext_pos - file_start_pos - 1;

  filename = "";
  if (file_start_pos < full_path.length()) { 
    if(ext_pos > file_start_pos && file_start_pos <= full_path.length()) { 
      filename = full_path.substr(file_start_pos + 1, filename_size);
    } 
  }
}

void FileIO::read_profile_vectors(const string path2vector, 
                                  vector<ProfileVector<int> > &profile_vectors) {

  if(path2vector.empty()) return;
  string line;
  ifstream in_vector_file(path2vector.c_str(), ifstream::in);
  string sBit; int iBit;
  size_t counter = 0;
  if(in_vector_file.is_open()) {
    IProfileVector profile_vector;
    vector<int> freq_profiles;
    while(getline(in_vector_file, line)) {
      stringstream ss(line);
      while(!ss.eof()) {
        counter++;
        if(counter++ < 2)
          ss>>sBit;
        else { 
          ss>>iBit;
          freq_profiles.push_back(iBit);
        }
      }
      profile_vector.set_profile_vector(freq_profiles);
      profile_vectors.push_back(profile_vector);
      freq_profiles.clear();
    }
  }
  in_vector_file.close();
}

size_t FileIO::get_vector_size(const string v_filepath) {
  string line;
  ifstream in_file(v_filepath.c_str(), ifstream::in);
  size_t counter = 0;
  if(in_file.is_open()) {
    while(getline(in_file, line)) {
      istringstream vector_ss(line); 
      string tmp;
      while(vector_ss>>tmp) counter++;
      break;
    } 
  }
  in_file.close();
  return counter;
}

void FileIO::read_vectors(const string path2vector, const string path2coord, 
                          vector<ProfileVector<int> > &profile_vectors) {

  if(path2vector.empty() || path2coord.empty()) return;
  string vector_line, coord_line;

  size_t v_size = get_vector_size(path2vector);
  size_t c_size = get_vector_size(path2coord); 
  // cout<<"v_size "<<v_size<<" c_size "<<c_size<<endl;

  ifstream in_coord_file(path2coord.c_str(), ifstream::in);
  ifstream in_vector_file(path2vector.c_str(), ifstream::in);

  if(in_vector_file.is_open() && in_coord_file.is_open()) {
    vector<int> freq_profiles;
    vector<float> ca_coords;
    string svElement, scElement;

    while(getline(in_vector_file, vector_line) && getline(in_coord_file, coord_line)) {
      IProfileVector profile_vector;
      size_t counter = 0;
      istringstream vector_ss(vector_line), coord_ss(coord_line);
      // cout<<"Vector line "<<vector_line<<endl; cout<<"Coordinate  "<<coord_line<<endl;
      while(coord_ss || vector_ss) { 
        coord_ss>>scElement; vector_ss>>svElement;  
        counter++;
        if(counter > 2){
          if(counter <= v_size) freq_profiles.push_back(atoi(svElement.c_str())); 
          if(counter <= c_size)  ca_coords.push_back(atof(scElement.c_str()));
        } 
      }
      profile_vector.set_profile_vector(freq_profiles);
      profile_vector.set_ca_coord(ca_coords);
      profile_vectors.push_back(profile_vector);
      ca_coords.clear(); freq_profiles.clear();
    }
    in_vector_file.close();
    in_coord_file.close();
  }
}
