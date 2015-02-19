#include<SubModel.hh>

#include<iostream>
#include<cstring>
#include<cstdlib>
#include<cmath>
#include<SubModelIdentity.hh>

using namespace std;

SubModel::SubModel() {
  verbose_= false;  
}

SubModel::SubModel(const string line) {
  string element;
  bool isFirst = true;
  vector<float> coords;
  stringstream ssline(line);
  while(ssline>>element) {
    if(isFirst) {
      submodel_id_ = element;
      isFirst = false;

      char *converted_char = new char[element.size() + 1];
      copy(element.begin(), element.end(), converted_char);
      converted_char[element.size()] = '\0';
      char *delimited_char = strtok(converted_char, "-");
      size_t parser_count = 0;
      while( delimited_char != NULL) {
        if(parser_count == 1) model_name_    = delimited_char; 
        if(parser_count == 2) chain_id_      = delimited_char; 
        if(parser_count == 3) start_rsd_pos_ = atoi(delimited_char);
        if(parser_count == 4) gap_           = atoi(delimited_char); 
        delimited_char = strtok(NULL, "-");
      }
      delete[] converted_char;
    }
    else {
      coords.push_back(atof(element.c_str()));
    }
  }
  if(!coords.empty()) 
    prof_vect_.set_ca_coord(coords);
}

SubModel::SubModel(string model_name, string chain_id, size_t start_rsd_pos) {
  model_name_ = model_name;
  chain_id_   = chain_id;
  start_rsd_pos_ = start_rsd_pos;
  s_sec_struct_  = "";
  s_subsequence_ = "";
}

SubModel::SubModel(string model_name, string chain_id, size_t start_rsd_pos, size_t gap) {
  model_name_ = model_name;
  chain_id_   = chain_id;
  start_rsd_pos_ = start_rsd_pos;
  gap_ = gap;
  s_sec_struct_  = "";
  s_subsequence_ = "";
}

SubModel::SubModel(string model_name, string chain_id, string frag_seq, 
                   string frag_ss, size_t start_rsd_pos) {
  model_name_ = model_name;
  chain_id_   = chain_id;
  start_rsd_pos_ = start_rsd_pos;
  s_sec_struct_ = frag_ss;
  s_subsequence_ = frag_seq;
  // sec_struct_    = encode_string(frag_ss, SubModel::EncodeWhat.sec_struct);
  // subsequence_   = encode_string(frag_seq, SubModel::EncodeWhat.sequence);
  // cout<<"ROS::"<<frag_seq<<" "<<subsequence_<<endl;
  // cout<<"ROS::"<<frag_ss<<" "<<sec_struct_<<endl;
  verbose_ = false;  
}

unsigned int SubModel::encode_string(const string input_string) {
  if(input_string.empty()==0) return 0;
  unsigned int encoded_value = 0;

  map<string, short> lookup_map;

  for(size_t i = 0; i < input_string.length(); i++) {
    string key = input_string.substr(i, 1).c_str(); 
    unsigned int value = lookup_map[key];
    encoded_value = encoded_value | value; 
    encoded_value = encoded_value<<3;
  }
  return encoded_value; 
}

float SubModel::get_carmsd(vector<float> xyzs) {
  ProfileVector<int> prof_vect;
  prof_vect.set_ca_coord(xyzs);
  vector<float> xyz_coords = prof_vect.get_ca_coord(); 
  return prof_vect_.ca_rmsd(xyz_coords);
}

float SubModel::compute_carmsd(SubModel submodel) {
  vector<float> xyz_coords = submodel.get_prof_vect().get_ca_coord(); 
  return prof_vect_.ca_rmsd(xyz_coords);
}

float SubModel::compute_cosine(SubModel submodel, const size_t frag_size, 
                               const float bin_size, const float min_dist, 
                               const size_t vector_size) {
  vector<float> ca_xyz;
  if(prof_vect_.size() <= 0) {
    vector<float> ca_xyz = prof_vect_.get_ca_coord();
    compute_frequency_profile_vector(ca_xyz, frag_size, bin_size, min_dist, vector_size);
  }
  ca_xyz = submodel.get_prof_vect().get_ca_coord();
  submodel.compute_frequency_profile_vector(ca_xyz, frag_size, bin_size, min_dist, vector_size);
  ProfileVector<int> profile_vector = submodel.get_prof_vect(); 
  return prof_vect_.compute_cosine(profile_vector.get_profile_vector());
}

void SubModel::set_distance_vector(vector<float> xyz, const size_t fragment_size) {
  vector<float> distance_vector; 
  for(size_t i=0; i < fragment_size; i++) {
    float x1=xyz[i*3], y1=xyz[i*3+1], z1=xyz[i*3+2];
    for(size_t j=i+2; j < fragment_size; j++) {
      float x2=xyz[j*3], y2=xyz[j*3+1], z2=xyz[j*3+2];
      float distance = sqrt(pow((x2-x1), 2) + pow((y2-y1), 2) + pow((z2-z1), 2)); 
      distance_vector.push_back(distance);
    }
  }
  fprof_vect_.set_profile_vector(distance_vector);
}

float SubModel::compute_distance_vector_cosine(SubModel submodel, const size_t frag_size) {
  vector<float> ca_xyz;
  if(fprof_vect_.size() <= 0) {
    fprof_vect_.set_ca_coord(prof_vect_.get_ca_coord());
    // ca_xyz = fprof_vect_.get_ca_coord();
    ca_xyz = fprof_vect_.get_coords();
    set_distance_vector(ca_xyz, frag_size); 

    // cout<<"Reference vector: "<<get_submodel_id()<<": ";
    // fprof_vect_.verbose();
    // cout<<"Reference Matrix: "<<endl;
    // get_matrix(ca_xyz, frag_size); 
    ca_xyz.clear();
  }

  // ca_xyz = submodel.get_prof_vect().get_ca_coord();
  ca_xyz = submodel.get_prof_vect().get_coords();
  // cout<<"Size "<<ca_xyz.size()<<endl;
  submodel.set_distance_vector(ca_xyz, frag_size);

  // cout<<"Model vector: "<<submodel.get_submodel_id()<<": ";
  // submodel.get_fprof_vect().verbose();

  // cout<<"Model Matrix: "<<endl;
  // submodel.get_matrix(ca_xyz, frag_size); 

  ProfileVector<float> fprof_vect = submodel.get_fprof_vect(); 
  // return fprof_vect_.compute_cosine(fprof_vect.get_profile_vector());
  return fprof_vect_.pearsonCC(fprof_vect.get_profile_vector());
}

float SubModel::compute_QScore(SubModel submodel, const size_t frag_size) {
  vector<float> ca_xyz;
  if(fprof_vect_.size() <= 0) {
    fprof_vect_.set_ca_coord(prof_vect_.get_ca_coord());
    // ca_xyz = fprof_vect_.get_ca_coord();
    ca_xyz = fprof_vect_.get_coords();
    set_distance_vector(ca_xyz, frag_size); 
    // cout<<"Reference vector: "<<get_submodel_id()<<": ";
    // fprof_vect_.verbose();
    // cout<<"Reference Matrix: "<<endl;
    // get_matrix(ca_xyz, frag_size); 
    ca_xyz.clear();
  }

  // ca_xyz = submodel.get_prof_vect().get_ca_coord();
  ca_xyz = submodel.get_prof_vect().get_coords();
  submodel.set_distance_vector(ca_xyz, frag_size);

  // cout<<"Model vector: "<<submodel.get_submodel_id()<<": ";
  // submodel.get_fprof_vect().verbose();

  // cout<<"Model Matrix: "<<endl;
  // submodel.get_matrix(ca_xyz, frag_size); 

  ProfileVector<float> fprof_vect = submodel.get_fprof_vect(); 
  // return fprof_vect_.compute_cosine(fprof_vect.get_profile_vector());
  // return fprof_vect_.pearsonCC(fprof_vect.get_profile_vector());
  return fprof_vect_.QScore(fprof_vect.get_profile_vector());
}


void SubModel::get_matrix(vector<float> xyz, const size_t fragment_size) {
  float matrix[fragment_size][fragment_size];
  for(size_t i=0; i < fragment_size; i++) 
    for(size_t j=0; j < fragment_size; j++) 
      matrix[i][j] = 0.00; 

  for(size_t i=0; i < fragment_size; i++) {
    float x1=xyz[i*3], y1=xyz[i*3+1], z1=xyz[i*3+2];
    for(size_t j=i+1; j < fragment_size; j++) {
      float x2=xyz[j*3], y2=xyz[j*3+1], z2=xyz[j*3+2];
      float dist = sqrt(pow((x2-x1), 2) + pow((y2-y1), 2) + pow((z2-z1), 2)); 
      matrix[i][j] = dist; 
    }
  }

  // cout<<endl<<"Matrix:"<<endl;
  // just for debugging
  for(size_t i=0; i < fragment_size; i++) {
    for(size_t j=0; j < fragment_size; j++) {
      cout.setf(ios::fixed, ios::floatfield);
      cout.width(7);
      cout.precision(2);
      cout<<matrix[i][j]<<"  ";
    }
    cout<<endl;
  }
  cout<<endl;
}

void SubModel::compute_frequency_profile_vector( vector<float> xyz, 
                      const size_t frag_size, const float bin_size, 
                      const float min_dist, const size_t vector_size) {
  vector<int> prof_freq_vector; 
  for(size_t i = 0; i < vector_size; i++)
    prof_freq_vector.push_back(0);

  float  matrix[frag_size][frag_size]; // just for debugging
  for(size_t i=0; i < frag_size; i++) {
    float x1=xyz[i*3], y1=xyz[i*3+1], z1=xyz[i*3+2];
    for(size_t j=i+1; j < frag_size; j++) {
      float x2=xyz[j*3], y2=xyz[j*3+1], z2=xyz[j*3+2];
      float dist = sqrt(pow((x2-x1), 2) + pow((y2-y1), 2) + pow((z2-z1), 2)); 
      int index  = int(ceil((dist - min_dist) / bin_size));
      // cout<<"ROS: rsd "<<i + 1<<" "<<j + 1<<" index "<<index<<" dist "<<dist<<" ";
      // cout<<"ROS: coord "<<x1<<" "<<y1<<" "<<z1<<" "<<x2<<" "<<y2<<" "<<z2<<endl;
      if(vector_size < size_t(index)) {
        cout<<"warning: index is out of max set "<<index<<" dist "<<dist<<endl;
        index = (int)vector_size;
      }
      ++prof_freq_vector[index];
      matrix[i][j] = dist; // just for debugging
    }
  }
  // cout<<endl<<"Matrix:"<<endl;
  // // just for debugging
  // for(size_t i=0; i < frag_size; i++) {
  //   for(size_t j=0; j < frag_size; j++) {
  //     cout<<matrix[i][j]<<",";
  //   }
  //   cout<<endl;
  // }
  // cout<<endl;

  prof_vect_.set_profile_vector(prof_freq_vector);
  // if(verbose_) prof_vect_.verbose();
}

void SubModel::set_ca_coord(vector<float> ca_coord) {
  prof_vect_.set_ca_coord(ca_coord);
}


SubModelIdentity SubModel::get_sub_model_identity() {
  unsigned sub_ss = 5; //ROS:: write this function soon
  unsigned sub_seq = 5;
  return SubModelIdentity(model_name_, chain_id_, start_rsd_pos_, 
                          sub_seq_length_, sub_ss, sub_seq);
}

string SubModel::get_submodel_id() {
  if(submodel_id_.empty()) { 
    if((model_name_.empty())) return "";
    string start_rsd_pos = int2string(start_rsd_pos_, 5); 
    string gap = int2string(gap_, 3); 
    string id = model_name_ + "_" + start_rsd_pos + "_" + gap; 
    return id;
  }
  else 
    return submodel_id_;
}

void SubModel::std_out_coords(ofstream &std_out) {
  prof_vect_.std_out_coords(std_out);
}

void SubModel::set_gap_length(size_t gap) {
  gap_ = gap;
}

size_t SubModel::get_gap_length() {
  return gap_;
}

