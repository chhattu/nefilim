// Copyright (C)  2012 Zhang Initative Research Unit
// Distance.cc -description
// written by Rojan Shrestha 

#include<vector>
#include<time.h>
#include<stdlib.h>
#include<iostream>
#include<time.h>
#include<fstream>
#include<sstream>
#include<string>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include<Model.hh>
#include<FileIO.hh>
#include<SubModel.hh>
#include<Distance.hh>

using namespace std;

Distance::Distance() {
  refSubModels_.clear(); 
}

map<unsigned int, size_t> Distance::generate_random_numbers(const size_t num_of_records, 
                                                            const size_t num_of_refs) {
  map<unsigned int, size_t> random_numbers;
  srand(time(NULL));  
  size_t i = 1;
  while(i<=num_of_refs) {
    unsigned int random_number = rand() % num_of_records + 1; 
    if(random_number > 0) {
      random_numbers[random_number]=1;
      i++;
    }
  }
  return random_numbers;
}

void Distance::set_reference_submodels(const string pdb_path, size_t frag_size) {
  typedef core::pose::Pose Pose;
  vector<SubModel> submodels;
  map<string, Pose> poses;
  FileIO fileio;
  fileio.read_pose_no_stout(pdb_path, poses); 
  map<string, Pose>::iterator mit=poses.begin();
  Model model(pdb_path, mit->second);
  model.filter_similar_intra_submodels(mit->second, submodels, frag_size, 0.5); 
  if(refSubModels_.size() > 0) 
    refSubModels_.clear();
  refSubModels_ = submodels;
}

void Distance::set_reference_submodels(const string filepath) {
  try {
    ifstream input_file(filepath.c_str(), ifstream::in);
    if(input_file.is_open()) {
      string line;
      while(getline(input_file, line)) 
        refSubModels_.push_back(SubModel(line));
      input_file.close();
    }
    else {
      cout<<"Fatal: references are not found."<<endl;
      cout<<__FILE__<<" line number "<<__LINE__<<endl;
      exit(1);
    }
  }
  catch(exception &e) { 
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number "<<__LINE__<<endl;
    exit(1);
  }
}

inline unsigned int Distance::count_num_of_fragments(const string filepath) {
  try {
    ifstream input_file(filepath.c_str(), ifstream::in);
    unsigned int counter = 0;
    if(input_file.is_open()) {
      string line;
      while(getline(input_file, line))
        counter++;
      input_file.close();
    }
    return counter; 
  }
  catch(exception &e) { 
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number "<<__LINE__<<endl;
    return 0;
    exit(1);
  }
}

void Distance::get_reference_submodels(const size_t num_of_refs, const string filepath) {
  const unsigned int num_of_records = count_num_of_fragments(filepath); 
  map<unsigned int, size_t> rand_nums = generate_random_numbers(num_of_records, num_of_refs);
  try {
    ifstream input_file(filepath.c_str(), ifstream::in);
    if(input_file.is_open()) {
      refSubModels_.reserve(rand_nums.size());
      string line;
      unsigned int counter = 1;
      map<unsigned int, size_t>::iterator mit = rand_nums.begin();
      while(getline(input_file, line)) {
        if(counter++ == mit->first) {
          refSubModels_.push_back(SubModel(line));
          mit++;
        }
      }
    }
    input_file.close();
  }
  catch(exception &e) { 
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number "<<__LINE__<<endl;
    exit(1);
  }
}

void Distance::dump_reference_submodels(const char* filepath) {
  try {
    ofstream output_file(filepath, ifstream::out);
    vector<SubModel>::iterator vit = refSubModels_.begin();
    while( vit != refSubModels_.end()) {
      output_file<<vit->get_submodel_id(); 
      vit->std_out_coords(output_file); 
      output_file<<endl;
      vit++;
    }
    output_file.close();
  }
  catch(exception &e) {
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number "<<__LINE__<<endl;
    exit(1);
  }
}

void Distance::compute_and_dump_carmsd(const string input_path) {
  vector<SubModel>::iterator it = refSubModels_.begin(); 
  for(; it != refSubModels_.end(); it++) {
    try {
      ifstream input_stream(input_path.c_str(), ifstream::in);
      string submodel_id = it->get_submodel_id(); 
      size_t found_place = submodel_id.find(".pdb"); 
      string filename; 
      if(found_place > 0 && found_place <= submodel_id.size()) { 
        filename = submodel_id.replace(found_place, 4, ".out");
        // filename = it->get_submodel_id() + ".out"; 
      }
      ofstream output_stream(filename.c_str());
      string line; 
      if(input_stream.is_open()) { 
        while(getline(input_stream, line)) {
          SubModel submodel(line);
          output_stream<<submodel.get_submodel_id(); 
          output_stream.setf(ios::fixed, ios::floatfield);
          output_stream.width(7);
          output_stream.precision(2);
          output_stream<<it->compute_carmsd(submodel);
          output_stream<<endl;
        }
      }
      input_stream.close();
      output_stream.close();
    }
    catch(exception &e) {
      cout<<e.what()<<endl;
      cout<<__FILE__<<" line number: "<<__LINE__<<endl;
      exit(1);
    }
  }
}

void Distance::compute_and_dump_carmsd(const string input_path, const string output_path) {
  try {
    ifstream input_stream(input_path.c_str(), ifstream::in);
    ofstream output_stream(output_path.c_str());
    string line; 
    if(input_stream.is_open()) { 
      while(getline(input_stream, line)) {
        SubModel submodel(line);
        output_stream<<submodel.get_submodel_id(); 

        for(vector<SubModel>::iterator it = refSubModels_.begin(); 
                                       it != refSubModels_.end(); 
                                       it++) {
          output_stream.setf(ios::fixed, ios::floatfield);
          output_stream.width(7);
          output_stream.precision(2);
          output_stream<<it->compute_carmsd(submodel);
        }
        output_stream<<endl;
      }
    }
    input_stream.close();
    output_stream.close();
  }
  catch(exception &e) {
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number: "<<__LINE__<<endl;
    exit(1);
  }
}

void Distance::read_reference_submodels(const char* ref_filepath) {
  try {
    ifstream input_stream(ref_filepath, ifstream::in);
    if(input_stream.is_open()) { 
      string line;
      while(getline(input_stream, line)) 
        refSubModels_.push_back(SubModel(line));
      input_stream.close();
    }
    else {
      cout<<"Fatal: file is not found."<<endl;
      cout<<__FILE__<<" line number: "<<__LINE__<<endl;
      exit(1);
    }
  }
  catch(exception &e) {
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number: "<<__LINE__<<endl;
    exit(1);
  }
}

map<float, vector<string> > Distance::read_distance_to_ref_submodel(string ref_submodel_id) {
  map<float, vector<string> > distances_to_reference;
  try {
    const string ref_file_path = ref_submodel_id + ".out";
    cout<<"Ref "<<ref_file_path <<endl;
    ifstream input_stream(ref_file_path.c_str(), ifstream::in);
    if(input_stream.is_open()) { 
      string line;
      while(getline(input_stream, line)) { 
        stringstream str_stream(line); 
        string sub_model_id;
        float distance = 0.0;
        str_stream>>sub_model_id>>distance;
        insert_vector_in_map(sub_model_id, distance, distances_to_reference); 
      }
      input_stream.close();
    }
    else {
      cout<<"Fatal: reference distance file is not found."<<endl;
      exit(1);
    }
  }
  catch(exception &e) {
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number: "<<__LINE__<<endl;
    exit(1);
  }
  return distances_to_reference;
}

void Distance::read_distance_to_ref_submodel(const char* ref_file_path, 
                                    map<float, vector<string> > &distances_to_reference) {
  try {
    ifstream input_stream(ref_file_path, ifstream::in);
    if(input_stream.is_open()) { 
      string line;
      float predistance = 0.0;
      bool isFirst = true;
      vector<string> submodel_ids;
      submodel_ids.reserve(10000);
      while(getline(input_stream, line)) { 
        stringstream str_stream(line); 
        string submodel_id;
        float distance = 0.0;
        if(isFirst) {
          predistance = distance;
          isFirst      = false;
        }
        str_stream>>submodel_id>>distance; 
        if(predistance == distance) 
          submodel_ids.push_back(submodel_id);
        else {
          distances_to_reference[predistance] = submodel_ids;
          predistance = distance;
          submodel_ids.clear();
          submodel_ids.reserve(10000);
          submodel_ids.push_back(submodel_id);
        }
      }
    }
    input_stream.close();
  }
  catch(exception &e) {
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number: "<<__LINE__<<endl;
    exit(1);
  }
}

void Distance::read_distance_to_ref_submodel(const char* ref_file_path, 
                              map<float, vector<string> > &distances_to_reference,
                              float lower_bound, float upper_bound, const float threshold) {
  try {
    ifstream input_stream(ref_file_path, ifstream::in);
    if(input_stream.is_open()) { 
      string line;
      float predistance = 0.0;
      bool isFirst = true;
      vector<string> submodel_ids; submodel_ids.reserve(10000);
      lower_bound -= threshold;
      upper_bound += threshold;
      // cout<<"Lower and upper with threshold: "
      //     <<lower_bound<<" "<<upper_bound<<" "<<threshold<<endl;

      while(getline(input_stream, line)) { 
        stringstream str_stream(line); 
        string submodel_id;
        float distance = 0.0;
        str_stream>>submodel_id>>distance; 
        if(isFirst) {
          predistance = distance;
          isFirst     = false;
        }
        if(distance > upper_bound) break;
        if(distance < lower_bound) continue;
        if(predistance == distance) 
          submodel_ids.push_back(submodel_id);
        else {
          distances_to_reference[predistance] = submodel_ids;
          predistance = distance;
          submodel_ids.clear();
          submodel_ids.reserve(10000);
          submodel_ids.push_back(submodel_id);
        }
      }
      cout<<"# query hits: "<<distances_to_reference.size()<<endl;
      input_stream.close();
    }
    else {
      cout<<ref_file_path<<" is not found."<<endl;
      cout<<__FILE__<<" line number: "<<__LINE__<<endl;
      exit(1);
    }
  }
  catch(exception &e) {
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number: "<<__LINE__<<endl;
    exit(1);
  }
}

inline void Distance::insert_vector_in_map(string submodel_id, float distance, 
                                 map<float, vector<string> > &distances) {
  pair<map<float, vector<string> >::iterator, bool> it;
  vector<string> submodel_ids;
  submodel_ids.push_back(submodel_id);
  it = distances.insert(pair<float, vector<string> >(distance, submodel_ids));
  if(!it.second) {
    submodel_ids = distances[distance];
    submodel_ids.push_back(submodel_id);
    distances[distance] = submodel_ids;
  }
  else { 
    distances[distance] = submodel_ids;
  }
  // debugging
  // cout<<"Size "<<distances.size()<<endl;
  // map<float, vector<string> >::iterator mit = distances.begin();
  // while(mit != distances.end()) {
  //   vector<string> submodel_ids; 
  //   cout<<"Fragments are: ";
  //   cout<<mit->first<<" ";
  //   submodel_ids = mit->second;
  //   vector<string>::iterator vit = submodel_ids.begin();
  //   while(vit != submodel_ids.end()) {
  //     cout<<(*vit)<<" ";
  //     vit++;
  //   }
  //   mit++;
  //   cout<<endl;
  // }
} 

inline void Distance::sort_by_map_value(map<string, float> scores, vector<t_pair> & t_vector) {
  if(!t_vector.empty()) t_vector.clear();
  if(scores.empty()) {
    t_vector.clear();
    return;
  }
  vector<t_pair> s_vector(scores.begin(), scores.end());
  sort(s_vector.begin(), s_vector.end(), comp());
  t_vector = s_vector;
}


inline void  Distance::get_bounds(vector<t_pair> sorted_scores, float &lower_bound, float &upper_bound) { 
  lower_bound = 0; upper_bound = 0; 

  vector<t_pair>::iterator first_mit = sorted_scores.begin();
  vector<t_pair>::iterator last_mit  = sorted_scores.end();
  last_mit--;   

  lower_bound = first_mit->second; 
  upper_bound = last_mit->second;  
}

void  Distance::search(map<float, vector<string> > database, 
                       vector<t_pair> sorted_scores, 
                       map<string, map<string, float> > &found_fragments, 
                       const float threshold) {
  found_fragments.clear();

  size_t database_size = 0;
  map<float, vector<string> >::iterator rit = database.begin();
  for(;rit != database.end(); rit++) 
    database_size += rit->second.size();
  cout<<"# fragments in database      : "<<database_size<<endl;
  cout<<"# fragments to be search (#) : "<<sorted_scores.size()<<endl;

  map<float, vector<string> >::iterator current_mvit = database.begin();
  vector<t_pair>::iterator mit = sorted_scores.begin(); 
  map<float, vector<string> >::iterator first_mvit, last_mvit;

  while(mit != sorted_scores.end()) {
    first_mvit = database.lower_bound(mit->second - threshold);
    last_mvit  = database.upper_bound(mit->second + threshold);

    map<string, float> inner_map;
    current_mvit = first_mvit; 
    while( current_mvit != last_mvit) {
      vector<string> submodel_ids = current_mvit->second;
      vector<string>::iterator vit = submodel_ids.begin(); 
      for(; vit != submodel_ids.end(); vit++) 
        inner_map[(*vit)] = current_mvit->first;
      current_mvit++;
    }
    found_fragments[mit->first] = inner_map;
    mit++;
  }
}

inline void Distance::convert_vector_submodels_to_map(SubModel ref_submodel, 
                                                      vector<SubModel> submodels, 
                                                      map<string, float> &to_be_searched){
  float ca_rmsd = 0;
  for(vector<SubModel>::iterator vit = submodels.begin(); vit != submodels.end(); vit++) {
    ca_rmsd = ref_submodel.compute_carmsd((*vit));
    to_be_searched[vit->get_submodel_id()] = ca_rmsd;
    cout<<"Distance to reference:  "<<ca_rmsd<<endl; 
  }
}

void Distance::compute_distance_to_references_n_search(vector<SubModel> submodels, 
                                          map<string, map<string, float> > &found_fragments,
                                          const float threshold) {
  vector<SubModel>::iterator it = refSubModels_.begin(); 
  float upper_bound = 0, lower_bound = 0;

  while(it != refSubModels_.end()) {
    vector<t_pair> sorted_scores;
    map<string, float> to_be_searched; 
    cout<<"Reference models:       "<<it->get_submodel_id().c_str()<<endl;
    convert_vector_submodels_to_map((*it), submodels, to_be_searched); 
    sort_by_map_value(to_be_searched, sorted_scores);
    get_bounds(sorted_scores, lower_bound, upper_bound); 

    clock_t start = clock(), diff;
    cout<<"Reading reference file: "<<it->get_submodel_id()<<endl;
    cout<<"Lower and upper bounds are: "<<lower_bound<<" "<<upper_bound;
    cout<<endl;

    map<float, vector<string> > distances_to_reference;  
    string submodel_path = it->get_submodel_id() + ".out";
    // cout<<"Path to submodel "<<submodel_path<<endl;
    read_distance_to_ref_submodel(submodel_path.c_str(), distances_to_reference, 
                                  lower_bound, upper_bound, threshold); 

    clock_t end = clock();
    diff = end - start;
    cout<<"CPU time to read the file: "<<diff/(double)CLOCKS_PER_SEC<<" seconds."<<endl;
    start = clock();
    map<string, map<string, float> > inter_found_fragments; 

    search(distances_to_reference, sorted_scores, inter_found_fragments, threshold); 

    //just for debugging...
    // show_found_fragments(inter_found_fragments); 

    end = clock();
    diff = end - start;
    cout<<"CPU time taken to search "<<diff / (double)CLOCKS_PER_SEC <<" seconds."<<endl;
    start = clock();
    found_fragments = intersection(found_fragments, inter_found_fragments); 
    end = clock();
    diff = end - start;
    cout<<"CPU time taken to make set "<<diff / (double)CLOCKS_PER_SEC <<" seconds."<<endl;
    it++;
    cout<<endl;
  }
}

void Distance::compute_distance_to_references(vector<SubModel> submodels, 
                                              map<string, float>  &carmsd_to_refs) {
  vector<SubModel>::iterator it = refSubModels_.begin(); 
  while(it != refSubModels_.end()) {
    map<string, float> distances; 
    vector<SubModel>::iterator vit = submodels.begin();
    for(; vit != submodels.end(); vit++) 
      distances[vit->get_submodel_id()] = it->compute_carmsd((*vit));
    it++;
  }
}

map<string, map<string, float> > inline Distance::intersection(
                                  map<string, map<string, float> > lv_maps, 
                                  map<string, map<string, float> > rv_maps) {

  map<string, map<string, float> > intersections;
  map<string, map<string, float> >::iterator vmit1, vmit2;  
  if(lv_maps.size() == 0) {
    intersections = rv_maps;
    return intersections;
  }
  vmit1 = lv_maps.begin();
  vmit2 = rv_maps.begin();
  while(vmit1 != lv_maps.end() && vmit2 != rv_maps.end()) {
    intersections[vmit1->first] = (map_intersection(vmit1->second, vmit2->second)); 
    vmit1++; vmit2++;
  } 
  return intersections;
}

map<string, float> inline Distance::map_intersection(map<string, float> l_map, 
                                                     map<string, float> r_map) {
  map<string, float>::iterator mit1, mit2;
  mit1 = l_map.begin(); 
  mit2 = r_map.begin();
  map<string, float> intersections;
  while(mit1 != l_map.end() &&  mit2 != r_map.end()) {
    if(mit1->first == mit2->first) {
      intersections[mit1->first] = mit1->second; 
      mit1++; mit2++;
    }
    else {
      if(mit1->first > mit2->first) 
        mit2++;
      else 
        mit1++;
    }
  }
  return intersections;
}

void Distance::show_found_fragments(map<string, map<string, float> > found_fragments) {
  if(found_fragments.empty()) {
    cout<<"Warning: not found!!!"<<endl;
    return;
  }

  map<string, map<string, float> >::iterator mmit = found_fragments.begin();
  while(mmit != found_fragments.end()) {
    map<string, float> found_submodels = mmit->second;
    map<string, float>::iterator mit = found_submodels.begin();
    cout<<"Query submodel: "<<mmit->first<<endl;
    if(!found_submodels.empty()) 
      cout<<"Query hits: "<<found_submodels.size()<<endl;
    while(mit != found_submodels.end()) {
      cout<<fixed;
      cout<<setprecision (2);
      cout<<"  "<<mit->first<<endl;
      mit++;
    }
    mmit++;
  }
}

void Distance::get_selected_submodels(s_map2 found_elements, const string input_path,  
                                      map<string, SubModel> &submodels) {
  submodels.clear();
  vector<string> submodel_ids;
  // map<string, int> m_submodel_ids;
  s_map2::iterator smit2 = found_elements.begin();
  for(; smit2 != found_elements.end(); smit2++) {
    // cout<<"Reference model "<<smit2->first<<endl;
    map<string, float> selected_submodels = smit2->second;
    map<string, float>::iterator mit = selected_submodels.begin();
    for(; mit != selected_submodels.end(); mit++) { 
      // m_submodel_ids[mit->first] = 1;
      submodel_ids.push_back(mit->first); 
    }
  }
  sort(submodel_ids.begin(), submodel_ids.end());
  // map<string,int>::iterator mit = m_submodel_ids.begin();
  // for(;mit != m_submodel_ids.end(); mit++) {
  //   cout<<"Sorted: "<<mit->first<<endl;
  // }

  ifstream input_stream(input_path.c_str(), ifstream::in);
  try {
    string line; 
    if(input_stream.is_open()) { 
      vector<string>::iterator vit = submodel_ids.begin();
      while(getline(input_stream, line)) {
        if(vit==submodel_ids.end()) break;
        size_t found=line.find((*vit));
        if(found >= 0 && found <= line.size()) {
          SubModel submodel(line);
          submodels[(*vit)] = submodel;
          vit++;
        }
      }
    }
    else {
      cout<<"No coordinate file."<<endl;
      cout<<__FILE__<<" line number: "<<__LINE__<<endl;
      exit(1);
    }
  }
  catch(exception e) {
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number: "<<__LINE__<<endl;
    exit(1);
  }
}

void Distance::compute_and_dump_cosine(vector<SubModel> ref_submodels,  
                                       map<string, SubModel> submodels, 
                                       s_map2 found_elements, 
                                       const string output_filename) {
  try {
    ofstream output_stream(output_filename.c_str());
    vector<SubModel>::iterator vit = ref_submodels.begin();
    for(; vit != ref_submodels.end(); vit++) {
      SubModel ref_submodel = (*vit); 
      map<string, float> submodel_ids  = found_elements[vit->get_submodel_id()];
      map<string, float>::iterator mit = submodel_ids.begin();
      for(; mit != submodel_ids.end(); mit++) { 
        SubModel submodel = submodels[mit->first];
        // cout<<"Submodel id "<<mit->first<<" ; "<<submodel.get_submodel_id()<<endl;
        // output_stream<<ref_submodel.get_submodel_id()<<" ";
        output_stream<<submodel.get_submodel_id()<<" ";
        output_stream.setf(ios::fixed, ios::floatfield);
        output_stream.width(7); output_stream.precision(2);
        output_stream<<ref_submodel.compute_distance_vector_cosine(submodel, 9); 
        output_stream.setf(ios::fixed, ios::floatfield);
        output_stream.width(7); output_stream.precision(2);
        output_stream<<ref_submodel.compute_cosine(submodel, 9, 0.95, 0, 20);
        output_stream.setf(ios::fixed, ios::floatfield);
        output_stream.width(7); output_stream.precision(2);
        output_stream<<ref_submodel.compute_QScore(submodel, 9);
        output_stream<<endl;
      }
    } 
    output_stream.close();
  }
  catch(exception &e) {
    cout<<e.what()<<endl;
    cout<<__FILE__<<" line number: "<<__LINE__<<endl;
    exit(1);
  }
}

