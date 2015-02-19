#ifndef SUBMODEL_HEADER_HH 
#define SUBMODEL_HEADER_HH 

#include<vector>
#include<sstream>
#include<iostream>
#include<iomanip>

#include<ProfileVector.cc>
#include<SubModelIdentity.hh>

using namespace std;

class SubModel {

  private:
    string model_name_;
    string chain_id_;
    int start_rsd_pos_;
    size_t sub_seq_length_; 
    size_t gap_;

    string submodel_id_;

    unsigned int sec_struct_;
    unsigned int subsequence_;

    string s_sec_struct_;
    string s_subsequence_;

    ProfileVector<int> prof_vect_;
    ProfileVector<float> fprof_vect_;
    ProfileVector<long> encoded_prof_vect_;

    bool verbose_;// = false;  
  
  private:
    //TODO:
    /*
     * make look up table using char
     * */
    map<string, short> aa_lookup_map;
    map<string, short> ss_lookup_map;

    inline void set_lookup_map() {
      ss_lookup_map["H"] = 0; ss_lookup_map["E"] = 1; ss_lookup_map["L"] = 2;
      // Nonpolar, aliphatic
      aa_lookup_map["G"] = 0; aa_lookup_map["A"] = 0; aa_lookup_map["P"] = 0; 
      aa_lookup_map["V"] = 0; aa_lookup_map["L"] = 0; aa_lookup_map["I"] = 0; 
      aa_lookup_map["M"] = 0; 
      //Aromatic
      aa_lookup_map["F"] = 1; aa_lookup_map["Y"] = 1; aa_lookup_map["W"] = 1; 
      //Polar and unchared
      aa_lookup_map["S"] = 2; aa_lookup_map["T"] = 2; aa_lookup_map["C"] = 2; 
      aa_lookup_map["N"] = 2; aa_lookup_map["Q"] = 2; 
      //Positively charged
      aa_lookup_map["K"] = 3; aa_lookup_map["H"] = 3; aa_lookup_map["R"] = 3; 
      //Negatively charged
      aa_lookup_map["D"] = 4; aa_lookup_map["E"] = 4; 
    }

    enum EncodeWhat { sec_struct, sequence }; 

    EncodeWhat encode_what_;

    unsigned int encode_string(const string input_string);

    string inline int2string(int input_int) {
      ostringstream bridge;
      bridge << input_int;
      return bridge.str();
    }

    string inline int2string(int input_int, int length) {
      char prev;
      ostringstream bridge;
      bridge.width(length);
      prev = bridge.fill('0');
      bridge << input_int;
      bridge.fill(prev);
      return bridge.str();
    }

  public:
    SubModel();

    SubModel(const string line);

    SubModel(string model_name, string chain_id, size_t start_rsd_pos); 

    SubModel(string model_name, string chain_id, size_t start_rsd_pos, size_t gap); 

    SubModel(string model_name, string chain_id, string frag_seq, 
             string frag_ss, size_t start_rsd_pos);

    void compute_frequency_profile_vector(vector<float> xyz, const size_t frag_size, 
          const float bin_size, const float min_dist, const size_t prof_vect_size); 

    SubModelIdentity get_sub_model_identity();

    ProfileVector<int> get_prof_vect() {
      return prof_vect_;
    } 

    ProfileVector<float> get_fprof_vect() {
      return fprof_vect_;
    } 

    ProfileVector<long> get_encoded_prof_vect() {
      return encoded_prof_vect_;
    } 

    vector<float> get_center_of_mass(vector<float> frag_xyzs);

    float get_carmsd(vector<float> xyzs); 

    float compute_carmsd(SubModel submodel); 

    void set_distance_vector(vector<float> xyz, const size_t fragment_size); 

    float compute_distance_vector_cosine(SubModel submodel, const size_t frag_size); 

    float compute_QScore(SubModel submodel, const size_t frag_size); 

    float compute_cosine(SubModel submodel, const size_t frag_size, const float bin_size, 
                         const float min_dist, const size_t vector_size); 

    void get_distance_vector(vector<float> xyz, 
                             const size_t fragment_size, 
                             vector<float> & distance_vector); 

    void get_matrix(vector<float> xyz, const size_t fragment_size);

    void set_ca_coord(vector<float> ca_coord);

    unsigned int get_sec_struct() {
      return sec_struct_;
    }

    unsigned int get_sub_sequence() {
      return subsequence_;
    }

    string get_s_sec_struct() {
      return s_sec_struct_;
    }

    string get_s_sub_sequence() {
      return s_subsequence_;
    }

    void verbose() {
      cout<<"Model name "<<model_name_<<" rsd_pos "
          <<start_rsd_pos_<< " length "<<start_rsd_pos_ <<endl;
    }

    float compute_local_carmsd(vector<float> xyz);

    string get_submodel_id(); 
    
    void std_out_coords(ofstream &std_out); 

    void set_gap_length(size_t gap);

    size_t get_gap_length();

};

#endif
