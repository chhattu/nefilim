// Copyright (C)  2012 Zhang Initative Research Unit
// Distance.hh -description
// written by Rojan Shrestha 

#ifndef  DISTANCE_HEADER_HH
#define  DISTANCE_HEADER_HH

#include<map>
#include<vector>
#include<iostream>

using namespace std;

class Distance {

  private:
    typedef std::pair<string, float> t_pair;
    typedef std::map<string, map<string, float> > s_map2;

    vector<SubModel> refSubModels_;

    map<unsigned int, size_t> generate_random_numbers(const size_t num_of_records, 
                                                      const size_t num_of_refs); 

    struct key_compare {
      template <typename Pair>
      bool operator() (Pair const &lhs, Pair const &rhs) const {
        return lhs.first == rhs.first;
      }
    };

    template <typename Map>
    bool map_compare (Map const &lhs, Map const &rhs) {
      return std::equal(lhs.begin(), lhs.end(), rhs.begin(), key_compare()); 
    }

    inline void sort_by_map_value(map<string, float> scores, vector<t_pair> & t_vector); 

    struct comp{ 
      bool operator()(const t_pair &pair1, const t_pair &pair2) {
        return pair1.second < pair2.second;
      }
    }; 

    inline void convert_vector_submodels_to_map(SubModel ref_submodel, 
                                                vector<SubModel> submodels, 
                                                map<string, float> &to_be_searched);
    // struct comp_string{ 
    //   bool operator()(string &s1, const string &s2) {
    //     string::const_iterator p1 = s1.begin();
    //     string::const_iterator p2 = s2.begin();

    //     while(p1 != s1.end() && p2 != s2.end()) {
    //       if(toupper(*p1) != toupper(*p2)) return (toupper(*p1) < toupper(*p2)) ? -1: 1;
    //       ++p1;
    //       ++p2;
    //     }
    //     return(s1.size() == s2.size()) ? 0 : (s.size() < s2.size())
    //   }
    // }; 


  public:
    Distance();

    inline unsigned int count_num_of_fragments(const string filepath); 

    void get_reference_submodels(const size_t num_of_refs, const string filepath); 

    void set_reference_submodels(const string pdb_path, size_t frag_size); 

    void set_reference_submodels(const string filepath); 

    void get_submodels(const size_t num_of_refs, const string filepath); 

    void get_submodels(const size_t num_of_refs, const string filepath, 
                       const unsigned int num_of_records); 

    void dump_reference_submodels(const char* filepath); 

    void compute_and_dump_carmsd(const string input_path); 

    void compute_and_dump_carmsd(const string input_path, const string output_path); 

    // void compute_and_dump_cosine(map<string, SubModel> submodels, 
    //                              s_map2 found_elements, const string output_filename); 

    void compute_and_dump_cosine(vector<SubModel> ref_submodels, map<string, SubModel> submodels, 
                                 s_map2 found_elements, const string output_filename); 

    void get_selected_submodels(s_map2 found_elements, const string input_path,  
                                map<string, SubModel> &submodels); 

    void read_reference_submodels(const char* filepath);

    map<float, vector<string> > read_distance_to_ref_submodel(string ref_file_path); 

    void read_distance_to_ref_submodel(const char* ref_file_path, 
                                       map<float, vector<string> > &distances_to_reference); 

    void read_distance_to_ref_submodel(const char* ref_file_path, 
                                       map<float, vector<string> > &distances_to_reference,
                                       float lower_bound, float upper_bound, const float threshold);

    void search(map<float, vector<string> > database, 
                vector<t_pair> to_be_searched, 
                map<string, map<string, float> > &found_fragments,
                const float threshold); 

    inline void insert_vector_in_map(string sub_model_id, float distance, 
                                 map<float, vector<string> > &distances); 

    void compute_distance_to_references(vector<SubModel> submodels, 
                                        map<string, float>  &carmsd_to_refs); 

    void compute_distance_to_references_n_search(vector<SubModel> submodels, 
                                                 map<string, map<string, float> > &found_fragments,
                                                 const float threshold); 


    inline map<string, map<string, float> > intersection(map<string, map<string, float> > lv_maps, 
                                                    map<string, map<string, float> > rv_maps); 

    inline map<string, float> map_intersection(map<string, float> l_map, 
                                               map<string, float> r_map); 

    void show_found_fragments(map<string, map<string, float> > found_fragments); 

    inline void  get_bounds(vector<t_pair> sorted_scores, float &lower_bound, float &upper_bound); 
};

#endif

