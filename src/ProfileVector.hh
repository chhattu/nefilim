#ifndef PROFILEVECTOR_HEADER_HH 
#define PROFILEVECTOR_HEADER_HH 

#include<math.h>
#include<map>
#include<vector>
#include<iostream>
#include<fstream>

using namespace std;

template<class T> class ProfileVector {

  private:
    vector<T> freq_prof_vector_;
    vector<float> coords_; //This is not changed one to dump
    vector<float> ca_coord_; //Debug

  private:
    T mean(ProfileVector<T> model_freq_prof_vector);

    T standard_deviation(ProfileVector<T> model_freq_prof_vector);

  public:
    
    ProfileVector<T>(); 

    ProfileVector<T>(vector<T> freq_prof_vector); 

    void set_ca_coord(vector<float> ca_coord);  

    vector<float> get_ca_coord();  

    vector<float> get_coords(); 

    void set_profile_vector(vector<T> freq_prof_vector);

    vector<T> get_profile_vector();

    size_t size() {
      return freq_prof_vector_.size();
    }

    void verbose(); 

    void verbose(ofstream &c_out); 

    void std_out_ca_coord(ofstream &std_out);

    void std_out_coords(ofstream &std_out); 

    float compute_cosine(ProfileVector<T> model_freq_prof_vector); 

    float compute_length();

    float compute_euclidean_distance(ProfileVector<T> model_freq_prof_vector); 

    float pearsonCC(ProfileVector<T> model_freq_prof_vector); 

    float QScore(ProfileVector<T> model_freq_prof_vector); 

    bool operator==(ProfileVector<T> model_freq_prof_vector) {
      float cosine = compute_cosine(model_freq_prof_vector);
      float ref_dist = compute_length();
      float model_dist = model_freq_prof_vector.compute_length();
      return (cosine == 0) && (model_dist == ref_dist);
    }

    T sum();

    int operator[](int i) const {
      return freq_prof_vector_[i];
    }

    bool empty() {
      return freq_prof_vector_.empty(); 
    }

    float ca_rmsd(vector<float> coord1); 

    void show_ca_coord() {
      int count = 1;
      for(vector<float>::iterator vit = ca_coord_.begin(); vit != ca_coord_.end(); vit++) {
        cout<<(*vit)<<" ";
        if(count % 3 ==0)
          cout<<endl;
        count++;
      }
    }
      
};

#endif
