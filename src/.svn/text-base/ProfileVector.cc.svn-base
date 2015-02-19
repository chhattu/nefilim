#include<iostream>
#include<fstream>
#include<map>
#include<vector>
#include<assert.h>

#include<ProfileVector.hh>
#include<rmsd.hh>

using namespace std;

template<class T> 
ProfileVector<T>::ProfileVector() {
}

template<class T> 
ProfileVector<T>::ProfileVector(vector<T> freq_prof_vector){ 
  freq_prof_vector_ = freq_prof_vector;
}

template<class T> 
void ProfileVector<T>::set_profile_vector(vector<T> freq_prof_vector) {
  freq_prof_vector_ = freq_prof_vector;
}

template<class T> 
vector<T> ProfileVector<T>::get_profile_vector() {
  return freq_prof_vector_;
}

template<class T>
void ProfileVector<T>::set_ca_coord(vector<float> ca_coord) {
  coords_ = ca_coord;
  const size_t mLen = ca_coord.size() / 3;
  float cx=0.0, cy=0.0, cz=0.0;
  int k3 = 0;
  for(size_t i = 0; i < mLen; i++) {
    cx += ca_coord[k3];
    cy += ca_coord[k3 + 1];
    cz += ca_coord[k3 + 2];
    k3 += 3;
  }
  cx/=mLen; cy/=mLen; cz/=mLen; k3 = 0; 
  for(size_t i = 0; i < mLen; i++) {
    ca_coord[k3]     -= cx;
    ca_coord[k3 + 1] -= cy;
    ca_coord[k3 + 2] -= cz;
    k3+=3;
  }
  ca_coord_ = ca_coord;
}

template<class T>
vector<float> ProfileVector<T>::get_ca_coord() {
  return ca_coord_; 
}

template<class T>
vector<float> ProfileVector<T>::get_coords() {
  return coords_; 
}


template<class T>
T ProfileVector<T>::sum(){
  T sum = 0;
  for(size_t i =0; i < freq_prof_vector_.size(); i++) {
    sum += freq_prof_vector_[i];  
  } 
  return sum;
}

template<class T>
void ProfileVector<T>::verbose() {
  if(freq_prof_vector_.empty()) return; 
  for(size_t i = 0; i < freq_prof_vector_.size(); i++) {
    cout<<freq_prof_vector_[i]<<" ";
  }
  cout<<endl;
}

template<class T>
void ProfileVector<T>::verbose(ofstream &c_out) {
  if(freq_prof_vector_.empty()) return; 
  for(size_t i = 0; i < freq_prof_vector_.size(); i++) {
    c_out<<freq_prof_vector_[i]<<" ";
  }
}

template<class T>
void ProfileVector<T>::std_out_ca_coord(ofstream &std_out) {
  if (ca_coord_.empty()) return;
  for(size_t i = 0; i < ca_coord_.size(); i++) {
    std_out.width(9); std_out.precision(3); std_out<<ca_coord_[i];
  }
}

template<class T>
void ProfileVector<T>::std_out_coords(ofstream &std_out) {
  if (coords_.empty()) return;
  for(size_t i = 0; i < coords_.size(); i++) {
    std_out.setf(ios::fixed, ios::floatfield);
    std_out.width(9); 
    std_out.precision(3); 
    std_out<<coords_[i];
  }
}

template<class T>
float ProfileVector<T>::compute_cosine(ProfileVector<T> model_freq_prof_vector) {
  assert(freq_prof_vector_.size() == model_freq_prof_vector.size());
  float sum_of_product = 0, sum_of_sq_ref = 0.0, sum_of_sq_model = 0.0;
  for(size_t index = 0; index < freq_prof_vector_.size(); index++) {
    sum_of_product += float(freq_prof_vector_[index]*model_freq_prof_vector[index]); 
    sum_of_sq_ref += pow(float(freq_prof_vector_[index]), 2);
    sum_of_sq_model += pow(float(model_freq_prof_vector[index]), 2);
  }
  if(sum_of_sq_ref>0 && sum_of_sq_model>0) 
    return float(sum_of_product/(sqrt(sum_of_sq_ref)*sqrt(sum_of_sq_model)));
  else
    return -1.0;
} 

template<class T>
float ProfileVector<T>::compute_length() {
  assert(freq_prof_vector_.size() > 0);
  float sum_of_diff = 0.0;
  for(size_t i = 0; i < freq_prof_vector_.size(); i++) 
    sum_of_diff += pow(freq_prof_vector_[i], 2); 
  return sqrt(sum_of_diff);
}

template<class T>
float ProfileVector<T>::compute_euclidean_distance(ProfileVector<T> model_freq_prof_vector) {
  assert(freq_prof_vector_.size() == model_freq_prof_vector.size());
  float sum_of_diff = 0, diff = 0;
  for(size_t i = 0; i < freq_prof_vector_.size(); i++) {
    diff = (freq_prof_vector_[i] - model_freq_prof_vector[i]);
    sum_of_diff += (diff*diff);
  }
  return sqrt(sum_of_diff);
}

template<class T>
float ProfileVector<T>::ca_rmsd(vector<float> coord2) {
  vector<float> coord1 = ca_coord_;
  assert(coord1.size() == coord2.size());

  const size_t mLen_ = coord1.size()/3;
  double** mPos1_ = new double*[3];
  double** mPos2_ = new double*[3];
  for(size_t k = 0; k < 3; ++k) {
    mPos1_[k] = new double[mLen_];
    mPos2_[k] = new double[mLen_];
  }
  
  int k3;
  double rmsd;
  for (size_t k = 0; k < mLen_; ++k) {
    k3 = k*3;
    mPos1_[0][k] = coord1[k3];
    mPos1_[1][k] = coord1[k3 + 1];
    mPos1_[2][k] = coord1[k3 + 2];
    mPos2_[0][k] = coord2[k3];
    mPos2_[1][k] = coord2[k3 + 1];
    mPos2_[2][k] = coord2[k3 + 2];
  }
  fast_rmsd(mPos1_, mPos2_, mLen_, &rmsd);

  for (size_t k = 0; k < 3; ++k) {
    delete [] mPos1_[k];
    delete [] mPos2_[k];
  }
  delete [] mPos1_;
  delete [] mPos2_;
  return (float)rmsd;
}

// I will do Numerical Recipe
template<class T>
T ProfileVector<T>::mean(ProfileVector<T> model_freq_prof_vector) {
  assert(model_freq_prof_vector.size() > 0);
  T N = model_freq_prof_vector.size();
  T sum = 0.0;
  for(size_t i = 0; i < (size_t)N; i++) {
    sum += model_freq_prof_vector[i];
  }
  return sum/N; 
}

template<class T>
T ProfileVector<T>::standard_deviation(ProfileVector<T> model_freq_prof_vector) {
  assert(model_freq_prof_vector.size() > 0);
  T xbar = mean(model_freq_prof_vector);
  T N = model_freq_prof_vector.size();
  T diff_sq = 0.0, diff = 0.0;
  for(size_t i = 0; i < (size_t)N; i++) {
    diff = (model_freq_prof_vector[i] - xbar);
    diff_sq += (diff * diff);
  }
  return sqrt(diff_sq/N);   
}

template<class T>
float ProfileVector<T>::pearsonCC(ProfileVector<T> model_freq_prof_vector) {
  ProfileVector<T> local_freq_prof_vector = freq_prof_vector_;
  assert(local_freq_prof_vector.size() == model_freq_prof_vector.size());
  const float TINY = 1.0e-20; 
  float N    = model_freq_prof_vector.size();
  float xbar = mean(local_freq_prof_vector);
  float ybar = mean(model_freq_prof_vector);

  float xt = 0.0, yt = 0.0;
  float sxx = 0.0, syy = 0.0, sxy = 0.0;
  for(size_t j = 0; j < N; j++) {
    xt = local_freq_prof_vector[j] - xbar; 
    yt = model_freq_prof_vector[j] - ybar; 
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  float r = (float) sxy / (sqrt(sxx*syy) + TINY);
  return r;
}

template<class T>
float ProfileVector<T>::QScore(ProfileVector<T> model_freq_prof_vector) {
  ProfileVector<T> local_freq_prof_vector = freq_prof_vector_;
  assert(local_freq_prof_vector.size() == model_freq_prof_vector.size());
  float N = model_freq_prof_vector.size();
  float diff = 0.0;
  for(size_t j = 0; j < N; j++) 
    diff += abs(local_freq_prof_vector[j]-model_freq_prof_vector[j]);
  return diff/N;
}

