#include<EncodeDecode.hh>

#include<iostream>
#include<cmath>

using namespace std;

EncodeDecode::EncodeDecode() {
}

EncodeDecode::EncodeDecode(const unsigned int frag_size, const unsigned int prof_vect_size) {
  const unsigned int bit_in_byte = 8;
  max_possible_value_ = get_max_possible_value(frag_size);
  shift_bits_ = (max_possible_value_ / bit_in_byte); 
  unsigned long size_of_long; 
  size_of_long = sizeof(size_of_long);

  reduced_vector_size_ = (unsigned int)floor(size_of_long / shift_bits_);
  reduction_factor_    = (unsigned int)ceil(prof_vect_size / shift_bits_); 
  if(reduction_factor_ > 8) reduction_factor_ = 8;
}

unsigned int inline EncodeDecode::get_max_possible_value(const unsigned int frag_size) {
  unsigned int value = (frag_size * (frag_size-1)) / 2;
  return value - (frag_size - 1);
}

void EncodeDecode::encode(const vector<int> input_vector, vector<long> &output_vector) {
  output_vector.clear();
  output_vector.reserve(reduced_vector_size_);
  unsigned long reduced_value = 0; 
  size_t counter = 0;
  cout<<"New vector size="<<reduced_vector_size_<<", bit shifts="<<shift_bits_<<", ";
  cout<<"reduction factor="<<reduction_factor_<<endl;
  for(vector<int>::const_iterator uit = input_vector.begin(); 
                                  uit != input_vector.end(); uit++) {
    reduced_value = reduced_value<<shift_bits_;
    reduced_value = reduced_value | (*uit); 
    if(++counter % reduction_factor_ == 0) {
      output_vector.push_back(reduced_value);
      reduced_value = 0;
    }
  }
  if(reduced_value > 0) 
    output_vector.push_back(reduced_value);
  cout<<reduced_value<<endl;
}

void EncodeDecode::set_frag_vect_sizes(const unsigned int frag_size, const unsigned int prof_vect_size) {
  const unsigned int bit_in_byte = 8;
  max_possible_value_ = get_max_possible_value(frag_size);
  shift_bits_ = (max_possible_value_ / bit_in_byte); 
  unsigned long size_of_long; 
  size_of_long = sizeof(size_of_long);

  reduced_vector_size_ = (unsigned int) floor(size_of_long / shift_bits_);
  reduction_factor_ = (unsigned int) ceil(prof_vect_size / shift_bits_); 
  if(reduction_factor_ > 8) reduction_factor_ = 8;
}

