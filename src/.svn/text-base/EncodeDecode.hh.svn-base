#ifndef ENCODEDECODE_HEADER_HH
#define ENCODEDECODE_HEADER_HH

#include<vector>

using namespace std;

class EncodeDecode {

  private:
    unsigned int max_possible_value_; 
    unsigned int shift_bits_;
    unsigned int reduced_vector_size_;
    unsigned int reduction_factor_; 
    

  private:
    unsigned int inline get_max_possible_value(const unsigned int frag_size);

  public:
    EncodeDecode();

    EncodeDecode(const unsigned int frag_size, const unsigned int prof_vect_size);

    void encode(const vector<int> input_vector, vector<long> &output_vector);

    void decode(const vector<int> input_vector, vector<long> &output_vector);

    void set_frag_vect_sizes(const unsigned int frag_size, const unsigned int prof_vect_size); 

};

#endif
