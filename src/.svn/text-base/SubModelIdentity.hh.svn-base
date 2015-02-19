#ifndef SUBMODELIDENTITY_HH
#define SUBMODELIDENTITY_HH

#include<string>
#include<iostream>

using namespace std;


class SubModelIdentity {
  private:
    string model_name_;
    string chain_id_;
    int start_rsd_pos_;
    size_t sub_struct_size_; 
    unsigned sub_ss_;
    unsigned sub_sequence_;

  public:
    SubModelIdentity();

    SubModelIdentity(string, string, size_t, size_t, unsigned, unsigned);

    string get_model_name() {
      return model_name_; 
    }

    string get_chain_id() {
      return chain_id_;
    }

    int get_start_rsd_pos() {
      return start_rsd_pos_;
    }

    size_t get_size() {
      return sub_struct_size_;
    }

    unsigned get_sub_sec_struct() {
      return sub_ss_;
    }

    unsigned get_sub_seq() {
      return sub_sequence_;
    }

    void verbose() {
      cout<<"PDB "<<model_name_<<" chain "<<chain_id_
          <<"rsd_pos "<<start_rsd_pos_<<" end_pos "<<start_rsd_pos_ + sub_struct_size_
          <<endl;
    }

};

#endif

