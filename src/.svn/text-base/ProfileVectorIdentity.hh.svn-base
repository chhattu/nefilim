#ifndef PROFILE_VECTOR_IDENTITY_HEADER_HH 
#define PROFILE_VECTOR_IDENTITY_HEADER_HH 

#include<iostream>
#include<vector>

#include<SubModelIdentity.hh>

using namespace std;

class ProfileVectorIdentity {
  private:
    float cosine_;
    float vector_length_;
    float distance_;
    float ca_rmsd_;

    SubModelIdentity sub_model_identity_;

  public: 
    ProfileVectorIdentity() {
    }

    ProfileVectorIdentity(float cosine, float length, float dist, SubModelIdentity sub_model_identity) {
      cosine_ = cosine;
      vector_length_ = length;
      distance_ = dist;
      sub_model_identity_ = sub_model_identity;
    }
    
    bool operator>(ProfileVectorIdentity &prof_vect_identity) const {
      return cosine_ > prof_vect_identity.cosine_;
    } 

    bool operator<(ProfileVectorIdentity &prof_vect_identity) const {
      return cosine_ < prof_vect_identity.cosine_;
    }

    bool operator>=(ProfileVectorIdentity &prof_vect_identity) const {
      return cosine_ >= prof_vect_identity.cosine_;
    } 

    bool operator<=(ProfileVectorIdentity &prof_vect_identity) const {
      return cosine_ < prof_vect_identity.cosine_;
    }

    bool operator==(ProfileVectorIdentity &prof_vect_identity) const {
      return cosine_ == prof_vect_identity.cosine_;
    }

    void update(SubModelIdentity sub_model_identity) {
      sub_model_identity =  sub_model_identity;
    }

    void set_ca_rmsd(float ca_rmsd) {
      ca_rmsd_ = ca_rmsd;
    } 

    void verbosity() {
      cout<<"Cosine "<<cosine_<<" length "<<vector_length_
          <<" distance "<<distance_<<" rmsd "<<ca_rmsd_
          <<endl;
      // for(vector<SubModelIdentity>::iterator vit = sub_model_identites_.begin();
      //     vit != sub_model_identites_.end(); vit++) 
      //   vit->verbose(); 
    }
};

#endif
