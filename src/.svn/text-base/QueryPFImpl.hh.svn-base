#ifndef QUERY_PROFILE_VECTOR_IMPL_HEADER_HH
#define QUERY_PROFILE_VECTOR_IMPL_HEADER_HH

#include<vector>


#include<ProfileVector.hh>
#include<FileIO.hh>

using namespace std;

class QueryPFImpl {

  private:
    typedef vector<ProfileVector<int> > IProfileVectors;
    typedef ProfileVector<int>  IProfileVector;
    IProfileVectors dprofile_vectors_;
    IProfileVectors qprofile_vectors_;

  public: 
    QueryPFImpl();

    bool set_query_vector(const string qv_filepath, const string qc_filepath);

    bool set_database_vector(const string dv_filepath, const string dc_filepath);

    void query();

};

#endif 
