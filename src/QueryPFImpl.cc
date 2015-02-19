#include<QueryPFImpl.hh>

#include<iostream>
#include<time.h>

using namespace std;


QueryPFImpl::QueryPFImpl() {

}

bool QueryPFImpl::set_query_vector(const string qv_filepath,const string qc_filepath) {
  if(qv_filepath.empty() || qc_filepath.empty()) return false;
  FileIO fileio;
  fileio.read_vectors(qv_filepath, qc_filepath, qprofile_vectors_);
  cout<<"QP profile vector "<<qprofile_vectors_.size()<<endl;
  return true;
}


bool QueryPFImpl::set_database_vector(const string dv_filepath,const string dc_filepath) {
  if(dv_filepath.empty() || dc_filepath.empty()) return false;
  FileIO fileio;  
  fileio.read_vectors(dv_filepath, dc_filepath, dprofile_vectors_);
  cout<<"QP profile vector "<<dprofile_vectors_.size()<<endl;
  return true;
}


void QueryPFImpl::query() {
  if(dprofile_vectors_.empty() || qprofile_vectors_.empty()) return ;
  const unsigned int total_qvector = qprofile_vectors_.size();
  const unsigned int total_dvector = dprofile_vectors_.size();
  time_t start, end;
  for(size_t i = 0; i < total_qvector; i++) {
    // cout<<"vector index "<<i+1<<endl;
    time(&start);
    for(size_t j = 0; j < total_dvector; j++) {
      //cout<<qprofile_vectors_[i].size()<<" "<<dprofile_vectors_[j].size()<<endl;
      float dist = qprofile_vectors_[i].compute_euclidean_distance(dprofile_vectors_[j]);
      float cosine = qprofile_vectors_[i].compute_cosine(dprofile_vectors_[j]);
      float ca_rmsd = qprofile_vectors_[i].ca_rmsd(dprofile_vectors_[j].get_ca_coord());
      float corr = qprofile_vectors_[i].pearsonCC(dprofile_vectors_[j]);
      cout<<"CA_RMSD "<<ca_rmsd<<" distance "<<dist<<" cosine "<<cosine<<" pearson cc "<<corr<<endl;
    }
    cout<<endl;
    time(&end);
    cout<<"Elapsed : "<<difftime(end, start)<<" seconds."<<endl;
  }
}


