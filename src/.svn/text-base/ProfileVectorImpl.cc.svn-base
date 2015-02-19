#include<fstream>

#include<ProfileVectorImpl.hh>
#include<FileIO.hh>
#include<Model.hh>

ProfileVectorImpl::ProfileVectorImpl() {
}

ProfileVectorImpl::ProfileVectorImpl(const string path2model_list) {
  path2model_list_ = path2model_list;
}

ProfileVectorImpl::ProfileVectorImpl(const string path2model_list, const size_t frag_size, 
                                     const float bin_size, const bool verbose) {
  path2model_list_ = path2model_list;
  frag_size_ = frag_size;
  bin_size_  = bin_size;
  min_dist_  = 0; 
  max_dist_  = 30.0; 
  vector_size_ = size_t(ceil((max_dist_-min_dist_) /bin_size_));
  verbose_   = verbose;

  encode_decode_.set_frag_vect_sizes(frag_size, vector_size_);
  if(verbose) {
    cout<<"Initialization: "<<endl;
    cout<<"Fragment size: "<<frag_size_<<endl;
    cout<<"Bin size: "<<bin_size_<<endl;
    cout<<"# dim: "<<vector_size_<<endl;
  }
}

void ProfileVectorImpl::compute_freq_vector_profiles() {
   FileIO fileio; 
   map<string, Pose> poses;
   fileio.read_poses_no_stout(path2model_list_, poses);
   if(poses.empty()) return;
   vector<SubModel> submodels;
   for(map<string, Pose>::iterator mit = poses.begin(); mit != poses.end(); mit++) {
     cout<<"PDB Models: "<<mit->first <<endl;
     Model model(mit->first, mit->second, verbose_);
     model.generate_sub_models(mit->second, submodels, frag_size_, bin_size_, min_dist_, vector_size_);
     if(submodels.empty()) continue; 
     vector<SubModel>::iterator vit = submodels.begin(); 
     for(;vit != submodels.end(); vit++) 
       encode_prof_vect(*vit); 
   }
   if(verbose_) {
      vector<ProfileVectorIdentity>::iterator  vit; 
      for(vit = profile_vectors_identities_.begin(); vit != profile_vectors_identities_.end(); vit++) {
        vit->verbosity();
     }
   }
}

void ProfileVectorImpl::std_out_coords(ofstream &fstream, vector<SubModel> submodels) {
  vector<SubModel>::iterator vit = submodels.begin(); 
  for(; vit != submodels.end(); vit++)  {
    string submodel_id = vit->get_submodel_id(); 
    if(submodel_id.empty()) continue;
    fstream<<submodel_id; 
    vit->std_out_coords(fstream); 
    fstream<<endl;
  }
}

void ProfileVectorImpl::filter_intra_submodels(float threshold) {
   FileIO fileio; 
   if(path2model_list_.empty()) return;
   vector<string> path2scores;
   fileio.read_model_path(path2model_list_, path2scores);
   ofstream outstream;
   outstream.open("coordinates.out", ios::out);

   if(path2scores.empty()) return; 
   vector<string>::iterator vit = path2scores.begin();
   for(; vit != path2scores.end(); vit++) {
     string model_id;
     Pose pose;
     fileio.read_pose_no_stout((*vit), model_id,  pose); 
     vector<SubModel> submodels;
     Model model(model_id, pose, verbose_);
     cout<<"Model: "<<model_id<<endl;
     model.filter_similar_intra_submodels(pose, submodels, frag_size_, threshold); 
     std_out_coords(outstream, submodels); 
     submodels.clear();
   }
   outstream.close();
}

void ProfileVectorImpl::filter_intra_submodels(vector<SubModel> &submodels, 
                                               const string pdb_path, float threshold) {
   FileIO fileio; 
   map<string, Pose> poses;
   fileio.read_pose_no_stout(pdb_path, poses);
   if(poses.empty()) return;

   map<string, Pose>::iterator mit = poses.begin(); 
   for(; mit != poses.end(); mit++) {
     cout<<"ROS Model:: "<<mit->first<<endl;
     Model model(mit->first, mit->second, verbose_);
     model.filter_similar_intra_submodels(mit->second, submodels, frag_size_, threshold); 
   }
}

void ProfileVectorImpl::compute_freq_vector_profiles(bool debug_mode, bool overlapping) {
   FileIO fileio; 
   map<string, Pose> poses;
   fileio.read_poses_no_stout(path2model_list_, poses);
   if(poses.empty()) return;

   ofstream c_file, ca_coord_file;
   c_file.open("profile_vectors.dst",ios::out);
   ca_coord_file.open("ca_coordinates.lst",ios::out);
   vector<SubModel> submodels;
   static size_t counter = 0; 

   for(map<string, Pose>::iterator mit = poses.begin(); mit != poses.end(); mit++) {
     Model model(mit->first, mit->second, verbose_);
     if(overlapping)
       model.generate_overlapping_sub_models(mit->second, submodels, frag_size_, 
                                             bin_size_, min_dist_, vector_size_);
     else
       model.generate_sub_models(mit->second, submodels, frag_size_, bin_size_, 
                                 min_dist_, vector_size_);

     if(submodels.empty()) continue; 
     for(vector<SubModel>::iterator vit = submodels.begin(); vit != submodels.end(); vit++)  {
       c_file.width(10); ca_coord_file.width(10);
       c_file<<++counter<<" "<<mit->first<<" ";
       ca_coord_file<<counter<<" "<<mit->first<<" ";
       IProfileVector model_profile = (*vit).get_prof_vect();
       if(debug_mode) {
         float distance = ref_prof_vect_.compute_euclidean_distance(model_profile);
         float ca_rmsd  = ref_prof_vect_.ca_rmsd(model_profile.get_ca_coord()); 
         c_file<<" "<<distance<<" "<<ca_rmsd<<" ";
       }
       model_profile.verbose(c_file);
       model_profile.std_out_ca_coord(ca_coord_file);
       c_file<<endl; ca_coord_file<<endl;
     }
   }
   c_file.close();
   ca_coord_file.close();

   if(verbose_) {
      vector<ProfileVectorIdentity>::iterator  vit = profile_vectors_identities_.begin(); 
      for(; vit != profile_vectors_identities_.end(); vit++) 
        vit->verbosity();
   }
}

void ProfileVectorImpl::encode_prof_vect(SubModel submodel) {
  static int count;
  IProfileVector model_profile = submodel.get_prof_vect();
  vector<long> encoded_prof_vects;
  cout<<"Initial vector: "; model_profile.verbose(); 
  float distance = ref_prof_vect_.compute_euclidean_distance(model_profile);
  float ca_rmsd  = ref_prof_vect_.ca_rmsd(model_profile.get_ca_coord()); 
  cout<<"count:"<<count<<" Eucliean Dist: "<<distance<<"  CA-RMSD "<<ca_rmsd<<endl; 
}

void ProfileVectorImpl::update_prof_vect_identity(SubModel submodel) {
  IProfileVector model_prof_vect       = submodel.get_prof_vect();
  SubModelIdentity sub_model_identity = submodel.get_sub_model_identity();
  
  float length   = model_prof_vect.compute_length(); 
  float cosine   = ref_prof_vect_.compute_cosine(model_prof_vect);
  float distance = ref_prof_vect_.compute_euclidean_distance(model_prof_vect);
  float ca_rmsd  = ref_prof_vect_.ca_rmsd(model_prof_vect.get_ca_coord()); 

  //update 
  ProfileVectorIdentity pro_vect_identity(cosine, length, distance, sub_model_identity);
  pro_vect_identity.set_ca_rmsd(ca_rmsd);
  profile_vectors_identities_.push_back(pro_vect_identity);
}

void ProfileVectorImpl::set_verbosity(const bool verbose) {
  verbose_ = verbose;
}


bool ProfileVectorImpl::set_ref_prof_vect(const string path2pdb) {
  FileIO fileio;
  string model_id;
  Pose pose;
  if(!fileio.read_pose_no_stout(path2pdb, model_id, pose)) return false;
  cout<<"Reference pdb: "<<model_id<<endl;
  cout<<"Sequenec  "<<pose.sequence()<<endl;
  cout<<"Secondary "<<pose.secstruct()<<endl;

  vector<SubModel> submodels;
  Model model(model_id, pose, verbose_);
  model.generate_overlapping_sub_models(pose, submodels, frag_size_, bin_size_, min_dist_, vector_size_);
  if(submodels.empty()) return false;
  ref_prof_vect_ = submodels[0].get_prof_vect();
  for(size_t i = 0; i< submodels.size(); i++) {
    submodels[i].verbose();
    cout<<submodels[i].get_s_sec_struct()<<" "<<submodels[i].get_s_sub_sequence()<<endl; 
  }

  cout<<"Ref vector: "; ref_prof_vect_.verbose(); cout<<endl;
  return true;
}
