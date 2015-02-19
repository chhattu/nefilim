#include<ModelEnergy.hh>

class core::scoring::Energies;

ModelEnergy::ModelEnergy() {}

ModelEnergy::ModelEnergy(Rosetta rosetta) {
  rosetta_ = rosetta;
}

void ModelEnergy::compute_rosetta_energy(Pose pose, vector<float> &rsd_energies) { 
                                        // vector<float> rsd_wt_energies) {
  ScoreFunctionOP score_fxn = rosetta_.get_cen_score_fxn();
  float total_energy = (*score_fxn)(pose);
  cout<<"Total energy "<<total_energy<<endl;
  // float total_energy = 0;
  size_t total_rsd   = pose.total_residue();
  cout<<"Residue length " <<total_rsd<<endl;
  Energies energies  = pose.energies();
  EnergyMap scorefxn_weights = energies.weights();

  // vector<float> rsd_energies; 
  rsd_energies.reserve(total_rsd);
  for(size_t rsd_pos = 1; rsd_pos <= total_rsd; rsd_pos++) {
    float rsd_energy = energies.residue_total_energy(rsd_pos);
    rsd_energies.push_back(rsd_energy);
  //   // EnergyMap const &emap  = energies.residue_total_energies(rsd_count);
  //   // for(EnergyMap::const_iterator it = emap.begin(), it_end = emap.end(); it != it_end; ++it) {
  //   //   ScoreType const scoretype = ScoreType( it - emap.begin() + 1 );
  //   //   if ( scorefxn_weights_[scoretype] != 0 ) {
  //   //   float weighted_score = scorefxn_weights_[scoretype] * (*it);
  //   //   total_summed_energy += weighted_score;
  //   //   }
  //   // } 
  }
}

