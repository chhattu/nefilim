#ifndef MODEL_ENERGY_HEADER_HH 
#define MODEL_ENERGY_HEADER_HH 

#include<iostream>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>



#include <init_rosetta.hh>


typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
typedef core::scoring::ScoreType ScoreType;
typedef core::scoring::EnergyMap EnergyMap;
typedef core::scoring::Energies Energies;

typedef core::pose::Pose Pose;

using namespace std;

class ModelEnergy {
  
  private:
    string model_name;
    float total_ros_energy_;
    float individual_rsd_ros_energy_;
    Rosetta rosetta_;


  public:
    ModelEnergy();

    ModelEnergy(Rosetta rosetta);

    void compute_rosetta_energy(Pose pose, vector<float> &scores);
    
};

#endif
