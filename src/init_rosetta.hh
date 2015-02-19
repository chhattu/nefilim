#ifndef INIT_HEADER_HH
#define INIT_HEADER_HH

//rosetta headers
#include <protocols/init.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <core/types.hh>

class Rosetta {
  private:
    int init_;
    // core::fragment::FragSetOP fragset_small_;
    // core::fragment::FragSetOP fragset_large_;
    core::scoring::ScoreFunctionOP fullatom_scorefxn_;
    core::scoring::ScoreFunctionOP centroid_scorefxn_;

  
  public:
    Rosetta(); //just for object status 

    Rosetta(int init); //just for object status 

    Rosetta(int argc, char* argv[]); 

    core::scoring::ScoreFunctionOP get_cen_score_fxn(); 

    core::scoring::ScoreFunctionOP get_fullatom_score_fxn(); 

    int get_init();  
};

#endif //end of INIT_HEADER_HH
