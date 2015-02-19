#include <init_rosetta.hh>

Rosetta::Rosetta() {
}

Rosetta::Rosetta(int init) {
  init_ = init; 
}

Rosetta::Rosetta(int argc, char* argv[]) {
  protocols::init(argc, argv); 
  fullatom_scorefxn_ = core::scoring::getScoreFunction(true);
  centroid_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("score5");
  init_ = 1;
}

core::scoring::ScoreFunctionOP Rosetta::get_cen_score_fxn() {
  return centroid_scorefxn_;
}

core::scoring::ScoreFunctionOP Rosetta::get_fullatom_score_fxn() {
  return fullatom_scorefxn_;
}

int Rosetta::get_init() {
  return init_;
}
