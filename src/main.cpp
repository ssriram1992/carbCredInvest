#include "bondad.h"
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

#define NUM_SCEN 1

int main() {
  /*

######
#     #    ##     #####    ##
#     #   #  #      #     #  #
#     #  #    #     #    #    #
#     #  ######     #    ######
#     #  #    #     #    #    #
######   #    #     #    #    #

   */

  ccI::LeadPar Chil_Par(0,  // Import limit
                        0,  // Export limit
                        80, // Min reqd consum
                        150 // Carbon credit
  );
  ccI::LeadPar Arge_Par(0,  // Import limit
                        0,  // Export limit
                        80, // Min reqd consum
                        0   // Carbon credit
  );

  ccI::FollPar<NUM_SCEN> Chi_Foll(
      {0.8, 4, 5},  // Quad cost
      {10, 12, 12}, // Lin Costs
      {             // std::vector<double>{32743, 11956, 9729},
       std::vector<double>{32743, 11956, 9729}}, // Capacities
      {0, 0, 0},                                 // Emission cost
      {"Coal", "NG", "Green"}                    // Names
  );

  ccI::FollPar<NUM_SCEN> Arg_Foll(
      {0.6, 0.8, 5.5}, // Quad cost
      {10, 12, 12},    // Lin Costs
      {                // std::vector<double>{32743*3, 11956*3.5, 9729*3},
       std::vector<double>{32743 * 3, 11956.0 * 3.5, 9729 * 3}}, // Capacities
      {0, 0, 0},              // Emission cost
      {"Coal", "NG", "Green"} // Names
  );

  std::array<ccI::DemPar, NUM_SCEN> Chil_Dem = {// ccI::DemPar{90e3, 0.95},
                                                ccI::DemPar{90e3, 0.95}};

  std::array<ccI::DemPar, NUM_SCEN> Arge_Dem = {// ccI::DemPar{160e3, 0.4},
                                                ccI::DemPar{160e3, 0.4}};

  ccI::LeadAllPar<NUM_SCEN> Chile(3,                    // Number of followers
                                  std::string("Chile"), // Name
                                  Chi_Foll,             // Follower details
                                  Chil_Dem,             // Demand details
                                  Chil_Par,             // Leader Params
                                  // {0.5, 0.5}            // Probabilities
                                  {1});

  ccI::LeadAllPar<NUM_SCEN> Argentina(3, // Number of followers
                                      std::string("Argentina"), // Name
                                      Arg_Foll, // Follower details
                                      Arge_Dem, // Demand details
                                      Arge_Par, // Leader Params
                                                // {0.5, 0.5} // Probabilities
                                      {1});

  /*

  #     #
  ##   ##   ####   #####   ######  #
  # # # #  #    #  #    #  #       #
  #  #  #  #    #  #    #  #####   #
  #     #  #    #  #    #  #       #
  #     #  #    #  #    #  #       #
  #     #   ####   #####   ######  ######

   */
  boost::log::core::get()->set_filter(boost::log::trivial::severity >=
                                      boost::log::trivial::info);
  GRBEnv env;
  ccI::EPEC<NUM_SCEN> epec(&env);
  epec.addCountry(Chile).addCountry(Argentina);
  epec.finalize();
  epec.setAlgorithm(Game::EPECalgorithm::innerApproximation);
  // epec.setAlgorithm(Game::EPECalgorithm::combinatorialPNE);
  epec.setAggressiveness(50);
  epec.setAddPolyMethod(Game::EPECAddPolyMethod::random);
  epec.setNumThreads(3);
  epec.findNashEq();
  epec.writeSolution(1, "dat/ChileArgentinaCarbMix");

  return 0;
}
