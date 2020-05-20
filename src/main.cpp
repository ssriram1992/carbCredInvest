#include "carbCredInv.h"
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

template <class t, class u>
std::map<t, u> makeMap(const std::vector<t> keys, const std::vector<u> values) {
  std::map<t, u> mm{};
  if (keys.size() != values.size()) {
    BOOST_LOG_TRIVIAL(fatal)
        << "makeMap(): Potential error. Inconsistent sizes for keys and values";
    return mm;
  }
  const auto sz = keys.size();
  for (unsigned int ii = 0; ii < sz; ++ii)
    mm[keys.at(ii)] = values.at(ii);

  return mm;
}

template <typename t> void print(const std::vector<t> &v) {
  for (auto const &x : v)
    std::cout << x << "\t";
  std::cout << "\n";
}

template <typename t, unsigned long int N>
void print(const std::array<t, N> &v) {
  for (auto const &x : v)
    std::cout << x << "\t";
  std::cout << "\n";
}

template <typename t, typename u> void print(const std::map<t, u> &mm) {
  for (auto const &x : mm)
    std::cout << x.first << "\t-\t" << x.second << "\n";
}

int main(int argc, char *argv[]) {
  // First argument 1 for full enumeration
  using std::array;
  using std::cout;
  using std::map;
  using std::pair;
  using std::string;
  using std::vector;
  using linQuad = pair<double, double>;
  /*

######
#     #    ##     #####    ##
#     #   #  #      #     #  #
#     #  #    #     #    #    #
#     #  ######     #    ######
#     #  #    #     #    #    #
######   #    #     #    #    #

   */

  // Template parameters
  constexpr unsigned int NUM_CLEAN{2};
  constexpr unsigned int NUM_DIRTY{1};
  constexpr unsigned int NUM_SCEN{2};

  // Number of clean and dirty Energy
  const vector<string> dirtyEnergy{"gas"};
  const vector<string> cleanEnergy{"solar", "wind"};
  const vector<string> energy = [](const vector<string> v1,
                                   const vector<string> v2) {
    vector<string> e{v1};
    for (const auto z : v2)
      e.push_back(z);
    return e;
  }(dirtyEnergy, cleanEnergy);

  map<string, double> s1 = {{string("wind"), 0.4}, {string("solar"), 0.6}};
  map<string, double> s2 = {{string("wind"), 0.6}, {string("solar"), 0.4}};
  const array<map<string, double>, NUM_SCEN> capFac1 = {s1, s2};

  // s1["wind"] = 0.65;
  // s1["solar"] = 0.45;
  // s2["wind"] = 0.75;
  // s2["solar"] = 0.35;

  // map<string, double> s1 = {{string("wind"), 0.65}, {string("solar"), 0.45}};
  // map<string, double> s2 = {{string("wind"), 0.75}, {string("solar"), 0.35}};

  const array<map<string, double>, NUM_SCEN> capFac2 = {s1, s2};

  // Constant costs across producers
  map<string, linQuad> prodCost;
  map<string, linQuad> invCost;
  map<string, double> emitCost;
  // zERO production costs for clean energy
  for (const auto z : cleanEnergy) {
    prodCost[z] = {0, 0};
    emitCost[z] = 0;
  }
  // prodCost["coal"] = {12, 0.1};
  // emitCost["coal"] = 2;
  emitCost["gas"] = 1;
  prodCost["gas"] = {6, 1.5};
  invCost["wind"] = {15, 1};
  invCost["solar"] = {25, 4};

  // Country 1, first Follower
  cci::FollPar<NUM_SCEN> c1F1(0, string("c1F1"));
  c1F1.name = string{"F1"};
  c1F1.productionCosts = prodCost;
  c1F1.investmentCosts = invCost;
  c1F1.emissionCosts = emitCost;
  c1F1.renewCapAdjust = capFac1;
  c1F1.capacities =
      makeMap<string, double>(energy, vector<double>{// 5,
                                                     200, 20, 20});

  // Country 1, second follower
  cci::FollPar<NUM_SCEN> c1F2(0, string("c1F2"));
  c1F2.name = string{"F2"};
  c1F2.productionCosts = prodCost;
  c1F2.investmentCosts = invCost;
  c1F2.emissionCosts = emitCost;
  c1F2.renewCapAdjust = capFac1;
  c1F2.capacities =
      makeMap<string, double>(energy, vector<double>{// 5,
                                                     200, 20, 20});

  // Country 2, first Follower
  cci::FollPar<NUM_SCEN> c2F1(0, string("c2F1"));
  c2F1.name = string{"F3"};
  c2F1.productionCosts = prodCost;
  c2F1.investmentCosts = invCost;
  c2F1.emissionCosts = emitCost;
  c2F1.renewCapAdjust = capFac2;
  c2F1.capacities =
      makeMap<string, double>(energy, vector<double>{// 4,
                                                     200, 20, 20});
  // c2F1.investmentCosts["wind"].first = 1;
  c2F1.investmentCosts["solar"].first = 1;

  // Country 2, second follower
  cci::FollPar<NUM_SCEN> c2F2(0, string("c2F2"));
  c2F2.name = string{"F4"};
  c2F2.productionCosts = prodCost;
  c2F2.investmentCosts = invCost;
  c2F2.emissionCosts = emitCost;
  c2F2.renewCapAdjust = capFac2;
  // c2F2.investmentCosts["wind"].first = 8;
  c2F2.investmentCosts["solar"].first = 3;
  c2F2.capacities =
      makeMap<string, double>(energy, vector<double>{// 2,
                                                     200, 20, 20});

  // Country 1 Leader Par
  cci::LeadPar Country1Par(0, // Min reqd consum
                           0, // Carbon credit
                           0, // emitVals
                           0, 0, makeMap(cleanEnergy, vector<double>{0, 0}),
                           makeMap(cleanEnergy, vector<double>{0, 0}),
                           300, // prodnVal
                           10   // Tax
  );

  cci::LeadPar Country2Par(0, // Min reqd consum
                           0, // Carbon credit
                           0, // emitVals
                           0, 0, makeMap(cleanEnergy, vector<double>{0, 0}),
                           makeMap(cleanEnergy, vector<double>{0, 0}),
                           100, // prodnVal
                           5    // Tax
  );

  linQuad DP1s1 = {800, 1};
  linQuad DP1s2 = {1200, 1};
  linQuad DP2s1 = {1200, 1};
  linQuad DP2s2 = {800, 1};
  array<linQuad, NUM_SCEN> DP1 = {DP1s1, DP1s2};
  array<linQuad, NUM_SCEN> DP2 = {DP2s1, DP2s2};

  cci::LeadAllPar<NUM_SCEN> c1(
      2, "c1", vector<cci::FollPar<NUM_SCEN>>{c1F1, c1F2}, Country1Par, DP1,
      array<double, NUM_SCEN>{1.0 / NUM_SCEN, 1.0 / NUM_SCEN}
      // Probability of the scenarios
  );
  cci::LeadAllPar<NUM_SCEN> c2(
      2, "c2", vector<cci::FollPar<NUM_SCEN>>{c2F1, c2F2}, Country2Par, DP2,
      array<double, NUM_SCEN>{1.0 / NUM_SCEN, 1.0 / NUM_SCEN}
      // Probability of the scenarios
  );

  cout << "\ncleanEnergy\n";
  print(cleanEnergy);
  cout << "\ndirtyEnergy\n";
  print(dirtyEnergy);
  cout << "\nenergy\n";
  print(energy);
  cout << "\nprodCost\n";
  print(prodCost);
  cout << "\ninvCost\n";
  print(invCost);
  cout << "\nemitCost\n";
  print(emitCost);

  cci::commonData supplyData(1, 0.5);

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
                                      boost::log::trivial::debug);
  GRBEnv env;
  try {
    cci::EPEC<NUM_DIRTY, NUM_CLEAN, NUM_SCEN> epec(&env, dirtyEnergy,
                                                   cleanEnergy, supplyData);

    epec.addCountry(c1).addCountry(c2);
    epec.finalize();
    epec.setAlgorithm(Game::EPECalgorithm::innerApproximation);
    // epec.setPureNE(true);
    epec.setNumThreads(2);
    std::string temp("");
    if (argc >= 2)
      temp = std::string(argv[1]);
    if (temp == "1")
      epec.setAlgorithm(Game::EPECalgorithm::fullEnumeration);
    cout << "Now finding Nash Equilibrium";
    std::string temp2 = R"(
   #
  # #    #        ####    ####   #####      #     #####  #    #  #    #
 #   #   #       #    #  #    #  #    #     #       #    #    #  ##  ##
#     #  #       #       #    #  #    #     #       #    ######  # ## #
#######  #       #  ###  #    #  #####      #       #    #    #  #    #
#     #  #       #    #  #    #  #   #      #       #    #    #  #    #
#     #  ######   ####    ####   #    #     #       #    #    #  #    #

			)";
    cout << temp2;

    // epec.setAggressiveness(1);
    // epec.setTimeLimit(1000);
    epec.setAddPolyMethod(Game::EPECAddPolyMethod::reverse_sequential);
    epec.findNashEq();
    epec.writeSolution(1, "dat/Sol");
    epec.writeLcpModel("dat/lcpmodel.lp");
    epec.writeLcpModel("dat/lcpmodel.sol");
    epec.appendSolution4XL("dat/solLog", 35008, true);
  } catch (const string &s) {
    cerr << s << "\nOops\n";
    throw;
  } catch (const GRBException &e) {
    cerr << "GRBException: " << e.getMessage() << '\n';
  }
  cout << R"( 
######
#     #   ####   #    #  ######
#     #  #    #  ##   #  #
#     #  #    #  # #  #  #####
#     #  #    #  #  # #  #
#     #  #    #  #   ##  #
######    ####   #    #  ######
)";

  return 0;
}
