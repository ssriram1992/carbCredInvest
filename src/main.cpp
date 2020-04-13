#include "carbCredInv.h"
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

template <class t, class u>
std::map<t, u> makeMap(const std::vector<t> keys, const std::vector<u> values) {
  std::map<t, u> mm{};
  if (keys.size() != values.size())
    return mm;
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

int main() {
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
  constexpr unsigned int NUM_DIRTY{2};
  constexpr unsigned int NUM_SCEN{1};

  // Number of clean and dirty Energy
  const vector<string> dirtyEnergy{"coal", "gas"};
  const vector<string> cleanEnergy{"wind", "solar"};
  const vector<string> energy = [](const vector<string> v1,
                                   const vector<string> v2) {
    vector<string> e{v1};
    for (const auto z : v2)
      e.push_back(z);
    return e;
  }(dirtyEnergy, cleanEnergy);

  map<string, double> s1 = {{string("wind"), 0.7}, {string("solar"), 0.4}};
  map<string, double> s2 = {{string("wind"), 0.8}, {string("solar"), 0.3}};
  const array<map<string, double>, NUM_SCEN> capFac1 = {s1};

  s1["wind"] = 0.65;
  s1["solar"] = 0.45;
  s2["wind"] = 0.75;
  s2["solar"] = 0.35;

  // map<string, double> s1 = {{string("wind"), 0.65}, {string("solar"), 0.45}};
  // map<string, double> s2 = {{string("wind"), 0.75}, {string("solar"), 0.35}};

  const array<map<string, double>, NUM_SCEN> capFac2 = {s1};

  // Constant costs across producers
  map<string, linQuad> prodCost;
  map<string, linQuad> invCost;
  map<string, double> emitCost;
  // zERO production costs for clean energy
  for (const auto z : cleanEnergy) {
    prodCost[z] = {0.1, 0};
    emitCost[z] = 0;
  }
  prodCost["coal"] = {12, 0.1};
  emitCost["coal"] = 2;
  emitCost["gas"] = 1;
  prodCost["gas"] = {15, 0.5};
  invCost["wind"] = {1, 0};
  invCost["solar"] = {1, 0};

  // Country 1, first Follower
  cci::FollPar<NUM_SCEN> c1F1(15, string("c1F1"));
  c1F1.name = string{"c1F1"};
  c1F1.productionCosts = prodCost;
  c1F1.investmentCosts = invCost;
  c1F1.emissionCosts = emitCost;
  c1F1.renewCapAdjust = capFac1;
  c1F1.capacities =
      makeMap<string, double>(energy, vector<double>{50, 30, 10, 0});

  // Country 1, second follower
  cci::FollPar<NUM_SCEN> c1F2(10, string("c1F2"));
  c1F2.name = string{"c1F2"};
  c1F2.productionCosts = prodCost;
  c1F2.investmentCosts = invCost;
  c1F2.emissionCosts = emitCost;
  c1F2.renewCapAdjust = capFac1;
  c1F2.capacities =
      makeMap<string, double>(energy, vector<double>{5, 2, 4, 0});

  // Country 2, first Follower
  cci::FollPar<NUM_SCEN> c2F1(15, string("c2F1"));
  c2F1.name = string{"c2F1"};
  c2F1.productionCosts = prodCost;
  c2F1.investmentCosts = invCost;
  c2F1.emissionCosts = emitCost;
  c2F1.renewCapAdjust = capFac2;
  c2F1.capacities =
      makeMap<string, double>(energy, vector<double>{4, 3, 5, 0});
  c2F1.investmentCosts["wind"].first = 1;
  c2F1.investmentCosts["solar"].first = 1;

  // Country 2, second follower
  cci::FollPar<NUM_SCEN> c2F2(10, string("c2F2"));
  c2F2.name = string{"c2F2"};
  c2F2.productionCosts = prodCost;
  c2F2.investmentCosts = invCost;
  c2F2.emissionCosts = emitCost;
  c2F2.renewCapAdjust = capFac2;
  c2F2.investmentCosts["wind"].first = 8;
  c2F2.investmentCosts["solar"].first = 3;
  c2F2.capacities =
      makeMap<string, double>(energy, vector<double>{2, 5, 10, 0});

  // Country 1 Leader Par
  cci::LeadPar Country1Par(-1, // Import limit
                           -1, // Export limit
                           12, // Min reqd consum
                           20, // Carbon credit
                           10, // emitVals
                           1, 2, makeMap(cleanEnergy, vector<double>{100, 100}),
                           makeMap(cleanEnergy, vector<double>{0, 0}));

  cci::LeadPar Country2Par(-1, // Import limit
                           -1, // Export limit
                           13, // Min reqd consum
                           22, // Carbon credit
                           10, // emitVals
                           1, 2, makeMap(cleanEnergy, vector<double>{100, 100}),
                           makeMap(cleanEnergy, vector<double>{0, 0}));

  linQuad DP1s1 = {1000, 1};
  // linQuad DP1s2 = {1200, 0.95};
  linQuad DP2s1 = {800, 1.05};
  // linQuad DP2s2 = {900, 1};
  array<linQuad, NUM_SCEN> DP1 = {DP1s1};
  array<linQuad, NUM_SCEN> DP2 = {DP2s1};

  cci::LeadAllPar<NUM_SCEN> c1(1, "c1", vector<cci::FollPar<NUM_SCEN>>{c1F2},
                               Country1Par, DP1, array<double, NUM_SCEN>{1});
  cci::LeadAllPar<NUM_SCEN> c2(1, "c2", vector<cci::FollPar<NUM_SCEN>>{c2F1},
                               Country2Par, DP2, array<double, NUM_SCEN>{1});

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

  /*
#     #
##   ##   ####   #####   ######  #
# # # #  #    #  #    #  #       #
#  #  #  #    #  #    #  #####   #
#     #  #    #  #    #  #       #
#     #  #    #  #    #  #       #
#     #   ####   #####   ######  ######

*/

  GRBEnv env;
  cci::EPEC<NUM_DIRTY, NUM_CLEAN, NUM_SCEN> epec(&env, dirtyEnergy,
                                                 cleanEnergy);
  epec.addCountry(c1).addCountry(c2);
  epec.finalize();
  epec.setAlgorithm(Game::EPECalgorithm::innerApproximation);
	epec.setPureNE(true);
  // epec.setAlgorithm(Game::EPECalgorithm::fullEnumeration);
  cout << "Now finding Nash Equilibrium";
  std::string temp = R"(
   #
  # #    #        ####    ####   #####      #     #####  #    #  #    #
 #   #   #       #    #  #    #  #    #     #       #    #    #  ##  ##
#     #  #       #       #    #  #    #     #       #    ######  # ## #
#######  #       #  ###  #    #  #####      #       #    #    #  #    #
#     #  #       #    #  #    #  #   #      #       #    #    #  #    #
#     #  ######   ####    ####   #    #     #       #    #    #  #    #

			)";
  cout << temp;
  boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::debug);

  // epec.setAggressiveness(1);
	epec.setTimeLimit(1000);
  epec.findNashEq();
  epec.writeSolution(1, "dat/Sol");

  return 0;
}
