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

void solveProb(GRBEnv *env, const cci::probData pd, bool findPne) {
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
  constexpr unsigned int NUM_DIRTY{2};
  constexpr unsigned int NUM_SCEN{1};

  // Number of clean and dirty Energy
  const vector<string> dirtyEnergy{"coal", "gas"};
  const vector<string> cleanEnergy{"solar", "wind"};
  const vector<string> energy = [](const vector<string> v1,
                                   const vector<string> v2) {
    vector<string> e{v1};
    for (const auto z : v2)
      e.push_back(z);
    return e;
  }(dirtyEnergy, cleanEnergy);

  map<string, double> s1 = {{string("wind"), 0.348}, {string("solar"), 0.245}};
  // map<string, double> s2 = {{string("wind"), 0.6}, {string("solar"), 0.4}};
  const array<map<string, double>, NUM_SCEN> capFac1 = {s1};
  // const array<map<string, double>, NUM_SCEN> capFac2 = {s1};

  // Constant costs across producers
  map<string, linQuad> prodCost;
  map<string, linQuad> invCost;
  map<string, double> emitCost;

  prodCost["coal"] = {pd.prodCoalLin, pd.prodCoalQuad};
  prodCost["gas"] = {pd.prodGasLin, pd.prodGasQuad};
  prodCost["solar"] = {pd.prodSolLin, pd.prodSolQuad};
  prodCost["wind"] = {pd.prodWindLin, pd.prodWindQuad};

  emitCost["coal"] = pd.emitCoal;
  emitCost["gas"] = pd.emitGas;
  emitCost["solar"] = pd.emitSol;
  emitCost["wind"] = pd.emitWind;

  invCost["solar"] = {pd.invSolLin, pd.invSolQuad};
  invCost["wind"] = {pd.invWindLin, pd.invWindQuad};

  // Country 1, first Follower
  cci::FollPar<NUM_SCEN> c1F1(0, string("c1F1"));
  c1F1.name = string{"F1"};
  c1F1.productionCosts = prodCost;
  c1F1.investmentCosts = invCost;
  c1F1.emissionCosts = emitCost;
  c1F1.renewCapAdjust = capFac1;
  c1F1.capacities = makeMap<string, double>(
      energy, vector<double>{pd.c1f1CapCoal, pd.c1f1CapGas, pd.c1f1CapSol,
                             pd.c1f1CapWind});

  // Country 1, second follower
  cci::FollPar<NUM_SCEN> c1F2(0, string("c1F2"));
  c1F2.name = string{"F2"};
  c1F2.productionCosts = prodCost;
  c1F2.investmentCosts = invCost;
  c1F2.emissionCosts = emitCost;
  c1F2.renewCapAdjust = capFac1;
  c1F2.capacities = makeMap<string, double>(
      energy, vector<double>{pd.c1f2CapCoal, pd.c1f2CapGas, pd.c1f2CapSol,
                             pd.c1f2CapWind});

  // Country 2, first Follower
  cci::FollPar<NUM_SCEN> c2F1(0, string("c2F1"));
  c2F1.name = string{"F3"};
  c2F1.productionCosts = prodCost;
  c2F1.investmentCosts = invCost;
  c2F1.emissionCosts = emitCost;
  c2F1.renewCapAdjust = capFac1;
  c2F1.capacities = makeMap<string, double>(
      energy, vector<double>{pd.c2f1CapCoal, pd.c2f1CapGas, pd.c2f1CapSol,
                             pd.c2f1CapWind});
  // c2F1.investmentCosts["wind"].first = 1;
  // c2F1.investmentCosts["solar"].first = 1;

  // Country 2, second follower
  cci::FollPar<NUM_SCEN> c2F2(0, string("c2F2"));
  c2F2.name = string{"F4"};
  c2F2.productionCosts = prodCost;
  c2F2.investmentCosts = invCost;
  c2F2.emissionCosts = emitCost;
  c2F2.renewCapAdjust = capFac1;
  // c2F2.investmentCosts["wind"].first = 8;
  // c2F2.investmentCosts["solar"].first = 3;
  c2F2.capacities = makeMap<string, double>(
      energy, vector<double>{pd.c2f2CapCoal, pd.c2f2CapGas, pd.c2f2CapSol,
                             pd.c2f2CapWind});

  // Country 1 Leader Par
  cci::LeadPar Country1Par(
      0,                                          // Min reqd consum
      0,                                          // Carbon credit
      pd.c1emitCost,                              // emitVals
      0, 0,                                       // emitQuad, emitCross
      makeMap(cleanEnergy, vector<double>{0, 0}), // cleanInvVal
      makeMap(cleanEnergy, vector<double>{0, 0}), // cleanInvCrossVal
      pd.c1Pp,                                    // prodnVal
      pd.c1tax                                    // Tax
  );

  cci::LeadPar Country2Par(
      0,                                          // Min reqd consum
      0,                                          // Carbon credit
      pd.c2emitCost,                              // emitVals
      0, 0,                                       // emitQuad, emitCross
      makeMap(cleanEnergy, vector<double>{0, 0}), // cleanInvVal
      makeMap(cleanEnergy, vector<double>{0, 0}), // cleanInvCrossVal
      pd.c2Pp,                                    // prodnVal
      pd.c2tax                                    // Tax
  );

  Country1Par.follGenNash = pd.c1GenNash;
  Country2Par.follGenNash = pd.c2GenNash;

  Country1Par.expCleanMix = pd.c1CleanMix;
  Country2Par.expCleanMix = pd.c2CleanMix;

	Country1Par.rps = pd.c1rps;
	Country2Par.rps = pd.c2rps;

	BOOST_LOG_TRIVIAL(debug) << pd.c1rps<<" " <<pd.c2rps <<" RPS values << ";


  c1F1.carbonLimitFrac = pd.c1CreditSplit; // 0.5 * 2;
  c1F2.carbonLimitFrac = pd.c1CreditSplit; // 0.5 * 2;

  c2F1.carbonLimitFrac = pd.c1CreditSplit; // 0.5 * 2;
  c1F2.carbonLimitFrac = pd.c1CreditSplit; // 0.5 * 2;

  // linQuad DP1s1 = {3690, 0.6};
  linQuad DP1s1 = {pd.c1demInt, pd.c1demSl};
  // linQuad DP1s2 = {1200, 1};
  linQuad DP2s1 = {pd.c2demInt, pd.c2demSl};
  // linQuad DP2s2 = {800, 1};
  array<linQuad, NUM_SCEN> DP1 = {DP1s1};
  array<linQuad, NUM_SCEN> DP2 = {DP2s1};

  cci::LeadAllPar<NUM_SCEN> c1(
      2, "c1", vector<cci::FollPar<NUM_SCEN>>{c1F1, c1F2}, Country1Par, DP1,
      array<double, NUM_SCEN>{1.0 / NUM_SCEN} // Probability of the scenarios
  );
  cci::LeadAllPar<NUM_SCEN> c2(
      2, "c2", vector<cci::FollPar<NUM_SCEN>>{c2F1, c2F2}, Country2Par, DP2,
      array<double, NUM_SCEN>{1.0 / NUM_SCEN} // Probability of the scenarios
  );

  cci::commonData supplyData(pd.suIn, pd.suSl);

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
  try {
    // if (probId%10 == 0)
    std::string pref = pd.prefix;
    auto probId = pd.probId;
    cci::EPEC<NUM_DIRTY, NUM_CLEAN, NUM_SCEN> epec(env, dirtyEnergy,
                                                   cleanEnergy, supplyData);

    epec.addCountry(c1).addCountry(c2);
    epec.finalize();
    epec.setAlgorithm(Game::EPECalgorithm::innerApproximation);
    epec.setPureNE(findPne);
    epec.setNumThreads(2);
    // epec.setAggressiveness(1);
    epec.setTimeLimit(60);
    epec.setAddPolyMethod(Game::EPECAddPolyMethod::reverse_sequential);
    epec.findNashEq();

    { // Printing probability position
      for (int ii = 0; ii < 2; ii++) {
        cout << "Country: " << ii << '\n';
        unsigned int n_lead = epec.getNPoly_Lead(ii);
        cout << "Number of polyhedra: " << n_lead << '\n';
        for (unsigned int jj = 0; jj < n_lead; jj++)
          cout << "Leader " << ii << " poly " << jj
               << " probability posn: " << epec.getPosition_Probab(ii, jj)
               << " with a value " << epec.getVal_Probab(ii, jj) << '\n';
      }
      std::string prefix = pref + std::to_string(probId);
      epec.writeMNEposns("dat/" + prefix + "MNE", 1);
    }
    // epec.writeSolution(1, "dat/Sol" + std::to_string(probId));
    // epec.writeLcpModel("dat/lcpmodel_" + std::to_string(probId) + ".lp");
    // epec.writeLcpModel("dat/lcpmodel_" + std::to_string(probId) + ".sol");

    {
      std::string prefix = pref + std::to_string(probId);
      epec.appendSolution4XL("dat/" + prefix + "solLog", probId, false);
      // epec.writeSolution(1, "dat/" + prefix + "Sol" +
      // std::to_string(probId));
      epec.writeLcpModel("dat/" + prefix + "lcpmodel.lp");
      if (epec.getStatistics().status == Game::EPECsolveStatus::timeLimit) {
        std::ofstream el;
        el.open("dat/" + prefix + "timeLimit.log", ios::app);
        cout << "Time Limit reached \n";
        el.close();
        epec.writeSolution(1, "dat/" + prefix + "Sol");
        epec.writeLcpModel("dat/" + prefix + "lcpmodel.lp");
        // epec.writeLcpModel("dat/lcpmodel_" + std::to_string(probId) +
        // ".sol");
      }
      if (epec.getStatistics().status != Game::EPECsolveStatus::nashEqFound) {
        std::ofstream el;
        el.open("dat/" + prefix + "errorLog.log", ios::app);
        el.close();
      } else
        cout << "Successful solve\n";
    }

  } catch (...) {
    BOOST_LOG_TRIVIAL(fatal) << "ERRROR!!!!";
  }
}

int main(int argc, char *argv[]) {
  GRBEnv e;
  std::string filename;
  if (argc > 1)
    filename = std::string(argv[1]);
  else {
    std::cout << R"(
No input file given. Usage:
carbCredInv_PNE <filename.dat> <findPureNE=0> <reWriteInputFile=0>
Enter 1 to continue with default file. 0 to abort.
)";
    int choice;
    std::cin >> choice;
    if (choice == 1)
      filename = std::string("src/pd.dat");
    else
      return 0;
  }

  bool findPne{false};
  if (argc > 2)
    if (std::string(argv[2]) == std::string("1"))
      findPne = true;

  const cci::probData pd = cci::readProbData(filename);
  if (argc > 3)
    if (std::string(argv[3]) == std::string("1"))
      cci::writeProbData(filename, pd);

  solveProb(&e, pd, findPne);

  // solveProb(e, 0, 100, 100, 0, 0, 0, 1e-4, 0.5, -0.15, "BaseCase");
  /*
   solveProb(e, 3,     // problem Id
             100, 100, // coutnry production preference
             82, 82,   // Tax
             0, 1e-4,  // supply curve
             0.5,    // Proportion of carbon credits to be given to each
   follower 0.15, // In expectation % of production from clean sources. "Base");
                                                 */

  //   cout << R"(

  // #######
  // #     #  #    #  ######  #####
  // #     #  #    #  #       #    #
  // #     #  #    #  #####   #    #
  // #     #  #    #  #       #####
  // #     #   #  #   #       #   #
  // #######    ##    ######  #    #

  // )";
  return 0;
}
