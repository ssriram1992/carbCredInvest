#pragma once

struct commonData {
  double suppInt{100};
  double suppSlope{1};
  commonData(double sI, double sS) : suppInt{sI}, suppSlope{sS} {}
};

/// @brief Stores the parameters of the follower in a country model
template <unsigned int n_Scen> struct FollPar {
  /// For Eaxh prodType, linear and quadratic cost of production
  std::map<std::string, std::pair<double, double>> productionCosts;
  /// Production capacity of each prodType
  std::map<std::string, double> capacities;
  /// Renewable capacity adjustment (capacity factor) for each scenario
  std::array<std::map<std::string, double>, n_Scen> renewCapAdjust;
  ///	Linear and QUadratic Investment costs for renewables
  std::map<std::string, std::pair<double, double>> investmentCosts;
  /// Emission costs for unit quantity of the fuel. Emission costs feature only
  /// on the leader's problem
  std::map<std::string, double> emissionCosts;

  double carbonCreditInit = 0;
  double carbonLimitFrac = 0.5; // Number between 0 to 1
  std::string name = {};        ///< Optional Names for the Followers.

  FollPar(double cCI = 0, std::string name = "")
      : carbonCreditInit{cCI}, name{name} {}
};

/// @brief Stores the parameters of the leader in a country model
struct LeadPar {
  double consum_limit =
      0; ///< Government does not want the price to exceed this limit
  double carbCreditInit = 0;

  // Investment terms in the objective
  std::map<std::string, double>
      cleanInvVal{}; ///< Linear coefficient in the objective for each type of
                     ///< domestic investment
  std::map<std::string, double>
      cleanInvCrossVal{}; ///< Coefficient in the objective for  products of
                          ///< domestic investments and investment by others
  // Never considering square of domestic investment as it causes loss of
  // convexity

  // Emission terms in the objective
  double emissionVal; ///< Linear coefficient in the objective for domestic
                      ///< emission
  double emissionValQuad{0}; ///< For the product, should we also consider
                             ///< square of domestic emissions
  double emissionCrossVal;   ///< Coefficient in the objective for  products of
                             ///< domestic emission and emission by others
  double prodnVal{100};

  bool follGenNash{
      false}; ///< Whether the followers have the constraint that the total
              ///< carbon consumed by all the followers together is at most the
              ///< credits that their corresponding leader has.

  double taxCarbon =
      -1; ///< -1 implies tax is a variable. Any positive value fixes tax there

  LeadPar(double consum_limit = 0, double carbCred = 0, double emitVal = 0,
          double emitQuad = 0, double emitCross = 0,
          std::map<std::string, double> clInvVal = {},
          std::map<std::string, double> clInvCros = {}, double prodnVal = 100,
          double taxCarbon = 0)
      : consum_limit{consum_limit}, carbCreditInit{carbCred},
        cleanInvVal{clInvVal}, cleanInvCrossVal{clInvCros},
        emissionVal{emitVal}, emissionValQuad{emitQuad},
        emissionCrossVal{emitCross}, prodnVal{prodnVal}, taxCarbon{taxCarbon} {}
};

/// @brief Stores the parameters of a country model
template <unsigned int n_Scen> struct LeadAllPar {
  unsigned int n_followers; ///< Number of followers in the country
  std::string name;         ///< Country Name
  std::vector<cci::FollPar<n_Scen>> FollowerParam =
      {};                        ///< A vector to hold Follower Parameters
  cci::LeadPar LeaderParam = {}; ///< A struct to hold Leader Parameters
  std::array<std::pair<double, double>, n_Scen> DemandParam = {
      {}}; ///< To hold Demand Parameters - intercept and slope
  std::array<double, n_Scen> scenProb = {{}};
  LeadAllPar(unsigned int n_followers, std::string name,
             std::vector<cci::FollPar<n_Scen>> FP, cci::LeadPar LP,
             std::array<std::pair<double, double>, n_Scen> DP,
             std::array<double, n_Scen> sP = {})
      : n_followers{n_followers}, name{name}, FollowerParam{FP},
        LeaderParam{LP}, DemandParam{DP}, scenProb{sP} {}
};

/// @brief Stores a single Instance
template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
struct EPECInstance {
  std::vector<cci::LeadAllPar<n_Scen>> Countries = {}; ///< LeadAllPar vector

  EPECInstance(std::string filename) {
    this->load(filename);
  } ///< Constructor from instance file
  EPECInstance(std::vector<cci::LeadAllPar<n_Scen>> Countries_)
      : Countries{Countries_} {}
  ///< Constructor from instance objects

  void load(std::string filename);
  ///< Reads the EPECInstance from a file

  void save(std::string filename);
  ///< Writes the EPECInstance from a file
};

enum class LeaderVars {
  Followers,
  CarbImp,
  CarbTax,
  CarbFollLim,
  TotInv,
  TotEmission,
  DualVar,
  ConvHullDummy,
  End
};
