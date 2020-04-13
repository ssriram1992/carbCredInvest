#pragma once

/**
 * @file src/carbCredInv.h Using EPECSolve to solve problems arising in
 * international energy markets with climate-conscious
 * countries.
 */

#include <armadillo>
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <epecsolve.h>
#include <gurobi_c++.h>
#include <iostream>
#include <memory>
#include <set>

namespace cci {

// const std::vector<std::string> dirtyEnergy {"coal", "gas"};
// const std::vector<std::string> cleanEnergy {"wind", "solar"};

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
  std::string name = {}; ///< Optional Names for the Followers.

  FollPar(double cCI = 0, std::string name = "")
      : carbonCreditInit{cCI}, name{name} {}
};

/// @brief Stores the parameters of the leader in a country model
struct LeadPar {
  double import_limit = -1; ///< Maximum net import in the country. If no limit,
                            ///< set the value as -1;
  double export_limit = -1; ///< Maximum net export in the country. If no limit,
                            ///< set the value as -1;
  double consum_limit =
      0; ///< Government does not want the price to exceed this limit
  double carbCreditInit = 100;

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

  LeadPar(double imp_lim = -1, double exp_lim = -1, double consum_limit = 0,
          double carbCred = 0, double emitVal = 0, double emitQuad = 0,
          double emitCross = 0, std::map<std::string, double> clInvVal = {},
          std::map<std::string, double> clInvCros = {})
      : import_limit{imp_lim}, export_limit{exp_lim},
        consum_limit{consum_limit}, carbCreditInit{carbCred},
        cleanInvVal{clInvVal}, cleanInvCrossVal{clInvCros},
        emissionVal{emitVal}, emissionValQuad{emitQuad}, emissionCrossVal{
                                                             emitCross} {}
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
  CarbExp,
  CarbImp,
  CarbPrice,
  TotInv,
  TotEmission,
  DualVar,
  ConvHullDummy,
  End
};

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
std::ostream &operator<<(std::ostream &ost, const FollPar<n_Scen> P);
std::ostream &operator<<(std::ostream &ost, const std::pair<double, double> P);

std::ostream &operator<<(std::ostream &ost, const LeadPar P);

template <unsigned int n_Scen>
std::ostream &operator<<(std::ostream &ost, const LeadAllPar<n_Scen> P);

std::ostream &operator<<(std::ostream &ost, const LeaderVars l);

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
std::ostream &operator<<(std::ostream &ost,
                         EPECInstance<n_Dirty, n_Clean, n_Scen> I);

using LeadLocs = std::map<LeaderVars, unsigned int>;

void increaseVal(LeadLocs &L, const LeaderVars start, const unsigned int val,
                 const bool startnext = true);
void decreaseVal(LeadLocs &L, const LeaderVars start, const unsigned int val,
                 const bool startnext = true);

void init(LeadLocs &L);

LeaderVars operator+(cci::LeaderVars a, int b);
/*

  ####   #         ##     ####    ####
 #    #  #        #  #   #       #
 #       #       #    #   ####    ####
 #       #       ######       #       #
 #    #  #       #    #  #    #  #    #
  ####   ######  #    #   ####    ####


 #####   #####    ####    #####   ####    #####   #   #  #####   ######
 #    #  #    #  #    #     #    #    #     #      # #   #    #  #
 #    #  #    #  #    #     #    #    #     #       #    #    #  #####
 #####   #####   #    #     #    #    #     #       #    #####   #
 #       #   #   #    #     #    #    #     #       #    #       #
 #       #    #   ####      #     ####      #       #    #       ######

*/

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
class EPEC : public Game::EPEC {
  // Mandatory virtuals
private:
  void make_obj_leader(const unsigned int i, Game::QP_objective &QP_obj);
  virtual void updateLocs() override {
    for (unsigned int i = 0; i < this->getNcountries(); ++i) {
      LeadLocs &Loc = this->Locations.at(i);
      cci::decreaseVal(Loc, cci::LeaderVars::ConvHullDummy,
                       Loc[cci::LeaderVars::ConvHullDummy + 1] -
                           Loc[cci::LeaderVars::ConvHullDummy]);
      cci::increaseVal(Loc, cci::LeaderVars::ConvHullDummy,
                       this->convexHullVariables.at(i));
    }
  };
  virtual void postfinalize() override{};
  virtual void initializeSoln(arma::vec &init_x) const;
  // override;

public:
  // Rest
  static constexpr double social_cost_of_carbon{1000};

private:
  const std::vector<std::string> dirtyEnergy;
  const std::vector<std::string> cleanEnergy;
  /// cci::energy should contain everything in dirtyEnergy and cleanEnergy
  std::vector<std::string> energy;

  mutable std::map<std::string, double> dataMap{};

  static constexpr unsigned int LL_MC_count =
      0; // No market clearing for follower

  static constexpr unsigned int FollInv{0};
  static constexpr unsigned int FollProdDirty{n_Clean};
  static constexpr unsigned int FollProdClean{FollProdDirty + n_Dirty * n_Scen};
  static constexpr unsigned int FollCarbBuy{n_Clean +
                                            (n_Clean + n_Dirty) * n_Scen};
  static constexpr unsigned int FollCarbSel{FollCarbBuy + 1};
  static constexpr unsigned int FollEnd{FollCarbSel + 1};

  static constexpr unsigned int FollVarCount{FollEnd};

  std::vector<LeadAllPar<n_Scen>> AllLeadPars =
      {}; ///< The parameters of each leader in the EPEC game
  std::vector<std::shared_ptr<Game::QP_Param>> MC_QP =
      {}; ///< The QP corresponding to the market clearing condition of each
          ///< player
  std::vector<LeadLocs> Locations =
      {}; ///< Location of variables for each country

  std::map<std::string, unsigned int> name2nos = {};
  std::vector<arma::sp_mat>
      LeadConses{};                   ///< Stores each leader's constraint LHS
  std::vector<arma::vec> LeadRHSes{}; ///< Stores each leader's constraint RHS

  // Super low level
  /// Checks that the parameter given to add a country is valid. Does not have
  /// obvious errors
  bool ParamValid(const LeadAllPar<n_Scen> &Param) const;

  /// Makes the lower level quadratic program object for each follower.
  void make_LL_QP(const LeadAllPar<n_Scen> &Params, const unsigned int follower,
                  Game::QP_Param *Foll, const LeadLocs &Loc) const noexcept;

  /// Makes the leader constraint matrix and RHS
  void make_LL_LeadCons(arma::sp_mat &LeadCons, arma::vec &LeadRHS,
                        const LeadAllPar<n_Scen> &Param,
                        const cci::LeadLocs &Loc = {},
                        const unsigned int import_lim_cons = 1,
                        const unsigned int export_lim_cons = 1) const noexcept;

  void make_MC_cons(arma::sp_mat &MCLHS, arma::vec &MCRHS) const override;

  void WriteCountry(const unsigned int i, const std::string filename,
                    const arma::vec x, const bool append = true) const;

  void WriteCountryMCprice(const std::string filename, const arma::vec x,
                           const bool append = false) const;

  void WriteFollower(const unsigned int i, const unsigned int j,
                     const std::string filename, const arma::vec x) const;

public: // Attributes
  // double timeLimit = {-1}; ///< Controls the timeLimit (s) for findNashEq

  EPEC() = delete;

  EPEC(GRBEnv *env, const std::vector<std::string> dirtyEnergy,
       const std::vector<std::string> cleanEnergy)
      : Game::EPEC(env), dirtyEnergy{dirtyEnergy}, cleanEnergy{cleanEnergy} {
    if (dirtyEnergy.size() != n_Dirty || cleanEnergy.size() != n_Clean) {
      throw std::string(
          "Error in EPEC(): Illegal size of dirtyEnergy/cleanEnergy");
    }
    this->energy.clear();
    this->energy.insert(energy.end(), dirtyEnergy.cbegin(), dirtyEnergy.cend());
    this->energy.insert(energy.end(), cleanEnergy.cbegin(), cleanEnergy.cend());
  }

  ///@brief %cci a Standard Nash-Cournot game within a country
  EPEC &addCountry(
      /// The Parameter structure for the leader
      LeadAllPar<n_Scen> Params);

  unsigned int getPosition(const unsigned int countryCount,
                           const LeaderVars var = LeaderVars::Followers) const;

  unsigned int getPosition(const std::string countryCount,
                           const LeaderVars var = LeaderVars::Followers) const;

  EPEC &unlock();

  std::unique_ptr<GRBModel> Respond(const std::string name,
                                    const arma::vec &x) const;
  // Data access methods
  Game::NashGame *get_LowerLevelNash(const unsigned int i) const;

  Game::LCP *playCountry(std::vector<Game::LCP *> countries);

  // Writing model files
  void write(const std::string filename, const unsigned int i,
             bool append = true) const;

  void write(const std::string filename, bool append = true) const;

  void readSolutionJSON(const std::string filename);

  void writeSolutionJSON(std::string filename, const arma::vec x,
                         const arma::vec z) const;

  void writeSolution(const int writeLevel, std::string filename) const;

  ///@brief Get the current EPECInstance loaded
  const EPECInstance<n_Dirty, n_Clean, n_Scen> getInstance() const {
    return EPECInstance<n_Dirty, n_Clean, n_Scen>(this->AllLeadPars);
  }
};

} // namespace cci

// Gurobi functions
std::string to_string(const GRBVar &var);

std::string to_string(const GRBConstr &cons, const GRBModel &model);

// ostream functions
namespace cci {
enum class prn { label, val };

std::ostream &operator<<(std::ostream &ost, cci::prn l);
} // namespace cci

template <unsigned int n_Scen>
cci::FollPar<n_Scen> operator+(const cci::FollPar<n_Scen> &F1,
                               const cci::FollPar<n_Scen> &F2);

/*


######
#     #  ######  ######  #    #
#     #  #       #       ##   #
#     #  #####   #####   # #  #
#     #  #       #       #  # #   ###
#     #  #       #       #   ##   ###
######   ######  #       #    #   ###


*/

/* Examples */
template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::make_obj_leader(
    const unsigned int
        i, ///< The location of the country whose objective is to be made
    Game::QP_objective
        &QP_obj ///< The object where the objective parameters are to be stored.
    )
/**
 * Makes the objective function of each country.
 */
{
  const unsigned int nEPECvars = this->getnVarinEPEC();
  const unsigned int nThisCountryvars =
      this->Locations.at(i).at(cci::LeaderVars::End);
  const LeadAllPar<n_Scen> &Params = this->AllLeadPars.at(i);
  const LeadLocs &Loc = this->Locations.at(i);

	BOOST_LOG_TRIVIAL(debug) << "In cci::EPEC::make_obj_leader";
  QP_obj.Q.zeros(nThisCountryvars, nThisCountryvars);
  QP_obj.c.zeros(nThisCountryvars);
  QP_obj.C.zeros(nThisCountryvars, nEPECvars - nThisCountryvars);

  // Total investment locally
  for (unsigned int ii = 0; ii < n_Clean; ++ii) {
    QP_obj.c.at(Loc.at(LeaderVars::TotInv) + ii) =
        -1 * Params.LeaderParam.cleanInvVal.at(cleanEnergy.at(ii));
    // Negative sign because investment has to be maximized
  }
  // Total investment cross terms
  for (unsigned int ii = 0; ii < n_Clean; ++ii) {
    for (unsigned int cc = 0; cc < this->getNcountries(); ++cc) {
      QP_obj.C(Loc.at(LeaderVars::TotInv) + ii,
               this->getPosition(this->getNcountries() - 1,
                                 cci::LeaderVars::TotInv + ii) -
                   nThisCountryvars) =
          Params.LeaderParam.cleanInvCrossVal.at(cleanEnergy.at(ii));
    }
  }

  // Emissions local
  QP_obj.Q(Loc.at(LeaderVars::TotEmission), Loc.at(LeaderVars::TotEmission)) =
      Params.LeaderParam.emissionValQuad;
  QP_obj.c(Loc.at(LeaderVars::TotEmission)) = Params.LeaderParam.emissionVal;
  // Emission cross terms
  for (unsigned int cc = 0; cc < this->getNcountries(); ++cc) {
    QP_obj.C(Loc.at(LeaderVars::TotEmission),
             this->getPosition(this->getNcountries() - 1,
                               cci::LeaderVars::TotEmission) -
                 nThisCountryvars) = Params.LeaderParam.emissionCrossVal;
  }

  // Carbon credit trade term
  QP_obj.C.at(
      Loc.at(LeaderVars::CarbExp),
      this->getPosition(this->getNcountries() - 1, cci::LeaderVars::End) -
          nThisCountryvars) = 1;
  QP_obj.C.at(
      Loc.at(LeaderVars::CarbImp),
      this->getPosition(this->getNcountries() - 1, cci::LeaderVars::End) -
          nThisCountryvars) = -1;
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
cci::EPEC<n_Dirty, n_Clean, n_Scen> &
cci::EPEC<n_Dirty, n_Clean, n_Scen>::addCountry(
    cci::LeadAllPar<n_Scen> Params) {

  LeadLocs Loc;
  cci::init(Loc);

  // Allocate so much space for each of these types of variables
  // For follower's carbon purchase from government.
  cci::increaseVal(Loc, LeaderVars::Followers,
                   Params.n_followers * this->FollVarCount);
  cci::increaseVal(Loc, LeaderVars::CarbExp, 1);
  cci::increaseVal(Loc, LeaderVars::CarbImp, 1);
  cci::increaseVal(Loc, LeaderVars::CarbPrice, 1);
  cci::increaseVal(Loc, LeaderVars::TotInv, n_Clean);
  cci::increaseVal(Loc, LeaderVars::TotEmission, 1);

  const unsigned int LeadVars =
      Loc.at(cci::LeaderVars::End) - this->FollVarCount * Params.n_followers;
  // Should be 4+n_Clean - one for CarbExp, one for CarbImp, one for CarbPrice,
  // and n_Clean each for TotInv and TotEmission

  // Leader Constraints
  const short int import_lim_cons = [&Params]() {
    if (Params.LeaderParam.import_limit >= 0)
      return 1u;
    else
      return 0u;
  }();
  const short int export_lim_cons = [&Params]() {
    if (Params.LeaderParam.export_limit >= 0)
      return 1u;
    else
      return 0u;
  }();

  arma::sp_mat LeadCons(import_lim_cons +     // Import limit constraint
                            export_lim_cons + // Export limit constraint
                            n_Scen +          // Min consumption
                            2*n_Clean +       // Investment summing
                            2 +               // Emission summing
                            1,                // Carbon credits >=0
                        Loc[cci::LeaderVars::End] - this->LL_MC_count);
  arma::vec LeadRHS(LeadCons.n_rows, arma::fill::zeros);

  std::vector<std::shared_ptr<Game::QP_Param>> FollowersVec{};
  // Create the QP_Param* for each follower
  try {
    for (unsigned int follower = 0; follower < Params.n_followers; follower++) {
      auto Foll = std::make_shared<Game::QP_Param>(this->env);
      this->make_LL_QP(Params, follower, Foll.get(), Loc);
      FollowersVec.push_back(Foll);
    }
    // Make Leader Constraints
    this->make_LL_LeadCons(LeadCons, LeadRHS, Params, Loc, import_lim_cons,
                           export_lim_cons);
  } catch (const char *e) {
    cerr << e << '\n';
    throw;
  } catch (std::string e) {
    cerr << "String in cci::EPEC<n_Dirty,n_Clean,n_Scen>::addCountry : " << e
         << '\n';
    throw;
  } catch (GRBException &e) {
    cerr << "GRBException in cci::EPEC<n_Dirty,n_Clean,n_Scen>::addCountry : "
         << e.getErrorCode() << ": " << e.getMessage() << '\n';
    throw;
  } catch (std::exception &e) {
    cerr << "Exception in cci::EPEC<n_Dirty,n_Clean,n_Scen>::addCountry : "
         << e.what() << '\n';
    throw;
  }

  // Lower level Market clearing constraints - None
  arma::sp_mat MC(this->LL_MC_count,
                  LeadVars + FollVarCount * Params.n_followers);
  arma::vec MCRHS(this->LL_MC_count, arma::fill::zeros);

  // Convert the country QP to a NashGame
  auto N = std::make_shared<Game::NashGame>(this->env, FollowersVec, MC, MCRHS,
                                            LeadVars - this->LL_MC_count,
                                            LeadCons, LeadRHS);
  this->name2nos[Params.name] = this->countries_LL.size();
  this->countries_LL.push_back(N);
  cci::increaseVal(
      Loc, cci::LeaderVars::DualVar,
      N->getNduals()); // N->getNduals() will sum the number of constraints in
                       // each lower level QP and provide the sum. Indeed, this
                       // is the number of dual variables for the lower level.
  this->Locations.push_back(Loc);

  this->EPEC::LocEnds.push_back(&this->Locations.back().at(LeaderVars::End));
  this->EPEC::convexHullVariables.push_back(0);

  this->LeadConses.push_back(N->RewriteLeadCons()); // Not mandatory!
  this->AllLeadPars.push_back(Params);
  this->Game::EPEC::n_MCVar = 1;
  return *this;
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::make_LL_QP(
    const cci::LeadAllPar<n_Scen> &Params, ///< The Parameters object
    const unsigned int follower, ///< Which follower's QP has to be made?
    Game::QP_Param
        *Foll, ///< Non-owning pointer to the Follower QP_Param object
    const cci::LeadLocs
        &Loc ///< LeadLocs object for accessing different leader locations.
) const noexcept
/**
 * @brief Makes Lower Level Quadratic Programs
 * @details Sets the constraints and objective for the lower level problem
 * (i.e., the follower)
 */
{
  const unsigned int n_Foll = Params.n_followers;
  const unsigned int LeadVars =
      Loc.at(cci::LeaderVars::End) - this->FollVarCount * n_Foll;
  const unsigned int otherVars = LeadVars + (n_Foll - 1) * FollVarCount;

  arma::sp_mat Q(FollVarCount, FollVarCount);
  arma::sp_mat C(FollVarCount, otherVars);

  // n_Dirty*n_Scen for dirty energy
  // infrastructural limits
  // n_Clean*n_Scen for clean energy
  // infrastructural limits
  // n_Scen for carbon credit limits
  constexpr unsigned int n_Const = (n_Dirty + n_Clean) * n_Scen + n_Scen;

  const auto &Follparam = Params.FollowerParam.at(follower);

  arma::sp_mat A(n_Const, otherVars);
  arma::sp_mat B(n_Const, FollVarCount);
  arma::vec c(FollVarCount), b(n_Const);
  c.fill(0);
  b.fill(0);
  A.zeros();
  B.zeros();
  C.zeros();
  b.zeros();
  Q.zeros();
  c.zeros();

  // OBJECTIVE
  // Investment cost
  for (unsigned int ii = 0; ii < n_Clean; ii++) {
    Q(FollInv + ii, FollInv + ii) =
        Follparam.investmentCosts.at(cleanEnergy.at(ii)).second;
    c(FollInv + ii) = Follparam.investmentCosts.at(cleanEnergy.at(ii)).first;
  }
  // Dirty energy production cost and demand term
  for (unsigned int ii = 0; ii < n_Dirty; ii++) {
    for (unsigned int scen = 0; scen < n_Scen; ++scen) {
      // position of cost
      const auto pos{FollProdDirty + n_Dirty * scen + ii};
      // quadratic part
      Q(pos, pos) = Params.scenProb.at(scen) *
                    (Follparam.productionCosts.at(dirtyEnergy.at(ii)).second +
                     2 * Params.DemandParam.at(scen).second);
      // linear part
      c(pos) = Params.scenProb.at(scen) *
               (Follparam.productionCosts.at(dirtyEnergy.at(ii)).first -
                Params.DemandParam.at(scen).first);
    }
    // Energy Demand - off diagonal terms
    for (unsigned int jj = 0; jj < ii; jj++) {
      for (unsigned int scen = 0; scen < n_Scen; ++scen) {
        const auto pos1{FollProdDirty + n_Dirty * scen + ii};
        const auto pos2{FollProdDirty + n_Dirty * scen + jj};
        Q(pos1, pos2) =
            Params.scenProb.at(scen) * (2 * Params.DemandParam.at(scen).second);
        Q(pos2, pos1) = Q(pos1, pos2);
      }
    }
  }
  // Clean energy Production cost and demand
  for (unsigned int ii = 0; ii < n_Clean; ii++) {
    for (unsigned int scen = 0; scen < n_Scen; ++scen) {
      // position of cost
      const auto pos{FollProdClean + n_Clean * scen + ii};
      // quadratic part
      Q(pos, pos) = Params.scenProb.at(scen) *
                    (Follparam.productionCosts.at(cleanEnergy.at(ii)).second +
                     2 * Params.DemandParam.at(scen).second);
      // linear part
      c(pos) = Params.scenProb.at(scen) *
               (Follparam.productionCosts.at(cleanEnergy.at(ii)).first -
                Params.DemandParam.at(scen).first);
    }
    // Energy Demand - off diagonal terms
    for (unsigned int jj = 0; jj < ii; jj++) {
      for (unsigned int scen = 0; scen < n_Scen; ++scen) {
        const auto pos1{FollProdClean + n_Clean * scen + ii};
        const auto pos2{FollProdClean + n_Clean * scen + jj};
        Q(pos1, pos2) =
            Params.scenProb.at(scen) * (2 * Params.DemandParam.at(scen).second);
        Q(pos2, pos1) = Q(pos1, pos2);
      }
    }
  }
  // Energy Demand - off diagonal cross terms
  for (unsigned int ii = 0; ii < n_Clean; ii++) {
    for (unsigned int jj = 0; jj < ii; jj++) {
      for (unsigned int scen = 0; scen < n_Scen; ++scen) {
        const auto pos1{FollProdDirty + n_Dirty * scen + ii};
        const auto pos2{FollProdClean + n_Clean * scen + jj};
        Q(pos1, pos2) =
            Params.scenProb.at(scen) * (2 * Params.DemandParam.at(scen).second);
        Q(pos2, pos1) = Q(pos1, pos2);
      }
    }
  }

  // C - the crossterms
  // Carbon price by government times carbon purchased
  C(FollCarbBuy, Loc.at(LeaderVars::CarbPrice) - FollVarCount) = 1;
  C(FollCarbSel, Loc.at(LeaderVars::CarbPrice) - FollVarCount) = -1;
  // C - the crossterms
  for (unsigned int scen = 0; scen < n_Scen; scen++) {
    for (unsigned int ii = 0; ii < n_Clean + n_Dirty; ii++) {
      for (unsigned int jj = 0; jj < n_Clean + n_Dirty; jj++) {
        // Lambda that gives the position of the
        // current production variable
        const auto pos = [scen](const unsigned int i_nos) {
          if (i_nos < n_Dirty)
            return FollProdDirty + n_Dirty * scen + i_nos;
          else
            return FollProdClean + n_Clean * scen + i_nos;
        };
        const auto pos1 = pos(ii);
        const auto pos2 = pos(jj);
        for (unsigned int kk = 0; kk < n_Foll - 1; ++kk) {
          C(pos1, kk * FollVarCount + pos2) =
              Params.scenProb.at(scen) *
              (2 * Params.DemandParam.at(scen).second);
        }
      }
    }
  }

  // OBJECTIVE Described

  // CONSTRAINTS
  const auto &infCap = Follparam.capacities;
  const auto &renCapAdj = Follparam.renewCapAdjust;
  unsigned int constrCount{0};
  for (unsigned int scen = 0; scen < n_Scen; ++scen) {
    // Infrastructural limit
    for (unsigned int ii = 0; ii < n_Dirty; ++ii) {
      B(constrCount, FollProdDirty + scen * n_Dirty + ii) = 1;
      b(constrCount) = infCap.at(dirtyEnergy.at(ii));
      // Next constraint
      constrCount++;
    }
    for (unsigned int ii = 0; ii < n_Clean; ++ii) {
      const auto energ = cleanEnergy.at(ii);
      B(constrCount, FollProdClean + scen * n_Clean + ii) = 1;
      b(constrCount) = infCap.at(energ) * renCapAdj.at(scen).at(energ);
      // Next constraint
      constrCount++;
    }
    // Carbon credits limit Constraint
    B(constrCount, FollCarbBuy) = -1;
    B(constrCount, FollCarbSel) = +1;
    for (unsigned int ii = 0; ii < n_Dirty; ++ii)
      B(constrCount, FollProdDirty + scen * n_Dirty + ii) =
          Follparam.emissionCosts.at(dirtyEnergy.at(ii));
    for (unsigned int ii = 0; ii < n_Clean; ++ii)
      B(constrCount, FollProdClean + scen * n_Clean + ii) =
          Follparam.emissionCosts.at(cleanEnergy.at(ii));
    b(constrCount) = Follparam.carbonCreditInit;
  }
  // CONSTRAINTS Described

  Foll->set(Q, C, A, B, c, b);
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::make_LL_LeadCons(
    arma::sp_mat
        &LeadCons,      ///< The LHS matrix of leader constraints (for output)
    arma::vec &LeadRHS, ///< RHS vector for leader constraints (for output)
    const LeadAllPar<n_Scen> &Params,   ///< All country specific parameters
    const cci::LeadLocs &Loc,           ///< Location of variables
    const unsigned int import_lim_cons, ///< Does a constraint on import limit
                                        ///< exist or no limit?
    const unsigned int export_lim_cons  ///< Does a constraint on export limit
                                        ///< exist or no limit?
) const noexcept
/**
 * Makes the leader level constraints for a country.
 */
{
  // LeadCons and LeadRHS already have the right size
  LeadCons.zeros();
  LeadRHS.zeros();
  unsigned int constrCount{0};
  // Import limit constraints
  if (import_lim_cons) {
    LeadCons(constrCount, Loc.at(LeaderVars::CarbExp)) = -1;
    LeadCons(constrCount, Loc.at(LeaderVars::CarbImp)) = 1;
    LeadRHS(constrCount) = Params.LeaderParam.import_limit;
    constrCount++;
  }
  // export limit constraints
  if (export_lim_cons) {
    LeadCons(constrCount, Loc.at(LeaderVars::CarbExp)) = 1;
    LeadCons(constrCount, Loc.at(LeaderVars::CarbImp)) = -1;
    LeadRHS(constrCount) = Params.LeaderParam.export_limit;
    constrCount++;
  }
  // Scenario constraints
  for (unsigned int scen = 0; scen < n_Scen; scen++) {
    // Minimum Consumption  constraint
    // Each producer production
    for (unsigned int ff = 0; ff < Params.n_followers; ff++) {
      // Each producer's dirty production
      for (unsigned int ii = 0; ii < n_Dirty; ++ii)
        LeadCons(constrCount,
                 FollVarCount * ff + FollProdDirty + scen * n_Dirty + ii) = -1;
      // Each producer's clean production
      for (unsigned int ii = 0; ii < n_Clean; ++ii)
        LeadCons(constrCount,
                 FollVarCount * ff + FollProdClean + scen * n_Clean + ii) = -1;
    }
    LeadRHS(constrCount) = -std::max(0.0, Params.LeaderParam.consum_limit);
    constrCount++;
  }
  // Investment summing

  for (unsigned int ii = 0; ii < n_Clean; ++ii) {
    for (unsigned int ff = 0; ff < Params.n_followers; ++ff) {
      LeadCons(constrCount, FollVarCount * ff + FollInv + ii) = 1;
      LeadCons(constrCount + 1, FollVarCount * ff + FollInv + ii) = -1;
    }
    LeadCons(constrCount, Loc.at(LeaderVars::TotInv) + ii) = -1;
    LeadCons(constrCount + 1, Loc.at(LeaderVars::TotInv) + ii) = 1;
  }
  constrCount += 2;
  // Expected Emission summing
  for (unsigned int scen = 0; scen < n_Scen; ++scen) {
    const auto probab =
        Params.scenProb.at(scen); // Probability of the current scenario

    for (unsigned int ff = 0; ff < Params.n_followers; ++ff) {
      // Dirty Producers
      for (unsigned int ii = 0; ii < n_Dirty; ++ii) {
        const auto &pt = dirtyEnergy.at(ii);
        const auto emitCost = Params.FollowerParam.at(ff).emissionCosts.at(pt);
        LeadCons(constrCount, FollVarCount * ff + FollProdDirty +
                                  n_Dirty * scen + ii) = probab * emitCost;
        LeadCons(constrCount + 1,
                 FollVarCount * ff + FollProdDirty + n_Dirty * scen + ii) =
            -1 * probab * emitCost;
      }
      // Clean Producers
      for (unsigned int ii = 0; ii < n_Clean; ++ii) {
        const auto &pt = cleanEnergy.at(ii); // player type
        const auto emitCost = Params.FollowerParam.at(ff).emissionCosts.at(pt);
        LeadCons(constrCount + 1, FollVarCount * ff + FollProdClean +
                                      n_Clean * scen + ii) = probab * emitCost;
        LeadCons(constrCount + 1,
                 FollVarCount * ff + FollProdClean + n_Clean * scen + ii) =
            -1 * probab * emitCost;
      }
      LeadCons(constrCount, Loc.at(LeaderVars::TotEmission)) = -1;
      LeadCons(constrCount + 1, Loc.at(LeaderVars::TotEmission)) = 1;
    }
  }
  constrCount += 2;
  // Carbon credit positivity constraint
  LeadCons(constrCount, Loc.at(LeaderVars::CarbExp)) = 1;
  LeadCons(constrCount, Loc.at(LeaderVars::CarbImp)) = -1;
  for (unsigned int ff = 0; ff < Params.n_followers; ff++) {
    LeadCons(constrCount, FollVarCount * ff + FollCarbBuy) = 1;
    LeadCons(constrCount, FollVarCount * ff + FollCarbSel) = -1;
  }
  LeadRHS(constrCount) = Params.LeaderParam.carbCreditInit;
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::make_MC_cons(arma::sp_mat &MCLHS,
                                                       arma::vec &MCRHS) const
/** @brief Returns leader's Market clearing constraints in matrix form
 * @details
 */
{
  if (!this->finalized)
    throw std::string(
        "Error in cci::EPEC<n_Dirty,n_Clean,n_Scen>::make_MC_cons: This "
        "function can be "
        "run only AFTER calling finalize()");
  // Output matrices
  MCRHS.zeros(1);
  MCLHS.zeros(1, this->getnVarinEPEC());
  if (this->getNcountries() > 1) {
    for (unsigned int i = 0; i < this->getNcountries(); ++i) {
      // The MC constraint for each leader country's carbon trade market
      MCLHS(0, this->getPosition(i, LeaderVars::CarbExp)) = 1;
      MCLHS(0, this->getPosition(i, LeaderVars::CarbImp)) = -1;
    }
  }
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
unsigned int cci::EPEC<n_Dirty, n_Clean, n_Scen>::getPosition(
    const unsigned int countryCount, const cci::LeaderVars var) const
/**
 * @brief Gets position of a variable in a country.
 */
{
  if (countryCount >= this->getNcountries())
    throw std::string(
        "Error in cci::EPEC<n_Dirty,n_Clean,n_Scen>::getPosition: Bad Country "
        "Count");
  return this->LeaderLocations.at(countryCount) +
         this->Locations.at(countryCount).at(var);
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
unsigned int cci::EPEC<n_Dirty, n_Clean, n_Scen>::getPosition(
    const std::string countryName, const cci::LeaderVars var) const
/**
 * @brief Gets position of a variable in a country given the country name and
 * the variable.
 */
{
  return this->getPosition(name2nos.at(countryName), var);
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
Game::NashGame *cci::EPEC<n_Dirty, n_Clean, n_Scen>::get_LowerLevelNash(
    const unsigned int i) const
/**
 * @brief Returns a non-owning pointer to the @p i -th country's lower level
 * NashGame
 */
{
  return this->countries_LL.at(i).get();
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
cci::EPEC<n_Dirty, n_Clean, n_Scen> &
cci::EPEC<n_Dirty, n_Clean, n_Scen>::unlock()
/**
 * @brief Unlocks an EPEC model
 * @details A finalized model cannot be edited unless it is unlocked first.
 * @internal EPEC::finalize() performs "finalizing" acts on an object.
 * @warning Exclusively for debugging purposes for developers. Don't call this
 * function, unless you know what you are doing.
 */
{
  this->finalized = false;
  return *this;
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
std::unique_ptr<GRBModel>
cci::EPEC<n_Dirty, n_Clean, n_Scen>::Respond(const std::string name,
                                             const arma::vec &x) const {
  return this->Game::EPEC::Respond(this->name2nos.at(name), x);
}
//

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::write(const std::string filename,
                                                const unsigned int i,
                                                bool append) const {
  std::ofstream file;
  file.open(filename, append ? ios::app : ios::out);
  const LeadAllPar<n_Scen> &Params = this->AllLeadPars.at(i);
  file << "**************************************************\n";
  file << "COUNTRY: " << Params.name << '\n';
  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << Params;
  file << "**************************************************\n\n\n\n\n";
  file.close();
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::write(const std::string filename,
                                                bool append) const {
  if (append) {
    std::ofstream file;
    file.open(filename, ios::app);
    file << "\n\n\n\n\n";
    file << "##################################################\n";
    file << "############### COUNTRY PARAMETERS ###############\n";
    file << "##################################################\n";
  }
  for (unsigned int i = 0; i < this->getNcountries(); ++i)
    this->write(filename, i, (append || i));
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::writeSolution(
    const int writeLevel, std::string filename) const {
  /**
   * @brief Writes the computed Nash Equilibrium in the EPEC instance
   * @p writeLevel is an integer representing the write configuration. 0: only
   * Json solution; 1: only human readable solution; 2:both
   */
  if (this->Stats.status == Game::EPECsolveStatus::nashEqFound) {
    if (writeLevel == 1 || writeLevel == 2) {
      this->WriteCountryMCprice(filename + ".txt", this->sol_x, false);
      for (unsigned int ell = 0; ell < this->getNcountries(); ++ell)
        this->WriteCountry(ell, filename + ".txt", this->sol_x, true);
      this->write(filename + ".txt", true);
    }
    if (writeLevel == 2 || writeLevel == 0)
      // this->writeSolutionJSON(filename, this->sol_x, this->sol_z);
      ;
  } else {
    cerr << "Error in cci::EPEC<n_Dirty,n_Clean,n_Scen>::writeSolution: no "
            "solution to write."
         << '\n';
  }
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::WriteCountryMCprice(
    const std::string filename, const arma::vec x, const bool append) const {
  using cci::prn;
  std::ofstream file;
  file.open(filename, append ? ios::app : ios::out);
  file << "**************************************************\n";
  file << "INTERNATIONAL PRICES\n";
  file << "**************************************************\n";
  file << prn::label << "Carbon Price: " << prn::val
       << x.at(this->getPosition(this->getNcountries() - 1, LeaderVars::End));

  file << "**************************************************\n";

  file << "\n\n\n";
  // FILE OPERATIONS END
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::WriteCountry(
    const unsigned int i, const std::string filename, const arma::vec x,
    const bool append) const {
  // if (!lcp) return;
  // const LeadLocs& Loc = this->Locations.at(i);
  using cci::prn;
  std::ofstream file;
  file.open(filename, append ? ios::app : ios::out);
  // FILE OPERATIONS START
  const LeadAllPar<n_Scen> &Params = this->AllLeadPars.at(i);
  file << "**************************************************\n";
  file << "COUNTRY: " << Params.name << '\n';
  file << "**************************************************\n\n";
  // Carbon price and trades
  double carbDomPrice = x.at(this->getPosition(i, cci::LeaderVars::CarbPrice));
  double carbImp = x.at(this->getPosition(i, cci::LeaderVars::CarbImp));
  double carbExp = x.at(this->getPosition(i, cci::LeaderVars::CarbExp));
  double carbInit = Params.LeaderParam.carbCreditInit;

  file << "Carbon credit details\n";
  file << prn::label << "Initial credit: " << prn::val << carbInit << "\n";
  file << prn::label << "Net Export: " << prn::val << carbExp - carbImp << "\n";
  file << prn::label << "Domestic price: " << prn::val << carbDomPrice << "\n";

  // Follower productions
  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << "FOLLOWER DETAILS:\n";
  file.close();
  for (unsigned int j = 0; j < Params.n_followers; ++j)
    this->WriteFollower(i, j, filename, x);
  file.open(filename, ios::app);

  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << "AGGREGATED DETAILS:\n";

  double Prod{0}, Price{0};
  file << "\n";
  for (unsigned int scen = 0; scen < n_Scen; scen++) {
    double prob = Params.scenProb.at(scen);
    file << "Scenario " << scen + 1 << " with a probability " << prob << '\n';
    double prod{0};
    for (unsigned int j = 0; j < Params.n_followers; j++)
      prod += dataMap.at("prod_" + std::to_string(i) + "_" + std::to_string(j) +
                         "_" + std::to_string(scen));
    file << prn::label << " Net production: " << prn::val << prod << "\n";
    double price = Params.DemandParam.at(scen).first -
                   Params.DemandParam.at(scen).second* prod;
    file << prn::label << " Domestic energy price: " << prn::val << price
         << "\n";
    Prod += (prod * prob);
    Price += (price * prob);
    file << '\n';
  }

  file << "\nEXPECTED VALUES\n";
  file << prn::label << "NET PRODUCTION: " << prn::val << Prod << "\n";
  file << prn::label << "DOMESTIC ENERGY PRICE: " << prn::val << Price << "\n";
  //
  // file << "\n\n\n";
  // FILE OPERATIONS END
}
//
template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::WriteFollower(
    const unsigned int i, ///< Leader number
    const unsigned int j, ///< Follower number
    const std::string filename, const arma::vec x) const {
  using cci::prn;
  std::ofstream file;
  file.open(filename, ios::app);

  const LeadAllPar<n_Scen> &Params = this->AllLeadPars.at(i);
  //
  std::string name;
  try {
    name = Params.name + " --- " + Params.FollowerParam.at(j).name;
  } catch (...) {
    name = "Follower " + std::to_string(j) + " of leader " + std::to_string(i);
  }
  file << "\n" << name << "\n\n";

  const unsigned int foll_loc{this->getPosition(i, LeaderVars::Followers) +
                              j * FollVarCount};
  file << prn::label << "CarbCred purchased: " << prn::val
       << x.at(foll_loc + FollCarbBuy) << '\n';
  file << prn::label << "CarbCred sold: " << prn::val
       << x.at(foll_loc + FollCarbSel) << '\n';
  double carbNet = x.at(foll_loc + FollCarbBuy) - x.at(foll_loc + FollCarbSel);
  file << prn::label << (carbNet > 0 ? "Net purchased: " : "Net sold: ")
       << prn::val << std::abs(carbNet) << '\n';
  file << "Clean Investments\n";
  for (unsigned int ii = 0; ii < n_Clean; ++ii) {
    file << prn::label << "Investments in " + this->cleanEnergy.at(ii) + ": "
         << prn::val << x.at(foll_loc + FollInv + ii) << '\n';
  }

  double Prod{0};
  // , Trade{0};
  for (unsigned int scen = 0; scen < n_Scen; scen++) {
    double prob = Params.scenProb.at(scen);
    double prod{0};
    file << "\nScenario " << scen + 1 << " with a probability " << prob << '\n';
    file << "\nDirty means of production quantities: \n";
    unsigned int loc{foll_loc + FollProdDirty + scen * n_Dirty};
    for (unsigned int ii = 0; ii < n_Dirty; ++ii) {
      prod += x.at(loc + ii);
      file << prn::label << this->dirtyEnergy.at(ii) + ": " << prn::val
           << x.at(loc + ii) << '\n';
    }

    file << "\nClean means of production: \n";
    loc = {foll_loc + FollProdClean + scen * n_Clean};
    for (unsigned int ii = 0; ii < n_Clean; ++ii) {
      prod += x.at(loc + ii);
      file << prn::label << this->cleanEnergy.at(ii) + ": " << prn::val
           << x.at(loc + ii) << '\n';
    }

    Prod += prod * prob;
    dataMap.insert(std::pair<std::string, double>(
        "prod_" + std::to_string(i) + " " + std::to_string(j) + "_" +
            std::to_string(scen),
        Prod));
    file << prn::label << " Quantity produced: " << prn::val << prod << "\n";
    // double trade = x.at(foll_loc + FollCbuyScen + scen) -
    //                x.at(foll_loc + FollCselScen + scen);
    // Trade += trade * prob;
    // file << prn::label
    //      << ((trade > 0) ? " CarbCred purchased in 2nd stage mkt: "
    //                      : " CarbCred sold in 2nd stage mkt:")
    //      << prn::val << std::abs(trade) << "\n";
  }

  file << "\nEXPECTED VALUES\n";
  file << prn::label << " QUANTITY PRODUCED: " << prn::val << Prod << "\n";
  //
  dataMap.insert(std::pair<std::string, double>(
      "expProd_" + std::to_string(i) + " " + std::to_string(j), Prod));
  // file << prn::label
  //      << ((Trade > 0) ? " CARBCRED PURCHASED IN 2ND STAGE MKT: "
  //                      : " CARBCRED SOLD IN 2ND STAGE MKT:")
  //      << prn::val << std::abs(Trade) << "\n";
  file << '\n';

  file.close();
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
bool cci::EPEC<n_Dirty, n_Clean, n_Scen>::ParamValid(
    const LeadAllPar<n_Scen> &Params ///< Object whose validity is to be tested
) const
/**
 * @brief Checks the Validity of cci::LeadAllPar object
 */
{
  std::cout << Params.name << " parameters not checked! Warning!!!\n";
  return true;
}

template <unsigned int n_Scen>
std::ostream &cci::operator<<(std::ostream &ost,
                              const cci::LeadAllPar<n_Scen> P) {
  ost << "\n\n";
  ost << "***************************"
      << "\n";
  ost << "Leader Complete Description"
      << "\n";
  ost << "***************************"
      << "\n"
      << "\n";
  ost << cci::prn::label << "Number of followers"
      << ":" << cci::prn::val << P.n_followers << "\n "
      << "\n";
  ost << '\n' << P.LeaderParam << '\n';
  // ost << "Follower Parameters: " << '\n';
  // ost << "********************" << '\n';
  // for (unsigned int ff = 0; ff < P.n_followers; ff++)
  //   ost << P.FollowerParam.at(ff) << '\n';
  ost << "\n Random scenario probabilities: ";
  for (unsigned int i = 0; i < n_Scen; ++i) {
    ost << "Scenario " << i + 1
        << " happens with probability: " << P.scenProb.at(i) << '\n';
  }
  ost << "\n";
  // ost 	<< P.DemandParam << "\n";
  ost << "\n";
  return ost;
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
std::ostream &cci::operator<<(std::ostream &ost, const cci::FollPar<n_Scen> P) {
  using cci::prn;

  ost << "Follower: " << P.name << '\n';
  for (const auto &ccc : P.capacities) {
    const auto prodName = ccc.first;
    ost << "\t *****" << prodName << "*****\n";
    const auto prodCost = P.productionCosts.at(prodName);
    ost << "\t " << prn::label << "Prodn cost: " << prn::val << prodCost.first
        << "q+0.5" << prodCost.second << "q*q\n";
    ost << "\t" << prn::label << "Capacity: " << prn::val << ccc.second << '\n';
    if (P.investmentCosts.find(prodName) != P.investmentCosts.end()) {
      const auto invc = P.investmentCosts.at(prodName);
      ost << "\t" << prn::label << "Investment cost:" << prn::val << invc.first
          << "y + 0.5" << invc.second << "y*y\n";
    }
    ost << "\t" << prn::label << "Emission cost: " << prn::val
        << P.emissionCosts.at(prodName) << '\n';
  }
  ost << '\n';
  return ost;
}

// template <int n_Scen>
// std::ostream &operator<<(std::ostream &ost,
// const std::array<cci::DemPar, n_Scen> P) {
// ost << "Demand Parameters: " << '\n';
// ost << "******************" << '\n';
// for (unsigned int i = 0; i < n_Scen; ++i)
// ost << "\tScenario " << i + 1 << ": P\t\t =\t\t " << P[i].alpha << "\t-\t"
// << P[i].beta << cci::prn::label << " x Q\n";
// return ost;
// }

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::initializeSoln(
    arma::vec &init_x) const {
  double chiToArg_Carb{20};

  init_x.at(this->getPosition(0, cci::LeaderVars::CarbExp)) = chiToArg_Carb;
  init_x.at(this->getPosition(1, cci::LeaderVars::CarbExp)) = -chiToArg_Carb;
  init_x.at(this->getPosition(0, cci::LeaderVars::CarbImp)) = -chiToArg_Carb;
  init_x.at(this->getPosition(1, cci::LeaderVars::CarbImp)) = chiToArg_Carb;

  init_x.at(this->getPosition(1, cci::LeaderVars::End)) = -2000;
  init_x.at(this->getPosition(1, cci::LeaderVars::End) + 1) = -1000;

  // init_x.at(this->getPosition(0,
}
