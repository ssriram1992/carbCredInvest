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
  /// Renewable capacity adjustment for each scenario
  std::array<double, n_Scen> renewCapAdjust;
  ///	Linear and QUadratic Investment costs for renewables
  std::map<std::string, std::pair<double, double>> investmentCosts;
  /// Emission costs for unit quantity of the fuel. Emission costs feature only
  /// on the leader's problem
  std::map<std::string, double> emissionCosts;

  double carbonCreditInit = 0;
  std::string name = {}; ///< Optional Names for the Followers.
  FollPar() = default;
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

  LeadPar(double imp_lim = -1, double exp_lim = -1, double consum_limit = 0,
          double carbCred = 0)
      : import_limit{imp_lim}, export_limit{exp_lim},
        consum_limit{consum_limit}, carbCreditInit{carbCred} {}
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
  DualVar,
  ConvHullDummy,
  End
};

template <unsigned int n_Scen>
std::ostream &operator<<(std::ostream &ost, const FollPar<n_Scen> P);

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

/********************/
// Class definition
/********************/
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
  const std::vector<const std::string> dirtyEnergy;
  const std::vector<const std::string> cleanEnergy;
  /// cci::energy should contain everything in dirtyEnergy and cleanEnergy
  const std::vector<std::string> energy;
  static constexpr unsigned int LL_MC_count =
      0; // No market clearing for follower

  static constexpr unsigned int FollInv{0};
  static constexpr unsigned int FollProdDirty{n_Clean};
  static constexpr unsigned int FollProdClean{FollProdDirty+n_Dirty*n_Scen};
  static constexpr unsigned int FollCarbBuy{n_Clean + (n_Clean+n_Dirty)*n_Scen};
  static constexpr unsigned int FollCarbSel{FollCarbBuy+1};
  static constexpr unsigned int FollEnd{FollCarbSel+1};

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
                        const unsigned int export_lim_cons = 1,
                        const unsigned int consum_lim_cons = 1) const noexcept;

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

  EPEC(GRBEnv *env, const std::vector<const std::string> dirtyEnergy,
       const std::vector<const std::string> cleanEnergy)
      : Game::EPEC(env), dirtyEnergy{dirtyEnergy}, cleanEnergy{cleanEnergy} {
    if (dirtyEnergy.size() != n_Dirty || cleanEnergy.size() != n_Clean) {
      throw std::string(
          "Error in EPEC(): Illegal size of dirtyEnergy/cleanEnergy");
    }
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

  QP_obj.Q.zeros(nThisCountryvars, nThisCountryvars);
  QP_obj.c.zeros(nThisCountryvars);
  QP_obj.C.zeros(nThisCountryvars, nEPECvars - nThisCountryvars);

  // Total emission
  for (unsigned int j = 0; j < Params.n_followers; j++) {
    for (unsigned int ell = 0; ell < n_Scen; ell++) {
      QP_obj.c.at(FollVarCount * j + FollProdScen + ell) =
          Params.scenProb.at(ell) * Params.FollowerParam.emission_costs.at(j) *
          social_cost_of_carbon;
    }
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
  // Energy international trade term
  for (unsigned int j = 0; j < n_Scen; ++j) {
    QP_obj.C.at(
        Loc.at(LeaderVars::EnergExpScen) + j,
        this->getPosition(this->getNcountries() - 1, cci::LeaderVars::End) -
            nThisCountryvars + 1 + j) = Params.scenProb.at(j);
    QP_obj.C.at(
        Loc.at(LeaderVars::EnergImpScen) + j,
        this->getPosition(this->getNcountries() - 1, cci::LeaderVars::End) -
            nThisCountryvars + 1 + j) = -Params.scenProb.at(j);
  }
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
cci::EPEC<n_Dirty, n_Clean, n_Scen> &
cci::EPEC<n_Dirty, n_Clean, n_Scen>::addCountry(
    cci::LeadAllPar<n_Scen> Params) {
  // const unsigned int addnlLeadVars) {
  if (this->finalized)
    throw std::string(
        "Error in cci::EPEC<n_Dirty,n_Clean,n_Scen>::addCountry: EPEC object "
        "finalized. Call "
        "EPEC::unlock() to unlock this object first and then edit.");

  bool noError = false;
  try {
    noError = this->ParamValid(Params);
  } catch (const char *e) {
    cerr << "Error in cci::EPEC<n_Dirty,n_Clean,n_Scen>::addCountry: " << e
         << '\n';
  } catch (std::string e) {
    cerr << "String: Error in cci::EPEC<n_Dirty,n_Clean,n_Scen>::addCountry: "
         << e << '\n';
  } catch (std::exception &e) {
    cerr
        << "Exception: Error in cci::EPEC<n_Dirty,n_Clean,n_Scen>::addCountry: "
        << e.what() << '\n';
  }

  if (!noError)
    return *this;

  constexpr unsigned int LeadVars = 4;
  // One each for Total Carbon imported, total carbon exported, Carbon price set
  // and total renewable investment in the country

  LeadLocs Loc;
  cci::init(Loc);

  // Allocate so much space for each of these types of variables
  // For follower's carbon purchase from government.
  cci::increaseVal(Loc, LeaderVars::Followers,
                   Params.n_followers * this->FollVarCount);
  cci::increaseVal(Loc, LeaderVars::CarbExp, 1);
  cci::increaseVal(Loc, LeaderVars::CarbImp, 1);
  cci::increaseVal(Loc, LeaderVars::CarbPrice, 1);
  cci::increaseVal(Loc, LeaderVars::TotInv, 1);

  // Leader Constraints
  constexpr short int consum_lim_cons{n_Scen};
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
                            consum_lim_cons + // Min consumption
                            2 +               // Investment summing
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
                           export_lim_cons, consum_lim_cons);
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

  // std::cout << "Nx, Ny, Ncons\n";
  // for(const auto &pl:FollowersVec)
  // {
  // std::cout<<pl->getNx()<<" "<<pl->getNy()<<" "<<pl->getA().n_rows<<"\n";
  // }
  // std::cout<<"MC: "<<MC.n_rows<<" "<<MC.n_cols<<"\n";
  // std::cout<< "MCRHS: "<<MCRHS.n_rows<<"\n";
  // std::cout <<" LeadCons: "<<LeadCons.n_rows<<" "<<LeadCons.n_cols<<"\n";
  // std::cout<<"LeaderVars: "<<LeadVars<<" " <<this->FollVarCount<<"\n";

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
  this->Game::EPEC::n_MCVar = 1 + n_Scen;
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
    Q(FollInv+ii,FollInv+ ii) = Follparam.investmentCosts.at(cleanEnergy.at(ii)).second;
    c(FollInv+ii) = Follparam.investmentCosts.at(cleanEnergy.at(ii)).first;
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
  // Clean energy Production cost
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
        const auto pos1{FolProdDirty + n_Dirty * scen + ii};
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
	for(unsigned int scen=0; scen<n_Scen;

  for (unsigned int i = 0; i < n_Scen; ++i) {
    // Quantity produced by other folowers times quantity produced by self terms
    // - for each scenario
    for (unsigned int j = 0; j < Params.n_followers - 1; ++j) {
      C(FollProdScen + i, FollVarCount * j + FollProdScen + i) =
          Params.DemandParam.at(i).beta * Params.scenProb.at(i);
    }
  }
  // OBJECTIVE Described

  // CONSTRAINTS
  for (unsigned int i = 0; i < n_Scen; ++i) {
    // Infrastructural capacity constraint
    B(i, FollProdScen + i) = 1;
    b(i) = Params.FollowerParam.capacities.at(i).at(follower);
    // Carbon credit limit constraint
    B(n_Scen + i, FollCarb) = -1;
    B(n_Scen + i, FollProdScen + i) =
        Params.FollowerParam.emission_costs.at(follower);
    B(n_Scen + i, FollCbuyScen + i) = -1;
    B(n_Scen + i, FollCselScen + i) = 1;
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
    const unsigned int export_lim_cons, ///< Does a constraint on export limit
                                        ///< exist or no limit?
    const unsigned int consum_lim_cons  ///< Does a constraint on consumption
                                        ///< limit exist or no limit?
) const noexcept
/**
 * Makes the leader level constraints for a country.
 */
{
  // LeadCons and LeadRHS already have the right size
  LeadCons.zeros();
  LeadRHS.zeros();
  // Scenario constraints
  for (unsigned int i = 0; i < n_Scen; i++) {
    // Import limit constraints
    if (import_lim_cons) {
      LeadCons(i, Loc.at(LeaderVars::EnergExpScen) - this->LL_MC_count + i) =
          -1;
      LeadCons(i, Loc.at(LeaderVars::EnergImpScen) - this->LL_MC_count + i) = 1;
      LeadRHS(i) = Params.LeaderParam.import_limit;
    }
    // export limit constraints
    if (export_lim_cons) {
      LeadCons(import_lim_cons + i,
               Loc.at(LeaderVars::EnergExpScen) - this->LL_MC_count + i) = 1;
      LeadCons(import_lim_cons + i,
               Loc.at(LeaderVars::EnergImpScen) - this->LL_MC_count + i) = -1;
      LeadRHS(import_lim_cons + i) = Params.LeaderParam.export_limit;
    }
    // Minimum Consumption  constraint
    LeadCons(import_lim_cons + export_lim_cons + i,
             Loc.at(LeaderVars::EnergExpScen) - this->LL_MC_count + i) = 1;
    LeadCons(import_lim_cons + export_lim_cons + i,
             Loc.at(LeaderVars::EnergImpScen) - this->LL_MC_count + i) = -1;
    for (unsigned int j = 0; j < Params.n_followers; j++)
      LeadCons(import_lim_cons + export_lim_cons + i,
               FollVarCount * j + FollProdScen + i) = -1;

    LeadRHS(import_lim_cons + export_lim_cons + i) =
        -std::max(0.0, Params.LeaderParam.consum_limit);
  }
  // Carbon credit positivity constraint
  const unsigned int lastCons =
      import_lim_cons + export_lim_cons + consum_lim_cons;
  LeadCons(lastCons, Loc.at(LeaderVars::CarbExp) - this->LL_MC_count) = 1;
  LeadCons(lastCons, Loc.at(LeaderVars::CarbImp) - this->LL_MC_count) = -1;
  for (unsigned int j = 0; j < Params.n_followers; j++)
    LeadCons(lastCons, FollVarCount * j + FollCarb) = 1;
  LeadRHS(lastCons) = Params.LeaderParam.carbCredit_init;
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
  MCRHS.zeros(1 + n_Scen);
  MCLHS.zeros(1 + n_Scen, this->getnVarinEPEC());
  if (this->getNcountries() > 1) {
    for (unsigned int i = 0; i < this->getNcountries(); ++i) {
      // The MC constraint for each leader country's carbon trade market
      MCLHS(0, this->getPosition(i, LeaderVars::CarbExp)) = 1;
      MCLHS(0, this->getPosition(i, LeaderVars::CarbImp)) = -1;
      // The MC constraint for each leader country's energy trade market,
      // scenario-wise
      for (unsigned int j = 0; j < n_Scen; ++j) {
        MCLHS(1 + j, this->getPosition(i, LeaderVars::EnergExpScen) + j) = 1;
        MCLHS(1 + j, this->getPosition(i, LeaderVars::EnergImpScen) + j) = -1;
      }
    }
  }
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
bool cci::EPEC<n_Dirty, n_Clean, n_Scen>::dataCheck(
    const bool
        chkAllLeadPars, ///< Checks if
                        ///< cci::EPEC<n_Dirty,n_Clean,n_Scen>::AllLeadPars
                        ///< has size @p n
    const bool
        chkcountries_LL, ///< Checks if
                         ///< cci::EPEC<n_Dirty,n_Clean,n_Scen>::countries_LL
                         ///< has
    ///< size @p n
    const bool chkMC_QP, ///< Checks if cci::EPEC<n_Dirty,n_Clean,n_Scen>::MC_QP
                         ///< has size @p n
    const bool chkLeadConses, ///< Checks if
                              ///< cci::EPEC<n_Dirty,n_Clean,n_Scen>::LeadConses
                              ///< has size @p n
    const bool chkLeadRHSes,  ///< Checks if
                              ///< cci::EPEC<n_Dirty,n_Clean,n_Scen>::LeadRHSes
                              ///< has size @p n
    const bool chkLocations,  ///< Checks if
                              ///< cci::EPEC<n_Dirty,n_Clean,n_Scen>::Locations
                              ///< has size @p n
    const bool
        chkLeaderLocations, ///< Checks if
                            ///< cci::EPEC<n_Dirty,n_Clean,n_Scen>::LeaderLocations
                            ///< has
    ///< size @p n and cci::EPEC<n_Dirty,n_Clean,n_Scen>::nVarinEPEC is set
    const bool
        chkLeadObjec ///< Checks if cci::EPEC<n_Dirty,n_Clean,n_Scen>::LeadObjec
                     ///< has size @p n
) const
/**
 * Checks the data in cci::EPEC object, based on checking flags, @p n is the
 * number of countries in the cci::EPEC object.
 */
{
  if (!chkAllLeadPars && AllLeadPars.size() != this->getNcountries())
    return false;
  if (!chkcountries_LL && countries_LL.size() != this->getNcountries())
    return false;
  if (!chkMC_QP && MC_QP.size() != this->getNcountries())
    return false;
  if (!chkLeadConses && LeadConses.size() != this->getNcountries())
    return false;
  if (!chkLeadRHSes && LeadRHSes.size() != this->getNcountries())
    return false;
  if (!chkLocations && Locations.size() != this->getNcountries())
    return false;
  if (!chkLeaderLocations && LeaderLocations.size() != this->getNcountries())
    return false;
  if (!chkLeaderLocations && this->getnVarinEPEC() == 0)
    return false;
  if (!chkLeadObjec && LeadObjec.size() != this->getNcountries())
    return false;
  return true;
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
  const LeadAllPar<n_Scen> &Params = this->AllLeadPars.at(0);
  double Price{0};
  for (unsigned int scen = 0; scen < n_Scen; scen++) {
    double prob = Params.scenProb.at(scen);
    file << "\nScenario " << scen + 1 << " with a probability " << prob << '\n';
    double price{x.at(
        this->getPosition(this->getNcountries() - 1, LeaderVars::End) + scen)};
    file << prn::label << "Energy Price: " << prn::val << price;
    Price += price * prob;
  }
  file << prn::label << "\nEXPECTED ENERGY PRICE: " << prn::val << Price
       << "\n";
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
  double carbInit = Params.LeaderParam.carbCredit_init;

  file << "Carbon credit details\n";
  file << prn::label << "Initial credit: " << prn::val << carbInit << "\n";
  ;
  file << prn::label << "Net Export: " << prn::val << carbExp - carbImp << "\n";
  ;
  file << prn::label << "Domestic price: " << prn::val << carbDomPrice << "\n";
  ;
  double Expo{0}, Prod{0}, Price{0}, C2p{0};
  file << "\n";
  for (unsigned int scen = 0; scen < n_Scen; scen++) {
    double prob = Params.scenProb.at(scen);
    file << "Scenario " << scen + 1 << " with a probability " << prob << '\n';
    double prod{0};
    for (unsigned int j = 0; j < Params.n_followers; j++)
      prod += x.at(this->getPosition(i, LeaderVars::Followers) +
                   j * FollVarCount + FollProdScen + scen);
    file << prn::label << " Net production: " << prn::val << prod << "\n";
    double expo = x.at(this->getPosition(i, LeaderVars::EnergExpScen) + scen);
    double impo = x.at(this->getPosition(i, LeaderVars::EnergImpScen) + scen);
    file << prn::label << " Net energy export: " << prn::val << (expo - impo)
         << "\n";
    double price = Params.DemandParam.at(scen).alpha -
                   Params.DemandParam.at(scen).beta * (prod + impo - expo);
    file << prn::label << " Domestic energy price: " << prn::val << price
         << "\n";
    double c2p = x.at(this->getPosition(i, LeaderVars::Foll2Cpr) + scen);
    file << prn::label << " Carbon 2ndary mkt price: " << prn::val << c2p
         << "\n";
    Prod += (prod * prob);
    Expo += (expo - impo) * prob;
    Price += (price * prob);
    C2p += (c2p * prob);
    file << '\n';
  }
  file << "\nEXPECTED VALUES\n";
  file << prn::label << "NET PRODUCTION: " << prn::val << Prod << "\n";
  file << prn::label << "NET ENERGY EXPORT: " << prn::val << Expo << "\n";
  file << prn::label << "DOMESTIC CONSUMPTION: " << prn::val << Prod - Expo
       << "\n";
  file << prn::label << "DOMESTIC ENERGY PRICE: " << prn::val << Price << "\n";
  file << prn::label << "CARBON 2nd PRICE: " << prn::val << C2p << "\n";

  // Follower productions
  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << "FOLLOWER DETAILS:\n";
  file.close();
  for (unsigned int j = 0; j < Params.n_followers; ++j)
    this->WriteFollower(i, j, filename, x);
  //
  // file << "\n\n\n";
  // FILE OPERATIONS END
}
//
template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::WriteFollower(
    const unsigned int i, const unsigned int j, const std::string filename,
    const arma::vec x) const {
  using cci::prn;
  std::ofstream file;
  file.open(filename, ios::app);

  const LeadAllPar<n_Scen> &Params = this->AllLeadPars.at(i);
  //
  std::string name;
  try {
    name = Params.name + " --- " + Params.FollowerParam.names.at(j);
  } catch (...) {
    name = "Follower " + std::to_string(j) + " of leader " + std::to_string(i);
  }
  file << "\n" << name << "\n\n";

  const unsigned int foll_loc{this->getPosition(i, LeaderVars::Followers) +
                              j * FollVarCount};
  file << prn::label << "CarbCred purchased: " << prn::val
       << x.at(foll_loc + FollCarb) << '\n';
  double Prod{0}, Trade{0};
  for (unsigned int scen = 0; scen < n_Scen; scen++) {
    double prob = Params.scenProb.at(scen);
    file << "\nScenario " << scen + 1 << " with a probability " << prob << '\n';
    double prod = x.at(foll_loc + FollProdScen + scen);
    Prod += prod * prob;
    file << prn::label << " Quantity produced: " << prn::val << prod << "\n";
    double trade = x.at(foll_loc + FollCbuyScen + scen) -
                   x.at(foll_loc + FollCselScen + scen);
    Trade += trade * prob;
    file << prn::label
         << ((trade > 0) ? " CarbCred purchased in 2nd stage mkt: "
                         : " CarbCred sold in 2nd stage mkt:")
         << prn::val << std::abs(trade) << "\n";
  }

  file << "\nEXPECTED VALUES\n";
  file << prn::label << " QUANTITY PRODUCED: " << prn::val << Prod << "\n";
  file << prn::label
       << ((Trade > 0) ? " CARBCRED PURCHASED IN 2ND STAGE MKT: "
                       : " CARBCRED SOLD IN 2ND STAGE MKT:")
       << prn::val << std::abs(Trade) << "\n";
  file << '\n';

  file.close();
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
bool cci::EPEC<n_Dirty, n_Clean, n_Scen>::ParamValid(
    const LeadAllPar<n_Scen> &Params ///< Object whose validity is to be tested
) const
/**
 * @brief Checks the Validity of cci::LeadAllPar object
 * @details Checks the following:
 * 	-	Size of FollowerParam.costs_lin, FollowerParam.costs_quad,
 * FollowerParam.capacities, FollowerParam.emission_costs are all equal to @p
 * Params.n_followers -	@p DemandParam.alpha and @p DemandParam.beta are greater
 * than zero -	@p name is not empty -	@p name does not match with the name of
 * any other existing countries in the EPEC object.
 */
{
  if (Params.n_followers == 0)
    throw std::string("Error in EPEC::ParamValid(). 0 Followers?");
  if (Params.FollowerParam.costs_lin.size() != Params.n_followers ||
      Params.FollowerParam.costs_quad.size() != Params.n_followers ||
      Params.FollowerParam.emission_costs.size() != Params.n_followers)
    throw std::string("Error in EPEC::ParamValid(). Size Mismatch");
  for (unsigned int i = 0; i < n_Scen; ++i) {
    if (Params.DemandParam[i].alpha <= 0 || Params.DemandParam[i].beta <= 0)
      throw std::string(
          "Error in EPEC::ParamValid(). Invalid demand curve params");
    if (Params.FollowerParam.capacities[i].size() != Params.n_followers)
      throw std::string(
          "Error in EPEC::ParamValid(). Size Mismatch in capacity");
  }
  // Country should have a name!
  if (Params.name == "")
    throw std::string("Error in EPEC::ParamValid(). Country name empty");
  // Country should have a unique name
  for (const auto &p : this->AllLeadPars)
    if (Params.name.compare(p.name) == 0) // i.e., if the strings are same
      throw std::string("Error in EPEC::ParamValid(). Country name repetition");
  return true;
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
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
  ost << '\n' << P.LeaderParam << '\n' << P.FollowerParam << '\n';
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
  ost << "Follower Parameters: " << '\n';
  ost << "********************" << '\n';
  ost << cci::prn::label << "Linear Costs"
      << ":\t";
  for (auto a : P.costs_lin)
    ost << cci::prn::val << a;
  ost << '\n'
      << cci::prn::label << "Quadratic costs"
      << ":\t";
  for (auto a : P.costs_quad)
    ost << cci::prn::val << a;
  ost << '\n'
      << cci::prn::label << "Production capacities"
      << ":\n";
  // for (unsigned int i = 0; i < n_Scen; i++) {
  // for (auto sc : P.capacities.at(i)) {
  // ost << "\t Scenario " << i + 1 << ": (";
  // for (auto a : sc)
  // ost << cci::prn::val
  // << (a < 0 ? std::numeric_limits<double>::infinity() : a);
  // ost << ")\n ";
  // }
  // }
  ost << '\n'
      << cci::prn::label << "Carbon emission"
      << ":\t";
  for (auto a : P.emission_costs)
    ost << cci::prn::val << a;
  ost << '\n';
  return ost;
}

template <int n_Scen>
std::ostream &operator<<(std::ostream &ost,
                         const std::array<cci::DemPar, n_Scen> P) {
  ost << "Demand Parameters: " << '\n';
  ost << "******************" << '\n';
  for (unsigned int i = 0; i < n_Scen; ++i)
    ost << "\tScenario " << i + 1 << ": P\t\t =\t\t " << P[i].alpha << "\t-\t"
        << P[i].beta << cci::prn::label << " x Q\n";
  return ost;
}

template <unsigned int n_Scen>
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
