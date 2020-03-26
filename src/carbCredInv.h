#pragma once

/**
 * @file src/bienestar.h Using EPECSolve to solve problems arising in
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
#include <tuple>

namespace cci {

	template <unsigned int N> using gamsSet = const std::array<const std::string, N>;
	template <unsigned int N> using gamsNumSet = const std::array<const unsigned int, N>;
template<class t> using params = std::map<t, double>; 
template <class t, class u> using params2 = std::map < std::tuple<t,u> , double>;
template <class t, class u, class v> using params3 = std::map < std::tuple<t,u,v> , double>;

	gamsSet<3> countries {"c1", "c2", "c3"};
	gamsSet<3> producers {"p1", "p2", "p3"};
	gamsSet<3> prodtypes {"c1", "c2", "c3"};


/// @brief Stores the parameters of the follower in a country model
template <unsigned int num_scen> struct FollPar {
  std::vector<double> costs_quad =
      {}; ///< Quadratic coefficient of i-th follower's cost. Size of this
          ///< std::vector should be equal to n_followers
  std::vector<double> costs_lin =
      {}; ///< Linear  coefficient of i-th follower's cost. Size of this
          ///< std::vector should be equal to n_followers
  std::array<std::vector<double>, num_scen>
      capacities; ///< Production capacity of each follower. Size of this
                  ///< std::vector should be equal to n_followers
  std::vector<double> emission_costs =
      {}; ///< Emission costs for unit quantity of the fuel. Emission costs
          ///< feature only on the leader's problem
  std::vector<std::string> names = {}; ///< Optional Names for the Followers.
  FollPar(std::vector<double> costs_quad_, std::vector<double> costs_lin_,
          std::array<std::vector<double>, num_scen> capacities_,
          std::vector<double> emission_costs_ = {},
          std::vector<std::string> names_ = {})
      : costs_quad{costs_quad_}, costs_lin{costs_lin_},
        emission_costs{emission_costs_}, names{names_} {
    for (unsigned int i = 0; i < num_scen; ++i)
      capacities[i] = capacities_[i];
  }
};

/// @brief Stores the parameters of the demand curve in a country model
struct DemPar {
  double alpha = 100; ///< Intercept of the demand curve. Written as: Price =
                      ///< alpha - beta*(Total quantity in domestic market)
  double beta = 2; ///< Slope of the demand curve. Written as: Price = alpha -
                   ///< beta*(Total quantity in domestic market)
  DemPar(double alpha = 100, double beta = 2) : alpha{alpha}, beta{beta} {};
};

/// @brief Stores the parameters of the leader in a country model
struct LeadPar {
  double import_limit = -1; ///< Maximum net import in the country. If no limit,
                            ///< set the value as -1;
  double export_limit = -1; ///< Maximum net export in the country. If no limit,
                            ///< set the value as -1;
  double consum_limit =
      0; ///< Government does not want the price to exceed this limit
  double carbCredit_init = 100;

  LeadPar(double imp_lim = -1, double exp_lim = -1, double consum_limit = 0,
          double carbCred = 0)
      : import_limit{imp_lim}, export_limit{exp_lim},
        consum_limit{consum_limit}, carbCredit_init{carbCred} {}
};

/// @brief Stores the parameters of a country model
template <unsigned int num_scen> struct LeadAllPar {
  unsigned int n_followers; ///< Number of followers in the country
  std::string name;         ///< Country Name
  cci::FollPar<num_scen> FollowerParam =
      {};                         ///< A struct to hold Follower Parameters
  cci::LeadPar LeaderParam = {}; ///< A struct to hold Leader Parameters
  std::array<cci::DemPar, num_scen> DemandParam = {
      {}}; ///< A struct to hold Demand Parameters
  std::array<double, num_scen> scenProb = {{}};
  LeadAllPar(unsigned int n_foll, std::string name, cci::FollPar<num_scen> FP,
             std::array<cci::DemPar, num_scen> DP, cci::LeadPar LP,
             std::array<double, num_scen> probab)
      : n_followers{n_foll}, name{name}, FollowerParam{FP}, LeaderParam{LP} {
    for (unsigned int i = 0; i < num_scen; ++i) {
      DemandParam[i] = DP[i];
      scenProb[i] = probab[i];
    }
  }
};

/// @brief Stores a single Instance
template <unsigned int num_scen> struct EPECInstance {
  std::vector<cci::LeadAllPar<num_scen>> Countries = {}; ///< LeadAllPar vector

  EPECInstance(std::string filename) {
    this->load(filename);
  } ///< Constructor from instance file
  EPECInstance(std::vector<cci::LeadAllPar<num_scen>> Countries_)
      : Countries{Countries_} {}
  ///< Constructor from instance objects

  void load(std::string filename);
  ///< Reads the EPECInstance from a file

  void save(std::string filename);
  ///< Writes the EPECInstance from a file
};

enum class LeaderVars {
  Followers,
  Foll2Cpr,
  CarbPrice,
  CarbExp,
  CarbImp,
  EnergExpScen,
  EnergImpScen,
  DualVar,
  ConvHullDummy,
  End
};

template <unsigned int num_scen>
std::ostream &operator<<(std::ostream &ost, const FollPar<num_scen> P);

std::ostream &operator<<(std::ostream &ost, const DemPar P);

template <int num_scen>
std::ostream &operator<<(std::ostream &ost,
                         const std::array<DemPar, num_scen> P);

std::ostream &operator<<(std::ostream &ost, const LeadPar P);

template <unsigned int num_scen>
std::ostream &operator<<(std::ostream &ost, const LeadAllPar<num_scen> P);

std::ostream &operator<<(std::ostream &ost, const LeaderVars l);

template <unsigned int num_scen>
std::ostream &operator<<(std::ostream &ost, EPECInstance<num_scen> I);

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
template <unsigned int num_scen> class EPEC : public Game::EPEC {
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
  virtual void initializeSoln(arma::vec &init_x) const override;
  // override;

public:
  // Rest
  static constexpr double social_cost_of_carbon{1000};

private:
  static constexpr unsigned int LL_MC_count = num_scen;
  static constexpr unsigned int FollVarCount{1 + 3 * num_scen};

  static constexpr unsigned int FollCarb{0};
  static constexpr unsigned int FollProdScen{1};
  static constexpr unsigned int FollCbuyScen{1 + num_scen};
  static constexpr unsigned int FollCselScen{1 + 2 * num_scen};
  static constexpr unsigned int FollEnd{1 + 3 * num_scen};

  std::vector<LeadAllPar<num_scen>> AllLeadPars =
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

  bool dataCheck(const bool chkAllLeadPars = true,
                 const bool chkcountriesLL = true, const bool chkMC_QP = true,
                 const bool chkLeadConses = true,
                 const bool chkLeadRHSes = true, const bool chkLocations = true,
                 const bool chkLeaderLocations = true,
                 const bool chkLeadObjec = true) const;

  // Super low level
  /// Checks that the parameter given to add a country is valid. Does not have
  /// obvious errors
  bool ParamValid(const LeadAllPar<num_scen> &Param) const;

  /// Makes the lower level quadratic program object for each follower.
  void make_LL_QP(const LeadAllPar<num_scen> &Params,
                  const unsigned int follower, Game::QP_Param *Foll,
                  const LeadLocs &Loc) const noexcept;

  /// Makes the leader constraint matrix and RHS
  void make_LL_LeadCons(arma::sp_mat &LeadCons, arma::vec &LeadRHS,
                        const LeadAllPar<num_scen> &Param,
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

  EPEC(GRBEnv *env) : Game::EPEC(env) {}

  ///@brief %cci a Standard Nash-Cournot game within a country
  EPEC &addCountry(
      /// The Parameter structure for the leader
      LeadAllPar<num_scen> Params,
      /// Create columns with 0s in it. To handle additional dummy leader
      /// variables.
      const unsigned int addnlLeadVars = 0);

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
  const EPECInstance<num_scen> getInstance() const {
    return EPECInstance<num_scen>(this->AllLeadPars);
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

template <unsigned int num_scen>
cci::FollPar<num_scen> operator+(const cci::FollPar<num_scen> &F1,
                                  const cci::FollPar<num_scen> &F2);

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
template <unsigned int num_scen>
void cci::EPEC<num_scen>::make_obj_leader(
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
  const LeadAllPar<num_scen> &Params = this->AllLeadPars.at(i);
  const LeadLocs &Loc = this->Locations.at(i);

  QP_obj.Q.zeros(nThisCountryvars, nThisCountryvars);
  QP_obj.c.zeros(nThisCountryvars);
  QP_obj.C.zeros(nThisCountryvars, nEPECvars - nThisCountryvars);

  // Total emission
  for (unsigned int j = 0; j < Params.n_followers; j++) {
    for (unsigned int ell = 0; ell < num_scen; ell++) {
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
  for (unsigned int j = 0; j < num_scen; ++j) {
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

template <unsigned int num_scen>
cci::EPEC<num_scen> &
cci::EPEC<num_scen>::addCountry(cci::LeadAllPar<num_scen> Params,
                                 const unsigned int addnlLeadVars) {
  if (this->finalized)
    throw std::string(
        "Error in cci::EPEC<num_scen>::addCountry: EPEC object "
        "finalized. Call "
        "EPEC::unlock() to unlock this object first and then edit.");

  bool noError = false;
  try {
    noError = this->ParamValid(Params);
  } catch (const char *e) {
    cerr << "Error in cci::EPEC<num_scen>::addCountry: " << e << '\n';
  } catch (std::string e) {
    cerr << "String: Error in cci::EPEC<num_scen>::addCountry: " << e << '\n';
  } catch (std::exception &e) {
    cerr << "Exception: Error in cci::EPEC<num_scen>::addCountry: " << e.what()
         << '\n';
  }
  if (!noError)
    return *this;

  constexpr unsigned int LeadVars = num_scen + 3 + 2 * num_scen;
  // C2Price - one per scenario, CarbonPrice, CarbExp, CarbImp - 1 each,
  // EnergyExpScen, EnergyImpScen - one  per scenario

  LeadLocs Loc;
  cci::init(Loc);

  // Allocate so much space for each of these types of variables
  // For follower's carbon purchase from government.
  cci::increaseVal(Loc, LeaderVars::Followers, Params.n_followers);
  // Follower production in each scenario
  cci::increaseVal(Loc, LeaderVars::Followers, Params.n_followers * num_scen);
  // Follower Carbon selling in each scenario 2nd stage
  cci::increaseVal(Loc, LeaderVars::Followers, Params.n_followers * num_scen);
  // Follower carbon buying in each scenario 2nd stage
  cci::increaseVal(Loc, LeaderVars::Followers, Params.n_followers * num_scen);

  cci::increaseVal(Loc, LeaderVars::Foll2Cpr, num_scen);
  cci::increaseVal(Loc, LeaderVars::CarbPrice, 1);
  cci::increaseVal(Loc, LeaderVars::CarbExp, 1);
  cci::increaseVal(Loc, LeaderVars::CarbImp, 1);
  cci::increaseVal(Loc, LeaderVars::EnergExpScen, num_scen);
  cci::increaseVal(Loc, LeaderVars::EnergImpScen, num_scen);

  // Leader Constraints
  constexpr short int consum_lim_cons{num_scen};
  const short int import_lim_cons = [&Params]() {
    if (Params.LeaderParam.import_limit >= 0)
      return num_scen;
    else
      return 0u;
  }();
  const short int export_lim_cons = [&Params]() {
    if (Params.LeaderParam.export_limit >= 0)
      return num_scen;
    else
      return 0u;
  }();

  arma::sp_mat LeadCons(import_lim_cons +     // Import limit constraint
                            export_lim_cons + // Export limit constraint
                            consum_lim_cons + // Min consumption
                            1,                // Carbon credits >=0
                        Loc[cci::LeaderVars::End] - this->LL_MC_count);
  arma::vec LeadRHS(import_lim_cons + export_lim_cons + consum_lim_cons + 1,
                    arma::fill::zeros);

  std::vector<std::shared_ptr<Game::QP_Param>> Players{};
  // Create the QP_Param* for each follower
  try {
    for (unsigned int follower = 0; follower < Params.n_followers; follower++) {
      auto Foll = std::make_shared<Game::QP_Param>(this->env);
      this->make_LL_QP(Params, follower, Foll.get(), Loc);
      Players.push_back(Foll);
    }
    // Make Leader Constraints
    this->make_LL_LeadCons(LeadCons, LeadRHS, Params, Loc, import_lim_cons,
                           export_lim_cons, consum_lim_cons);
  } catch (const char *e) {
    cerr << e << '\n';
    throw;
  } catch (std::string e) {
    cerr << "String in cci::EPEC<num_scen>::addCountry : " << e << '\n';
    throw;
  } catch (GRBException &e) {
    cerr << "GRBException in cci::EPEC<num_scen>::addCountry : "
         << e.getErrorCode() << ": " << e.getMessage() << '\n';
    throw;
  } catch (std::exception &e) {
    cerr << "Exception in cci::EPEC<num_scen>::addCountry : " << e.what()
         << '\n';
    throw;
  }

  // Lower level Market clearing constraints - One for each scenario - 2nd stage
  // carbon clearing
  arma::sp_mat MC(num_scen, LeadVars + FollVarCount * Params.n_followers);
  arma::vec MCRHS(num_scen, arma::fill::zeros);

  for (unsigned int i = 0; i < num_scen; ++i) {
    for (unsigned int j = 0; j < Params.n_followers; j++) {
      MC(i, j * FollVarCount + cci::EPEC<num_scen>::FollCbuyScen + i) = 1;
      MC(i, j * FollVarCount + cci::EPEC<num_scen>::FollCselScen + i) = -1;
    }
  }

  // std::cout << "Nx, Ny, Ncons\n";
  // for(const auto &pl:Players)
  // {
  // std::cout<<pl->getNx()<<" "<<pl->getNy()<<" "<<pl->getA().n_rows<<"\n";
  // }
  // std::cout<<"MC: "<<MC.n_rows<<" "<<MC.n_cols<<"\n";
  // std::cout<< "MCRHS: "<<MCRHS.n_rows<<"\n";
  // std::cout <<" LeadCons: "<<LeadCons.n_rows<<" "<<LeadCons.n_cols<<"\n";
  // std::cout<<"LeaderVars: "<<LeadVars<<" " <<this->FollVarCount<<"\n";

  // Convert the country QP to a NashGame
  auto N = std::make_shared<Game::NashGame>(this->env, Players, MC, MCRHS,
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
  this->Game::EPEC::n_MCVar = 1 + num_scen;
  return *this;
}

template <unsigned int num_scen>
void cci::EPEC<num_scen>::make_LL_QP(
    const cci::LeadAllPar<num_scen> &Params, ///< The Parameters object
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
      Loc.at(cci::LeaderVars::End) - FollVarCount * n_Foll;
  const unsigned int otherVars = LeadVars + (n_Foll - 1) * FollVarCount;

  arma::sp_mat Q(FollVarCount, FollVarCount);
  arma::sp_mat C(FollVarCount, otherVars);
  // Two constraints for each scenario. One saying that you should be less than
  // the infrastructural capacity. Another saying that you should produce less
  // than the carbon credit availability.
  constexpr unsigned int n_Const = 2 * num_scen;
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
  for (unsigned int i = 0; i < num_scen; ++i) {
    Q(FollProdScen + i, FollProdScen + i) =
        (Params.FollowerParam.costs_quad.at(follower) +
         2 * Params.DemandParam.at(i).beta) *
        Params.scenProb.at(i);
    c(FollProdScen + i) = (Params.FollowerParam.costs_lin.at(follower) -
                           Params.DemandParam.at(i).alpha) *
                          Params.scenProb.at(i);
  }
  // C - the crossterms
  // Carbon price by government times carbon purchased
  C(FollCarb, Loc.at(LeaderVars::CarbPrice) - FollVarCount) = 1;
  // C - the crossterms

  for (unsigned int i = 0; i < num_scen; ++i) {
    // Quantity produced by other folowers times quantity produced by self terms
    // - for each scenario
    for (unsigned int j = 0; j < Params.n_followers - 1; ++j) {
      C(FollProdScen + i, FollVarCount * j + FollProdScen + i) =
          Params.DemandParam.at(i).beta * Params.scenProb.at(i);
    }
    // Country Import/Export terms in demand curve
    C(FollProdScen + i, Loc.at(LeaderVars::EnergImpScen) + i - FollVarCount) =
        Params.DemandParam.at(i).beta * Params.scenProb.at(i);
    C(FollProdScen + i, Loc.at(LeaderVars::EnergExpScen) + i - FollVarCount) =
        -Params.DemandParam.at(i).beta * Params.scenProb.at(i);
    // C(FollProdScen + i, Loc.at(LeaderVars::EnergImpScen) -this->LL_MC_count +
    // i - FollVarCount) = Params.DemandParam.at(i).beta *
    // Params.scenProb.at(i);
    // C(FollProdScen + i, Loc.at(LeaderVars::EnergExpScen) -this->LL_MC_count +
    // i - FollVarCount) = -Params.DemandParam.at(i).beta *
    // Params.scenProb.at(i);
    // 2nd Stage Domestic Carbon trade cost.
    C(FollCbuyScen + i, Loc.at(LeaderVars::Foll2Cpr) + i - FollVarCount) =
        Params.scenProb.at(i);
    C(FollCselScen + i, Loc.at(LeaderVars::Foll2Cpr) + i - FollVarCount) =
        -Params.scenProb.at(i);
    // C(FollCbuyScen + i, Loc.at(LeaderVars::Foll2Cpr) -this->LL_MC_count + i -
    // FollVarCount) = Params.scenProb.at(i);
    // C(FollCselScen + i, Loc.at(LeaderVars::Foll2Cpr) -this->LL_MC_count + i -
    // FollVarCount) = -Params.scenProb.at(i);
  }
  // OBJECTIVE Described

  // CONSTRAINTS
  for (unsigned int i = 0; i < num_scen; ++i) {
    // Infrastructural capacity constraint
    B(i, FollProdScen + i) = 1;
    b(i) = Params.FollowerParam.capacities.at(i).at(follower);
    // Carbon credit limit constraint
    B(num_scen + i, FollCarb) = -1;
    B(num_scen + i, FollProdScen + i) =
        Params.FollowerParam.emission_costs.at(follower);
    B(num_scen + i, FollCbuyScen + i) = -1;
    B(num_scen + i, FollCselScen + i) = 1;
  }
  // CONSTRAINTS Described

  Foll->set(Q, C, A, B, c, b);
}

template <unsigned int num_scen>
void cci::EPEC<num_scen>::make_LL_LeadCons(
    arma::sp_mat
        &LeadCons,      ///< The LHS matrix of leader constraints (for output)
    arma::vec &LeadRHS, ///< RHS vector for leader constraints (for output)
    const LeadAllPar<num_scen> &Params, ///< All country specific parameters
    const cci::LeadLocs &Loc,          ///< Location of variables
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
  for (unsigned int i = 0; i < num_scen; i++) {
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

template <unsigned int num_scen>
void cci::EPEC<num_scen>::make_MC_cons(arma::sp_mat &MCLHS,
                                        arma::vec &MCRHS) const
/** @brief Returns leader's Market clearing constraints in matrix form
 * @details
 */
{
  if (!this->finalized)
    throw std::string(
        "Error in cci::EPEC<num_scen>::make_MC_cons: This function can be "
        "run only AFTER calling finalize()");
  // Output matrices
  MCRHS.zeros(1 + num_scen);
  MCLHS.zeros(1 + num_scen, this->getnVarinEPEC());
  if (this->getNcountries() > 1) {
    for (unsigned int i = 0; i < this->getNcountries(); ++i) {
      // The MC constraint for each leader country's carbon trade market
      MCLHS(0, this->getPosition(i, LeaderVars::CarbExp)) = 1;
      MCLHS(0, this->getPosition(i, LeaderVars::CarbImp)) = -1;
      // The MC constraint for each leader country's energy trade market,
      // scenario-wise
      for (unsigned int j = 0; j < num_scen; ++j) {
        MCLHS(1 + j, this->getPosition(i, LeaderVars::EnergExpScen) + j) = 1;
        MCLHS(1 + j, this->getPosition(i, LeaderVars::EnergImpScen) + j) = -1;
      }
    }
  }
}

template <unsigned int num_scen>
bool cci::EPEC<num_scen>::dataCheck(
    const bool chkAllLeadPars, ///< Checks if cci::EPEC<num_scen>::AllLeadPars
                               ///< has size @p n
    const bool
        chkcountries_LL, ///< Checks if cci::EPEC<num_scen>::countries_LL has
    ///< size @p n
    const bool
        chkMC_QP, ///< Checks if cci::EPEC<num_scen>::MC_QP has size @p n
    const bool chkLeadConses, ///< Checks if cci::EPEC<num_scen>::LeadConses
                              ///< has size @p n
    const bool chkLeadRHSes,  ///< Checks if cci::EPEC<num_scen>::LeadRHSes has
                              ///< size @p n
    const bool chkLocations,  ///< Checks if cci::EPEC<num_scen>::Locations has
                              ///< size @p n
    const bool
        chkLeaderLocations, ///< Checks if cci::EPEC<num_scen>::LeaderLocations
                            ///< has
    ///< size @p n and cci::EPEC<num_scen>::nVarinEPEC is set
    const bool chkLeadObjec ///< Checks if cci::EPEC<num_scen>::LeadObjec has
                            ///< size @p n
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

template <unsigned int num_scen>
unsigned int cci::EPEC<num_scen>::getPosition(const unsigned int countryCount,
                                               const cci::LeaderVars var) const
/**
 * @brief Gets position of a variable in a country.
 */
{
  if (countryCount >= this->getNcountries())
    throw std::string(
        "Error in cci::EPEC<num_scen>::getPosition: Bad Country Count");
  return this->LeaderLocations.at(countryCount) +
         this->Locations.at(countryCount).at(var);
}

template <unsigned int num_scen>
unsigned int cci::EPEC<num_scen>::getPosition(const std::string countryName,
                                               const cci::LeaderVars var) const
/**
 * @brief Gets position of a variable in a country given the country name and
 * the variable.
 */
{
  return this->getPosition(name2nos.at(countryName), var);
}

template <unsigned int num_scen>
Game::NashGame *
cci::EPEC<num_scen>::get_LowerLevelNash(const unsigned int i) const
/**
 * @brief Returns a non-owning pointer to the @p i -th country's lower level
 * NashGame
 */
{
  return this->countries_LL.at(i).get();
}

template <unsigned int num_scen>
cci::EPEC<num_scen> &cci::EPEC<num_scen>::unlock()
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

template <unsigned int num_scen>
std::unique_ptr<GRBModel>
cci::EPEC<num_scen>::Respond(const std::string name,
                              const arma::vec &x) const {
  return this->Game::EPEC::Respond(this->name2nos.at(name), x);
}
//

template <unsigned int num_scen>
void cci::EPEC<num_scen>::write(const std::string filename,
                                 const unsigned int i, bool append) const {
  std::ofstream file;
  file.open(filename, append ? ios::app : ios::out);
  const LeadAllPar<num_scen> &Params = this->AllLeadPars.at(i);
  file << "**************************************************\n";
  file << "COUNTRY: " << Params.name << '\n';
  file << "- - - - - - - - - - - - - - - - - - - - - - - - - \n";
  file << Params;
  file << "**************************************************\n\n\n\n\n";
  file.close();
}

template <unsigned int num_scen>
void cci::EPEC<num_scen>::write(const std::string filename,
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

template <unsigned int num_scen>
void cci::EPEC<num_scen>::writeSolution(const int writeLevel,
                                         std::string filename) const {
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
    cerr
        << "Error in cci::EPEC<num_scen>::writeSolution: no solution to write."
        << '\n';
  }
}

template <unsigned int num_scen>
void cci::EPEC<num_scen>::WriteCountryMCprice(const std::string filename,
                                               const arma::vec x,
                                               const bool append) const {
  using cci::prn;
  std::ofstream file;
  file.open(filename, append ? ios::app : ios::out);
  file << "**************************************************\n";
  file << "INTERNATIONAL PRICES\n";
  file << "**************************************************\n";
  file << prn::label << "Carbon Price: " << prn::val
       << x.at(this->getPosition(this->getNcountries() - 1, LeaderVars::End));
  const LeadAllPar<num_scen> &Params = this->AllLeadPars.at(0);
  double Price{0};
  for (unsigned int scen = 0; scen < num_scen; scen++) {
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

template <unsigned int num_scen>
void cci::EPEC<num_scen>::WriteCountry(const unsigned int i,
                                        const std::string filename,
                                        const arma::vec x,
                                        const bool append) const {
  // if (!lcp) return;
  // const LeadLocs& Loc = this->Locations.at(i);
  using cci::prn;
  std::ofstream file;
  file.open(filename, append ? ios::app : ios::out);
  // FILE OPERATIONS START
  const LeadAllPar<num_scen> &Params = this->AllLeadPars.at(i);
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
  for (unsigned int scen = 0; scen < num_scen; scen++) {
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
template <unsigned int num_scen>
void cci::EPEC<num_scen>::WriteFollower(const unsigned int i,
                                         const unsigned int j,
                                         const std::string filename,
                                         const arma::vec x) const {
  using cci::prn;
  std::ofstream file;
  file.open(filename, ios::app);

  const LeadAllPar<num_scen> &Params = this->AllLeadPars.at(i);
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
  for (unsigned int scen = 0; scen < num_scen; scen++) {
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

template <unsigned int num_scen>
bool cci::EPEC<num_scen>::ParamValid(
    const LeadAllPar<num_scen>
        &Params ///< Object whose validity is to be tested
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
  for (unsigned int i = 0; i < num_scen; ++i) {
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

template <unsigned int num_scen>
std::ostream &cci::operator<<(std::ostream &ost,
                               const cci::LeadAllPar<num_scen> P) {
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
  for (unsigned int i = 0; i < num_scen; ++i) {
    ost << "Scenario " << i + 1
        << " happens with probability: " << P.scenProb.at(i) << '\n';
  }
  ost << "\n";
  // ost 	<< P.DemandParam << "\n";
  ost << "\n";
  return ost;
}

template <unsigned int num_scen>
std::ostream &cci::operator<<(std::ostream &ost,
                               const cci::FollPar<num_scen> P) {
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
  // for (unsigned int i = 0; i < num_scen; i++) {
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

template <int num_scen>
std::ostream &operator<<(std::ostream &ost,
                         const std::array<cci::DemPar, num_scen> P) {
  ost << "Demand Parameters: " << '\n';
  ost << "******************" << '\n';
  for (unsigned int i = 0; i < num_scen; ++i)
    ost << "\tScenario " << i + 1 << ": P\t\t =\t\t " << P[i].alpha << "\t-\t"
        << P[i].beta << cci::prn::label << " x Q\n";
  return ost;
}

template <unsigned int num_scen>
void cci::EPEC<num_scen>::initializeSoln(arma::vec &init_x) const {
  double chiToArg_Carb{20};

  init_x.at(this->getPosition(0, cci::LeaderVars::CarbExp)) = chiToArg_Carb;
  init_x.at(this->getPosition(1, cci::LeaderVars::CarbExp)) = -chiToArg_Carb;
  init_x.at(this->getPosition(0, cci::LeaderVars::CarbImp)) = -chiToArg_Carb;
  init_x.at(this->getPosition(1, cci::LeaderVars::CarbImp)) = chiToArg_Carb;

  init_x.at(this->getPosition(1, cci::LeaderVars::End)) = -2000;
  init_x.at(this->getPosition(1, cci::LeaderVars::End) + 1) = -1000;

  // init_x.at(this->getPosition(0,
}
