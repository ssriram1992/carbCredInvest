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

#include "bones/dataLoad.h"
#include "bones/helpers.h"
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
  static constexpr double social_cost_of_carbon{0};

private:
  const std::vector<std::string> dirtyEnergy;
  const std::vector<std::string> cleanEnergy;

  const commonData comDat;
  /// cci::energy should contain everything in dirtyEnergy and cleanEnergy
  std::vector<std::string> energy;

  mutable std::map<std::string, double> dataMap{};

  static constexpr unsigned int LL_MC_count =
      0; // No market clearing for follower

  static constexpr unsigned int FollInv{0};
  static constexpr unsigned int FollProdDirty{n_Clean};
  static constexpr unsigned int FollProdClean{FollProdDirty + n_Dirty * n_Scen};
  static constexpr unsigned int FollCarbBuy{
      FollProdClean +
      n_Clean * n_Scen}; // =  n_Clean +(n_Clean + n_Dirty) * n_Scen
  static constexpr unsigned int FollEnd{FollCarbBuy + 1};

  static constexpr unsigned int FollVarCount{
      FollEnd}; // n_Clean +(n_Clean + n_Dirty) * n_Scen + 1

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
                        const cci::LeadLocs &Loc = {}) const noexcept;

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
       const std::vector<std::string> cleanEnergy, const commonData cD)
      : Game::EPEC(env), dirtyEnergy{dirtyEnergy},
        cleanEnergy{cleanEnergy}, comDat{cD} {
    if (dirtyEnergy.size() != n_Dirty || cleanEnergy.size() != n_Clean) {
      throw std::string(
          "Error in EPEC(): Illegal size of dirtyEnergy/cleanEnergy");
    }
    this->energy.clear();
    this->energy.insert(energy.end(), dirtyEnergy.cbegin(), dirtyEnergy.cend());
    this->energy.insert(energy.end(), cleanEnergy.cbegin(), cleanEnergy.cend());
  }

  ///@brief add a Standard Nash-Cournot game within a country
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

  void appendSolution4XL(const std::string filename, int problemID,
                         bool append = true) const;

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

#include "meat/Helpers.cpp"
#include "meat/MainModel.cpp"
#include "meat/Write.cpp"

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

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::initializeSoln(
    arma::vec &init_x) const {
  double chiToArg_Carb{20};

  // init_x.at(this->getPosition(0, cci::LeaderVars::CarbExp)) = chiToArg_Carb;
  // init_x.at(this->getPosition(1, cci::LeaderVars::CarbExp)) = -chiToArg_Carb;
  init_x.at(this->getPosition(0, cci::LeaderVars::CarbImp)) = -chiToArg_Carb;
  init_x.at(this->getPosition(1, cci::LeaderVars::CarbImp)) = chiToArg_Carb;

  // init_x.at(this->getPosition(1, cci::LeaderVars::End)) = -2000;
  // init_x.at(this->getPosition(1, cci::LeaderVars::End) + 1) = -1000;

  // init_x.at(this->getPosition(0,
}
