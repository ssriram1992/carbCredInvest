#pragma once
#include "../bones/dataLoad.h"
#include "../carbCredInv.h"

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
      if (cc == i) continue;
      QP_obj.C(Loc.at(LeaderVars::TotInv) + ii,
               this->getPosition(cc, cci::LeaderVars::TotInv) + ii -
                   (cc > i ? nThisCountryvars : 0)) =
          Params.LeaderParam.cleanInvCrossVal.at(cleanEnergy.at(ii));
    }
  }

  // Emissions local
  QP_obj.Q(Loc.at(LeaderVars::TotEmission), Loc.at(LeaderVars::TotEmission)) =
      Params.LeaderParam.emissionValQuad;
  QP_obj.c(Loc.at(LeaderVars::TotEmission)) = Params.LeaderParam.emissionVal;
  // Emission cross terms
  for (unsigned int cc = 0; cc < this->getNcountries(); ++cc) {
    if (cc == i) continue;
    QP_obj.C(Loc.at(LeaderVars::TotEmission),
             this->getPosition(cc, cci::LeaderVars::TotEmission) -
                 (cc > i ? nThisCountryvars : 0)) =
        Params.LeaderParam.emissionCrossVal;
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

  // Q term for leader's income from selling carbon credits to followers
  // const unsigned int carbPricePos = Loc.at(LeaderVars::CarbPrice);
  // for (unsigned int ff = 0; ff < Params.n_followers; ff++) {
  //   const unsigned int ffCarbBuy = FollVarCount * ff + FollCarbBuy;
  //   const unsigned int ffCarbSel = FollVarCount * ff + FollCarbSel;
  //   QP_obj.Q.at(carbPricePos, ffCarbBuy) = 1;
  //   QP_obj.Q.at(ffCarbBuy, carbPricePos) = 1;
  //   QP_obj.Q.at(carbPricePos, ffCarbSel) = -1;
  //   QP_obj.Q.at(ffCarbSel, carbPricePos) = -1;
  // }

  // Nonconvex term
  QP_obj.c.at(Loc.at(cci::LeaderVars::NonConv)) = 1;

  // // Social cost of carbon
  // QP_obj.c.at(Loc.at(LeaderVars::CarbExp)) = social_cost_of_carbon;
  // QP_obj.c.at(Loc.at(LeaderVars::CarbImp)) = -social_cost_of_carbon;

  // for (unsigned int ff = 0; ff < Params.n_followers; ff++) {
  //   QP_obj.c.at(FollVarCount * ff + FollCarbBuy) = social_cost_of_carbon;
  //   QP_obj.c.at(FollVarCount * ff + FollCarbSel) = -social_cost_of_carbon;
  // }
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
  cci::increaseVal(Loc, LeaderVars::CarbBuy, 1);
  cci::increaseVal(Loc, LeaderVars::NonConv, 1);
  cci::increaseVal(Loc, LeaderVars::TotInv, n_Clean);
  cci::increaseVal(Loc, LeaderVars::TotEmission, 1);

  const unsigned int LeadVars =
      Loc.at(cci::LeaderVars::End) - this->FollVarCount * Params.n_followers;
  // Should be 6+n_Clean - one for CarbExp, one for CarbImp, one for CarbPrice,
  // one for CarbBuy, one for NonConv and n_Clean each for TotInv and
  // TotEmission

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
                            2 * n_Clean +     // Investment summing
                            2 +               // Emission summing
                            2 +               // CarbBuy summing
                            1 +               // CarbPrice Up Bnd
                            1 +               // McCormick
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
  for (unsigned int ii = 0; ii < n_Dirty; ii++) {
    for (unsigned int jj = 0; jj < n_Clean; jj++) {
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
            return FollProdClean + n_Clean * scen + (i_nos-n_Dirty);
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
      // for (const auto &dictkv:infCap)
      // {
      // BOOST_LOG_TRIVIAL(trace) <<" Dict val:
      // "<<dictkv.first<<"---"<<dictkv.second;
      // }
      b(constrCount) = infCap.at(dirtyEnergy.at(ii));
      // Next constraint
      constrCount++;
    }
    for (unsigned int ii = 0; ii < n_Clean; ++ii) {
      const auto energ = cleanEnergy.at(ii);
      B(constrCount, FollProdClean + scen * n_Clean + ii) = 1;
      B(constrCount, FollInv + ii) = -1 * renCapAdj.at(scen).at(energ);
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
  BOOST_LOG_TRIVIAL(trace) << "Inside cci::EPEC::make_LL_QP:end";
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
    constrCount += 2;
  }
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
        LeadCons(constrCount, FollVarCount * ff + FollProdClean +
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
  LeadCons(constrCount, Loc.at(LeaderVars::CarbBuy)) = -1;

  LeadRHS(constrCount) = Params.LeaderParam.carbCreditInit;
  constrCount++;
  // CarbBuy summing
  LeadCons(constrCount, Loc.at(LeaderVars::CarbBuy)) = 1;
  LeadCons(constrCount + 1, Loc.at(LeaderVars::CarbBuy)) = -1;
  for (unsigned int ff = 0; ff < Params.n_followers; ff++) {
    LeadCons(constrCount, FollVarCount * ff + FollCarbBuy) = -1;
    LeadCons(constrCount, FollVarCount * ff + FollCarbSel) = 1;
    LeadCons(constrCount + 1, FollVarCount * ff + FollCarbBuy) = 1;
    LeadCons(constrCount + 1, FollVarCount * ff + FollCarbSel) = -1;
  }
  LeadRHS(constrCount) = 0;
  LeadRHS(constrCount + 1) = 0;
  constrCount += 2;
  // Upper bound for leader Carbon Price
  LeadCons(constrCount, Loc.at(LeaderVars::CarbPrice)) = 1;
  LeadRHS(constrCount) = Params.LeaderParam.maxCarbPrice;
  // McCormick - NonConvex >= maxCarbPrice * CarbBuy + maxCarbBuy * carbPrice -
  // maxCarbBuy*maxCarbPrice
  const double maxCarbBuy = [&Params]() {
    double ans;
    for (unsigned int ff = 0; ff < Params.n_followers; ff++)
      ans += Params.FollowerParam.at(ff).carbonCreditInit;
    if (Params.LeaderParam.carbCreditInit)
      BOOST_LOG_TRIVIAL(fatal)
          << "Error in cci::EPEC::make_LL_LeadCons: Leader has non-zero "
             "carbon-credits to start with!";
    return ans;
  }();
  LeadCons(constrCount, Loc.at(LeaderVars::NonConv)) = -1;
  LeadCons(constrCount, Loc.at(LeaderVars::CarbBuy)) =
      Params.LeaderParam.maxCarbPrice;
  LeadCons(constrCount, Loc.at(LeaderVars::CarbPrice)) = maxCarbBuy;
  LeadRHS(constrCount) = maxCarbBuy * Params.LeaderParam.maxCarbPrice;
  constrCount++;
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
