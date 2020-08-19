#pragma once
#include "../bones/dataLoad.h"
#include "../carbCredInv.h"
#include <algorithm>

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

  BOOST_LOG_TRIVIAL(trace) << "In cci::EPEC::make_obj_leader";
  QP_obj.Q.zeros(nThisCountryvars, nThisCountryvars);
  QP_obj.c.zeros(nThisCountryvars);
  QP_obj.C.zeros(nThisCountryvars, nEPECvars - nThisCountryvars);

  // Total investment locally
  for (unsigned int ii = 0; ii < n_Clean; ++ii) {
    QP_obj.c.at(Loc.at(LeaderVars::TotInv) + ii) =
        -1 * Params.LeaderParam.cleanInvVal.at(cleanEnergy.at(ii));
    // Negative sign because investment has to be maximized
  }
  BOOST_LOG_TRIVIAL(trace)
      << "In cci::EPEC::make_obj_leader: Local Investments over";

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
  BOOST_LOG_TRIVIAL(trace)
      << "In cci::EPEC::make_obj_leader: Investments cross over";
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
  BOOST_LOG_TRIVIAL(trace) << "In cci::EPEC::make_obj_leader: Emission over";
  // Production value
  for (unsigned int ff = 0; ff < Params.n_followers; ++ff) {
    for (unsigned int scen = 0; scen < n_Scen; ++scen) {
      const auto probab =
          Params.scenProb.at(scen); // Probability of the current scenario
      for (unsigned int ii = 0; ii < n_Dirty; ++ii)
        QP_obj.c(FollVarCount * ff + FollProdDirty + n_Dirty * scen + ii) +=
            -1 * probab * Params.LeaderParam.prodnVal;
      for (unsigned int ii = 0; ii < n_Clean; ++ii)
        QP_obj.c(FollVarCount * ff + FollProdClean + n_Clean * scen + ii) +=
            -1 * probab * Params.LeaderParam.prodnVal;
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "In cci::EPEC::make_obj_leader: Production over";
  // Carbon credit trade term
  for (unsigned int ff = 0; ff < Params.n_followers; ++ff) {
    QP_obj.c(FollVarCount * ff + FollCarbBuy) = Params.LeaderParam.taxCarbon;
  }
  QP_obj.c(Loc.at(LeaderVars::CarbImp)) = this->comDat.suppInt;
  QP_obj.Q(Loc.at(LeaderVars::CarbImp), Loc.at(LeaderVars::CarbImp)) =
      this->comDat.suppSlope * 2;
  for (unsigned int cc = 0; cc < this->getNcountries(); ++cc) {
    if (cc == i)
      continue;
    QP_obj.C(Loc.at(LeaderVars::CarbImp),
             this->getPosition(cc, cci::LeaderVars::CarbImp) -
                 (cc > i ? nThisCountryvars : 0)) = this->comDat.suppSlope;
  }
  BOOST_LOG_TRIVIAL(trace) << "In cci::EPEC::make_obj_leader: Carbon Tradeover";
  QP_obj.Q = QP_obj.Q / 10000;
  QP_obj.C = QP_obj.C / 10000;
  QP_obj.c = QP_obj.c / 10000;
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
  cci::increaseVal(Loc, LeaderVars::CarbImp, 1);
  cci::increaseVal(Loc, LeaderVars::CarbTax, 1);
  // cci::increaseVal(Loc, LeaderVars::CarbFollLim, Params.n_followers);
  cci::increaseVal(Loc, LeaderVars::TotInv, n_Clean);
  cci::increaseVal(Loc, LeaderVars::TotEmission, 1);

  const unsigned int LeadVars =
      Loc.at(cci::LeaderVars::End) - this->FollVarCount * Params.n_followers;
  // Should be 3+n_Clean - one for CarbImp,
  // one for CarbSell and n_Clean each for TotInv and
  // TotEmission

  // Leader Constraints
  arma::sp_mat LeadCons(
      (Params.LeaderParam.consum_limit > 0 ? 1 : 0) +    // Minimum consumption
          n_Clean +                                      // Investment summing
          2 +                                            // Emission summing
          (Params.LeaderParam.taxCarbon < 0 ? 0 : 2) +   // Fixing the tax?
          (Params.LeaderParam.expCleanMix > 0 ? 1 : 0) + // Clean mix constraint
          1,                                             // Carbon credits >=0
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
    this->make_LL_LeadCons(LeadCons, LeadRHS, Params, Loc);
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
  this->Game::EPEC::n_MCVar = 0;
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

  const auto &Follparam = Params.FollowerParam.at(follower);

  arma::sp_mat Q(FollVarCount, FollVarCount);
  arma::sp_mat C(FollVarCount, otherVars);
  // 1 for generalized Nash
  // n_Dirty*n_Scen for dirty energy
  // infrastructural limits
  // n_Clean*n_Scen for clean energy
  // infrastructural limits
  // n_Scen for emission limits
  // 1 for Carbon Credit Fractional Limit
  const unsigned int n_Const =
      (Params.LeaderParam.follGenNash ? 1 : 0) + // Generalized Nash constraint
      (n_Dirty + n_Clean) * n_Scen + n_Scen +
			(Params.LeaderParam.rps > 0? n_Scen:0) +		// Renewable portfolio standard constr
      // (Follparam.carbonLimitFrac >= 1 ? 0 : 1) + // CarbFollLim constraint
      1;

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

  // c - the crossterms
  // Carbon tax by government times carbon purchased
  // c(FollCarbBuy) = Params.LeaderParam.taxCarbon;
  C(FollCarbBuy, Loc.at(LeaderVars::CarbTax) - FollVarCount) = 1;
  // C(FollCarbSel, Loc.at(LeaderVars::CarbTax) - FollVarCount) = -1;

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

            return FollProdClean + n_Clean * scen + (i_nos - n_Dirty);

        };
        const auto pos1 = pos(ii);
        const auto pos2 = pos(jj);
        for (unsigned int kk = 0; kk < n_Foll - 1; ++kk) {
          C(pos1, kk * FollVarCount + pos2) =
              Params.scenProb.at(scen) * (Params.DemandParam.at(scen).second);
        }
      }
    }
  }

  // OBJECTIVE Described
  // CONSTRAINTS
  const auto &infCap = Follparam.capacities;
  const auto &renCapAdj = Follparam.renewCapAdjust;
  unsigned int constrCount{0};
  if (Params.LeaderParam.follGenNash) {
    //  General Nash constraint
    B(constrCount, FollCarbBuy) = 1;
    A(constrCount, Loc.at(cci::LeaderVars::CarbImp) - FollVarCount) = -1;
    for (unsigned int ff = 0; ff < Params.n_followers - 1; ++ff)

      A(constrCount, ff * FollVarCount + FollCarbBuy) =
          1; // Summing up Carb buy of all other followers
    b(constrCount) = Params.LeaderParam.carbCreditInit;
    constrCount++;
  }

  for (unsigned int scen = 0; scen < n_Scen; ++scen) {
    // RPS constraints
    if (Params.LeaderParam.rps > 0){
      const double rps = Params.LeaderParam.rps;
      for (unsigned int ii = 0; ii < n_Dirty; ++ii) {
        B(constrCount, FollProdDirty + n_Dirty*scen +ii) = rps;
      }
      for (unsigned int ii = 0; ii < n_Clean; ++ii) {
        B(constrCount, FollProdClean + n_Clean*scen +ii) = rps-1;
      }
      constrCount++;
    }
  
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
    // Carbon credits emission limit Constraint
    B(constrCount, FollCarbBuy) = -1;
    for (unsigned int ii = 0; ii < n_Dirty; ++ii)
      B(constrCount, FollProdDirty + scen * n_Dirty + ii) =
          Follparam.emissionCosts.at(dirtyEnergy.at(ii));
    for (unsigned int ii = 0; ii < n_Clean; ++ii)
      B(constrCount, FollProdClean + scen * n_Clean + ii) =
          Follparam.emissionCosts.at(cleanEnergy.at(ii));
    b(constrCount) = Follparam.carbonCreditInit;
    constrCount++;
  }

  // carbonLimitFrac constraint
  B(constrCount, FollCarbBuy) = 1;
  A(constrCount, Loc.at(cci::LeaderVars::CarbImp) - FollVarCount) =
      -Follparam.carbonLimitFrac;
  constrCount++;
  // CONSTRAINTS Described
  if (constrCount != n_Const)
    BOOST_LOG_TRIVIAL(error)
        << "Inside cci::EPEC::make_LL_QP: constrCount != n_Const: "
        << constrCount << " !=" << n_Const;

  Foll->set(Q, C, A, B, c, b);
  BOOST_LOG_TRIVIAL(trace) << "Inside cci::EPEC::make_LL_QP:end";
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPEC<n_Dirty, n_Clean, n_Scen>::make_LL_LeadCons(
    arma::sp_mat
        &LeadCons,      ///< The LHS matrix of leader constraints (for output)
    arma::vec &LeadRHS, ///< RHS vector for leader constraints (for output)
    const LeadAllPar<n_Scen> &Params, ///< All country specific parameters
    const cci::LeadLocs &Loc) const noexcept
/**
 * Makes the leader level constraints for a country.
 */
{
  // LeadCons and LeadRHS already have the right size
  LeadCons.zeros();
  LeadRHS.zeros();
  BOOST_LOG_TRIVIAL(trace) << "make_LL_LeadCons: Nos: constraints: "
                           << LeadRHS.n_rows;
  unsigned int constrCount{0};
  // consum_limit constraint
  if (Params.LeaderParam.consum_limit > 0) {
    for (unsigned int scen = 0; scen < n_Scen; ++scen) {
      const auto probab =
          Params.scenProb.at(scen); // Probability of the current scenario
      for (unsigned int ff = 0; ff < Params.n_followers; ++ff) {
        // Dirty Producers
        for (unsigned int ii = 0; ii < n_Dirty; ++ii) {
          LeadCons(constrCount, FollVarCount * ff + FollProdDirty +
                                    n_Dirty * scen + ii) = -probab;
        }
        // Clean Producers
        for (unsigned int ii = 0; ii < n_Clean; ++ii) {
          LeadCons(constrCount, FollVarCount * ff + FollProdClean +
                                    n_Clean * scen + ii) = -probab;
        }
      }
    }
    LeadRHS(constrCount) = -Params.LeaderParam.consum_limit;
    constrCount++;
  }
  BOOST_LOG_TRIVIAL(trace) << "make_LL_LeadCons: Consum limit over: "
                           << constrCount << " "
                           << Params.LeaderParam.consum_limit;

  // Investment summing
  for (unsigned int ii = 0; ii < n_Clean; ++ii) {
    for (unsigned int ff = 0; ff < Params.n_followers; ++ff) {
      LeadCons(constrCount, FollVarCount * ff + FollInv + ii) = -1;
    }

    LeadCons(constrCount, Loc.at(LeaderVars::TotInv) + ii) = 1;
    constrCount += 1;
  }

  BOOST_LOG_TRIVIAL(trace) << "make_LL_LeadCons: Investment summing over: "
                           << constrCount;

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
        // BOOST_LOG_TRIVIAL(trace) << "make_LL_LeadCons: Clean:" <<
        // FollVarCount * ff + FollProdClean +
        //                           n_Clean * scen + ii;
      }
    }
  }
  LeadCons(constrCount, Loc.at(LeaderVars::TotEmission)) = -1;
  LeadCons(constrCount + 1, Loc.at(LeaderVars::TotEmission)) = 1;

  constrCount += 2;
  BOOST_LOG_TRIVIAL(trace) << "make_LL_LeadCons: Exp Emission summing over: "
                           << constrCount;
  // CarbBuy summing
  // LeadCons(constrCount, Loc.at(LeaderVars::CarbBuy)) = 1;
  LeadCons(constrCount, Loc.at(LeaderVars::CarbImp)) = -1;
  for (unsigned int ff = 0; ff < Params.n_followers; ff++)
    LeadCons(constrCount, FollVarCount * ff + FollCarbBuy) = 1;
  LeadRHS(constrCount) = Params.LeaderParam.carbCreditInit;
  constrCount += 1;
  BOOST_LOG_TRIVIAL(trace) << "make_LL_LeadCons: CarbBuy over: " << constrCount;

  // Tax Constraints
  const double taxAmt = Params.LeaderParam.taxCarbon;
  if (taxAmt >= 0) {
    LeadCons(constrCount, Loc.at(LeaderVars::CarbTax)) = 1;
    LeadRHS(constrCount) = taxAmt;
    LeadCons(constrCount + 1, Loc.at(LeaderVars::CarbTax)) = -1;
    LeadRHS(constrCount + 1) = -taxAmt;
    constrCount += 2;
  }
  BOOST_LOG_TRIVIAL(trace) << "make_LL_LeadCons: Tax over: " << constrCount;
  // clean mix constraints.
  if (Params.LeaderParam.expCleanMix > 0) {
    const double fraction = Params.LeaderParam.expCleanMix;
    for (unsigned int scen = 0; scen < n_Scen; ++scen) {
      const auto probab =
          Params.scenProb.at(scen); // Probability of the current scenario

      for (unsigned int ff = 0; ff < Params.n_followers; ++ff) {
        // Dirty Producers
        for (unsigned int ii = 0; ii < n_Dirty; ++ii) {
          LeadCons(constrCount, FollVarCount * ff + FollProdDirty +
                                    n_Dirty * scen + ii) = probab * fraction;
        }
        // Clean Producers
        for (unsigned int ii = 0; ii < n_Clean; ++ii) {
          LeadCons(constrCount,
                   FollVarCount * ff + FollProdClean + n_Clean * scen + ii) =
              probab * (fraction - 1);
        }
      }
    }
    constrCount++;
  }
  BOOST_LOG_TRIVIAL(trace) << "make_LL_LeadCons: cleanMix over: " << constrCount
                           << " " << Params.LeaderParam.expCleanMix;

  if (constrCount != LeadCons.n_rows)
    BOOST_LOG_TRIVIAL(error)
        << "Inside cci::EPEC::make_LL_LeadCons: constrCount != n_Const: "
        << constrCount << " !=" << LeadCons.n_rows;
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
  MCRHS.zeros(0);
  MCLHS.zeros(0, this->getnVarinEPEC());
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
GRBQuadExpr
cci::EPEC<n_Dirty, n_Clean, n_Scen>::make_lcp_objective(GRBModel *model) {
  GRBQuadExpr expr{0};
  for (unsigned int cc = 0; cc < this->getNcountries(); ++cc) {
    unsigned int posn = this->getPosition(cc, LeaderVars::CarbImp);
    expr += model->getVarByName("x_" + std::to_string(posn));
  }
  return expr;
}
