#include "carbCredInv.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <gurobi_c++.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

using namespace std;
void cci::init(LeadLocs &L) {
  for (LeaderVars l = cci::LeaderVars::Followers; l != cci::LeaderVars::End;
       l = l + 1)
    L[l] = 0;
  L[cci::LeaderVars::End] = 0;
};

void cci::increaseVal(LeadLocs &L, const LeaderVars start,
                      const unsigned int val, const bool startnext)
/**
 * Should be called ONLY after initializing @p L by calling cci::init
 */
{
  LeaderVars start_rl = startnext ? start + 1 : start;
  for (LeaderVars l = start_rl; l != cci::LeaderVars::End; l = l + 1)
    L[l] += val;
  L[cci::LeaderVars::End] += val;
  // BOOST_LOG_TRIVIAL(error)<<"End location changed to:
  // "<<L[cci::LeaderVars::End];
}

void cci::decreaseVal(LeadLocs &L, const LeaderVars start,
                      const unsigned int val, const bool startnext)
/**
 * Should be called ONLY after initializing @p L by calling cci::init
 */
{
  LeaderVars start_rl = startnext ? start + 1 : start;
  for (LeaderVars l = start_rl; l != cci::LeaderVars::End; l = l + 1)
    L[l] -= val;
  L[cci::LeaderVars::End] -= val;
  // BOOST_LOG_TRIVIAL(error)<<"End location changed to:
  // "<<L[cci::LeaderVars::End];
}

cci::LeaderVars cci::operator+(cci::LeaderVars a, int b) {
  return static_cast<LeaderVars>(static_cast<int>(a) + b);
}

string to_string(const GRBConstr &cons, const GRBModel &model) {
  const GRBVar *vars = model.getVars();
  const int nVars = model.get(GRB_IntAttr_NumVars);
  ostringstream oss;
  oss << cons.get(GRB_StringAttr_ConstrName) << ":\t\t";
  constexpr double eps = 1e-5;
  // LHS
  for (int i = 0; i < nVars; ++i) {
    double coeff = model.getCoeff(cons, vars[i]);
    if (abs(coeff) > eps) {
      char sign = (coeff > eps) ? '+' : ' ';
      oss << sign << coeff << to_string(vars[i]) << "\t";
    }
  }
  // Inequality/Equality and RHS
  oss << cons.get(GRB_CharAttr_Sense) << "\t" << cons.get(GRB_DoubleAttr_RHS);
  return oss.str();
}

string to_string(const GRBVar &var) {
  string name = var.get(GRB_StringAttr_VarName);
  return name.empty() ? "unNamedvar" : name;
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
void cci::EPECInstance<n_Dirty, n_Clean, n_Scen>::save(string filename) {
  /**
   * @brief Writes the current EPEC instance to the standard JSON instance file
   * @p filename dictates the name of the JSON instance file
   * @p epec contains the @p EPECInstance object with the data
   */
  std::cout << filename;
}

ostream &cci::operator<<(ostream &ost, const cci::prn ll) {
  switch (ll) {
  case cci::prn::label:
    ost << std::left << std::setw(50);
    break;
  case cci::prn::val:
    ost << std::right << std::setprecision(2) << std::setw(16) << std::fixed;
    break;
  default:
    break;
  }
  return ost;
}

ostream &cci::operator<<(ostream &ost, const std::pair<double, double> P) {
  ost << "Demand Parameters: " << '\n';
  ost << "******************" << '\n';
  ost << "Price\t\t =\t\t " << P.first << "\t-\t" << P.second
      << "  x   Quantity" << '\n';
  return ost;
}

ostream &cci::operator<<(ostream &ost, const cci::LeadPar P) {
  ost << "Leader Parameters: " << '\n';
  ost << "******************" << '\n';
  ost << std::fixed;
  ost << cci::prn::label << "Initial carbon credit allotted"
      << ":" << cci::prn::val << P.carbCreditInit;
  ost << '\n';
  ost << '\n';
  ost << cci::prn::label << "Min reqd consumption"
      << ":" << cci::prn::val
      << (P.consum_limit < 0 ? std::numeric_limits<double>::infinity()
                             : P.consum_limit)
      << '\n';
  ost << cci::prn::label
      << "Weights for domestic energy production: " << cci::prn::val
      << P.prodnVal << '\n';
  ost << cci::prn::label << "Carbon Tax: " << cci::prn::val << P.taxCarbon
      << '\n';
  ost << cci::prn::label
      << "Emission penalties (linear, quadratic): " << cci::prn::val
      << P.emissionVal << ", " << P.emissionValQuad;

  ost << '\n';
  ost << "Investment Incentives: \n";
  for (const auto &cc : P.cleanInvVal)
    ost << cci::prn::label << "Investment " << cc.first << cci::prn::val
        << cc.second << '\n';

  return ost;
}

template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int num_scen>
ostream &
cci::operator<<(ostream &ost,
                const cci::EPECInstance<n_Dirty, n_Clean, num_scen> I) {
  ost << "EPEC Instance: " << '\n';
  ost << "******************" << '\n';
  for (auto a : I.Countries)
    ost << a << '\n';
  return ost;
}

ostream &cci::operator<<(ostream &ost, const cci::LeaderVars l) {
  switch (l) {
  case cci::LeaderVars::Followers:
    ost << "cci::LeaderVars::Followers";
    break;
  case cci::LeaderVars::CarbImp:
    ost << "cci::LeaderVars::CarbImp";
    break;
  case cci::LeaderVars::CarbTax:
    ost << "cci::LeaderVars::CarbTax";
    break;
  case cci::LeaderVars::TotInv:
    ost << "cci::LeaderVars::TotInv";
    break;
  case cci::LeaderVars::TotEmission:
    ost << "cci::LeaderVars::TotEmission";
    break;
  case cci::LeaderVars::DualVar:
    ost << "cci::LeaderVars::DualVar";
    break;
  case cci::LeaderVars::ConvHullDummy:
    ost << "cci::LeaderVars::ConvHullDummy";
    break;
  case cci::LeaderVars::End:
    ost << "cci::LeaderVars::End";
    break;
  default:
    cerr << "Incorrect argument to ostream& operator<<(ostream& ost, const "
            "LeaderVars l)";
  };
  return ost;
}

template <unsigned int num_scen>
cci::FollPar<num_scen> operator+(const cci::FollPar<num_scen> &F1,
                                 const cci::FollPar<num_scen> &F2) {
  std::vector<double> cq, cl, ec;
  std::array<std::vector<double>, num_scen> cap;
  std::vector<std::string> nm;

  cq.insert(cq.end(), F1.costs_quad.begin(), F1.costs_quad.end());
  cq.insert(cq.end(), F2.costs_quad.begin(), F2.costs_quad.end());

  cl.insert(cl.end(), F1.costs_lin.begin(), F1.costs_lin.end());
  cl.insert(cl.end(), F2.costs_lin.begin(), F2.costs_lin.end());

  for (unsigned int i = 0; i < num_scen; i++) {
    cap[i].insert(cap[i].end(), F1[i].capacities.begin(),
                  F1[i].capacities.end());
    cap[i].insert(cap[i].end(), F2[i].capacities.begin(),
                  F2[i].capacities.end());
  }

  ec.insert(ec.end(), F1.emission_costs.begin(), F1.emission_costs.end());
  ec.insert(ec.end(), F2.emission_costs.begin(), F2.emission_costs.end());

  nm.insert(nm.end(), F1.names.begin(), F1.names.end());
  nm.insert(nm.end(), F2.names.begin(), F2.names.end());

  return cci::FollPar<num_scen>(cq, cl, cap, ec, nm);
}

void cci::writeProbData(const std::string filename, const cci::probData &pd) {
  std::ofstream ff;
  ff.open(filename, std::ios::out);
  ff << "prefix  " << pd.prefix << '\n';
  ff << "probId  " << pd.probId << '\n';
  ff << "suIn  " << pd.suIn << '\n';
  ff << "suSl  " << pd.suSl << '\n';
  ff << "c1Pp  " << pd.c1Pp << '\n';
  ff << "c2Pp  " << pd.c2Pp << '\n';
  ff << "c1CreditSplit  " << pd.c1CreditSplit << '\n';
  ff << "c2CreditSplit  " << pd.c2CreditSplit << '\n';
  ff << "c1CleanMix  " << pd.c1CleanMix << '\n';
  ff << "c2CleanMix  " << pd.c2CleanMix << '\n';
  ff << "c1tax  " << pd.c1tax << '\n';
  ff << "c2tax  " << pd.c2tax << '\n';
  ff << "prodCoalLin  " << pd.prodCoalLin << '\n';
  ff << "prodGasLin  " << pd.prodGasLin << '\n';
  ff << "prodSolLin  " << pd.prodSolLin << '\n';
  ff << "prodWindLin  " << pd.prodWindLin << '\n';
  ff << "prodCoalQuad  " << pd.prodCoalQuad << '\n';
  ff << "prodGasQuad  " << pd.prodGasQuad << '\n';
  ff << "prodSolQuad  " << pd.prodSolQuad << '\n';
  ff << "prodWindQuad  " << pd.prodWindQuad << '\n';
  ff << "invSolLin  " << pd.invSolLin << '\n';
  ff << "invWindLin  " << pd.invWindLin << '\n';
  ff << "invSolQuad  " << pd.invSolQuad << '\n';
  ff << "invWindQuad  " << pd.invWindQuad << '\n';
  ff << "emitCoal  " << pd.emitCoal << '\n';
  ff << "emitGas  " << pd.emitGas << '\n';
  ff << "emitSol  " << pd.emitSol << '\n';
  ff << "emitWind  " << pd.emitWind << '\n';
  ff << "c1f1CapCoal  " << pd.c1f1CapCoal << '\n';
  ff << "c1f1CapGas  " << pd.c1f1CapGas << '\n';
  ff << "c1f1CapSol  " << pd.c1f1CapSol << '\n';
  ff << "c1f1CapWind  " << pd.c1f1CapWind << '\n';
  ff << "c2f1CapCoal  " << pd.c2f1CapCoal << '\n';
  ff << "c2f1CapGas  " << pd.c2f1CapGas << '\n';
  ff << "c2f1CapSol  " << pd.c2f1CapSol << '\n';
  ff << "c2f1CapWind  " << pd.c2f1CapWind << '\n';
  ff << "c1f2CapCoal  " << pd.c1f2CapCoal << '\n';
  ff << "c1f2CapGas  " << pd.c1f2CapGas << '\n';
  ff << "c1f2CapSol  " << pd.c1f2CapSol << '\n';
  ff << "c1f2CapWind  " << pd.c1f2CapWind << '\n';
  ff << "c2f2CapCoal  " << pd.c2f2CapCoal << '\n';
  ff << "c2f2CapGas  " << pd.c2f2CapGas << '\n';
  ff << "c2f2CapSol  " << pd.c2f2CapSol << '\n';
  ff << "c2f2CapWind  " << pd.c2f2CapWind << '\n';
  ff << "c1GenNash  " << pd.c1GenNash << '\n';
  ff << "c2GenNash  " << pd.c2GenNash << '\n';
  ff << "c1demInt  " << pd.c1demInt << '\n';
  ff << "c1demSl  " << pd.c1demSl << '\n';
  ff << "c2demInt  " << pd.c2demInt << '\n';
  ff << "c2demSl  " << pd.c2demSl << '\n';
  ff << "c1emitCost  " << pd.c1emitCost << '\n';
  ff << "c2emitCost  " << pd.c2emitCost << '\n';
  ff.close();
}

cci::probData cci::readProbData(const std::string filename) {
  cci::probData pd;
  std::ifstream ff;
  ff.open(filename, std::ios::in);
  auto readLine = [&pd](std::ifstream &ff, const std::string name) {
    if (name == "prefix")
      ff >> pd.prefix;
    if (name == "probId")
      ff >> pd.probId;
    if (name == "suIn")
      ff >> pd.suIn;
    if (name == "suSl")
      ff >> pd.suSl;
    if (name == "c1Pp")
      ff >> pd.c1Pp;
    if (name == "c2Pp")
      ff >> pd.c2Pp;
    if (name == "c1CreditSplit")
      ff >> pd.c1CreditSplit;
    if (name == "c2CreditSplit")
      ff >> pd.c2CreditSplit;
    if (name == "c1CleanMix")
      ff >> pd.c1CleanMix;
    if (name == "c2CleanMix")
      ff >> pd.c2CleanMix;
    if (name == "c1tax")
      ff >> pd.c1tax;
    if (name == "c2tax")
      ff >> pd.c2tax;
    if (name == "prodCoalLin")
      ff >> pd.prodCoalLin;
    if (name == "prodGasLin")
      ff >> pd.prodGasLin;
    if (name == "prodSolLin")
      ff >> pd.prodSolLin;
    if (name == "prodWindLin")
      ff >> pd.prodWindLin;
    if (name == "prodCoalQuad")
      ff >> pd.prodCoalQuad;
    if (name == "prodGasQuad")
      ff >> pd.prodGasQuad;
    if (name == "prodSolQuad")
      ff >> pd.prodSolQuad;
    if (name == "prodWindQuad")
      ff >> pd.prodWindQuad;
    if (name == "invSolLin")
      ff >> pd.invSolLin;
    if (name == "invWindLin")
      ff >> pd.invWindLin;
    if (name == "invSolQuad")
      ff >> pd.invSolQuad;
    if (name == "invWindQuad")
      ff >> pd.invWindQuad;
    if (name == "emitCoal")
      ff >> pd.emitCoal;
    if (name == "emitGas")
      ff >> pd.emitGas;
    if (name == "emitSol")
      ff >> pd.emitSol;
    if (name == "emitWind")
      ff >> pd.emitWind;
    if (name == "c1f1CapCoal")
      ff >> pd.c1f1CapCoal;
    if (name == "c1f1CapGas")
      ff >> pd.c1f1CapGas;
    if (name == "c1f1CapSol")
      ff >> pd.c1f1CapSol;
    if (name == "c1f1CapWind")
      ff >> pd.c1f1CapWind;
    if (name == "c2f1CapCoal")
      ff >> pd.c2f1CapCoal;
    if (name == "c2f1CapGas")
      ff >> pd.c2f1CapGas;
    if (name == "c2f1CapSol")
      ff >> pd.c2f1CapSol;
    if (name == "c2f1CapWind")
      ff >> pd.c2f1CapWind;
    if (name == "c1f2CapCoal")
      ff >> pd.c1f2CapCoal;
    if (name == "c1f2CapGas")
      ff >> pd.c1f2CapGas;
    if (name == "c1f2CapSol")
      ff >> pd.c1f2CapSol;
    if (name == "c1f2CapWind")
      ff >> pd.c1f2CapWind;
    if (name == "c2f2CapCoal")
      ff >> pd.c2f2CapCoal;
    if (name == "c2f2CapGas")
      ff >> pd.c2f2CapGas;
    if (name == "c2f2CapSol")
      ff >> pd.c2f2CapSol;
    if (name == "c2f2CapWind")
      ff >> pd.c2f2CapWind;
    if (name == "c1GenNash")
      ff >> pd.c1GenNash;
    if (name == "c2GenNash")
      ff >> pd.c2GenNash;
    if (name == "c1demInt")
      ff >> pd.c1demInt;
    if (name == "c1demSl")
      ff >> pd.c1demSl;
    if (name == "c2demInt")
      ff >> pd.c2demInt;
    if (name == "c2demSl")
      ff >> pd.c2demSl;
    if (name == "c1emitCost")
      ff >> pd.c1emitCost;
    if (name == "c2emitCost")
      ff >> pd.c2emitCost;
  };
  while (!ff.eof()) {
    std::string name;
    ff >> name;
    readLine(ff, name);
  }
  ff.close();
  return pd;
}
