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
                             : P.consum_limit);
  ost << '\n';
  ost << "Investment Incentives: \n";
  for (const auto &cc : P.cleanInvVal)
    ost << cci::prn::label << "Investment " << cci::prn::val << cc.first << ";"
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
