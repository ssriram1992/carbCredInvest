#include "bondad.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <gurobi_c++.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
// #include <rapidjson/document.h>
// #include <rapidjson/istreamwrapper.h>
// #include <rapidjson/prettywriter.h>
// #include <rapidjson/stringbuffer.h>
#include <vector>

// using namespace rapidjson;
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

//
// template <unsigned int num_scen>
// void cci::EPEC<num_scen>::writeSolutionJSON(string filename, const arma::vec
// x, const arma::vec z) const {
// /**
// * @brief Writes the computed Nash Equilibrium in the standard JSON solution
// * file
// * @p filename dictates the name of the .JSON solution file
// */
// StringBuffer s;
// PrettyWriter<StringBuffer> writer(s);
// writer.StartObject();
// writer.Key("Meta");
// writer.StartObject();
// writer.Key("isPureEquilibrium");
// writer.Bool(this->isPureStrategy());
// writer.Key("nCountries");
// writer.Uint(this->getNcountries());
// writer.Key("nFollowers");
// writer.StartArray();
// for (unsigned i = 0; i < this->getNcountries(); i++)
// writer.Uint(this->AllLeadPars.at(i).n_followers);
// writer.EndArray();
// writer.Key("Countries");
// writer.StartArray();
// for (unsigned i = 0; i < this->getNcountries(); i++) {
// writer.StartObject();
// writer.Key("FollowerStart");
// writer.Uint(this->getPosition(i, cci::LeaderVars::FollowerStart));
// writer.Key("NetImport");
// writer.Uint(this->getPosition(i, cci::LeaderVars::NetImport));
// writer.Key("NetExport");
// writer.Uint(this->getPosition(i, cci::LeaderVars::NetExport));
// writer.Key("CountryImport");
// writer.Uint(this->getPosition(i, cci::LeaderVars::CountryImport));
// writer.Key("Caps");
// writer.Uint(this->getPosition(i, cci::LeaderVars::Caps));
// writer.Key("Tax");
// writer.Uint(this->getPosition(i, cci::LeaderVars::Tax));
// if (this->AllLeadPars.at(i).LeaderParam.tax_revenue) {
// writer.Key("QuadraticTax");
// writer.Uint(this->getPosition(i, cci::LeaderVars::TaxQuad));
// }
// writer.Key("DualVar");
// writer.Uint(this->getPosition(i, cci::LeaderVars::DualVar));
// writer.Key("ConvHullDummy");
// writer.Uint(this->getPosition(i, cci::LeaderVars::ConvHullDummy));
// writer.Key("End");
// writer.Uint(this->getPosition(i, cci::LeaderVars::End));
// writer.Key("ShadowPrice");
// writer.Uint(
// this->getPosition(this->getNcountries() - 1, cci::LeaderVars::End) +
// i);
// writer.EndObject();
// }
// writer.EndArray();
// writer.EndObject();
// writer.Key("Solution");
// writer.StartObject();
// writer.Key("x");
// writer.StartArray();
// for (unsigned i = 0; i < x.size(); i++)
// writer.Double(x.at(i));
// writer.EndArray();
// writer.Key("z");
// writer.StartArray();
// for (unsigned i = 0; i < z.size(); i++)
// writer.Double(z.at(i));
// writer.EndArray();
// writer.EndObject();
// writer.EndObject();
// ofstream file(filename + ".json");
// file << s.GetString();
// }

// template <unsigned int num_scen>
// void cci::EPEC<num_scen>::readSolutionJSON(string filename) {
// /**
// * @brief Reads the solution file and load it in the current EPEC instance
// * **/
// ifstream ifs(filename + ".json");
// if (ifs.good()) {
// IStreamWrapper isw(ifs);
// Document d;
// try {
// d.ParseStream(isw);
// const Value &x = d["Solution"].GetObject()["x"];
// // const Value &z = d["Solution"].GetObject()["z"];
// arma::vec new_x;
// // arma::vec new_z;
// new_x.zeros(x.GetArray().Size());
// // new_z.zeros(z.GetArray().Size());
//
// for (SizeType i = 0; i < this->getnVarinEPEC(); i++)
// new_x.at(i) = x[i].GetDouble();
//
// // for (SizeType i = 0; i < this->getnVarinEPEC(); i++)
// // new_z.at(i) = z[i].GetDouble();
// ifs.close();
// this->warmstart(new_x);
// } catch (exception &e) {
// cerr << "Exception in cci::readSolutionJSON : cannot read instance "
// "file."
// << e.what() << '\n';
// throw;
// } catch (...) {
// cerr << "Exception in cci::readSolutionJSON: cannot read instance file."
// << '\n';
// throw;
// }
// } else {
// cerr << "Exception in cci::readSolutionJSON : file instance not found."
// << '\n';
// throw;
// }
// }

template <unsigned int num_scen>
void cci::EPECInstance<num_scen>::save(string filename) {
  /**
   * @brief Writes the current EPEC instance to the standard JSON instance file
   * @p filename dictates the name of the JSON instance file
   * @p epec contains the @p EPECInstance object with the data
   */
  // StringBuffer s;
  // PrettyWriter<StringBuffer> writer(s);
  // writer.StartObject();
  // writer.Key("nCountries");
  // writer.Uint(this->Countries.size());
  // writer.Key("Countries");
  // writer.StartArray();
  // for (unsigned i = 0; i < this->Countries.size(); i++) {
  // writer.StartObject();
  //
  // writer.Key("nFollowers");
  // writer.Uint(this->Countries.at(i).n_followers);
  //
  // writer.Key("Name");
  // string currName = this->Countries.at(i).name;
  // char nameArray[currName.length() + 1];
  // strcpy(nameArray, currName.c_str());
  // writer.String(nameArray);
  //
  // writer.Key("DemandParam");
  // writer.StartObject();
  // writer.Key("Alpha");
  // writer.Double(this->Countries.at(i).DemandParam.alpha);
  // writer.Key("Beta");
  // writer.Double(this->Countries.at(i).DemandParam.beta);
  // writer.EndObject();
  //
  // writer.Key("TransportationCosts");
  // writer.StartArray();
  // for (unsigned j = 0; j < this->Countries.size(); j++)
  // writer.Double(this->TransportationCosts(i, j));
  // writer.EndArray();
  //
  // writer.Key("LeaderParam");
  // writer.StartObject();
  // writer.Key("ImportLimit");
  // writer.Double(this->Countries.at(i).LeaderParam.import_limit);
  // writer.Key("ExportLimit");
  // writer.Double(this->Countries.at(i).LeaderParam.export_limit);
  // writer.Key("PriceLimit");
  // writer.Double(this->Countries.at(i).LeaderParam.price_limit);
  // writer.Key("TaxRevenue");
  // writer.Bool(this->Countries.at(i).LeaderParam.tax_revenue);
  // writer.Key("TaxationType");
  // switch (this->Countries.at(i).LeaderParam.tax_type) {
  // case cci::TaxType::StandardTax:
  // writer.Int(0);
  // break;
  // case cci::TaxType::SingleTax:
  // writer.Int(1);
  // break;
  // default:
  // writer.Int(2);
  // }
  // writer.EndObject();
  //
  // writer.Key("Followers");
  // writer.StartObject();
  //
  // writer.Key("Names");
  // writer.StartArray();
  // for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++) {
  // string currName = this->Countries.at(i).FollowerParam.names.at(j);
  // char nameArray[currName.length() + 1];
  // strcpy(nameArray, currName.c_str());
  // writer.String(nameArray);
  // }
  // writer.EndArray();
  //
  // writer.Key("Capacities");
  // writer.StartArray();
  // for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
  // writer.Double(this->Countries.at(i).FollowerParam.capacities.at(j));
  // writer.EndArray();
  //
  // writer.Key("LinearCosts");
  // writer.StartArray();
  // for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
  // writer.Double(this->Countries.at(i).FollowerParam.costs_lin.at(j));
  // writer.EndArray();
  //
  // writer.Key("QuadraticCosts");
  // writer.StartArray();
  // for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
  // writer.Double(this->Countries.at(i).FollowerParam.costs_quad.at(j));
  // writer.EndArray();
  //
  // writer.Key("EmissionCosts");
  // writer.StartArray();
  // for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
  // writer.Double(this->Countries.at(i).FollowerParam.emission_costs.at(j));
  // writer.EndArray();
  //
  // writer.Key("TaxCaps");
  // writer.StartArray();
  // for (unsigned j = 0; j < this->Countries.at(i).n_followers; j++)
  // writer.Double(this->Countries.at(i).FollowerParam.tax_caps.at(j));
  // writer.EndArray();
  //
  // writer.EndObject();
  //
  // writer.EndObject();
  // }
  // writer.EndArray();
  // writer.EndObject();
  // ofstream file(filename + ".json");
  // file << s.GetString();
  // file.close();
  std::cout << filename;
}

// template <unsigned int num_scen>
// void cci::EPECInstance<num_scen>::load(string filename) {
// /**
// * @brief Reads an instance file and return a vector of @p LeadAllPar that can
// * be fed to the EPEC class
// * @p filename dictates the name of the JSON instance file
// */
// ifstream ifs(filename + ".json");
// if (ifs.good()) {
// IStreamWrapper isw(ifs);
// Document d;
// try {
// d.ParseStream(isw);
// vector<cci::LeadAllPar> LAP = {};
// int nCountries = d["nCountries"].GetInt();
// arma::sp_mat TrCo;
// TrCo.zeros(nCountries, nCountries);
// for (int j = 0; j < nCountries; ++j) {
// const Value &c = d["Countries"].GetArray()[j].GetObject();
//
// cci::FollPar FP;
// const Value &cap = c["Followers"]["Capacities"];
// for (SizeType i = 0; i < cap.GetArray().Size(); i++) {
// FP.capacities.push_back(cap[i].GetDouble());
// }
// const Value &lc = c["Followers"]["LinearCosts"];
// for (SizeType i = 0; i < lc.GetArray().Size(); i++) {
// FP.costs_lin.push_back(lc[i].GetDouble());
// }
// const Value &qc = c["Followers"]["QuadraticCosts"];
// for (SizeType i = 0; i < qc.GetArray().Size(); i++) {
// FP.costs_quad.push_back(qc[i].GetDouble());
// }
// const Value &ec = c["Followers"]["EmissionCosts"];
// for (SizeType i = 0; i < ec.GetArray().Size(); i++) {
// FP.emission_costs.push_back(ec[i].GetDouble());
// }
// const Value &tc = c["Followers"]["TaxCaps"];
// for (SizeType i = 0; i < tc.GetArray().Size(); i++) {
// FP.tax_caps.push_back(tc[i].GetDouble());
// }
// const Value &nm = c["Followers"]["Names"];
// for (SizeType i = 0; i < nm.GetArray().Size(); i++) {
// FP.names.push_back(nm[i].GetString());
// }
// for (SizeType i = 0; i < c["TransportationCosts"].GetArray().Size();
// i++) {
// TrCo.at(j, i) = c["TransportationCosts"].GetArray()[i].GetDouble();
// }
// bool tax_revenue = false;
// if (c["LeaderParam"].HasMember("TaxRevenue")) {
// tax_revenue = c["LeaderParam"].GetObject()["TaxRevenue"].GetBool();
// }
// unsigned int tax_type = 0;
// if (c["LeaderParam"].HasMember("TaxationType")) {
// tax_type = c["LeaderParam"].GetObject()["TaxationType"].GetInt();
// }
// LAP.push_back(cci::LeadAllPar(
// FP.capacities.size(), c["Name"].GetString(), FP,
// {c["DemandParam"].GetObject()["Alpha"].GetDouble(),
// c["DemandParam"].GetObject()["Beta"].GetDouble()},
// {c["LeaderParam"].GetObject()["ImportLimit"].GetDouble(),
// c["LeaderParam"].GetObject()["ExportLimit"].GetDouble(),
// c["LeaderParam"].GetObject()["PriceLimit"].GetDouble(),
// tax_revenue, tax_type}));
// }
// ifs.close();
// this->Countries = LAP;
// this->TransportationCosts = TrCo;
// } catch (exception &e) {
// cerr << "Exception in cci::load : cannot read instance file." << e.what()
// << '\n';
// throw;
// } catch (...) {
// cerr << "Exception in cci::load : cannot read instance file." << '\n';
// throw;
// }
// } else {
// cerr << "Exception in cci::load : file instance not found." << '\n';
// throw;
// }
// }
ostream &cci::operator<<(ostream &ost, const cci::prn l) {
  switch (l) {
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

ostream &cci::operator<<(ostream &ost, const cci::DemPar P) {
  ost << "Demand Parameters: " << '\n';
  ost << "******************" << '\n';
  ost << "Price\t\t =\t\t " << P.alpha << "\t-\t" << P.beta << "  x   Quantity"
      << '\n';
  return ost;
}

ostream &cci::operator<<(ostream &ost, const cci::LeadPar P) {
  ost << "Leader Parameters: " << '\n';
  ost << "******************" << '\n';
  ost << std::fixed;
  ost << cci::prn::label << "Initial carbon credit allotted"
      << ":" << cci::prn::val << P.carbCredit_init;
  ost << '\n';
  ost << cci::prn::label << "Export Limit"
      << ":" << cci::prn::val
      << (P.export_limit < 0 ? std::numeric_limits<double>::infinity()
                             : P.export_limit);
  ost << '\n';
  ost << cci::prn::label << "Import Limit"
      << ":" << cci::prn::val
      << (P.import_limit < 0 ? std::numeric_limits<double>::infinity()
                             : P.import_limit);
  ost << '\n';
  ost << cci::prn::label << "Min reqd consumption"
      << ":" << cci::prn::val
      << (P.consum_limit < 0 ? std::numeric_limits<double>::infinity()
                             : P.consum_limit);
  ost << '\n';
  return ost;
}

template <unsigned int num_scen>
ostream &cci::operator<<(ostream &ost, const cci::EPECInstance<num_scen> I) {
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
  case cci::LeaderVars::Foll2Cpr:
    ost << "cci::LeaderVars::Foll2Cpr";
    break;
  case cci::LeaderVars::CarbPrice:
    ost << "cci::LeaderVars::CarbPrice";
    break;
  case cci::LeaderVars::CarbExp:
    ost << "cci::LeaderVars::CarbExp";
    break;
  case cci::LeaderVars::CarbImp:
    ost << "cci::LeaderVars::CarbImp";
    break;
  case cci::LeaderVars::EnergExpScen:
    ost << "cci::LeaderVars::EnergExpScen";
    break;
  case cci::LeaderVars::EnergImpScen:
    ost << "cci::LeaderVars::EnergImpScen";
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
