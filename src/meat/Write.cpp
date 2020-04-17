
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
  file << prn::label << "Carbon Price: " << prn::val << x.at(this->getPosition(this->getNcountries() - 1, LeaderVars::End))<<'\n';

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
  file << prn::label << "Net carbon Export: " << prn::val << carbExp - carbImp << "\n";
  file << prn::label << "Domestic carbon price: " << prn::val << carbDomPrice << "\n";

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
        "prod_" + std::to_string(i) + "_" + std::to_string(j) + "_" +
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
