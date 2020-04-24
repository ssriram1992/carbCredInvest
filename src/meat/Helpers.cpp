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
