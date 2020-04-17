template <unsigned int n_Dirty, unsigned int n_Clean, unsigned int n_Scen>
std::ostream &operator<<(std::ostream &ost, const FollPar<n_Scen> P);
std::ostream &operator<<(std::ostream &ost, const std::pair<double, double> P);

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
