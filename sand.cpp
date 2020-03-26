#include <array>
#include <iostream>
#include <map>
#include <tuple>

enum class fruits { apple, orange, grapes };
enum class flowers { rose, lotus, tulips };

template <class t> using params = std::map<t, double>;
template <class t, class u> using params2 = std::map<std::tuple<t, u>, double>;
template <class t, class u, class v>
using params3 = std::map<std::tuple<t, u, v>, double>;

template <long unsigned int N>
using gamsSet = const std::array<const std::string, N>;

// std::array<std::string, 3>
gamsSet<3> arr_Fruits{"apple", "orange", "grape"};

template <class t, unsigned int N>
std::map<t, double> initParam(const std::array<t, N> arr) {
  std::map<t, double> myMap;
  for (const auto &x : arr)
    myMap.at(x) = 0;
  return myMap;
}

int main() {
  using std::cout;
  params<fruits> taste;
  params<flowers> beauty;
  taste[fruits::orange] = -1;
  taste[fruits::grapes] = 0;
  beauty[flowers::rose] = 0;

  params<std::string> a_taste = initParam(arr_Fruits);
  for (const auto &x : a_taste)
    cout << x.first << " " << x.second << "\n";
  for (const auto &x : arr_Fruits) {
    a_taste[x] = 1;
  }

  params2<fruits, flowers> value;
  value[{fruits::orange, flowers::rose}] = 10;

  cout << value[{fruits::orange, flowers::rose}] << " "
       << value[{fruits::grapes, flowers::rose}] << "\n";
  cout << "Beauty: ";
  // for (auto x:flowers) cout << beauty[x] << " ";

  cout << "\n";
  for (const auto &x : a_taste)
    cout << x.first << " " << x.second << "\n";

  cout << "\n";

  return 0;
}
