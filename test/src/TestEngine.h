// STL
#include <tuple>
#include <vector>
#include <type_traits>

template<int N>
using IntType = std::integral_constant<int, N>;

//All method arguments have been consummed, we can call the test functor
template<
  template<typename...> class TestT,
  typename... Tn>
bool TestImp() {
  return TestT<Tn...>::test();
};

//Case where at least 1 tuple has been fully consummed : stub it with true value
template<
  template<typename...> class TestT,
  typename... Tn,
  typename... Vn>
bool TestImp( std::tuple<> const& un, Vn const&... vn ) {
  return true;
};

//Generic case: While we have a tuple (used a a type list here), we consume the
//first element of the tuple, and recursively call TestImp until all tuples
//have been consummed. In the later case, we can instanciate the test class with
//a full set of types
template<
  template<typename...> class TestT,
  typename... Tn,
  typename U0, typename... Un,
  typename... Vn>
bool TestImp( std::tuple<U0,Un...> const& un, Vn const&... vn ) {
  return TestImp<TestT,Tn...,U0>(vn...) &&
  TestImp<TestT,Tn...>(std::tuple<Un...>(),vn...);
};

