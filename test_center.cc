#include <iostream>
#include <range/v3/all.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

using namespace ranges;

using vertex_position = std::array<double, 3>;
// using vertex_position = std::tuple<double, double, double, double>;
// using vertex_position = std::tuple<std::size_t,std::size_t>;

struct tuple_sum {
  template <class T, std::size_t... Idx>
  T impl_(T a, T b, std::index_sequence<Idx...>) {
    return {std::get<Idx>(a) + std::get<Idx>(b)...};
  }

  template <class... Elems>
  std::tuple<Elems...> operator()(std::tuple<Elems...> a,
                                  std::tuple<Elems...> b) {
    return impl_(a, b, std::make_index_sequence<sizeof...(Elems)>());
  }

  template <class T, std::size_t D>
  std::array<T, D> operator()(std::array<T, D> a, std::array<T, D> b) {
    return impl_(a, b, std::make_index_sequence<D>());
  }
};

auto center(std::vector<vertex_position> const &v) {
  auto vec_sum = accumulate(v, vertex_position{}, tuple_sum{});
  return std::apply(
      [&v](auto... e) {
        return vertex_position{e / static_cast<double>(v.size())...};
      },
      vec_sum);
}

int main() {
  std::vector<vertex_position> v{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
                                 {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};

  std::apply([](auto... v) { ((std::cout << v << " "), ...) << std::endl; },
             center(v));
}
