#include <iostream>
#include <range/v3/all.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

using namespace ranges;

//   auto vertices = cell_ids |
//                   view::transform([&g](std::size_t i) { return g.c2v(i); }) |
//                   view::transform([](auto &v_of_c) { return v_of_c; });

int main() {
  std::vector<std::size_t> cell_ids = view::ints(0) | view::take(100);
  auto trans = cell_ids | view::transform([](auto i) {
                 return std::vector{i, i + 1};
               }) |
               view::transform([](auto i) { return i; });
}
