#include <iostream>
#include <map>
#include <range/v3/all.hpp>
#include <vector>

using namespace ranges;

template <typename T> struct access_ {
  std::size_t id;
  std::vector<T> const &data;

  T operator()() { return data[id]; }
};
template <class T> access_(std::size_t, std::vector<T> const &)->access_<T>;

template <class T0, class T1> struct plus_ {
  T0 t0;
  T1 t1;

  auto operator()() { return t0() + t1(); }
};
template <class T0, class T1> plus_(T0, T1)->plus_<T0, T1>;

auto neighs(std::map<std::size_t, std::vector<std::size_t>> const &c2c) {
  return view::transform([&c2c](std::size_t t) { return c2c.at(t); });
}

auto sum(std::vector<double> &data) {
  return view::transform([&data](auto e) {
    double tmp = 0;
    for (auto n : e) {
      tmp += data[n];
    }
    return tmp;
  });
}

int main() {
  std::map<std::size_t, std::vector<std::size_t>> c2c;
  c2c[0] = std::vector<std::size_t>{1, 5};
  c2c[1] = std::vector<std::size_t>{0, 2};
  c2c[2] = std::vector<std::size_t>{1, 3, 7};
  c2c[3] = std::vector<std::size_t>{2};
  c2c[4] = std::vector<std::size_t>{5};
  c2c[5] = std::vector<std::size_t>{4, 0, 6};
  c2c[6] = std::vector<std::size_t>{5, 7};
  c2c[7] = std::vector<std::size_t>{2, 6};

  std::vector<std::size_t> cells = view::ints(0) | view::take(c2c.size());

  std::vector<double> a(c2c.size(), 1.);
  a[2] = 2.;
  std::vector<double> b(c2c.size(), 3.);

  //   triangles | transform([](t) { data[t] - c *data[t] + t | cell_neigh(elem)
  //   })
  for (auto e : a) {
    std::cout << e << " ";
  }
  std::cout << std::endl;

  auto copy_stencil = cells | view::transform([&](std::size_t t) {
                        return plus_{access_{t, a}, access_{t, b}};
                      });

  for (auto c : copy_stencil) {
    std::cout << c() << " ";
  }
  std::cout << std::endl;

  auto sum_neighbours = cells | neighs(c2c) | sum(a);

  for (auto c : sum_neighbours) {
    std::cout << c << " ";
  }
  std::cout << std::endl;
}
