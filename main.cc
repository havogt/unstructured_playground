#include <array>
#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream> // TODO remove
#include <ostream>
#include <range/v3/all.hpp>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

using grid_size = std::tuple<std::size_t, std::size_t>;
using cell_position = std::tuple<std::size_t, std::size_t, bool>; // i, j, color
using edge_position =
    std::tuple<std::size_t, std::size_t, std::size_t>;        // i, j, color
using vertex_position = std::tuple<std::size_t, std::size_t>; // i, j
using edges_of_cell_t = std::array<std::size_t, 3>;
using vertices_of_cell_t = std::array<std::size_t, 3>;
using vertices_of_edge_t = std::array<std::size_t, 2>;

auto get_x(grid_size const &size) { return std::get<0>(size); }
auto get_y(grid_size const &size) { return std::get<1>(size); }

vertex_position from_vertex(std::size_t index, grid_size size) {
  std::size_t i = index % (get_x(size) + 1);
  std::size_t j = index / (get_x(size) + 1);
  return {i, j};
}

cell_position from_cell(std::size_t index, grid_size size) {
  bool lower_triangle = index % 2 == 0;
  std::size_t i = (index / 2) % get_x(size);
  std::size_t j = (index / 2 / get_x(size));
  return {i, j, lower_triangle};
}

// i, j might be a cell outside of the grid
edge_position from_edge(std::size_t index, grid_size size) {
  auto normal_size = 3 * get_x(size) * get_y(size);
  if (index < normal_size) {
    // these are the non-special edges (not right-most or bottom-most)
    return {}; // TODO
  } else if (index < normal_size + get_y(size)) {
    // right-most
    return {get_x(size), index - normal_size, 2};
  } else {
    assert(index < normal_size + get_y(size) + get_x(size) && "You messed up!");
    return {index - normal_size - get_y(size), get_y(size), 0};
  }
}

vertices_of_cell_t vertices_of_cell(std::size_t index, grid_size size) {
  auto [i, j, lower] = from_cell(index, size);
  auto top_left = (get_x(size) + 1) * j + i;
  auto bottom_left = (get_x(size) + 1) * (j + 1) + i;
  auto top_right = (get_x(size) + 1) * j + (i + 1);
  auto bottom_right = (get_x(size) + 1) * (j + 1) + (i + 1);
  if (lower)
    return {top_left, bottom_right, bottom_left};
  else
    return {top_left, top_right, bottom_right};
}

//  x x x x x (0)
//  x x
//  x   x
//  x     x
//  x       x
// (2)        (1)
// the right-most edges are labelled Nx*Ny*3...Nx*Ny*3*+Ny-1
// the lower-most edges are labelled Nx*Ny*3*+Ny...Nx*Ny*3*+Ny+Nx-1
std::array<std::size_t, 3> edges_of_quad(std::size_t i, std::size_t j,
                                         grid_size size) {
  auto base = 3 * get_x(size) * j + 3 * i;
  if (i == get_x(size)) {
    return {0, 0,
            get_x(size) * get_y(size) * 3 +
                j}; // only 2 element should be read TODO build safety
  } else if (j == get_y(size)) {
    return {get_x(size) * get_y(size) * 3 + get_y(size) + i, 0,
            0}; // only 0 element should be read TODO build safety
  } else
    return {base, base + 1, base + 2};
}

edges_of_cell_t edges_of_cell(std::size_t index, grid_size size) {
  auto [i, j, lower] = from_cell(index, size);
  if (lower) {
    return {edges_of_quad(i, j, size)[1], edges_of_quad(i, j + 1, size)[0],
            edges_of_quad(i, j, size)[2]};
  } else {
    return {edges_of_quad(i, j, size)[0], edges_of_quad(i + 1, j, size)[2],
            edges_of_quad(i, j, size)[1]};
  }
}

template <class T, class... Arg> using id = T;

struct tuple_sum {
  template <class T0, class T1, std::size_t... Idx>
  std::tuple<double, double> impl_(std::tuple<T0, T0> a, std::tuple<T1, T1> b,
                                   std::index_sequence<Idx...>) {
    return {std::get<Idx>(a) + std::get<Idx>(b)...};
  }

  template <class T0, class T1>
  std::tuple<double, double> operator()(std::tuple<T0, T0> a,
                                        std::tuple<T1, T1> b) {
    return impl_(a, b, std::make_index_sequence<2>());
  }
};

auto center(std::vector<vertex_position> const &v) {
  auto vec_sum =
      ranges::accumulate(v, std::tuple<double, double>{}, tuple_sum{});
  return std::apply(
      [&v](auto... e) {
        return std::tuple<double, double>{e / static_cast<double>(v.size())...};
      },
      vec_sum);
}

std::string vertices_to_vtk(std::size_t Nx, std::size_t Ny) {
  std::ostringstream output;
  output << "POINTS " << (Nx + 1) * (Ny + 1) << " float\n";
  for (std::size_t i = 0; i < Nx + 1; ++i)
    for (std::size_t j = 0; j < Ny + 1; ++j) {
      output << i << " " << j << " 0\n";
    }
  return output.str();
}

std::string cells_to_vtk(std::vector<vertices_of_cell_t> const &cells) {
  std::ostringstream output;
  output << "CELLS " << cells.size() << " " << cells.size() * 4 << "\n";
  for (auto const &v : cells) {
    output << "3 " << v[0] << " " << v[1] << " " << v[2] << "\n";
  }
  output << "CELL_TYPES " << cells.size() << "\n";
  for (auto const &v : cells) {
    output << "5\n";
  }
  return output.str();
}

namespace to_vtk_impl {
std::string make_data_header(std::string const &data_type,
                             std::size_t n_cells) {
  std::ostringstream output;
  output << "CELL_DATA " << n_cells << "\nSCALARS temperature " << data_type
         << " 1\nLOOKUP_TABLE default\n";
  return output.str();
}

std::string get_data_type(float) { return "float"; }

std::string get_data_type(double) { return "double"; }
} // namespace to_vtk_impl

template <typename T> std::string cell_data_to_vtk(std::vector<T> const &data) {
  std::ostringstream output;
  output << to_vtk_impl::make_data_header(to_vtk_impl::get_data_type(T{}),
                                          data.size());
  for (auto d : data) {
    output << d << "\n";
  }
  return output.str();
}

class grid {
private:
  std::vector<std::size_t> e2c_impl(std::size_t id) const {
    std::vector<std::size_t> cell_neighs;
    for (std::size_t cell = 0; cell < c2e_.size(); ++cell) {
      for (auto e : c2e_[cell]) {
        if (id == e)
          cell_neighs.push_back(cell);
      }
    }
    return cell_neighs;
  }

public:
  grid(std::size_t Nx, std::size_t Ny)
      : size_i(Nx), size_j(Ny), c2v_(2 * Nx * Ny), c2e_(2 * Nx * Ny),
        e2v_(3 * Nx * Ny) {
    std::size_t n_edges = 3 * Nx * Ny + Nx + Ny;
    std::size_t n_vertices = (Nx + 1) * (Ny + 1);

    for (std::size_t c = 0; c < n_cells(); ++c) {
      c2v_[c] = vertices_of_cell(c, {Nx, Ny});
      c2e_[c] = edges_of_cell(c, {Nx, Ny});
    }
    for (std::size_t e = 0; e < n_edges; ++e) {
      e2c_.push_back(e2c_impl(e));
    }
  }

  std::size_t n_cells() const { return 2 * size_i * size_j; }

  std::string to_vtk() const {
    std::string output =
        "# vtk DataFile Version 3.0\n2D scalar data\nASCII\nDATASET "
        "UNSTRUCTURED_GRID\n";
    output += vertices_to_vtk(size_i, size_j);
    output += cells_to_vtk(c2v_);
    return output;
  }

  vertices_of_cell_t c2v(std::size_t id) const { return c2v_[id]; }

  auto e2c(std::size_t id) const { return e2c_[id]; }

  auto c2c(std::size_t id) const {
    std::set<std::size_t> cell_neighs;
    for (auto edge : c2e_[id]) {
      for (auto cell : e2c(edge)) {
        if (cell != id)
          cell_neighs.emplace(cell);
      }
    }
    return cell_neighs;
  };

  grid_size size() const { return {size_i, size_j}; }

private:
  std::size_t size_i;
  std::size_t size_j;
  std::vector<vertices_of_cell_t> c2v_;
  std::vector<edges_of_cell_t> c2e_;
  std::vector<vertices_of_edge_t> e2v_;

  std::vector<std::vector<size_t>> e2c_;
};

auto first() {
  return ranges::view::transform([](auto e) { return std::get<0>(e); });
}
auto second() {
  return ranges::view::transform([](auto e) { return std::get<1>(e); });
}

template <typename T>
void init_cells_with_gauss(T width, grid const &g, std::vector<T> &data) {
  using namespace ranges;

  data.resize(g.n_cells());

  std::vector<std::size_t> cell_ids = view::ints(0) | view::take(g.n_cells());

  auto vertex_positions =
      cell_ids | view::transform([&g](std::size_t i) { return g.c2v(i); }) |
      view::transform([&g](auto v_of_c) {
        std::vector<vertex_position> res;
        for (auto v : v_of_c) {
          res.push_back(from_vertex(v, g.size()));
        }
        return center(res); // TODO why not possible with range?
      });

  auto cell_vertex_pos_tup = view::zip(cell_ids, vertex_positions);
  for (auto const &[cell_id, center] : cell_vertex_pos_tup) {
    data[cell_id] = exp(
        -width * (pow(std::get<0>(center) - double(get_x(g.size())) / 2.0, 2) +
                  pow(std::get<1>(center) - double(get_y(g.size())) / 2.0, 2)));
  }
}

void write_vtk(std::string filename, grid const &g,
               std::vector<double> const &data) {
  std::string output;
  output += g.to_vtk();
  output += cell_data_to_vtk(data);

  std::ofstream ofile;
  ofile.open(filename);
  ofile << output << "\n";
  ofile.close();
}

template <typename Vec> void diff_explicit(Vec &temperature, grid const &g) {
  using namespace ranges;
  std::vector<std::size_t> cell_ids = view::ints(0) | view::take(g.n_cells());

  double alpha = 0.5;

  Vec tmp(temperature.size());

  std::vector<decltype(g.c2c(0))> cashed_c2c;
  for (auto cell : cell_ids) {
    cashed_c2c.push_back(g.c2c(cell));
  }

  for (std::size_t timestep = 0; timestep < 100; ++timestep) {
    std::cout << timestep << std::endl;
    for (auto cell : cell_ids) {
      // auto neighs = g.c2c(cell);
      auto neighs = cashed_c2c[cell];
      double new_val = temperature[cell];
      new_val -= alpha * static_cast<double>(neighs.size()) * temperature[cell];
      for (auto n : neighs) {
        new_val += alpha * temperature[n];
      }
      tmp[cell] = new_val;
      if (cell % 1000 == 0)
        std::cout << cell << std::endl;
    }

    write_vtk("diff_" + std::to_string(timestep) + ".vtk", g, tmp);
    std::swap(temperature, tmp);
  }
}

int main() {
  std::size_t cell_id = 9;
  std::size_t Nx = 3;
  std::size_t Ny = 2;
  auto [i, j, lower] = from_cell(cell_id, {3, 2});
  std::cout << cell_id << ": " << i << ", " << j << ", lower=" << lower
            << std::endl;

  auto vertices = vertices_of_cell(cell_id, {Nx, Ny});
  std::cout << cell_id << ": (" << vertices[0] << "," << vertices[1] << ","
            << vertices[2] << ")" << std::endl;

  auto edges = edges_of_cell(5, {Nx, Ny});
  std::cout << "(" << edges[0] << "," << edges[1] << "," << edges[2] << ")"
            << std::endl;
  edges = edges_of_cell(6, {Nx, Ny});
  std::cout << "(" << edges[0] << "," << edges[1] << "," << edges[2] << ")"
            << std::endl;
  edges = edges_of_cell(11, {Nx, Ny});
  std::cout << "(" << edges[0] << "," << edges[1] << "," << edges[2] << ")"
            << std::endl;
  edges = edges_of_cell(10, {Nx, Ny});
  std::cout << "(" << edges[0] << "," << edges[1] << "," << edges[2] << ")"
            << std::endl;

  grid g_tmp(Nx, Ny);
  auto cells = g_tmp.e2c(14);
  for (auto c : cells) {
    std::cout << c << " ";
  }
  std::cout << "\n-----------" << std::endl;

  auto cells2 = g_tmp.c2c(2);
  for (auto c : cells2) {
    std::cout << c << " ";
  }
  std::cout << "\n-----------" << std::endl;

  grid g(100, 100);
  std::vector<double> data(g.n_cells(), 1);
  init_cells_with_gauss(0.01, g, data);

  diff_explicit(data, g);

  // write_vtk("tmp_file, data);
}
