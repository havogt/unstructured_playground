#include <array>
#include <cassert>
#include <iostream> // TODO remove
#include <tuple>
#include <vector>

using grid_size = std::tuple<std::size_t, std::size_t>;
using cell_position = std::tuple<std::size_t, std::size_t, bool>; // i, j, color
using edge_position =
    std::tuple<std::size_t, std::size_t, std::size_t>; // i, j, color
using edges_of_cell_t = std::array<std::size_t, 3>;
using vertices_of_cell_t = std::array<std::size_t, 3>;
using vertices_of_edge_t = std::array<std::size_t, 2>;

auto get_x(grid_size const &size) { return std::get<0>(size); }
auto get_y(grid_size const &size) { return std::get<1>(size); }

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
    return {}
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

class grid {
public:
  grid(std::size_t Nx, std::size_t Ny)
      : c2v(2 * Nx * Ny), c2e(2 * Nx * Ny), e2v(3 * Nx * Ny) {
    std::size_t n_cells = 2 * Nx * Ny;
    std::size_t n_edges = 3 * Nx * Ny + Nx + Ny;
    std::size_t n_vertices = (Nx + 1) * (Ny + 1);

    for (std::size_t c = 0; c < n_cells; ++c) {
      c2v[c] = vertices_of_cell(c, {Nx, Ny});
      c2e[c] = edges_of_cell(c, {Nx, Ny});
    }
    for (std::size_t e = 0; e < n_edges; ++e) {
    }
  }

private:
  std::vector<vertices_of_cell_t> c2v;
  std::vector<edges_of_cell_t> c2e;
  std::vector<vertices_of_edge_t> e2v;
};

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
}
