#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "tsp.hpp"

namespace {
constexpr double pi = 3.141592653589793;
} // namespace

/********************************************/
/**************** city class ****************/
/********************************************/

City::City(const std::string &name, const std::vector<double> &pos)
    : _name(name), _pos(pos) {}

const std::string &City::get_name() const { return _name; }
const std::vector<double> &City::get_pos() const { return _pos; }

double distance(const City &a, const City &b) {
  const std::vector<double> &a_pos = a.get_pos();
  const std::vector<double> &b_pos = b.get_pos();
  std::size_t ndim = a_pos.size();

  double dist = 0.0;
  for (std::size_t i = 0; i < ndim; i++) {
    double x = a_pos[i] - b_pos[i];
    dist += x * x;
  }

  return std::sqrt(dist);
}

std::vector<City> read_cities(const std::string &filename) {
  std::vector<City> cities;

  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Cannot open " + filename);
  }

  std::string line;
  std::getline(file, line); /* skip header */
  while (std::getline(file, line)) {
    std::stringstream ss(line);

    std::string name;
    std::vector<double> pos;

    std::getline(ss, name, ',');
    std::string field;
    while (std::getline(ss, field, ',')) {
      pos.emplace_back(std::stod(field));
    }

    cities.emplace_back(name, pos);
  }

  return cities;
}

std::vector<std::vector<double>>
distances_matrix(const std::vector<City> &cities) {
  std::size_t N = cities.size();
  std::vector<std::vector<double>> distances(N, std::vector<double>(N, 0.0));

  for (std::size_t i = 0; i < N; i++) {
    for (std::size_t j = i + 1; j < N; j++) {
      double d = distance(cities[i], cities[j]);
      distances[i][j] = d;
      distances[j][i] = d;
    }
  }

  return distances;
}

/********************************************/
/***************** tsp class ****************/
/********************************************/

TSP::TSP(const std::vector<City> &cities)
    : _cities(cities), _ncities(_cities.size()),
      _distances_matrix(distances_matrix(_cities)) {
  if (_ncities < 3) {
    throw std::runtime_error("Too few cities for the problem to be valid.");
  }
}
TSP::TSP(const std::string &filename) : TSP(read_cities(filename)) {}

std::size_t TSP::n() const { return _ncities; }
double TSP::d(std::size_t i, std::size_t j) const {
  return _distances_matrix[i][j];
}

const City &TSP::city(std::size_t i) const { return _cities[i]; }

std::vector<City> cities_in_a_square(std::size_t N, Random &rnd) {
  std::vector<City> cities{};
  for (std::size_t i = 0; i < N; i++) {
    std::string name = "city" + std::to_string(i);
    std::vector<double> pos{rnd.rannyu(), rnd.rannyu()};

    cities.emplace_back(name, pos);
  }

  return cities;
}

std::vector<City> cities_on_a_circle(std::size_t N, Random &rnd) {
  std::vector<City> cities{};
  for (std::size_t i = 0; i < N; i++) {
    std::string name = "city" + std::to_string(i);
    double theta = rnd.rannyu(0, 2 * pi);
    std::vector<double> pos{std::cos(theta), std::sin(theta)};

    cities.emplace_back(name, pos);
  }

  return cities;
}

double circle_optimum(const TSP &tsp) {
  std::size_t N = tsp.n();
  std::vector<std::size_t> indices(N);
  std::ranges::iota(indices.begin(), indices.end(), 0);

  auto angle = [&](std::size_t i) {
    const auto &pos = tsp.city(i).get_pos();
    return std::atan2(pos[1], pos[0]);
  };
  std::sort(indices.begin(), indices.end(),
            [&](std::size_t a, std::size_t b) { return angle(a) < angle(b); });

  double L = 0.0;
  for (std::size_t i = 0; i < N; i++) {
    L += tsp.d(indices[i], indices[(i + 1) % N]);
  }

  return L;
}

/********************************************/
/**************** route class ***************/
/********************************************/

Route::Route(std::vector<std::size_t> order, const TSP &tsp) : _order(order) {
  if (!validate_order(tsp.n(), _order)) {
    throw std::runtime_error(
        "Tried to build route with unvalid order of cities.");
  }
  _set_length(tsp);
}

void Route::_set_length(const TSP &tsp) {
  std::size_t ncities = _order.size();

  _length = 0.0;
  for (std::size_t i = 0; i < ncities - 1; i++) {
    _length += tsp.d(_order[i], _order[i + 1]);
  }
  _length += tsp.d(_order[ncities - 1], _order[0]);
}

double Route::length() const { return _length; }
bool Route::operator<(const Route &other) const {
  return _length < other._length;
}
const std::vector<std::size_t> &Route::order() const { return _order; }

void Route::_mutate_swap(const TSP &tsp, Random &rnd) {
  std::size_t N = _order.size();
  std::size_t i = rnd_index(1, N, rnd), j = rnd_index(1, N - 1, rnd);
  if (j >= i) {
    j++;
  }
  std::swap(_order[i], _order[j]);
  _set_length(tsp);
}

void Route::_mutate_shift(const TSP &tsp, Random &rnd) {
  std::size_t N = _order.size(), block_size = rnd_index(1, N - 1, rnd);
  std::size_t i_start = rnd_index(1, N - block_size, rnd),
              shift = rnd_index(1, N - block_size - i_start + 1, rnd);

  std::rotate(_order.begin() + i_start, _order.begin() + i_start + block_size,
              _order.begin() + i_start + block_size + shift);
  _set_length(tsp);
}

void Route::_mutate_block_swap(const TSP &tsp, Random &rnd) {
  std::size_t N = _order.size(), block_size = rnd_index(1, (N + 1) / 2, rnd);
  std::size_t i_start = rnd_index(1, N - 2 * block_size + 1, rnd),
              j_start =
                  rnd_index(i_start + block_size, N - block_size + 1, rnd);

  std::swap_ranges(_order.begin() + i_start,
                   _order.begin() + i_start + block_size,
                   _order.begin() + j_start);
  _set_length(tsp);
}

void Route::_mutate_invert(const TSP &tsp, Random &rnd) {
  std::size_t N = _order.size(), block_size = rnd_index(2, N, rnd);
  std::size_t i_start = rnd_index(1, N - block_size + 1, rnd);
  std::reverse(_order.begin() + i_start, _order.begin() + i_start + block_size);
  _set_length(tsp);
}

void Route::mutate(const TSP &tsp, Random &rnd) {
  if (rnd.rannyu() < _p_swap) {
    _mutate_swap(tsp, rnd);
  }
  if (rnd.rannyu() < _p_shift) {
    _mutate_shift(tsp, rnd);
  }
  if (rnd.rannyu() < _p_block_swap) {
    _mutate_block_swap(tsp, rnd);
  }
  if (rnd.rannyu() < _p_invert) {
    _mutate_invert(tsp, rnd);
  }
}

std::size_t rnd_index(std::size_t min, std::size_t max, Random &rnd) {
  return min + static_cast<std::size_t>((max - min) * rnd.rannyu());
}

std::vector<std::size_t> rnd_order(std::size_t N, Random &rnd) {
  std::vector<std::size_t> result;

  do {
    result = {0};
    std::vector<std::size_t> possible_indices;
    possible_indices.reserve(N - 1);
    for (std::size_t i = 1; i < N; i++) {
      possible_indices.emplace_back(i);
    }

    while (!possible_indices.empty()) {
      std::size_t i = rnd_index(0, possible_indices.size(), rnd);
      result.emplace_back(possible_indices[i]);

      std::swap(possible_indices[i], possible_indices.back());
      possible_indices.pop_back();
    }
  } while (!validate_order(N, result));

  return result;
}

bool validate_order(std::size_t N, const std::vector<std::size_t> &order) {
  if (order.size() != N || order[0] != 0) {
    return false;
  }

  std::vector<bool> checks(N, false);
  for (std::size_t i = 0; i < N; i++) {
    if (order[i] >= N || order[i] < 0) {
      return false;
    }
    checks[order[i]] = true;
  }
  for (auto check : checks) {
    if (!check) {
      return false;
    }
  }

  return true;
}

void sort(std::vector<Route> &routes) {
  std::sort(routes.begin(), routes.end());
}

double recompute_length(const Route &route, const TSP &tsp) {
  const auto &order = route.order();
  double L = 0.0;
  for (std::size_t i = 0; i < order.size(); i++)
    L += tsp.d(order[i], order[(i + 1) % order.size()]);
  return L;
}

/********************************************/
/************* population class *************/
/********************************************/

Population::Population(std::size_t N, const TSP &tsp, Random &rnd) {
  std::size_t ncities = tsp.n();
  for (std::size_t i = 0; i < N; i++) {
    _routes.emplace_back(rnd_order(ncities, rnd), tsp);
  }
  sort(_routes);
}

std::pair<Route, Route> Population::_crossover(std::size_t a, std::size_t b,
                                               const TSP &tsp,
                                               Random &rnd) const {
  if (rnd.rannyu() >= _p_crossover) {
    return {_routes[a], _routes[b]};
  }

  std::size_t N = tsp.n();
  std::size_t i_cut = rnd_index(2, N - 1, rnd);
  std::vector<std::size_t> order1 = _routes[a].order(),
                           order2 = _routes[b].order();

  std::vector<std::size_t> new_order1(order1.begin(), order1.begin() + i_cut),
      new_order2(order2.begin(), order2.begin() + i_cut);
  new_order1.reserve(N);
  new_order2.reserve(N);
  for (std::size_t i = 0; i < N; i++) {
    if (std::find(new_order1.begin(), new_order1.end(), order2[i]) ==
        new_order1.end()) {
      new_order1.emplace_back(order2[i]);
    }
    if (std::find(new_order2.begin(), new_order2.end(), order1[i]) ==
        new_order2.end()) {
      new_order2.emplace_back(order1[i]);
    }
  }

  return {Route(new_order1, tsp), Route(new_order2, tsp)};
}

std::pair<Route, Route> Population::_select_and_generate(const TSP &tsp,
                                                         Random &rnd) const {
  std::size_t N = _routes.size();

  double x = rnd.rannyu(0, N), y = rnd.rannyu(0, N - 1);
  if (y >= x) {
    y += 1;
  }
  std::size_t i = static_cast<std::size_t>(N *
                                           std::pow(x / N, _select_weight_pow)),
              j = static_cast<std::size_t>(N *
                                           std::pow(y / N, _select_weight_pow));

  std::pair<Route, Route> children = _crossover(i, j, tsp, rnd);
  children.first.mutate(tsp, rnd);
  children.second.mutate(tsp, rnd);
  return children;
}

const Route &Population::best() const { return _routes.front(); }

double Population::best_length() const { return _routes.front().length(); }

double Population::half_mean_length() const {
  std::size_t half_pop_size = _routes.size() / 2;
  double sum = 0.0;
  for (std::size_t i = 0; i < half_pop_size; i++) {
    sum += _routes[i].length();
  }
  return sum / half_pop_size;
}

void Population::check(const TSP &tsp) {
  for (auto route : _routes) {
    bool order_check = validate_order(tsp.n(), route.order());
    bool length_check =
        std::abs(recompute_length(route, tsp) - route.length()) < 1e-9;
    if (!(order_check && length_check)) {
      throw std::runtime_error("Population check failed");
    }
  }
}

void Population::next_generation(const TSP &tsp, Random &rnd) {
  std::vector<Route> new_routes;
  new_routes.reserve(_routes.size());

  while (new_routes.size() < _routes.size()) {
    auto children = _select_and_generate(tsp, rnd);

    new_routes.emplace_back(children.first);
    if (new_routes.size() < _routes.size()) {
      new_routes.emplace_back(children.second);
    }
  }

  sort(new_routes);
  _routes = std::move(new_routes);
}
