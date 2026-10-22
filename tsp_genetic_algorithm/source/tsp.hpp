#pragma once

#include <string>
#include <vector>

#include "random.hpp"

/********************************************/
/**************** city class ****************/
/********************************************/

class City {
private:
  std::string _name;
  std::vector<double> _pos{};

public:
  City(const std::string &, const std::vector<double> &);
  const std::string &get_name() const;
  const std::vector<double> &get_pos() const;
};

double distance(const City &, const City &);
std::vector<City> read_cities(const std::string &);
std::vector<std::vector<double>> distances_matrix(const std::vector<City> &);

/********************************************/
/***************** tsp class ****************/
/********************************************/

class TSP {
private:
  std::vector<City> _cities{};
  std::size_t _ncities = 0;
  std::vector<std::vector<double>> _distances_matrix{};

public:
  explicit TSP(const std::vector<City> &);
  explicit TSP(const std::string &);
  std::size_t n() const;
  const City &city(std::size_t) const;
  double d(std::size_t, std::size_t) const;
};

std::vector<City> cities_in_a_square(std::size_t, Random &);
std::vector<City> cities_on_a_circle(std::size_t, Random &);
double circle_optimum(const TSP &);

/********************************************/
/**************** route class ***************/
/********************************************/

class Route {
private:
  std::vector<std::size_t> _order{};
  double _length = 0.0;
  void _set_length(const TSP &);
  void _mutate_swap(const TSP &, Random &),
      _mutate_shift(const TSP &, Random &),
      _mutate_block_swap(const TSP &, Random &),
      _mutate_invert(const TSP &, Random &);
  double _p_swap = 0.05, _p_shift = 0.05, _p_block_swap = 0.05,
         _p_invert = 0.05;

public:
  Route(std::vector<std::size_t>, const TSP &);
  double length() const;
  bool operator<(const Route &) const;
  const std::vector<std::size_t> &order() const;
  void mutate(const TSP &, Random &);
};

/* random index in [min, max) */
std::size_t rnd_index(std::size_t min, std::size_t max, Random &rnd);
/* random order of indices [0, N), always starting from 0 */
std::vector<std::size_t> rnd_order(std::size_t N, Random &rnd);
bool validate_order(std::size_t, const std::vector<std::size_t> &);

void sort(std::vector<Route> &);
double recompute_length(const Route &, const TSP &);

/********************************************/
/************* population class *************/
/********************************************/

class Population {
private:
  std::vector<Route> _routes{};
  std::pair<Route, Route> _crossover(std::size_t, std::size_t, const TSP &,
                                     Random &) const;
  std::pair<Route, Route> _select_and_generate(const TSP &, Random &) const;
  double _select_weight_pow = 2.0, _p_crossover = 0.6;

public:
  Population(std::size_t, const TSP &, Random &rnd);
  const Route &best() const;
  double best_length() const;
  double half_mean_length() const;
  void check(const TSP &);
  void next_generation(const TSP &, Random &);
  void receive_migrant(const std::vector<std::size_t> &, const TSP &);
};
