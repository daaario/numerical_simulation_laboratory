#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include "tsp.hpp"

int main() {
  Random rnd;
  //  for (std::size_t N : {3, 4, 5, 6, 7, 10, 34}) {
  //    std::vector<City> cities;
  //    for (std::size_t i = 0; i < N; i++) { 
  //      cities.emplace_back("c" + std::to_string(i),
  //      std::vector<double>{rnd.rannyu(), rnd.rannyu()}); 
  //    }
  //    TSP tsp(cities);

  //    for (int t = 0; t < 20000; t++) {
  //      Route route(rnd_order(N, rnd), tsp);
  //      route.mutate(tsp, rnd);
  //      assert(validate_order(N, route.order()));
  //      assert(std::abs(recompute_length(route, tsp) - route.length()) < 1e-9);
  //    }

  //    Population pop(200, tsp, rnd);
  //    for (int g = 0; g < 50; g++) { 
  //      pop.next_generation(tsp, rnd); 
  //    }
  //    pop.check(tsp);
  //    std::cout << "N=" << N << " ok\n";
  //  }

  TSP tsp(cities_in_a_square(50, rnd));
  Population pop(500, tsp, rnd);

  std::ofstream hist("output/history.log");
  hist << "#gen best half_mean\n";
  hist << 0 << " " << pop.best_length() << " " << pop.half_mean_length() << std::endl;
  for (int g = 0; g < 1000; g++) {
    pop.next_generation(tsp, rnd);
    hist << g << " " << pop.best_length() << " " << pop.half_mean_length() << std::endl;
  }

  //  hist << "best:          " << pop.best_length() << "\n"
  //       << "optimum:       " << circle_optimum(tsp) << "\n"
  //       << "diff:          " << pop.best_length() - circle_optimum(tsp) << "\n"
  //       << "relative diff: " << (pop.best_length() - circle_optimum(tsp)) / circle_optimum(tsp) << std::endl;

  std::ofstream broute("output/best_route.log");
  const std::vector<std::size_t>& best_order = pop.best().order();

  for (auto i : best_order) {
    std::vector<double> pos = tsp.city(i).get_pos();
    for (std::size_t x = 0; x < pos.size() - 1; x++) {
      broute << pos[x] << ",";
    }
    broute << pos[pos.size() - 1] << std::endl;
  }
}
