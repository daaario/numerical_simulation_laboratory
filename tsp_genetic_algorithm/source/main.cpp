#include <fstream>

#include "tsp.hpp"
#include "random.hpp"

int main() {
  Random rnd("input/seed_in.dat");

  std::size_t ncities = 34;
  std::size_t nroutes = 500;
  std::size_t ngens = 1000;
  TSP tsp_square(cities_in_a_square(ncities, rnd));
  TSP tsp_circle(cities_on_a_circle(ncities, rnd));
  Population pop_circle(nroutes, tsp_circle, rnd);
  Population pop_square(nroutes, tsp_square, rnd);
  std::ofstream circle_hist("output/circle_history.log");
  std::ofstream square_hist("output/square_history.log");

  circle_hist << "GEN,BEST,HALF_AVG" << std::endl;
  square_hist << "GEN,BEST,HALF_AVG" << std::endl;
  for (std::size_t i = 0; i < ngens; i++) {
    pop_circle.next_generation(tsp_circle, rnd);
    pop_square.next_generation(tsp_square, rnd);

    circle_hist << i << "," << pop_circle.best_length() << "," << pop_circle.half_mean_length() << std::endl;
    square_hist << i << "," << pop_square.best_length() << "," << pop_square.half_mean_length() << std::endl;
  }

  std::ofstream circle_best("output/circle_best.log");
  const std::vector<std::size_t>& circle_best_order = pop_circle.best().order();
  for (auto i : circle_best_order) {
    std::vector<double> pos = tsp_circle.city(i).get_pos();
    for (std::size_t x = 0; x < pos.size() - 1; x++) {
      circle_best << pos[x] << ",";
    }
    circle_best << pos[pos.size() - 1] << std::endl;
  }

  std::ofstream square_best("output/square_best.log");
  const std::vector<std::size_t>& square_best_order = pop_square.best().order();
  for (auto i : square_best_order) {
    std::vector<double> pos = tsp_square.city(i).get_pos();
    for (std::size_t x = 0; x < pos.size() - 1; x++) {
      square_best << pos[x] << ",";
    }
    square_best << pos[pos.size() - 1] << std::endl;
  }

  return 0;
};
