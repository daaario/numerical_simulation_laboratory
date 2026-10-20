#include <iostream>

#include "tsp.hpp"
#include "random.hpp"

int main() {
  TSP tsp("input/capoluoghi.csv");
  Random rnd("input/seed_in.dat");
  Population pop(500, tsp, rnd);

  for (int i = 0; i < 200; i++) {
    pop.next_generation(tsp, rnd);
  }

  std::cout << "Work in progress..." << std::endl;

  return 0;
};
