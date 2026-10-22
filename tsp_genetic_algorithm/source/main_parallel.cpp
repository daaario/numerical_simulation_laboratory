#include "parallel.hpp"
#include "tsp.hpp"
#include <vector>

int main(int argc, char **argv) {
  MpiEnv mpi(argc, argv);
  Random rnd(mpi.rank());

  std::size_t ndim = 2, ncities = 34, nroutes = 500, ngens = 1000;
  std::vector<City> cities = generate_and_bcast_cities(ndim, ncities, rnd, mpi);

  TSP tsp(cities);
  Population pop(nroutes, tsp, rnd);

  for (std::size_t gen = 0; gen < ngens; gen++) {
    pop.next_generation(tsp, rnd);
  }
}
