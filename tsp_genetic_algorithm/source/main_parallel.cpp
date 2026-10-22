#include "parallel.hpp"
#include "tsp.hpp"
#include <fstream>
#include <vector>

int main(int argc, char **argv) {
  /* environment objects */
  MpiEnv mpi(argc, argv);
  Random rnd(mpi.rank());
  std::string output_dir = "output";
  std::string hist_filename =
      output_dir + "/history_rank_" + std::to_string(mpi.rank()) + ".log";
  std::string best_filename = output_dir + "/best_all_ranks.log";
  std::string cities_filename = "input/capoluoghi.csv";
  /* tsp objects */
  std::size_t nroutes = 500, ngens = 2000, nmigr = 40;
  // std::size_t ndim = 2, ncities = 34;
  // std::vector<City> cities = generate_and_bcast_cities(ndim, ncities, rnd,
  // mpi); TSP tsp(cities); /* randomly-generated cities */
  TSP tsp(cities_filename); /* cities read from a file */
  Population pop(nroutes, tsp, rnd);
  /* output objects */

  /* run optimization and print history */
  std::ofstream hist(hist_filename);
  hist << "GEN,BEST,BEST_HALF_AVG\n"
       << 0 << "," << pop.best_length() << "," << pop.half_mean_length()
       << "\n";

  for (std::size_t gen = 0; gen < ngens; gen++) {
    pop.next_generation(tsp, rnd);
    if ((gen + 1) % nmigr == 0) {
      migrate(pop, tsp, rnd, mpi);
    }

    hist << gen << "," << pop.best_length() << "," << pop.half_mean_length()
         << "\n";
  }

  /* find global best and print it */
  double local_best = pop.best_length();
  if (mpi.rank() == rank_of_best(local_best, mpi)) {
    std::ofstream best(best_filename);
    const std::vector<std::size_t> &best_order = pop.best().order();
    for (auto i : best_order) {
      std::vector<double> pos = tsp.city(i).get_pos();
      for (std::size_t idx = 0; idx < pos.size() - 1; idx++) {
        best << pos[idx] << ",";
      }
      best << pos.back() << "\n";
    }
  }
}
