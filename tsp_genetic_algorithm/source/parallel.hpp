#pragma once

#include "mpi.h"
#include "random.hpp"
#include "tsp.hpp"

class MpiEnv {
private:
  int _rank = 0, _size = 1;

public:
  MpiEnv(int &argc, char **&argv);
  ~MpiEnv();
  int rank() const, size() const;
};

std::vector<City> generate_and_bcast_cities(std::size_t, std::size_t, Random &,
                                            const MpiEnv &);
void checksum(const std::vector<City> &, const MpiEnv &);
int rank_of_best(double local_best, const MpiEnv &);
void migrate(Population &pop, const TSP &, Random &, const MpiEnv &);
