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
void checksum(const std::vector<City>&, MpiEnv&);
