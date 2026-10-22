#include <iostream>
#include "parallel.hpp"

MpiEnv::MpiEnv(int &argc, char **&argv) {
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &_size);
}

MpiEnv::~MpiEnv() { MPI_Finalize(); }

int MpiEnv::rank() const { return _rank; }

int MpiEnv::size() const { return _size; }

std::vector<City> generate_and_bcast_cities(std::size_t ndim,
                                            std::size_t ncities, Random &rnd,
                                            const MpiEnv &mpi) {
  std::vector<City> cities{};
  std::vector<double> flat_coord{};

  if (mpi.rank() == 0) {
    cities = cities_in_a_square(ncities, rnd);
    for (auto &c : cities) {
      for (auto x : c.get_pos()) {
        flat_coord.emplace_back(x);
      }
    }
  } else {
    flat_coord.resize(ndim * ncities);
  }
  MPI_Bcast(flat_coord.data(), ndim * ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (mpi.rank() == 0) {
    return cities;
  }

  for (std::size_t c = 0; c < ncities; c++) {
    std::string name = "city" + std::to_string(c);
    std::vector<double> pos{};
    for (std::size_t i = 0; i < ndim; i++) {
      pos.emplace_back(flat_coord[c * ndim + i]);
    }
    cities.emplace_back(name, pos);
  }

  return cities;
}

void checksum(const std::vector<City>& cities, MpiEnv& mpi) {
  double checksum = 0.0;
  for (auto& c : cities) {
    for (auto x : c.get_pos()) {
      checksum += x;
    }
  }

  double cmin, cmax;
  MPI_Reduce(&checksum, &cmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&checksum, &cmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (mpi.rank() == 0)
    std::cout << "cities checksum min/max: " << cmin << " " << cmax << std::endl;
}
