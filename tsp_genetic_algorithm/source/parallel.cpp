#include "parallel.hpp"
#include <iostream>

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

void checksum(const std::vector<City> &cities, const MpiEnv &mpi) {
  double checksum = 0.0;
  for (auto &c : cities) {
    for (auto x : c.get_pos()) {
      checksum += x;
    }
  }

  double cmin, cmax;
  MPI_Reduce(&checksum, &cmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&checksum, &cmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (mpi.rank() == 0)
    std::cout << "cities checksum min/max: " << cmin << " " << cmax
              << std::endl;
}

int rank_of_best(double local_best, const MpiEnv &mpi) {
  struct {
    double val;
    int rank;
  } in{local_best, mpi.rank()}, out; /* in-place definition of struct,    *
                                      * and declaration of two instances, *
                                      * only one of which is defined      */
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  return out.rank;
}

void migrate(Population &pop, const TSP &tsp, Random &rnd, const MpiEnv &mpi) {
  const int rank = mpi.rank(), size = mpi.size();
  if (size == 1) {
    return;
  }

  std::vector<int> pop_order(size);
  if (rank == 0) {
    std::vector<std::size_t> pop_order_sizet = rnd_order(size, rnd);
    pop_order = std::vector<int>(
        pop_order_sizet.begin(),
        pop_order_sizet.end()); /* range constructor: converts each element (ok
                                   for sending via MPI) */
  }

  MPI_Bcast(pop_order.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

  const int pos_in_pop_order =
      std::find(pop_order.begin(), pop_order.end(), rank) - pop_order.begin();
  const int destination_rank = pop_order[(pos_in_pop_order + 1) % size];
  const int source_rank = pop_order[(pos_in_pop_order - 1 + size) % size];

  const std::vector<std::size_t> &best_order = pop.best().order();
  std::vector<int> buffer(best_order.begin(), best_order.end());

  MPI_Sendrecv_replace(buffer.data(), buffer.size(), MPI_INT, destination_rank,
                       0, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  pop.receive_migrant(std::vector<std::size_t>(buffer.begin(), buffer.end()),
                      tsp);
}
