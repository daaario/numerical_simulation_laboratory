#include "mpi.h"
#include <iostream>

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int values[3];
  for (int i = 0; i < 3; i++) {
    if (rank == 0) {
      values[i] = i + 1;
      continue;
    }
    values[i] = 0;
  }

  std::cout << "Processo " << rank << ", pre-broadcast: (" << values[0] << ", " << values[1] << ", " << values[2] <<")" << std::endl;
  MPI_Bcast(values, 3, MPI_INTEGER, 0, MPI_COMM_WORLD);
  std::cout << "Processo " << rank << ", post-broadcast: (" << values[0] << ", " << values[1] << ", " << values[2] <<")" << std::endl;

  MPI_Finalize();
  return 0;
}
