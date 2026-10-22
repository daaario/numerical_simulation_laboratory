# Traveling Salesman Problem - Genetic algorithm
C++ code for solving the TSP for either a randomized collection of cities (generated on a unit circumference or inside a $[-1,1]^2$ square), or a list of cities read from *'input/capoluoghi.dat'*.
This is done using a GA with hardcoded parameters, with two different 

  - *'main.cpp'* compiles to *'run'* for solving it in a single process
  - *'main_parallel.cpp'* compiles to *'run_parallel'*, which makes use of the MPI library for parallelizing the process.

The makefile compiles everything using `mpicxx`. The parallel-computing version must be run using `mpiexec -n <number_of_parallel_processes> ./run_parallel`.
