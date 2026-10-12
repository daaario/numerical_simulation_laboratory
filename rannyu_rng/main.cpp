#include <iostream>
#include <vector>

#include "random.hpp"

int main(int argc, char* argv[]) {
  /* argv[1]: number of samples *
   * optional, argv[2]: name of distribution (gauss, exp, lorentz) *
   *    - if not provided: sample uniformly in [0, 1)
   * argv[3], argv[4], ... : parameters of the distribution*/

  std::vector<std::string> args(argv, argv+argc);
  if (args.size() < 2) {
    throw std::runtime_error("Error: at least one argument is required.");
  }

  Random rng;
  int N = std::stoi(args[1]);

  if (args.size() == 2) {
    for (int i = 0; i < N; i++) {
      std::cout << rng.rannyu() << std::endl;
    }
  } else if (args[2] == "gauss") {
    if (args.size() != 5) {
      throw std::runtime_error("Error: invalid number of arguments in sampling a Gaussian distribution.");
    }
    double mean = std::stof(args[3]);
    double sigma = std::stof(args[4]);

    for (int i = 0; i < N; i++) {
      std::cout << rng.gauss(mean, sigma) << std::endl;
    }
  } else if (args[2] == "exp") {
    if (args.size() != 4)
    {
      throw std::runtime_error("Error: invalid number of arguments in sampling a decreasing exponential distribution.");
    }
    double lambda = std::stof(args[3]);

    for (int i = 0; i < N; i++)
    {
      std::cout << rng.exponential(lambda) << std::endl;
    }
  } else if (args[2] == "lorentz") {
    if (args.size() != 5)
    {
      throw std::runtime_error("Error: invalid number of arguments in sampling a Cauchy-Lorentz distribution.");
    }
    double mean = std::stof(args[3]);
    double gamma = std::stof(args[4]);

    for (int i = 0; i < N; i++)
    {
      std::cout << rng.cauchy_lorentz(mean, gamma) << std::endl;
    }
  } else if (args[2] == "linear") {
    if (args.size() != 3)
    {
      throw std::runtime_error("Error: invalid number of arguments in sampling a decreasing linear distribution.");
    }
    for (int i = 0; i < N; i++)
    {
      std::cout << rng.linear_decreasing() << std::endl;
    }
  } else if (args[2] == "buffon") {
    if (args.size() != 6)
    {
      throw std::runtime_error("Error: invalid number of arguments in performing a series of throws.");
    }
    int N_blocks = std::stoi(args[1]);
    double L = std::stof(args[3]);
    double d = std::stof(args[4]);
    int N_throws_per_block = std::stoi(args[5]);

    for (int i = 0; i < N_blocks; i++) {
      rng.buffon(L, d, N_throws_per_block);
    }
  } else if (args[2] == "rw_lattice") {
    if (args.size() != 5)
    {
      throw std::runtime_error("Error: invalid number of arguments in performing a series of random walks.");
    }
    int N_reps = std::stoi(args[1]);
    double a = std::stof(args[3]);
    int N_steps = std::stoi(args[4]);

    for (int i = 0; i < N_reps; i++) {
      rng.rw_lattice(a, N_steps);
    }
  } else if (args[2] == "rw_continuum") {
    if (args.size() != 5)
    {
      throw std::runtime_error("Error: invalid number of arguments in performing a series of random walks.");
    }
    int N_reps = std::stoi(args[1]);
    double a = std::stof(args[3]);
    int N_steps = std::stoi(args[4]);

    for (int i = 0; i < N_reps; i++) {
      rng.rw_continuum(a, N_steps);
    }
  } else if (args[2] == "wiener") {
    if (args.size() != 5)
    {
      throw std::runtime_error("Error: invalid number of arguments in sampling a 1D Brownian motion.");
    }
    int N_reps = std::stoi(args[1]);
    double dt = std::stof(args[3]);
    int N_steps = std::stoi(args[4]);

    for (int i = 0; i < N_reps; i++) {
      rng.wiener(dt, N_steps);
    }
  } else {
    throw std::runtime_error("Error: user-provided probability distribution is unknown.");
  }

  return 0;
}
