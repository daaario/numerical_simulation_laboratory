#include "random.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

namespace {
constexpr std::array<int, 4> m = {2107, 4071, 1521, 502};
constexpr double twom12 = 1.0 / 4096.0; /* for efficiency */
constexpr double pi = 3.141592653589793;
} // namespace

Random::Random(std::string seed_input_filepath) {
  read_seed(seed_input_filepath);

  std::ifstream primes_stream("input/primes.dat");
  if (!primes_stream) {
    throw std::runtime_error(
        "Unable to open or find \'input/primes.dat\' file.");
  }

  std::pair<int, int> primes;
  primes_stream >> primes.first >> primes.second;

  set_seed(_seed, primes);
}

Random::Random() : Random("input/seed_in.dat") {}

void Random::read_seed(std::string seed_input_filepath) {
  std::ifstream seed_input_stream(seed_input_filepath);
  if (!seed_input_stream) {
    throw std::runtime_error("Unable to open or find \'" + seed_input_filepath +
                             "\' file.");
  }

  std::array<int, 4> seed;
  for (int i = 0; i < 4; i++) {
    seed_input_stream >> seed[3 - i];
  }

  for (auto s : seed) {
    if (s < 0 || s >= 4096) {
      throw std::runtime_error("Seed out of bounds in \'" +
                               seed_input_filepath + "\' file.");
    }
  }

  _seed = seed;
}

void Random::set_seed(const std::array<int, 4> new_seed,
                      const std::pair<int, int> primes) {
  _seed = new_seed;
  _n = {primes.second, primes.first, 0, 0};
}

void Random::log_last_seed() const {
  std::ofstream seed_log_stream("output/last_seed.log");
  if (!seed_log_stream) {
    throw std::runtime_error(
        "Unable to open or generate \'output/last_seed.log\' file.");
  }

  for (int i = 0; i < 4; i++) {
    seed_log_stream << " " << _seed[3 - i];
  }
  seed_log_stream << std::endl;
}

double Random::rannyu() {
  const std::array<int, 4> seed_old = _seed;
  std::array<int, 4> tmp{};
  double r = 0.;

  for (int i = 0; i < 4; i++) {
    /* tmp values from seed, m and n */
    tmp[i] = _n[i];
    for (int j = 0; j <= i; j++) {
      tmp[i] += seed_old[i - j] * m[j];
    }

    if (i > 0) {
      tmp[i] += tmp[i - 1] / 4096;
    }

    /* seed values from tmp */
    _seed[i] = tmp[i] % 4096;

    /* build random variable */
    r += _seed[i];
    r = r * twom12;
  }

  return r;
}

double Random::rannyu(double min, double max) {
  return min + (max - min) * rannyu();
}

double Random::gauss(double mean, double sigma) {
  double s = rannyu(), t = rannyu();
  /* normal distribution from 2 uniform numbers in [0,1) */
  double x = std::sqrt(-2. * std::log(1. - s)) * std::cos(2. * pi * t);
  return mean + x * sigma;
}

double Random::exponential(double lambda) {
  if (lambda <= 0.0) {
    throw std::runtime_error("Error: invalid parameters in sampling a "
                             "decreasing exponential distribution.");
  }
  return -std::log(1 - rannyu()) / lambda;
}

double Random::cauchy_lorentz(double mean, double gamma) {
  if (gamma <= 0.0) {
    throw std::runtime_error(
        "Error: invalid parameters in sampling a Cauchy-Lorentz distribution.");
  }
  return mean + gamma * std::tan(pi * (rannyu() - 0.5));
}

double Random::linear_decreasing() { return 1 - std::sqrt(1 - rannyu()); }

void Random::buffon(double L, double d, int N_throws) {
  int hits = 0;
  for (int i = 0; i < N_throws; i++) {
    double x = rannyu(), y = rannyu();
    double r2 = x * x + y * y;
    while (r2 > 1) {
      x = rannyu();
      y = rannyu();
      r2 = x * x + y * y;
    }
    double sin_theta = std::abs(y) / std::sqrt(r2);
    double y_center = rannyu() * d / 2;

    if (y_center <= sin_theta * L / 2) {
      hits++;
    }
  }
  std::cout << hits << std::endl;
  ;
}

void Random::rw_lattice(double a, int N_steps) {
  int x = 0, y = 0, z = 0;
  std::cout << 0 << " ";

  for (int i = 0; i < N_steps; i++) {
    double rand = 6*rannyu();
    if (rand < 1.0) {
      x++;
    } else if (rand < 2.0) {
      x--;
    } else if (rand < 3.0) {
      y++;
    } else if (rand < 4.0) {
      y--;
    } else if (rand < 5.0) {
      z++;
    } else {
      z--;
    }
    std::cout << std::abs(a) * std::sqrt(x*x + y*y + z*z) << " ";
  }
  std::cout << std::endl;
}

void Random::rw_continuum(double a, int N_steps) {
  double x = 0, y = 0, z = 0;
  std::cout << 0 << " ";

  for (int i = 0; i < N_steps; i++) {
    double theta = pi * (rannyu() - 0.5); // CHANGE!
    double phi = 2 * pi * rannyu();

    x += std::cos(theta) * std::cos(phi);
    y += std::cos(theta) * std::sin(phi);
    z += std::sin(theta);

    std::cout << std::abs(a) * std::sqrt(x*x + y*y + z*z) << " ";
  }
  std::cout << std::endl;
}

void Random::wiener(double dt, int N_steps) {
  double x = 0;
  std::cout << 0 << " ";

  for (int i = 0; i < N_steps; i++) {
    double dx = gauss(0, std::sqrt(dt));
    x += dx;
    std::cout << x << " ";
  }

  std::cout << std::endl;
}
