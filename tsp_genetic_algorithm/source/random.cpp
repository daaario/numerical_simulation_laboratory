#include "random.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
constexpr std::array<int, 4> m = {2107, 4071, 1521, 502};
constexpr double twom12 = 1.0 / 4096.0; /* for efficiency */
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

Random::Random(int s) {
  _seed = std::array<int, 4>{0, 0, 0, s};

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

Random::~Random() {
  try {
    log_last_seed();
  } catch (const std::exception &e) {
    std::cerr << "Warning: could not log last seed: " << e.what() << std::endl;
  }
}

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
