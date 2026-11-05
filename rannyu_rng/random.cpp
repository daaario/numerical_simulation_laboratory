#include "random.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

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
    double rand = 6 * rannyu();
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
    std::cout << std::abs(a) * std::sqrt(x * x + y * y + z * z) << " ";
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

    std::cout << std::abs(a) * std::sqrt(x * x + y * y + z * z) << " ";
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

void Random::hydrogen(std::array<int, 3> quant_num, int N, double delta) {
  auto hydrogen_prob_dist = [quant_num](std::array<double, 3> pos) {
    double r = 0.0;
    for (auto x : pos) {
      r += x * x;
    }
    r = std::sqrt(r);

    if (quant_num == std::array<int, 3>{1, 0, 0}) {
      return std::exp(-2 * r);
    } else if (quant_num == std::array<int, 3>{2, 1, 0}) {
      double z = pos[2];
      return z * z * std::exp(-r);
    } else {
      throw std::runtime_error(
          "State (n,l,m)=(" + std::to_string(quant_num[0]) + "," +
          std::to_string(quant_num[1]) + "," + std::to_string(quant_num[2]) +
          ") not supported.");
    }
  };

  auto pos_new = [this, delta](std::array<double, 3> &pos_old) {
    double theta = std::asin(2. * rannyu() - 1.);
    double phi = 2 * pi * rannyu();

    std::array<double, 3> pos{
        pos_old[0] + delta * std::cos(theta) * std::cos(phi),
        pos_old[1] + delta * std::cos(theta) * std::sin(phi),
        pos_old[2] + delta * std::sin(theta)};
    return pos;
  };

  std::array<double, 3> pos{0, 0, 1}, pos_old = pos;
  double prob = hydrogen_prob_dist(pos), prob_old = prob;
  int nacc = 0, natt = 0;
  for (int i = 0; i < N; i++) {
    do {
      pos = pos_new(pos_old);
      natt++;
    } while (rannyu() >= std::min(1.0, prob / prob_old));

    std::cout << pos[0] << "," << pos[1] << "," << pos[2] << "\n";
    nacc++;
    pos_old = pos;
    prob_old = prob;
  }
  std::cerr << "Acceptance ratio: " << double(nacc) / natt << std::endl;
}

void Random::hydrogen_nu(std::array<int, 3> quant_num, int N, double delta) {
  auto hydrogen_prob_dist = [quant_num](std::array<double, 3> pos) {
    double r = 0.0;
    for (auto x : pos) {
      r += x * x;
    }
    r = std::sqrt(r);

    if (quant_num == std::array<int, 3>{1, 0, 0}) {
      return std::exp(-2 * r);
    } else if (quant_num == std::array<int, 3>{2, 1, 0}) {
      double z = pos[2];
      return z * z * std::exp(-r);
    } else {
      throw std::runtime_error(
          "State (n,l,m)=(" + std::to_string(quant_num[0]) + "," +
          std::to_string(quant_num[1]) + "," + std::to_string(quant_num[2]) +
          ") not supported.");
    }
  };

  auto pos_new = [this, delta](std::array<double, 3> &pos_old) {
    std::array<double, 3> pos{pos_old[0] + gauss(0.0, delta),
                              pos_old[1] + gauss(0.0, delta),
                              pos_old[2] + gauss(0.0, delta)};
    return pos;
  };

  std::array<double, 3> pos{0, 0, 1}, pos_old = pos;
  double prob = hydrogen_prob_dist(pos_old), prob_old = prob;

  //  int nacc = 0, natt = 0;
  for (int i = 0; i < N; i++) {
    do {
      pos = pos_new(pos_old);
      prob = hydrogen_prob_dist(pos);
      //  natt++;
    } while (rannyu() >= std::min(1.0, prob / prob_old));

    std::cout << pos[0] << "," << pos[1] << "," << pos[2] << "\n";
    //  nacc++;
    pos_old = pos;
    prob_old = prob;
  }
  //  std::cerr << "Acceptance ratio: " << double(nacc) / natt << std::endl;
}

std::vector<double> Random::ex8_energy_block_avgs(double mu, double sigma,
                                                  int Nblk, int Nsmp) {
  /* integrand evaluation */
  auto integrand = [mu, sigma](double x) {
    double m2 = mu * mu, s2 = sigma * sigma, x2 = x * x;
    double s4 = s2 * s2, x4 = x2 * x2;
    return x4 - 2.5 * x2 - (x2 + m2 - s2) / (2 * s4) +
           mu * x * std::tanh(mu * x / s2) / s4;
  };

  /* probability density evaluation */
  auto prob = [mu, sigma](double x) {
    double exp_minus = std::exp(-(x - mu) * (x - mu) / (2 * sigma * sigma));
    double exp_plus = std::exp(-(x + mu) * (x + mu) / (2 * sigma * sigma));
    return (exp_plus + exp_minus) * (exp_plus + exp_minus);
  };

  /* sample the distributuion and compute block averages at the same time*/
  double delta_x = 2.5 * std::sqrt(sigma * sigma / 2 + mu * mu);
  auto x_new = [this, delta_x](double x_old) {
    return x_old + rannyu(-delta_x, delta_x);
  };

  double x = mu, x_old = x;
  double p = prob(x), p_old = p;
  //  int nacc = 0, natt = 0;
  std::vector<double> blk_avgs{};
  for (int block = 0; block < Nblk; block++) {
    double sum = 0.0;
    for (int smp = 0; smp < Nsmp; smp++) {
      /* generate new sample */
      do {
        x = x_new(x_old);
        p = prob(x);
        //  natt++;
      } while (rannyu() >= std::min(1.0, p / p_old));
      //  nacc++;
      x_old = x;
      p_old = p;

      /* add its contribution to the block average */
      sum += integrand(x);
    }
    blk_avgs.emplace_back(sum / Nsmp);
  }
  //  std::cerr << "Acceptance ratio: " << double(nacc) / natt << std::endl;

  return blk_avgs;
}

void Random::ex8_parameters_trajectory(double mu, double sigma, double T, int Nblk, int Nsmp) {
  std::ofstream log_stream("output/simulated_annealing_history.csv");
  log_stream << "TEMP,MU,SIGMA,E_EXP,E_ERR\n";

  /* generate new parameters */
  double delta = 0.1;
  auto new_parameters = [this, &delta](double &mu, double &sigma) {
    double theta = rannyu(0.0, 2 * pi), r = delta * std::sqrt(rannyu());
    mu = std::abs(mu + r * std::cos(theta));
    sigma = std::abs(sigma + r * std::sin(theta));
  };

  /* compute new energy */
  auto energy_with_error = [this, Nblk, Nsmp](double mu, double sigma) {
    std::vector blk_avgs = ex8_energy_block_avgs(mu, sigma, Nblk, Nsmp);
    double avg = 0.0;
    double avg2 = 0.0;
    for (auto E : blk_avgs) {
      avg += E;
      avg2 += E * E;
    }
    avg /= Nblk;
    avg2 /= Nblk;
    double err = std::sqrt(std::max(0.0, avg2 - avg * avg) / std::max(1, Nblk - 1));
    return std::pair<double, double>(avg, err);
  };

  /* compare new and old parameters results */
  auto accept_adjustment = [this, &T](double E, double E_old) {
    double prob_ratio = std::exp(-(E - E_old) / T);
    return rannyu() < std::min(1.0, prob_ratio);
  };

  /* cool down */
  int steps_per_temperature = 100;
  double alpha = 0.95; /*  TODO: decide */
  double T_min = 1e-2;
  auto cool_down = [alpha, &delta](double &T) { 
    T *= alpha; 
    delta *= std::sqrt(alpha);
  };

  /* SA LOOP */
  /* setup */
  auto [E, E_err] = energy_with_error(mu, sigma);
  log_stream << T << "," << mu << "," << sigma << "," << E << "," << E_err << std::endl;

  double E_old = E, E_err_old = E_err;
  double mu_old = mu, sigma_old = sigma;
  
  int current_step = 0;
  int acc_moves = 0, att_moves = 0;
  std::cerr << std::endl << "T = " << std::setprecision(1) << std::setw(5) << T 
            << "    mu = " << std::setprecision(2) << std::setw(5) << mu 
            << "    sigma = " << std::setprecision(2) << std::setw(5) << sigma 
            << "    E = " << std::setprecision(2) << std::setw(5) << E
            << std::endl;
  do {
    do {
      /* reset previous attempt */
      mu = mu_old, sigma = sigma_old;
      E = E_old, E_err = E_err_old;

      /* attemp new parameters */
      new_parameters(mu, sigma);
      auto [E_new, E_err_new] = energy_with_error(mu, sigma);
      E = E_new, E_err = E_err_new;

      /* update counter and print statusbar */
      att_moves++;
      std::cerr << "\racc: " << std::setw(5) << acc_moves << "\tatt: " << std::setw(5) << att_moves
                << "\tratio: " << double(acc_moves) / att_moves;
    } while (!accept_adjustment(E, E_old));
    /* confirm accepted move, log accepted parameters and their energy expectation */
    mu_old = mu, sigma_old = sigma;
    E_old = E, E_err_old = E_err;

    /* update counter and print statusbar */
    acc_moves++;
    std::cerr << "\racc: " << std::setw(5) << acc_moves << "\tatt: " << std::setw(5) << att_moves
              << "\tratio: " << double(acc_moves) / att_moves;

    /* if appropriate, lower temperature and zero counters */
    current_step++;
    if (current_step % steps_per_temperature == 0) {
      log_stream << T << "," << mu << "," << sigma << "," << E << "," << E_err << std::endl;

      cool_down(T);
      std::cerr << std::endl << "T = " << std::setprecision(2) << std::setw(5) << T 
                << "    mu = " << std::setprecision(2) << std::setw(5) << mu 
                << "    sigma = " << std::setprecision(2) << std::setw(5) << sigma 
                << "    E = " << std::setprecision(2) << std::setw(5) << E
                << std::endl;
      acc_moves = 0;
      att_moves = 0;
    }
  } while (T >= T_min);
}

std::vector<double> Random::ex8_sample_wavefunction(double mu, double sigma, int N) {
  double w_side = 1.0;                              // weight of each of the two Gaussians at +-mu
  double w_mid  = 2.0 * std::exp(-mu * mu / (sigma * sigma)); // weight of the Gaussian at 0
  double W = 2.0 * w_side + w_mid;

  auto sample = [this, w_side, W, mu, sigma]() {
    double u = rannyu() * W;
    double center;
    if (u < w_side) {
      center = +mu;
    }
    else if (u < 2.0 * w_side) {
      center = -mu;
    } else {
      center = 0.0;
    }
    return gauss(center, sigma / std::sqrt(2.0));
  };

  std::vector<double> samples{};
  for (int i = 0; i < N; i++) {
    samples.emplace_back(sample());
  }
  return samples;
}
