#include "system.hpp"

#include "particle.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>

/* ******************** */
/* CONSTRUCTION METHODS */
/* ******************** */

LJSystem::LJSystem() : System() {
  _read_lj_settings();
  _load_init_config();
  _define_measurables();
}

void LJSystem::_read_lj_settings() {
  std::ifstream settings_stream("input/lj_settings.dat");
  if (!settings_stream) {
    throw std::runtime_error(
        "Unable to open or find \'input/lj_settings.dat\' file.");
  }
  std::ofstream log_stream("output/logfile.log", std::ios::app);
  if (!log_stream) {
    throw std::runtime_error(
        "Unable to open or find \'output/logfile.log\' file.");
  }

  std::string line;
  while (std::getline(settings_stream, line)) {
    std::stringstream line_stream(line);
    std::string token;
    line_stream >> token;

    if (token == "ALGORITHM") {
      line_stream >> token;
      if (token == "VERLET") {
        log_stream << "ALGORITHM: VERLET" << std::endl;
        step = [this]() { _step_verlet(); };
        monitor.auto_temp_correction = true;
      } else if (token == "METROPOLIS") {
        log_stream << "ALGORITHM: METROPOLIS" << std::endl;
        step = [this]() { _step_metropolis(); };
        monitor.auto_mc_parameters = true;
      } else {
        throw std::runtime_error(
            "Error: unrecognized algorithm in \'input/lj_settings.dat\'.");
      }
    } else if (token == "NDIM") {
      line_stream >> _ndim;
      log_stream << "NDIM: " << _ndim << std::endl;
    } else if (token == "RHO") {
      line_stream >> _rho;
      _volume = _npart / _rho;
      _box_dim = std::pow(_volume, 1.0 / 3.0);
      log_stream << "RHO: " << _rho << "\n"
                 << "VOLUME: " << _volume << "\n"
                 << "BOX_DIM: " << _box_dim << std::endl;
    } else if (token == "RCUT") {
      line_stream >> _r_cut;
      log_stream << "RCUT:" << _r_cut << std::endl;
      if (_r_cut >= _box_dim / 2.0) {
        throw std::runtime_error("Invalid parameters: _r_cut (" +
                                 std::to_string(_r_cut) +
                                 ") must be smaller than half box dimension (" +
                                 std::to_string(_box_dim / 2.0) + ")!");
      }

      _vtail =
          8 * pi * _rho * _epsilon *
          (1 / (9 * std::pow(_r_cut, 9.0)) - 1 / (3 * std::pow(_r_cut, 3.0)));
      _ptail =
          32 * pi * _rho *
          (1 / (9 * std::pow(_r_cut, 9.0)) - 1 / (6 * std::pow(_r_cut, 3.0)));
    } else if (token == "NBINS_GOFR") {
      line_stream >> _n_bins_gofr;
      _bin_size_gofr = _box_dim / (2.0 * _n_bins_gofr);
      log_stream << "NBINS_GOFR: " << _n_bins_gofr << std::endl;
    } else if (token == "NBINS_POFV") {
      line_stream >> _n_bins_pofv;
      _v_max_pofv = 5.0 * std::sqrt(_temp / _m);
      _bin_size_pofv = _v_max_pofv / _n_bins_pofv;
      log_stream << "NBINS_POFV: " << _n_bins_pofv << std::endl;
    } else if (token == "DELTA") {
      line_stream >> _delta;
      log_stream << "DELTA: " << _delta << std::endl;
    } else if (token == "DT") {
      line_stream >> _dt;
      log_stream << "DT: " << _dt << std::endl;
    } else if (token == "ENDSETTINGS") {
      log_stream << "Lennard-Jones simulation settings reading complete."
                 << std::endl;
      break;
    }
  }
}

void LJSystem::_load_init_config() {
  _particles.set_size(_npart);

  if (!_restarted) {
    if (_ndim != 3) {
      throw std::runtime_error("FCC lattice is only defined for _ndim = 3.");
    }
    int ncells = static_cast<int>(std::round(std::cbrt(_npart / 4.0)));
    if (4 * ncells * ncells * ncells != _npart) {
      throw std::runtime_error(
          "FCC lattice requires NPART = 4 n^3 (32, 108, 256, 500, ...).");
    }

    double box_fill_perc = 1.0;
    const double a = box_fill_perc * _box_dim / ncells;
    const double basis[4][3] = {
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};

    std::vector<arma::vec> sites;
    sites.reserve(_npart);
    for (int ix = 0; ix < ncells; ix++) {
      for (int iy = 0; iy < ncells; iy++) {
        for (int iz = 0; iz < ncells; iz++) {
          for (const auto &b : basis) {
            arma::vec pos = {(ix + b[0]) * a - _box_dim / 2.0,
                             (iy + b[1]) * a - _box_dim / 2.0,
                             (iz + b[2]) * a - _box_dim / 2.0};
            sites.push_back(pos);
          }
        }
      }
    }

    double v_T = std::sqrt(_ndim * _temp / _m);
    arma::vec v_av(_ndim, arma::fill::zeros);

    std::size_t i = 0;
    for (auto &p : _particles) {
      p.set_x(sites[i++]);

      /* speed is v_T, along a random axis */
      arma::vec v(_ndim, arma::fill::zeros);
      int rand = int(_rng.rannyu() * 6.0);
      v(rand / 2) = (rand % 2 - 0.5) * 2.0 * v_T;
      p.set_v(v);
      v_av += v;
    }

    v_av /= _npart;
    double v_av2 = 0.0;
    for (auto &p : _particles) {
      p.set_v(p.get_v() - v_av);
      for (double v_component : p.get_v()) {
        v_av2 += v_component * v_component;
      }
    }
    v_av2 /= _npart;

    double v_factor = std::sqrt(_ndim * _temp / (v_av2 * _m));
    for (auto &p : _particles) {
      p.set_v(p.get_v() * v_factor);
    }
  } else {
    std::ifstream init_config_stream("output/last_config.log");
    if (!init_config_stream) {
      throw std::runtime_error(
          "Unable to open or find \'output/last_config.log\' file.");
    }

    std::ifstream init_param_stream("output/logfile.log");
    if (!init_param_stream) {
      throw std::runtime_error(
          "Unable to open or find \'output/logfile.log\' file.");
    }

    std::string line;
    for (auto &p : _particles) {
      std::getline(init_config_stream, line);
      std::stringstream line_stream(line);

      arma::vec x(_ndim), v(_ndim);
      for (int i = 0; i < _ndim; i++) {
        line_stream >> x(i);
      }
      for (int i = 0; i < _ndim; i++) {
        line_stream >> v(i);
      }

      p.set_x(x);
      p.set_v(v);
    }

    while (std::getline(init_param_stream, line)) {
      std::stringstream line_stream(line);
      std::string token;
      line_stream >> token;

      if (token == "TEMP_EQUIL:") {
        line_stream >> _temp;
        _beta = 1 / _temp;
      } else if (token == "NSTEPS_EQUIL:") {
        line_stream << _nsteps;
      } else if (token == "MC_DELTA:") {
        line_stream << _delta;
      }
    }
  }

  for (int i = 0; i < _npart; i++) {
    arma::vec x_old = _particles(i).get_x() - _dt * _particles(i).get_v();
    _particles(i).set_x_old(x_old);

    arma::vec a = _force(i) / _m;
    _particles(i).set_a(a);
  }
}

void LJSystem::_define_measurables() {
  /* scalar measurable quantities */
  _measurables.emplace_back("kinetic_energy_per_particle", /* computes K/n at
              each step, blocks estimate its average */
                            [](std::span<const double> args) {
                              double k = args[0];
                              return k;
                            });
  _measurables.emplace_back("potential_energy_per_particle", /* computes U/n at
              each step, blocks estimate its average */
                            [](std::span<const double> args) {
                              double u = args[1];
                              return u;
                            });
  _measurables.emplace_back("energy_per_particle", /* computes E/n at each
        step, blocks estimate its average */
                            [](std::span<const double> args) {
                              double k = args[0], u = args[1];
                              return k + u;
                            });
  _measurables.emplace_back("pressure", /* computes pressure at each step,
                                             block estimate is its average */
                            [this](std::span<const double> args) {
                              double k = args[0], virial = args[2];
                              return (2.0 / 3.0) * _rho * k +
                                     (_ptail * _npart + 4.0 * virial) / _volume;
                            });

  /* measurable distributions */
  _gofr.emplace(
      "gofr", _n_bins_gofr, [this](std::span<const double> distances) {
        std::vector<double> hist(_n_bins_gofr, 0.0);
        for (double r : distances) {
          int bin = static_cast<int>(r / _bin_size_gofr);
          if (bin >= 0 && bin < _n_bins_gofr) {
            double r_lo = bin * _bin_size_gofr, r_hi = r_lo + _bin_size_gofr;
            double shell_volume = 4 * std::numbers::pi *
                                  (std::pow(r_hi, 3) - std::pow(r_lo, 3)) / 3.0;
            hist[bin] += 2.0 * _volume / (shell_volume * _npart * (_npart - 1));
          }
        }
        return hist;
      });
  _pofv.emplace("pofv", _n_bins_pofv, [this](std::span<const double> speeds) {
    std::vector<double> hist(_n_bins_pofv, 0.0);
    for (double v : speeds) {
      int bin = static_cast<int>(v / _bin_size_pofv);
      if (bin >= 0 && bin < _n_bins_pofv) {
        hist[bin] += 1.0 / (_npart * _bin_size_pofv);
      }
    }
    return hist;
  });
}

/* ******************* */
/* BOUNDARY CONDITIONS */
/* ******************* */

arma::vec LJSystem::_pbc(arma::vec x) const {
  double half_box = _box_dim / 2;
  for (int i = 0; i < _ndim; i++) {
    x(i) -= std::floor((x(i) + half_box) / _box_dim) * _box_dim;
  }
  return x;
}

double LJSystem::_nearest_image_dist(arma::vec a, arma::vec b) const {
  return mod(_pbc(a - b));
}

/* ******************* */
/* ENERGIES AND FORCES */
/* ******************* */

double LJSystem::_kin_energy() const {
  double v_sum_sq = 0.0;
  for (auto &p : _particles) {
    double v = mod(p.get_v());
    v_sum_sq += v * v;
  }
  return 0.5 * v_sum_sq * _m;
}

double LJSystem::_pot_energy(int i, int j) const {
  double r = _nearest_image_dist(_particles(i).get_x(), _particles(j).get_x());
  if (i == j || r >= _r_cut) {
    return 0.0;
  } else {
    double sr6 = std::pow(_sigma / r, 6.0);
    double sr12 = sr6 * sr6;
    return 4 * _epsilon * (sr12 - sr6);
  }
}

double LJSystem::_pot_energy(int i) const {
  double U = 0.0;
  for (int j = 0; j < _npart; j++) {
    U += _pot_energy(i, j);
  }
  return U;
}

double LJSystem::_pot_energy() const {
  double U = 0.0;
  for (int i = 0; i < _npart; i++) {
    for (int j = i + 1; j < _npart; j++) {
      U += _pot_energy(i, j);
    }
  }
  return U;
}

arma::vec LJSystem::_force(int i) const {
  if (i < 0 || i >= _npart) {
    throw std::runtime_error("Error: computing total force on a "
                             "non-existent particle at index " +
                             std::to_string(i) + ".");
  }

  arma::vec F(_ndim, arma::fill::zeros);
  for (int j = 0; j < _npart; j++) {
    if (j == i) {
      continue;
    }
    arma::vec r_versor = _pbc(_particles(i).get_x() - _particles(j).get_x());
    double r = mod(r_versor);
    r_versor /= r;
    if (r >= _r_cut) {
      continue;
    }
    double sr6 = std::pow(_sigma / r, 6.0);
    double sr12 = sr6 * sr6;
    F += (24 * _epsilon * (2 * sr12 - sr6) / r) * r_versor;
  }

  return F;
}

double LJSystem::measure_temp() const {
  if (!monitor.auto_temp_correction) {
    return _temp;
  }
  return 2.0 * _kin_energy() / (k_B * _ndim * (_npart - 1));
};

/* *********************** */
/* SIMULATION LOOP METHODS */
/* *********************** */

void LJSystem::_step_verlet() {
  for (int i = 0; i < _npart; i++) {
    arma::vec x_new = _particles(i).get_x() + _dt * _particles(i).get_v() +
                      0.5 * _dt * _dt * _particles(i).get_a();
    x_new = _pbc(x_new);
    _particles(i).set_x(x_new);
  }
  arma::field<arma::vec> a_new(_npart);
  for (int i = 0; i < _npart; i++) {
    a_new(i) = _force(i) / _m;
  }
  for (int i = 0; i < _npart; i++) {
    arma::vec v_new =
        _particles(i).get_v() + 0.5 * (_particles(i).get_a() + a_new(i)) * _dt;
    _particles(i).set_v(v_new);
    _particles(i).set_a(a_new(i));
  }
}

void LJSystem::_step_metropolis() {
  /* TODO: add randomization of velocities? */
  for (int i = 0; i < _npart; i++) {
    int rnd_index = static_cast<int>(_rng.rannyu() * _npart);
    if (rnd_index < 0 || rnd_index >= _npart) {
      throw std::runtime_error(
          "A particle index has been generated out of bounds:" +
          std::to_string(rnd_index) + ".");
    }

    _particles(rnd_index).set_x_old(_particles(rnd_index).get_x());
    double E_old = _pot_energy(rnd_index);

    /* Trial transition distribution: uniform in [-_delta, _delta]^3 *
     * box for positions, VELOCITIES DO NOT GET CHANGED              */
    arma::vec delta_x(_ndim);
    for (int d = 0; d < _ndim; d++) {
      delta_x(d) = _rng.rannyu(-1.0, 1.0) * _delta;
    }
    arma::vec new_x = _pbc(_particles(rnd_index).get_x() + delta_x);
    _particles(rnd_index).set_x(new_x);

    /* Acceptance probability: non-normalized Boltzmann distribution  */
    double E_new = _pot_energy(rnd_index);
    double acc_prob = std::exp(-_beta * (E_new - E_old));

    if (_rng.rannyu() < acc_prob) {
      _particles(rnd_index).set_x_old(_particles(rnd_index).get_x());
      monitor.naccepted_up();
    } else {
      _particles(rnd_index).set_x(_particles(rnd_index).get_x_old());
    }

    monitor.nattempted_up();
  }
}

void LJSystem::measure() {
  double k = 0.0, u = 0.0, virial = 0.0;
  std::vector<double> distances, speeds;

  for (int i = 0; i < _npart; i++) {
    double v = mod(_particles(i).get_v());
    speeds.emplace_back(v);
    k += 0.5 * _m * v * v;

    for (int j = i + 1; j < _npart; j++) {
      double r =
          _nearest_image_dist(_particles(i).get_x(), _particles(j).get_x());
      distances.emplace_back(r);
      if (r < _r_cut) {
        double sr6 = std::pow(_sigma / r, 6.0);
        double sr12 = sr6 * sr6;
        u += 4.0 * _epsilon * (sr12 - sr6);
        virial += 4.0 * _epsilon * (sr12 - 0.5 * sr6);
      }
    }
  }

  k /= _npart;
  u = u / _npart + _vtail;

  std::array<double, 3> args{k, u, virial};
  for (auto &measurable : _measurables) {
    measurable.measure(args);
  }
  _gofr->measure(distances);
  _pofv->measure(speeds);
}

void LJSystem::adjust_parameters() {
  if (monitor.auto_mc_parameters && _delta <= _box_dim / 4.0) {
    double factor = monitor.get_acc_ratio() / 0.5;
    _delta *= std::max(std::min(factor, 2.0), 0.5);
  }
  if (monitor.auto_temp_correction) {
    double factor = std::sqrt(monitor.temp_target / monitor.temp);
    for (auto &p : _particles) {
      p.set_v(factor * p.get_v());
    }
    _temp = monitor.temp;
    _beta = 1.0 / _temp;
  }
}

void LJSystem::measure_for_monitoring() {
  double u = _pot_energy() / _npart;
  double T = measure_temp();
  monitor.save_measure(u, T);
}

void LJSystem::end_block(int block) {
  System::end_block(block);

  _gofr->end_block(block, _nsteps);
  _pofv->end_block(block, _nsteps);
}

void LJSystem::end_equilibration(double T) {
  set_temp(T);

  std::ofstream last_param_stream;
  last_param_stream.open("output/logfile.log", std::ios_base::app);
  last_param_stream << "TEMP_EQUIL: " << std::to_string(_temp) << "\n"
                    << "NSTEPS_EQUIL: " << std::to_string(_nsteps) << std::endl;
  if (!monitor.auto_equilibration) {
    return;
  }
  if (monitor.auto_mc_parameters) {
    last_param_stream << "MC_DELTA: " << std::to_string(_delta) << std::endl;
  }
}

/* **************** */
/* LOG AND CLEAN UP */
/* **************** */

void LJSystem::_log_last_config() const {
  std::ofstream last_config_stream("output/last_config.log");
  if (!last_config_stream) {
    std::cerr << "WARNING: Unable to log last configuration in "
                 "\'output/last_config.log\' file."
              << std::endl;
  }

  for (auto p : _particles) {
    for (int i = 0; i < _ndim; i++) {
      last_config_stream << std::to_string(p.get_x()(i)) << " ";
    }
    for (int i = 0; i < _ndim; i++) {
      last_config_stream << std::to_string(p.get_v()(i)) << " ";
    }
    last_config_stream << std::endl;
  }
}

void LJSystem::clean_up() {
  _gofr->finalize(_nblocks);
  _pofv->finalize(_nblocks);
  System::clean_up();
}
