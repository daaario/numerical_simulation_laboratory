#include "system.hpp"

#include "particle.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>

IsingSystem::IsingSystem() : System() {
  _read_ising_settings();
  _load_init_config();
  _define_measurables();
}

void IsingSystem::_read_ising_settings() {
  std::ifstream settings_stream("input/ising_settings.dat");
  if (!settings_stream) {
    throw std::runtime_error(
        "Unable to open or find \'input/ising_settings.dat\' file.");
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
      if (token == "GIBBS") {
        log_stream << "ALGORITHM: GIBBS" << std::endl;
        step = [this]() { _step_ising(); };
        _attempt_flip = [this](int k) { return _flip_gibbs(k); };
        monitor.auto_mc_parameters = false;
        monitor.auto_temp_correction = false; /* overrides settings file */
      } else if (token == "METROPOLIS") {
        log_stream << "ALGORITHM: METROPOLIS" << std::endl;
        step = [this]() { _step_ising(); };
        _attempt_flip = [this](int k) { return _flip_metropolis(k); };
        monitor.auto_mc_parameters = false;
        monitor.auto_temp_correction = false; /* overrides settings file */
      } else {
        throw std::runtime_error(
            "Error: unrecognized algorithm in \'input/ising_settings.dat\'.");
      }
    } else if (token == "J") {
      line_stream >> _J;
      log_stream << "J: " << _J << std::endl;
    } else if (token == "H") {
      line_stream >> _H;
      log_stream << "H: " << _H << std::endl;
    } else if (token == "ENDINPUT") {
      log_stream << "Ising simulation settings reading complete." << std::endl;
      break;
    }
  }
}

void IsingSystem::_load_init_config() {
  _particles.set_size(_npart);

  if (!_restarted) {
    for (int i = 0; i < _npart; i++) {
      if (_rng.rannyu() > 0.5) {
        _particles(i).flip_spin();
      }
    }
    return;
  }

  std::ifstream init_config_stream("output/last_config.log");
  if (!init_config_stream) {
    throw std::runtime_error(
        "Unable to open or find \'output/last_config.log\' file.");
  }

  std::string line;
  std::getline(init_config_stream, line);
  if (line != "MODEL: ISING") {
    throw std::runtime_error("Error: simulation model does not coincide with "
                             "'output/last_config.log' file.");
  }
  std::getline(init_config_stream, line);
  if (line != ("NPART: " + std::to_string(_npart))) {
    throw std::runtime_error("Particle number of last logged configuration "
                             "does not coincide with this simulation.");
  }

  int spin;
  for (int i = 0; i < _npart; i++) {
    init_config_stream >> spin;
    _particles(i).set_spin(spin);
  }
}

void IsingSystem::_define_measurables() {
  _measurables.emplace_back("energy_per_particle", /* compute energy
               per particle at each step, block estimate is its average */
                            [](std::span<const double> args) {
                              double e = args[0];
                              return e;
                            });
  _measurables.emplace_back("magnetization", /* compute magnetization at each
                 step, block estimate is its average */
                            [](std::span<const double> args) {
                              double m = args[1];
                              return m;
                            });
  _measurables.emplace_back(
      "specific_heat", /* compute energy at each step, block estimate is
its variance */
      [](std::span<const double> args) {
        double e = args[0];
        return e;
      },
      [this](double sum, double sum_sq, int nsteps) {
        double average = sum / nsteps, average_sq = sum_sq / nsteps;
        return this->_beta * this->_beta * (average_sq - average * average) *
               this->_npart;
      });
  _measurables.emplace_back(
      "susceptibility", /* computes magnetization at each step, block
estimate is its variance */
      [](std::span<const double> args) {
        double m = args[1];
        return m;
      },
      [this](double sum, double sum_sq, int nsteps) {
        double average = sum / nsteps, average_sq = sum_sq / nsteps;
        return this->_beta * (average_sq - average * average) * this->_npart;
      });
}

bool IsingSystem::_flip_metropolis(int k) {
  double delta_E = 2.0 * _particles(k).get_spin() *
                   (_J * (_particles(_pbc(k + 1)).get_spin() +
                          _particles(_pbc(k - 1)).get_spin()) +
                    _H);
  return _rng.rannyu() < exp(-_beta * delta_E);
}

bool IsingSystem::_flip_gibbs(int k) {
  double prob_spinup =
      1 / (1 + std::exp(-2 * _beta *
                        (_J * (_particles(_pbc(k - 1)).get_spin() +
                               _particles(_pbc(k + 1)).get_spin()) +
                         _H)));

  int target_spin;
  if (_rng.rannyu() < prob_spinup) {
    target_spin = 1;
  } else {
    target_spin = -1;
  }
  return _particles(k).get_spin() != target_spin;
}

void IsingSystem::_step_ising() {
  for (int i = 0; i < _npart; i++) {
    int chosen = static_cast<int>(_rng.rannyu() * _npart);
    if (chosen < 0 || chosen >= _npart) {
      throw std::runtime_error(
          "A particle index has been generated out of bounds:" +
          std::to_string(chosen) + ".");
    }

    if (_attempt_flip(chosen)) {
      _particles(chosen).flip_spin();
      monitor.naccepted_up();
    }
    monitor.nattempted_up();
  }
}

void IsingSystem::adjust_parameters() {
  if (monitor.auto_mc_parameters) {
    /*  TODO */
    return;
  }
  if (monitor.auto_temp_correction) {
    /* TODO */
    return;
  }
}

double IsingSystem::measure_temp() const {
  if (!monitor.auto_temp_correction) {
    return _temp;
  }
  /* TODO: define!! */
  return _temp;
}

void IsingSystem::measure_for_monitoring() {
  double spin_sum_nnpairs = 0.0, spin_sum = 0.0;
  for (int i = 0; i < _npart; i++) {
    spin_sum_nnpairs +=
        _particles(i).get_spin() * _particles(_pbc(i + 1)).get_spin();
    spin_sum += _particles(i).get_spin();
  }
  double E = (-_J * spin_sum_nnpairs - _H * spin_sum);

  double T = measure_temp();
  monitor.save_measure(E, T);
}

void IsingSystem::end_equilibration(double T) {
  set_temp(T);

  std::ofstream last_param_stream;
  last_param_stream.open("output/logfile.log", std::ios_base::app);
  last_param_stream << "TEMP_EQUIL: " << std::to_string(_temp) << "\n"
                    << "NSTEPS_EQUIL: " << std::to_string(_nsteps) << std::endl;
  if (!monitor.auto_equilibration) {
    return;
  }
  if (monitor.auto_mc_parameters) {
    /* ?? parameters to set in 'adjust_mc_parameters' */
  }
}

int IsingSystem::_pbc(int index) const {
  if (index >= _npart) {
    return index - _npart;
  } else if (index < 0) {
    return index + _npart;
  } else {
    return index;
  }
}

/* ***************** */
/* measuring methods */
/* ***************** */

void IsingSystem::measure() {
  double spin_sum_nnpairs = 0.0, spin_sum = 0.0;
  for (int i = 0; i < _npart; i++) {
    spin_sum_nnpairs +=
        _particles(i).get_spin() * _particles(_pbc(i + 1)).get_spin();
    spin_sum += _particles(i).get_spin();
  }
  double e = (-_J * spin_sum_nnpairs - _H * spin_sum) / _npart;
  double m = spin_sum / _npart;

  std::array<double, 2> args{e, m};
  for (auto &measurable : _measurables) {
    measurable.measure(args);
  }
}

void IsingSystem::_log_last_config() const {
  std::ofstream last_config_stream("output/last_config.log");
  if (!last_config_stream) {
    std::cerr << "WARNING: Unable to log last configuration in "
                 "\'output/last_config.log\' file."
              << std::endl;
  }
  last_config_stream << "MODEL: ISING" << std::endl;
  last_config_stream << "NPART: " << std::to_string(_npart) << std::endl;
  for (auto p : _particles) {
    last_config_stream << std::to_string(p.get_spin()) << " ";
  }
}
