#include "system.hpp"

#include "particle.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>

// 1. Use create_system() to create a unique_ptr to a System object: reading
// from 'input/sys_settings.dat' it chooses the right class
// 2. In creating the pointer it calls each derived class' constructor
// 3. Each derived class constructor first calls the base class' constructor
// System() to read common class members from 'input/sys_settings.dat'
// 4. Then class-specific members are read from specific files.

/* ************************ */
/*  simulation loop methods */
/* ************************ */

std::unique_ptr<System> create_system() {
  std::cout << "  ================ Initializing .... ==============="
            << std::endl;

  std::ifstream settings_stream("input/sys_settings.dat");
  if (!settings_stream) {
    throw std::runtime_error(
        "Unable to open or find \'input/sys_settings.dat\' file.");
  }
  std::ofstream log_stream("output/logfile.log");
  if (!log_stream) {
    throw std::runtime_error(
        "Unable to open or find \'output/logfile.log\' file.");
  }

  std::string token;
  while (settings_stream >> token) {
    if (token != "MODEL") {
      continue;
    }

    settings_stream >> token;
    if (token == "ISING") {
      log_stream << "Ising 1D simulation, system created." << std::endl;
      return std::make_unique<IsingSystem>();
    } else if (token == "LJ") {
      log_stream << "Lennard-Jones 3D simulation, system created." << std::endl;
      return std::make_unique<LJSystem>();
    } else {
      log_stream
          << "Error: unrecognized physical system in 'input/sys_settings.dat'."
          << std::endl;
      throw std::runtime_error(
          "Error: unrecognized physical system in 'input/sys_settings.dat'.");
    }
  }

  log_stream << "Error: no line starting with 'MODEL:' has been found in "
                "'input/sys_settings.dat'."
             << std::endl;
  throw std::runtime_error("Error: no line starting with 'MODEL:' has been "
                           "found in 'input/sys_settings.dat'.");
}

void equilibrate(System &sys) {
  int nsteps_equil = sys.get_nsteps_equil(),
      nblocks_equil = sys.get_nblocks_equil();

  double T = sys.get_temp();
  std::cout << "  === Equilibrating: initial temperature T = " << std::setw(3)
            << T << " ===" << std::endl;

  if (sys.monitor.auto_equilibration) {
    /* automatic equlibration */
    int block = 0;
    do {
      print_auto_equil(-1, block, sys);
      for (int i = 0; i < nsteps_equil; i++) {
        sys.step();
        sys.measure_for_monitoring(); /* computes scalar measure and appends
        it
                                         to sys.monitor._inst_measures */
        if ((i + 1) % 10 == 0) {
          print_auto_equil(i, block, sys);
        }
      }
      sys.monitor.end_block(
          nsteps_equil); /* updates: _m, _m_err,
                          * _acceptance_ratio, _temp, _temp_err */

      sys.adjust_parameters(); /* updates sys: _delta, velocities... */
      // if (block >= nblocks_equil - 1) {
      //   nsteps_equil = sys.monitor.get_correct_nsteps();
      //   sys.set_nsteps_equil(nsteps_equil);
      // }

      print_auto_equil(nsteps_equil - 1, block, sys);
      block++;
      if (block == 32) {
        break;
      } // TODO: tbd
    } while (!(sys.monitor.is_acc_ratio_ok() &&
               sys.monitor.is_at_equilibrium() && sys.monitor.is_temp_ok()));

    std::cout << std::endl;

    T = sys.monitor.temp;
    // int nsteps_corr = sys.monitor.get_correct_nsteps();
    // std::cout << "Block size set at " << nsteps_corr " steps." << std::endl;
    // sys.set_nsteps(nsteps_corr);
  } else {
    /* manual equilibration */
    print_manual_equil(0, nsteps_equil * nblocks_equil);

    std::vector<double> temps;
    for (int block = 0; block < nblocks_equil; block++) {
      temps.clear();
      for (int step = 0; step < nsteps_equil; step++) {
        sys.step();
        print_manual_equil(step, nsteps_equil);
        temps.emplace_back(sys.measure_temp());
      }
    }
    T = 0.0;
    for (auto temp : temps) {
      T += temp;
    }
    T /= nsteps_equil;
    std::cout << std::endl;
  }

  sys.end_equilibration(T);
}

void simulate(System &sys) {
  /* Main simulation loop */
  std::cout << "  ================= Simulating ... ================="
            << std::endl;

  int nsteps = sys.get_nsteps(), nblocks = sys.get_nblocks();
  for (int block = 0; block < nblocks; block++) {
    print_progress_bar(0, nsteps, block, nblocks);
    for (int step = 0; step < nsteps; step++) {
      sys.step();
      sys.measure();
      print_progress_bar(step, nsteps, block, nblocks);
    }

    sys.end_block(block);
  }
  std::cout << std::endl;
}

void clean_up(System &sys) {
  /* Clean up */
  sys.clean_up();
  std::cout << "  =================== Simulation finished. ==================\n"
            << std::endl;
}

/**************************/
/* printing progress bars */
/**************************/

void print_auto_equil(int step, int block, System &sys) {
  std::cout.precision(2);
  int nsteps = sys.get_nsteps_equil();

  int bar_width = 20;
  int bar_fill = std::round(double(step * bar_width) / nsteps);
  std::cout << "\r" << std::setw(9) << std::right
            << "block " + std::to_string(block + 1) << " [";
  for (int i = 0; i < bar_width; i++) {
    if (i < bar_fill) {
      std::cout << "=";
    } else {
      std::cout << " ";
    }
  }
  std::cout << "]  " << std::setw(std::to_string(nsteps).size()) << step + 1
            << "/" << nsteps << " steps";

  if (sys.monitor.auto_mc_parameters && block > 0) {
    double acc_ratio = sys.monitor.get_acc_ratio();
    int acc_perc = std::round(acc_ratio * 100);
    std::cout << "  |  acceptance ratio: " << std::right << std::setprecision(3)
              << acc_perc << "%";
  }

  if (sys.monitor.auto_temp_correction && block > 0) {
    double T = sys.monitor.temp, T_err = sys.monitor.temp_err;
    std::cout << "  |  T: " << std::right << std::setprecision(4) << T
              << " \u00B1 " << std::left << std::setprecision(2) << T_err;
  }

  if (block > 2) {
    double m = sys.monitor.get_m(), m_err = sys.monitor.get_m_err();
    std::cout << "  |  m: " << std::right << std::setprecision(4) << m
              << " \u00B1 " << std::left << std::setprecision(2) << m_err;
  }

  std::cout << "\t";
  std::cout.flush();
}

void print_manual_equil(int step, int nsteps) {
  int bar_width = 50;
  double perc_done = double(step) / (nsteps - 1);
  int bar_fill = std::round(perc_done * bar_width);

  std::cout << "\r[";
  for (int i = 0; i < bar_width; i++) {
    if (i < bar_fill) {
      std::cout << "=";
    } else {
      std::cout << " ";
    }
  }
  std::cout << "]  " << std::setw(3) << std::right
            << std::round(perc_done * 100) << "%  "
            << std::setw(std::to_string(nsteps).size()) << step + 1 << "/"
            << nsteps << " steps";
}

void print_progress_bar(int step, int nsteps, int block, int nblocks) {
  int blocks_bar_width = 20, steps_bar_width = 30;
  double perc_blocks =
             double(block * nsteps + step) / std::max(nblocks * nsteps - 1, 1),
         perc_steps = double(step) / std::max(nsteps - 1, 1);
  int blocks_bar_fill = std::round(perc_blocks * blocks_bar_width),
      steps_bar_fill = std::round(perc_steps * steps_bar_width);

  std::cout << "\r  block " << std::setw(2 * std::to_string(nblocks).size() + 1)
            << std::left
            << std::to_string(block + 1) + "/" + std::to_string(nblocks) << "  "
            << std::setw(3) << std::right << int(std::round(perc_blocks * 100))
            << "% [";
  for (int i = 0; i < blocks_bar_width; i++) {
    if (i < blocks_bar_fill) {
      std::cout << "=";
    } else {
      std::cout << " ";
    }
  }
  std::cout << "][";
  for (int i = 0; i < steps_bar_width; i++) {
    if (i < steps_bar_fill) {
      std::cout << "=";
    } else {
      std::cout << " ";
    }
  }
  std::cout << "]  " << std::setw(3) << std::right
            << int(std::round(perc_steps * 100)) << "%  "
            << "step " << std::setw(2 * std::to_string(nsteps).size() + 1)
            << std::to_string(step + 1) + "/" + std::to_string(nsteps);
}

/* ********************** */
/* initialization methods */
/* ********************** */

/* **** System class initialization **** */

System::System() {
  _read_sys_settings();

  if (_restarted == false) {
    _rng.read_seed("input/seed_in.dat");
  } else {
    _rng.read_seed("output/last_seed.log");
  }
}

void System::_read_sys_settings() {
  /* TODO: initialize output/acceptance.log file */

  std::ifstream settings_stream("input/sys_settings.dat");
  if (!settings_stream) {
    throw std::runtime_error(
        "Unable to open or find \'input/sys_settings.dat\' file.");
  }
  std::ofstream log_stream("output/logfile.log", std::ios::app);
  if (!log_stream) {
    throw std::runtime_error(
        "Unable to open or find \'output/logfile.log\' file.");
  }

  std::string line;

  bool auto_equilibration = false;
  double acceptance_tolerance, slope_stdev_precision, temp_stdev_precision;

  while (std::getline(settings_stream, line)) {
    std::stringstream line_stream(line);
    std::string token;
    line_stream >> token;

    if (token == "NPART") {
      line_stream >> _npart;
      log_stream << "NPART: " << _npart << std::endl;
    } else if (token == "TEMP") {
      line_stream >> _temp;
      _beta = 1.0 / (k_B * _temp);
      log_stream << "TEMPERATURE: " << _temp << std::endl;
    } else if (token == "NSTEPS") {
      line_stream >> _nsteps;
      log_stream << "NSTEPS: " << _nsteps << std::endl;
    } else if (token == "NBLOCKS") {
      line_stream >> _nblocks;
      log_stream << "NBLOCKS: " << _nblocks << std::endl;
    } else if (token == "RESTARTED") {
      line_stream >> token;
      if (token == "0") {
        _restarted = false;
      } else {
        _restarted = true;
      }
    } else if (token == "TEMP_STDEV_PREC") {
      line_stream >> temp_stdev_precision;
      log_stream << "TEMP_STDEV_PREC: " << temp_stdev_precision << std::endl;
    } else if (token == "AUTO_EQUILIBRATION") { /* only if new run */
      if (!_restarted) {
        line_stream >> token;
        if (token == "0") {
          auto_equilibration = false;
        } else {
          auto_equilibration = true;
        }
      }
    } else if (token == "NSTEPS_EQUIL") { /* only if new run */
      if (!_restarted) {
        line_stream >> _nsteps_equil;
        log_stream << "NSTEPS_EQUIL: " << _nsteps_equil << std::endl;
      }
    } else if (token == "NBLOCKS_EQUIL") { /* only if new run */
      if (!_restarted) {
        line_stream >> _nblocks_equil;
        log_stream << "NBLOCKS_EQUIL: " << _nblocks_equil << std::endl;
      }
    } else if (token == "ACC_RATIO_TOL") { /* only if new run and equilibration
                                               is automatic */
      if (!_restarted && auto_equilibration) {
        line_stream >> acceptance_tolerance;
        log_stream << "ACC_RATIO_TOL: " << acceptance_tolerance << std::endl;
      }
    } else if (token == "SLOPE_STDEV_PREC") { /* only if new run and
                                                equilibration is automatic */
      if (!_restarted && auto_equilibration) {
        line_stream >> slope_stdev_precision;
        log_stream << "SLOPE_STDEV_PREC: " << slope_stdev_precision
                   << std::endl;
      }
    } else if (token == "ENDINPUT") {
      break;
    }
  }

  monitor = Monitor(_nblocks_equil, acceptance_tolerance, auto_equilibration,
                    slope_stdev_precision, _temp, temp_stdev_precision);
  monitor.temp = _temp;
}

/* *********************** */
/* simulation loop methods */
/* *********************** */

int System::get_nsteps() const { return _nsteps; }
int System::get_nblocks() const { return _nblocks; }
int System::get_nsteps_equil() const { return _nsteps_equil; }
int System::get_nblocks_equil() const { return _nblocks_equil; }
void System::set_nsteps(int nsteps) { _nsteps = nsteps; }
void System::set_nsteps_equil(int nsteps_equil) {
  _nsteps_equil = nsteps_equil;
}

bool System::get_restarted() const { return _restarted; }

double System::get_temp() const { return _temp; }
void System::set_temp(double temp) {
  _temp = temp;
  _beta = 1 / temp;
}

/* ************************ */
/* finish & cleanup methods */
/* ************************ */

void System::end_block(int block) {
  for (auto &measurable : _measurables) {
    measurable.end_block(block, _nsteps);
  }
  /* TODO: log ratio accepted/attempted moves */
}

void System::clean_up() {
  _rng.log_last_seed();
  _log_last_config();
}
