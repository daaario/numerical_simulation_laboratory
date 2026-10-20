#include "measurable.hpp"

#include <fstream>
#include <iostream>

namespace {
double block_estimate_default(double sum, double /*sum_sq*/, double n_data) {
  return sum / n_data;
}
} // namespace

/*********************************************************/
/**************** scalar measurable class ****************/
/*********************************************************/
Measurable::Measurable(
    std::string name, std::function<double(std::span<const double>)> compute,
    std::function<double(double, double, int)> block_estimate)
    : _name(name), _output_stream("output/" + name + ".csv"),
      _compute(std::move(compute)), _block_estimate(std::move(block_estimate)) {
  if (!_output_stream) {
    throw std::runtime_error("Unable to open \'output/" + _name +
                             ".csv\' file.");
  }
  _output_stream << "BLOCK_N,BLOCK_AV,BLOCK_ERR,PROGR_AV,PROGR_ERR"
                 << std::endl;
}

Measurable::Measurable(std::string name,
                       std::function<double(std::span<const double>)> compute)
    : Measurable(name, compute, block_estimate_default) {}

void Measurable::measure(std::span<const double> args) {
  double X = _compute(args);
  _block_sum += X;
  _block_sum_sq += X * X;
}

void Measurable::end_block(int block_number, int nsteps) {
  double block_estimate = _block_estimate(_block_sum, _block_sum_sq, nsteps);
  double block_err = error(_block_sum, _block_sum_sq, nsteps);
  _global_av += block_estimate;
  _global_av_sq += block_estimate * block_estimate;

  _output_stream << std::to_string(block_number + 1) + "," +
                        std::to_string(block_estimate) + "," +
                        std::to_string(block_err) + "," +
                        std::to_string(_global_av / (block_number + 1)) + "," +
                        std::to_string(
                            error(_global_av, _global_av_sq, block_number + 1))
                 << std::endl;

  _block_sum = 0.0;
  _block_sum_sq = 0.0;
}

/****************************************************************/
/**************** measurable distributions class ****************/
/****************************************************************/
MeasurableHistogram::MeasurableHistogram(
    std::string name, int n_bins,
    std::function<std::vector<double>(std::span<const double>)> compute,
    std::function<double(double, double, int)> block_estimate)
    : _name(name), _n_bins(n_bins), _block_sum(n_bins, 0.0),
      _block_sum_sq(n_bins, 0.0), _global_av(n_bins, 0.0),
      _global_av_sq(n_bins, 0.0), _compute(std::move(compute)),
      _block_estimate(std::move(block_estimate)),
      _output_global_stream("output/" + name + "_global.csv"),
      _output_blocks_stream("output/" + name + "_blocks.csv") {
  if (!_output_global_stream) {
    throw std::runtime_error("Unable to open \'output/" + _name +
                             "_global.csv\' file.");
  }
  if (!_output_blocks_stream) {
    throw std::runtime_error("Unable to open \'output/" + _name +
                             "_blocks.csv\' file.");
  }
  _output_blocks_stream << "BLOCK,";
  for (int bin = 0; bin < _n_bins; bin++) {
    _output_blocks_stream << "BIN_" << bin + 1;
    if (bin < _n_bins - 1) {
      _output_blocks_stream << ",";
    }
  }
  _output_blocks_stream << std::endl;
  _output_global_stream << "GLOBAL_AV,ERROR" << std::endl;
}
MeasurableHistogram::MeasurableHistogram(
    std::string name, int n_bins,
    std::function<std::vector<double>(std::span<const double>)> compute)
    : MeasurableHistogram(name, n_bins, compute, block_estimate_default) {}

void MeasurableHistogram::measure(std::span<const double> args) {
  std::vector<double> hist = _compute(args);
  for (int bin = 0; bin < _n_bins; bin++) {
    _block_sum[bin] += hist[bin];
    _block_sum_sq[bin] += hist[bin] * hist[bin];
  }
}

void MeasurableHistogram::end_block(int block_number, int n_samples) {
  _output_blocks_stream << block_number + 1 << ",";
  for (int bin = 0; bin < _n_bins; bin++) {
    double bin_block_estimate =
        _block_estimate(_block_sum[bin], _block_sum_sq[bin], n_samples);
    _global_av[bin] += bin_block_estimate;
    _global_av_sq[bin] += bin_block_estimate * bin_block_estimate;

    if (bin < _n_bins - 1) {
      _output_blocks_stream << bin_block_estimate << ",";
    } else {
      _output_blocks_stream << bin_block_estimate << std::endl;
    }

    _block_sum[bin] = 0;
    _block_sum_sq[bin] = 0;
  }
}

void MeasurableHistogram::finalize(int n_blocks) {
  for (int bin = 0; bin < _n_bins; bin++) {
    double bin_global_av = _global_av[bin] / n_blocks;
    double err = error(_global_av[bin], _global_av_sq[bin], n_blocks);
    _output_global_stream << bin_global_av << "," << err << std::endl;
  }
}

/**********************************************/
/**************** common tools ****************/
/**********************************************/
double error(double sum, double sum_sq, int nblocks) {
  if (nblocks < 2) {
    return 0.0;
  }
  return std::sqrt(
      std::max(0.0, sum_sq / nblocks - std::pow(sum / nblocks, 2)) /
      (nblocks - 1));
}

/**************************************************/
/********** Thermalisation monitor class **********/
/**************************************************/

Monitor::Monitor(int blocks_window, double acceptance_tolerance,
                 bool auto_equil, double slope_stdev_precision,
                 double temp_trgt, double temp_stdev_precision)
    : _blocks_window(blocks_window),
      _acceptance_tolerance(acceptance_tolerance),
      _slope_stdev_precision(slope_stdev_precision),
      auto_equilibration(auto_equil), temp_target(temp_trgt),
      temp_stdev_precision(temp_stdev_precision) {}

/* delta corrections for LJ Metropolis calculations */

void Monitor::naccepted_up() { _naccepted++; }
void Monitor::nattempted_up() { _nattempted++; }
double Monitor::get_acc_ratio() const { return _acceptance_ratio; }
double Monitor::get_acc_target() const { return _acceptance_target; }

bool Monitor::is_acc_ratio_ok() const {
  return !auto_mc_parameters ||
         (std::abs(_acceptance_ratio - _acceptance_target) <
          _acceptance_tolerance);
}

/* measurable monitoring for
 * checking equilibrium */

void Monitor::save_measure(double value, double temp) {
  _inst_measures.emplace_back(value);
  _temps_block.emplace_back(temp);
}

std::pair<double, double> Monitor::_lin_reg_slope() const {
  int nblocks = _block_avgs.size();
  if (nblocks < 3) {
    return {0.0, 0.0};
  }

  double S = 0.0, S_x = 0.0, S_y = 0.0, S_xx = 0.0, S_xy = 0.0;
  for (int block = 0; block < nblocks; block++) {
    double w_i = 1 / std::max(_block_vars[block], 1e-8),
           y_i = _block_avgs[block];
    S += w_i;
    S_x += w_i * block;
    S_y += w_i * y_i;
    S_xx += w_i * block * block;
    S_xy += w_i * block * y_i;
  }
  double delta = S * S_xx - S_x * S_x;

  double m = (S * S_xy - S_x * S_y) / delta;
  double m_err = std::sqrt(S / delta);

  return {m, m_err};
}

double Monitor::get_m() const { return _m; }
double Monitor::get_m_err() const { return _m_err; }

bool Monitor::is_at_equilibrium() const {
  int nblocks = _block_avgs.size();
  if (nblocks < _blocks_window) {
    return false;
  }
  return !auto_equilibration ||
         (std::abs(_m / _m_err) <= _slope_stdev_precision);
}

/* equilibrium detection */

void Monitor::end_block(int nsteps) {
  /* acceptance ratio */
  if (auto_mc_parameters && _nattempted != 0) {
    _acceptance_ratio = static_cast<double>(_naccepted) / _nattempted;
  }
  _nattempted = 0, _naccepted = 0;

  /* stability of a measure */
  double block_avg = 0.0, block_avg2 = 0.0;
  double nmeasures = _inst_measures.size();
  if (nmeasures < nsteps) {
    throw std::runtime_error(
        "Error: \'Monitor::end_block()\' has been called while "
        "\'_inst_measures\' contains less than \'nsteps\' values.");
  }

  for (int i = 0; i < nsteps; i++) {
    double x = _inst_measures[nmeasures - 1 - i];
    block_avg += x;
    block_avg2 += x * x;
  }
  block_avg /= nsteps;
  block_avg2 /= nsteps;
  double block_var = block_avg2 - block_avg * block_avg;

  _block_avgs.emplace_back(block_avg);
  _block_vars.emplace_back(block_var);

  if (static_cast<int>(_block_avgs.size()) > _blocks_window) {
    _block_avgs.erase(_block_avgs.begin());
    _block_vars.erase(_block_vars.begin());
    _inst_measures.erase(_inst_measures.begin(),
                         _inst_measures.begin() + nsteps);
  }

  auto [m, m_err] = _lin_reg_slope();
  _m = m, _m_err = m_err;

  /* temperature correction */
  temp = 0.0, temp_err = 0.0;
  for (int i = 0; i < nsteps; i++) {
    double T = _temps_block[i];
    temp += T;
    temp_err += T * T;
  }
  temp /= nsteps, temp_err /= nsteps;
  temp_err = std::sqrt(temp_err - temp * temp);
  _temps_block.clear();
}

/* autocorrelation computation */

double Monitor::_autocorrelation_time() const {
  int n = _inst_measures.size();
  int k_max = n; /* TODO: to be determined! */
  double av = 0.0, av_sq = 0.0;
  for (int i = 0; i < n; i++) {
    double x = _inst_measures[i];
    av += x;
    av_sq += x * x;
  }
  av /= n;
  av_sq /= n;
  double var = av_sq - av * av;

  double tau = 0.5; /* k = 0 term */
  for (int k = 1; k < k_max; k++) {
    double C_k = 0.0;
    for (int i = 0; i < n - k; i++) {
      C_k += (_inst_measures[i] - av) * (_inst_measures[i + k] - av);
    }
    C_k /= var * (n - k);
    tau += C_k;
    if (k >= 5.0 * tau) {
      break;
    }
  }

  return 2.0 * tau;
}

int Monitor::get_correct_nsteps() const {
  return static_cast<int>(std::ceil(10.0 * _autocorrelation_time()));
}

/* temperature_correction */
bool Monitor::is_temp_ok() const {
  return !auto_temp_correction ||
         (std::abs(temp - temp_target) / temp_err <= temp_stdev_precision);
}
