#pragma once

#include <cmath>
#include <fstream>
#include <functional>
#include <span>
#include <string>

class Measurable {
private:
  std::string _name;
  std::ofstream _output_stream;
  std::function<double(std::span<const double>)> _compute;
  std::function<double(double, double, int)> _block_estimate;

  double _block_sum = 0.0, _block_sum_sq = 0.0;
  double _global_av = 0.0, _global_av_sq = 0.0;

public:
  Measurable(std::string, std::function<double(std::span<const double>)>,
             std::function<double(double, double, int)>);
  Measurable(std::string, std::function<double(std::span<const double>)>);
  void measure(std::span<const double>);
  void end_block(int, int);
};

class MeasurableHistogram {
private:
  std::string _name;
  int _n_bins;
  std::vector<double> _block_sum{}, _block_sum_sq{};
  std::vector<double> _global_av{}, _global_av_sq{};
  /* evaluation methods */
  std::function<std::vector<double>(std::span<const double>)> _compute;
  std::function<double(double, double, int)> _block_estimate;
  /* data streams */
  std::ofstream _output_global_stream;
  std::ofstream _output_blocks_stream;

public:
  MeasurableHistogram(
      std::string, int,
      std::function<std::vector<double>(std::span<const double>)>,
      std::function<double(double, double, int)>);
  MeasurableHistogram(
      std::string, int,
      std::function<std::vector<double>(std::span<const double>)>);
  void measure(std::span<const double>);
  void end_block(int, int);
  void finalize(int);
};

double error(double, double, int);

/**************************************************/
/********** Thermalisation monitor class **********/
/**************************************************/

class Monitor {
  /* TODO: save adjusted parameters in logfile to restart the simulation */
private:
  int _blocks_window = 16;
  /* acceptance ratio monitoring to adjust _delta */
  int _nattempted = 0, _naccepted = 0;
  double _acceptance_ratio = 0.0;
  double _acceptance_target = 0.5;
  double _acceptance_tolerance = 0.05;
  /* istantaneous measures of a selected observable
   * to detect equilibration */
  std::vector<double> _inst_measures;
  std::vector<double> _block_avgs, _block_vars;
  std::pair<double, double> _lin_reg_slope() const;
  double _m = 0.0, _m_err = 1.0;
  double _slope_stdev_precision = 2.0;
  /* temperature correction */
  std::vector<double> _temps_block;
  /* autocorrelation */
  double _autocorrelation_time() const;

public:
  Monitor() = default;
  Monitor(int, double, bool, double, double, double);
  /* acceptance ratio monitoring to adjust _delta */
  bool auto_mc_parameters = false;
  void naccepted_up(), nattempted_up();
  double get_acc_ratio() const;
  double get_acc_target() const;
  bool is_acc_ratio_ok() const;
  /* equilibration check */
  bool auto_equilibration = false;
  double get_m() const;
  double get_m_err() const;
  bool is_at_equilibrium() const;
  /* temperature correction */
  bool auto_temp_correction = false;
  double temp, temp_err = 0.0;
  double temp_target, temp_stdev_precision;
  bool is_temp_ok() const;

  void save_measure(double, double);
  void end_block(int);
  int get_correct_nsteps() const;
};
