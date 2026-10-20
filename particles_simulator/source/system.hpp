#pragma once

#include "measurable.hpp"
#include "random.hpp"
#include <armadillo>

namespace {
constexpr double pi = 3.141592653589793;
constexpr double k_B = 1.0;
} // namespace

class IsingParticle;
class PointParticle;
class Measurable;

class System {
protected:
  /* initialization methods */
  void _read_sys_settings();
  virtual void _load_init_config() = 0;
  /* simulation parameters */
  bool _restarted = false;
  int _nsteps, _nblocks, _nsteps_equil, _nblocks_equil;
  Random _rng;
  std::vector<Measurable> _measurables = {};
  /* physical parameters */
  int _npart;
  double _temp, _beta;
  /* finish and cleanup */
  virtual void _log_last_config() const = 0;

public:
  System();
  virtual ~System() = default;
  /* access methods */
  int get_nsteps() const, get_nblocks() const, get_nsteps_equil() const,
      get_nblocks_equil() const;
  bool get_restarted() const;
  double get_temp() const;
  /* automatic corrections */
  Monitor monitor;
  virtual double measure_temp() const = 0;
  void set_nsteps(int);
  void set_nsteps_equil(int);
  void set_temp(double);
  virtual void adjust_parameters() = 0;
  virtual void measure_for_monitoring() = 0;
  virtual void end_equilibration(double) = 0;
  /* main loop members */
  std::function<void()> step;
  virtual void measure() = 0;
  virtual void end_block(int);
  virtual void clean_up();
};

class IsingSystem : public System {
private:
  arma::field<IsingParticle> _particles;
  /* parameters */
  double _J, _H;
  /* general tools */
  int _pbc(int) const;
  /* simulation loop methods */
  std::function<bool(int)> _attempt_flip;
  bool _flip_metropolis(int);
  bool _flip_gibbs(int);
  void _step_ising();
  /* initialization */
  void _read_ising_settings();
  void _define_measurables();
  void _load_init_config() override;
  /* finish and cleanup */
  void _log_last_config() const override;

public:
  IsingSystem();

  double measure_temp() const override;
  void adjust_parameters() override;       // TODO: implement!!
  void measure_for_monitoring() override;  // TODO: implement!!
  void end_equilibration(double) override; // TODO: implement!!
  void measure() override;
};

class LJSystem : public System {
private:
  arma::field<PointParticle> _particles;
  /* parameters */
  int _ndim = 3;
  static constexpr double _epsilon = 1.0, _sigma = 1.0, _m = 1.0;
  double _box_dim, _volume, _rho;
  double _r_cut;
  double _ptail = 0.0, _vtail = 0.0;
  /* measurables */
  arma::vec _force(int) const;
  double _kin_energy() const;
  double _pot_energy(int, int) const;
  double _pot_energy(int) const;
  double _pot_energy() const;
  std::optional<MeasurableHistogram> _gofr, _pofv;
  int _n_bins_gofr, _n_bins_pofv;
  double _bin_size_gofr, _bin_size_pofv, _v_max_pofv;
  /* general tools */
  arma::vec _pbc(arma::vec) const;
  double _nearest_image_dist(arma::vec, arma::vec) const;

  /* simulation loop methods */
  double _dt;
  void _step_verlet();
  double _delta;
  void _step_metropolis();
  /* initialization */
  void _read_lj_settings();
  void _define_measurables();
  void _load_init_config() override;
  /* finish and cleanup */
  void _log_last_config() const override;

public:
  LJSystem();

  double measure_temp() const override;
  void adjust_parameters() override;
  void measure_for_monitoring() override;
  void end_equilibration(double) override;
  void measure() override;
  void end_block(int) override;
  void clean_up() override;
};

std::unique_ptr<System> create_system();
void equilibrate(System &);
void simulate(System &);
void clean_up(System &);

void print_progress_bar(int, int, int, int);
void print_auto_equil(int, int, System &);
void print_manual_equil(int, int);
