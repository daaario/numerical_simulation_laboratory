#include "particle.hpp"

#include <cassert>

IsingParticle::IsingParticle(int spin) : _spin(spin) {
  assert(spin == 1 || spin == -1);
}

int IsingParticle::get_spin() const { return _spin; }

void IsingParticle::set_spin(int spin) {
  assert(spin == 1 || spin == -1);
  _spin = spin;
}

void IsingParticle::flip_spin() { _spin *= -1; }

PointParticle::PointParticle(arma::vec x, arma::vec v) : _x(x), _v(v) {}

const arma::vec PointParticle::get_x() const { return _x; }
const arma::vec PointParticle::get_x_old() const { return _x_old; }
const arma::vec PointParticle::get_v() const { return _v; }
const arma::vec PointParticle::get_a() const { return _a; }
void PointParticle::set_x(arma::vec x) { _x = x; }
void PointParticle::set_x_old(arma::vec x_old) { _x_old = x_old; }
void PointParticle::set_v(arma::vec v) { _v = v; }
void PointParticle::set_a(arma::vec a) { _a = a; }

double mod(const arma::vec v) {
  double mod = 0.0;
  for (auto v_i : v) {
    mod += v_i * v_i;
  }
  return std::sqrt(mod);
}
