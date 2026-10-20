#pragma once

#include <armadillo>

class IsingParticle {
private:
  int _spin = 1;

public:
  IsingParticle() = default;
  explicit IsingParticle(int);

  int get_spin() const;
  void set_spin(int);

  void flip_spin();
};

class PointParticle {
private:
  arma::vec _x, _x_old, _v, _a;

public:
  PointParticle() = default;
  explicit PointParticle(arma::vec, arma::vec);

  const arma::vec get_x() const, get_x_old() const, get_v() const,
      get_a() const;
  void set_x(arma::vec), set_x_old(arma::vec), set_v(arma::vec),
      set_a(arma::vec);
};

double mod(const arma::vec);
