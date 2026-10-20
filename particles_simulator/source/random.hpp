#pragma once

#include <array>
#include <string>
#include <utility>

class Random {
private:
  std::array<int, 4> _seed;
  std::array<int, 4> _n;

public:
  Random(std::string);
  Random();
  void read_seed(std::string);
  void set_seed(const std::array<int, 4>, const std::pair<int, int>);
  void log_last_seed() const;

  double rannyu();
  double rannyu(double, double);
  double gauss(double, double);
  double exponential(double);
  double cauchy_lorentz(double, double);
};
