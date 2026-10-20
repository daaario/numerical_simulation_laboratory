#include "system.hpp"

int main() {
  std::unique_ptr<System> sys = create_system();

  equilibrate(*sys);
  simulate(*sys);
  clean_up(*sys);

  return 0;
}

/* functions are defined in 'system.cpp' */
