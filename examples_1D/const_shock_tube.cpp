#include "const_shock_tube.hpp"

const double dx = 0.001;
const double xmin = 0.0;
const double xmax = 1.0;
const int nx = int((xmax - xmin) / dx);
const double CFL = 0.7;
const double gamma = 5.0 / 3.0;

