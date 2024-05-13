#include <algorithm>
#include "minmod.hpp"


inline double minmod(double x, double y)
{
    double eps = 1e-20;

    int sign_x = (x > 0) - (x < 0);
    double abs_x = std::abs(x);

    return sign_x * std::max(std::min(abs_x, sign_x * y), eps);
}


