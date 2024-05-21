#include <vector>
#include <iostream>
#include "../../lib_IdealMHD_1D/const.hpp"
#include "../../lib_IdealMHD_1D/idealMHD_1D.hpp"

const double EPS = 1e-20;
const double dx = 0.001;
const double xmin = 0.0;
const double xmax = 1.0;
const int nx = int((xmax - xmin) / dx);
const double CFL = 0.7;
const double gamma_mhd = 5.0 / 3.0;
double dt = 0.0;
const int totalStep = 100;


int main()
{
    std::string directoryname = "results";
    std::string filenameWithoutStep = "shock_tube";

    double rho0, u0, v0, w0, bx0, by0, bz0, p0, e0;
    rho0 = u0 = v0 = w0 = bx0 = by0 = bz0 = p0 = 1.0;
    e0 = p0 / (gamma_mhd - 1.0)
      + 0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0)
      + 0.5 * (bx0 * bx0 + by0 * by0 + bz0 * bz0);

    std::vector<std::vector<double>> UInit(8, std::vector<double>(nx, 0.0));
    for (int i = 0; i < nx; i++) {
        UInit[0][i] = rho0;
        UInit[1][i] = u0;
        UInit[2][i] = v0;
        UInit[3][i] = w0;
        UInit[4][i] = bx0;
        UInit[5][i] = by0;
        UInit[6][i] = bz0;
        UInit[7][i] = e0;
    }


    IdealMHD1D idealMHD1D;

    idealMHD1D.initializeU(UInit);

    for (int step = 0; step < totalStep; step++) {
        if (step % 10 == 0) {
            idealMHD1D.save(directoryname, filenameWithoutStep, step);
        }
        idealMHD1D.oneStepRK2();
    }
}


