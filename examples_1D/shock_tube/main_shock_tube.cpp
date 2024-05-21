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
const int totalStep = 1;


int main()
{
    std::string directoryname = "results";
    std::string filenameWithoutStep = "shock_tube";

    IdealMHD1D idealMHD1D;
    std::vector<std::vector<double>> UInit(8, std::vector<double>(nx, 1.0));

    idealMHD1D.initializeU(UInit);

    for (int step = 0; step < totalStep; step++) {
        if (step % 10 == 0) {
            idealMHD1D.save(directoryname, filenameWithoutStep, step);
        }
        idealMHD1D.oneStepRK2();
    }
}


