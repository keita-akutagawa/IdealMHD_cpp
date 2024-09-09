#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../../lib_IdealMHD_2D/const.hpp"
#include "../../lib_IdealMHD_2D/idealMHD_2D_symmetric_y.hpp"

const double EPS = 1e-20;
const double PI = 3.141592653589793;
const double dtError = 1e100;

const double gamma_mhd = 5.0 / 3.0;
const double sheat_thickness = 1.0;
const double beta_upstream = 0.5 * 0.5;
const double rho0 = 1.0;
const double b0 = 1.0;
const double p0 = b0 * b0 / 2.0;

const double xmin = 0.0;
const double xmax = 1.0 * sheat_thickness;
const double dx = sheat_thickness / 10.0;
const int nx = int((xmax - xmin) / dx);
const double ymin = 0.0;
const double ymax = 200.0 * sheat_thickness;
const double dy = sheat_thickness / 10.0;
const int ny = int((ymax - ymin) / dy);

const double CFL = 0.7;
double dt = 0.0;
const int totalStep = 10000;
const int recordStep = 100;
double totalTime = 0.0;


int main()
{
    std::string directoryname = "results_current_sheet";
    std::string filenameWithoutStep = "current_sheet";
    std::ofstream logfile("results_current_sheet/log_current_sheet.txt");


    double rho, u, v, w, bx, by, bz, p, e;
    
    std::vector<std::vector<std::vector<double>>> UInit(8, std::vector(nx, std::vector<double>(ny, 0.0)));
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x, y;
            x = i * dx;
            y = j * dy;
            double yCenter = 0.5 * (ymax - ymin) + ymin;
            double VA = b0 / sqrt(rho0);

            rho = rho0 * (1.0 + pow(cosh((y - yCenter) / sheat_thickness), -2) / beta_upstream);
            u = 0.0;
            v = 1.0 * VA * 0.5 * (-tanh((j - 0.3 * ny) / 100) + 1.0)
              + 1.0 * VA * 0.5 * (-tanh((j - 0.7 * ny) / 100) - 1.0);
            w = 0.0;
            bx = b0 * tanh((y - yCenter) / sheat_thickness);
            by = 0.0;
            bz = 0.0;
            p = p0 * (beta_upstream + pow(cosh((y - yCenter) / sheat_thickness), -2));
            e = p / (gamma_mhd - 1.0)
               + 0.5 * rho * (u * u + v * v + w * w)
               + 0.5 * (bx * bx + by * by + bz * bz);

            UInit[0][i][j] = rho;
            UInit[1][i][j] = rho * u;
            UInit[2][i][j] = rho * v;
            UInit[3][i][j] = rho * w;
            UInit[4][i][j] = bx;
            UInit[5][i][j] = by;
            UInit[6][i][j] = bz;
            UInit[7][i][j] = e;
        }
    }


    IdealMHD2DSymmetricY idealMHD2D;

    idealMHD2D.initializeU(UInit);

    for (int step = 0; step < totalStep+1; step++) {
        if (step % recordStep == 0) {
            idealMHD2D.save(directoryname, filenameWithoutStep, step);
            logfile << std::to_string(step) << ","
                    << std::setprecision(4) << totalTime
                    << std::endl;
            std::cout << std::to_string(step) << " step done : total time is "
                      << std::setprecision(4) << totalTime
                      << std::endl;
        }
        
        idealMHD2D.oneStepRK2();

        if (idealMHD2D.checkCalculationIsCrashed()) {
            std::cout << "Calculation stopped! : " << step << " steps" << std::endl;
            return 0;
        }

        totalTime += dt;
    }
    
    return 0;
}


