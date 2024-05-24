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
const double shear_thickness = 1.0;
const double rr = 0.2;
const double br = 1.0;
const double beta = 2.0;
const double theta = PI / 2.0;
const double rho0 = 1.0;
const double b0 = 1.0;
const double p0 = beta * b0 * b0 / 2.0;
const double v0 = sqrt(b0 * b0 / rho0 + gamma_mhd * p0 / rho0);

const double xmin = 0.0;
const double xmax = 2.0 * PI * shear_thickness / 0.4;
const double xCenter = (xmax - xmin) / 2.0;
const double dx = shear_thickness / 8.0;
const int nx = int((xmax - xmin) / dx);
const double ymin = 0.0;
const double ymax = 2.0 * 10.0 * shear_thickness;
const double yCenter = (ymax - ymin) / 2.0;
const double dy = shear_thickness / 8.0;
const int ny = int((ymax - ymin) / dy);

const double CFL = 0.7;
double dt = 0.0;
const int totalStep = 10000;
double totalTime = 0.0;


int main()
{
    std::string directoryname = "results";
    std::string filenameWithoutStep = "KH_rr=0.2";
    std::ofstream logfile("log_KH_rr=0.2.txt");
    int recordStep = 100;


    double rho, u, v, w, bx, by, bz, p, e;
    
    std::vector<std::vector<std::vector<double>>> UInit(8, std::vector(nx, std::vector<double>(ny, 0.0)));
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double xPosition, yPosition;
            xPosition = i * dx - xCenter;
            yPosition = j * dy - yCenter;

            rho = rho0 / 2.0 * ((1.0 - rr) * tanh(yPosition / shear_thickness) + 1.0 + rr);
            u = -v0 / 2.0 * tanh(yPosition / shear_thickness);
            v = 0.02 * v0 * cos(2.0 * PI * xPosition / xmax) / pow(cosh(yPosition / shear_thickness), 2);
            w = 0.0;
            bx = b0 / 2.0 * ((1.0 - br) * tanh(yPosition / shear_thickness) + 1.0 + br) * cos(theta);
            by = 0.0;
            bz = b0 / 2.0 * ((1.0 - br) * tanh(yPosition / shear_thickness) + 1.0 + br) * sin(theta);
            p = beta * (bx * bx + by * by + bz * bz) / 2.0;
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


