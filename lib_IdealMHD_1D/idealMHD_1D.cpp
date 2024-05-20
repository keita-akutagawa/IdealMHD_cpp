#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "const.hpp"
#include "idealMHD_1D.hpp"


void IdealMHD1D::setU(
    const std::vector<std::vector<double>> UInit
)
{
    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            U[comp][i] = UInit[comp][i];
        }
    }
}


void IdealMHD1D::oneStepRK2()
{
    for (int comp = 0; comp < 8; comp++) {
        std::copy(U[comp].begin(), U[comp].end(), UBar[comp].begin());
    }
    calculateDt();
    fluxF = fluxSolver.getFluxF(U);

}


void IdealMHD1D::save(
    std::string directoryname, 
    std::string filenameWithoutStep, 
    int step
)
{
    std::string filename;
    filename = filenameWithoutStep + std::to_string(step);

    std::ofstream ofs(filename);
    ofs << std::fixed << std::setprecision(6);

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            ofs << U[comp][i] << ",";
        }
        ofs << std::endl;
    }
}


void IdealMHD1D::calculateDt()
{
    double rho, u, v, w, bx, by, bz, e, p, cs, ca;
    double maxSpeed;
    
    dt = 1e100; //十分大きくしておく
    for (int i = 0; i < nx; i++) {
        rho = U[0][i];
        u = U[1][i] / rho;
        v = U[2][i] / rho;
        w = U[3][i] / rho;
        bx = U[4][i];
        by = U[5][i];
        bz = U[6][i];
        e = U[7][i];
        p = (gamma_mhd - 1.0)
          * (e - 0.5 * rho * (u * u + v * v + w * w)
          - 0.5 * (bx * bx + by * by + bz * bz));
        
        cs = sqrt(gamma_mhd * p / rho);
        ca = sqrt((bx * bx + by * by + bz * bz) / rho);

        maxSpeed = std::abs(u) + sqrt(cs * cs + ca * ca);

        dt = std::min(dt, CFL * 1.0 / (maxSpeed + EPS));
    }
}

