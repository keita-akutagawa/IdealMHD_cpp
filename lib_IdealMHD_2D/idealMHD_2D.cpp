#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "const.hpp"
#include "idealMHD_2D.hpp"


IdealMHD2D::IdealMHD2D()
{
    U = std::vector(8, std::vector(nx, std::vector<double>(ny, 0.0)));
    UBar = std::vector(8, std::vector(nx, std::vector<double>(ny, 0.0)));
}


void IdealMHD2D::initializeU(
    const std::vector<std::vector<std::vector<double>>> UInit
)
{
    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                U[comp][i][j] = UInit[comp][i][j];
            }
        }
    }
}


void IdealMHD2D::oneStepRK2()
{
    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                UBar[comp][i][j] = U[comp][i][j];
            }
        }
    }

    calculateDt();

    fluxF = fluxSolver.getFluxF(U);

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 1; i < nx; i++) {
            UBar[comp][i] = U[comp][i]
                          - dt / dx * (fluxF.flux[comp][i] - fluxF.flux[comp][i-1]);
        }
        //周期境界条件
        UBar[comp][0] = U[comp][0] 
                      - dt / dx * (fluxF.flux[comp][0] - fluxF.flux[comp][nx-1]);
    }

    //これはどうにかすること。保守性が低い
    boundary.symmetricBoundary2nd(UBar);

    fluxF = fluxSolver.getFluxF(UBar);

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 1; i < nx; i++) {
            U[comp][i] = 0.5 * (U[comp][i] + UBar[comp][i]
                       - dt / dx * (fluxF.flux[comp][i] - fluxF.flux[comp][i-1]));
        }
        //周期境界条件
        U[comp][0] = 0.5 * (U[comp][0] + UBar[comp][0]
                   - dt / dx * (fluxF.flux[comp][0] - fluxF.flux[comp][nx-1]));
    }

    //これはどうにかすること。保守性が低い
    boundary.symmetricBoundary2nd(U);
}


void IdealMHD2D::save(
    std::string directoryname, 
    std::string filenameWithoutStep, 
    int step
)
{
    std::string filename;
    filename = directoryname + "/"
             + filenameWithoutStep + "_" + std::to_string(step)
             + ".txt";

    std::ofstream ofs(filename);
    ofs << std::fixed << std::setprecision(6);

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx-1; i++) {
            for (int j = 0; j < ny; j++) {
                ofs << U[comp][i][j] << ",";
            }
        }
        for (int j = 0; j < ny-1; j++) {
            ofs << U[comp][nx-1][j] << ",";
        }
        ofs << U[comp][nx-1][ny-1];
        ofs << std::endl;
    }
}


void IdealMHD2D::calculateDt()
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

        dt = std::min(dt, 1.0 / (maxSpeed / dx + EPS));
    }
    
    dt *= CFL;
}


// getter
std::vector<std::vector<std::vector<double>>> IdealMHD2D::getU()
{
    return U;
}

