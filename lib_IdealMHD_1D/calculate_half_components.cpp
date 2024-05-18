#include <algorithm>
#include <vector>
#include "const.hpp"
#include "calculate_half_components.hpp"


CalculateHalfComponents::CalculateHalfComponents()
{
    rho = std::vector<double>(nx, 0.0);
    u = std::vector<double>(nx, 0.0);
    v = std::vector<double>(nx, 0.0);
    w = std::vector<double>(nx, 0.0);
    bx = std::vector<double>(nx, 0.0);
    by = std::vector<double>(nx, 0.0);
    bz = std::vector<double>(nx, 0.0);
    p = std::vector<double>(nx, 0.0);

    rhoL = std::vector<double>(nx, 0.0);
    uL = std::vector<double>(nx, 0.0);
    vL = std::vector<double>(nx, 0.0);
    wL = std::vector<double>(nx, 0.0);
    bxL = std::vector<double>(nx, 0.0);
    byL = std::vector<double>(nx, 0.0);
    bzL = std::vector<double>(nx, 0.0);
    pL = std::vector<double>(nx, 0.0);

    rhoR = std::vector<double>(nx, 0.0);
    uR = std::vector<double>(nx, 0.0);
    vR = std::vector<double>(nx, 0.0);
    wR = std::vector<double>(nx, 0.0);
    bxR = std::vector<double>(nx, 0.0);
    byR = std::vector<double>(nx, 0.0);
    bzR = std::vector<double>(nx, 0.0);
    pR = std::vector<double>(nx, 0.0);

    tmpVector = std::vector<double>(nx, 0.0);
}


void CalculateHalfComponents::getPhysicalParameters(
    const std::vector<std::vector<double>> U
)
{
    for (int i = 0; i < nx; i++) {
        rho[i] = U[0][i];

        //U[1] = rho * u
        tmpVector[i] = U[1][i];
        u[i] = tmpVector[i] / rho[i];

        //U[2] = rho * v
        tmpVector[i] = U[2][i];
        v[i] = tmpVector[i] / rho[i];

        //U[3] = rho * w
        tmpVector[i] = U[3][i];
        w[i] = tmpVector[i] / rho[i];

        bx[i] = U[4][i];
        by[i] = U[5][i];
        bz[i] = U[6][i];

        //U[7] = e
        tmpVector[i] = U[7][i];
        p[i] = (gamma_mhd - 1.0)
             * (tmpVector[i]
             - 0.5 * rho[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i])
             - 0.5 * (bx[i] * bx[i] + by[i] * by[i] + bz[i] * bz[i])
             );
    }
}


void CalculateHalfComponents::calculateLeftComponents()
{ 
    muscl->getLeftComponent(rho, rhoL);
    muscl->getLeftComponent(u, uL);
    muscl->getLeftComponent(v, vL);
    muscl->getLeftComponent(w, wL);
    muscl->getLeftComponent(bx, bxL);
    muscl->getLeftComponent(by, byL);
    muscl->getLeftComponent(bz, bzL);
    muscl->getLeftComponent(p, pL);
}


void CalculateHalfComponents::calculateRightComponents()
{ 
    muscl->getRightComponent(rho, rhoR);
    muscl->getRightComponent(u, uR);
    muscl->getRightComponent(v, vR);
    muscl->getRightComponent(w, wR);
    muscl->getRightComponent(bx, bxR);
    muscl->getRightComponent(by, byR);
    muscl->getRightComponent(bz, bzR);
    muscl->getRightComponent(p, pR);
}

