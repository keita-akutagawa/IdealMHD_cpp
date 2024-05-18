#include <algorithm>
#include <vector>
#include "const.hpp"
#include "calculate_half_components.hpp"


CalculateHalfComponents::CalculateHalfComponents()
{
    componentsCenter.rho = std::vector<double>(nx, 0.0);
    componentsCenter.u = std::vector<double>(nx, 0.0);
    componentsCenter.v = std::vector<double>(nx, 0.0);
    componentsCenter.w = std::vector<double>(nx, 0.0);
    componentsCenter.bx = std::vector<double>(nx, 0.0);
    componentsCenter.by = std::vector<double>(nx, 0.0);
    componentsCenter.bz = std::vector<double>(nx, 0.0);
    componentsCenter.p = std::vector<double>(nx, 0.0);

    componentsLeft.rho = std::vector<double>(nx, 0.0);
    componentsLeft.u = std::vector<double>(nx, 0.0);
    componentsLeft.v = std::vector<double>(nx, 0.0);
    componentsLeft.w = std::vector<double>(nx, 0.0);
    componentsLeft.bx = std::vector<double>(nx, 0.0);
    componentsLeft.by = std::vector<double>(nx, 0.0);
    componentsLeft.bz = std::vector<double>(nx, 0.0);
    componentsLeft.p = std::vector<double>(nx, 0.0);

    componentsRight.rho = std::vector<double>(nx, 0.0);
    componentsRight.u = std::vector<double>(nx, 0.0);
    componentsRight.v = std::vector<double>(nx, 0.0);
    componentsRight.w = std::vector<double>(nx, 0.0);
    componentsRight.bx = std::vector<double>(nx, 0.0);
    componentsRight.by = std::vector<double>(nx, 0.0);
    componentsRight.bz = std::vector<double>(nx, 0.0);
    componentsRight.p = std::vector<double>(nx, 0.0);

    tmpVector = std::vector<double>(nx, 0.0);
}


void CalculateHalfComponents::getPhysicalParameters(
    const std::vector<std::vector<double>> U
)
{
    double rho, u, v, w, bx, by, bz, e, p;

    for (int i = 0; i < nx; i++) {
        rho = U[0][i];

        //U[1] = rho * u
        u = U[1][i] / rho;

        //U[2] = rho * v
        v = U[2][i] / rho;

        //U[3] = rho * w
        w = U[3][i] / rho;

        bx = U[4][i];
        by = U[5][i];
        bz = U[6][i];

        //U[7] = e
        e = U[7][i];
        p = (gamma_mhd - 1.0)
          * (e
          - 0.5 * rho * (u * u + v * v + w * w)
          - 0.5 * (bx * bx + by * by + bz * bz)
          );
        

        componentsCenter.rho[i] = rho;
        componentsCenter.u[i] = u;
        componentsCenter.v[i] = v;
        componentsCenter.w[i] = w;
        componentsCenter.bx[i] = bx;
        componentsCenter.by[i] = by;
        componentsCenter.bz[i] = bz;
        componentsCenter.p[i] = p;
    }
}


void CalculateHalfComponents::calculateLeftComponents()
{ 
    muscl->getLeftComponent(componentsCenter.rho, componentsLeft.rho);
    muscl->getLeftComponent(componentsCenter.u,   componentsLeft.u);
    muscl->getLeftComponent(componentsCenter.v,   componentsLeft.v);
    muscl->getLeftComponent(componentsCenter.w,   componentsLeft.w);
    muscl->getLeftComponent(componentsCenter.bx,  componentsLeft.bx);
    muscl->getLeftComponent(componentsCenter.by,  componentsLeft.by);
    muscl->getLeftComponent(componentsCenter.bz,  componentsLeft.bz);
    muscl->getLeftComponent(componentsCenter.p,   componentsLeft.p);
}


void CalculateHalfComponents::calculateRightComponents()
{ 
    muscl->getRightComponent(componentsCenter.rho, componentsRight.rho);
    muscl->getRightComponent(componentsCenter.u,   componentsRight.u);
    muscl->getRightComponent(componentsCenter.v,   componentsRight.v);
    muscl->getRightComponent(componentsCenter.w,   componentsRight.w);
    muscl->getRightComponent(componentsCenter.bx,  componentsRight.bx);
    muscl->getRightComponent(componentsCenter.by,  componentsRight.by);
    muscl->getRightComponent(componentsCenter.bz,  componentsRight.bz);
    muscl->getRightComponent(componentsCenter.p,   componentsRight.p);
}


Components CalculateHalfComponents::getCenterComponents()
{
    return componentsCenter;
}


Components CalculateHalfComponents::getLeftComponents()
{
    return componentsLeft;
}


Components CalculateHalfComponents::getRightComponents()
{
    return componentsRight;
}

