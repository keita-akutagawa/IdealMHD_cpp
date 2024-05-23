#include <algorithm>
#include <vector>
#include <iostream>
#include "const.hpp"
#include "calculate_half_components.hpp"


void CalculateHalfComponents::setPhysicalParameters(
    const std::vector<std::vector<double>>& U
)
{
    double rho, u, v, w, bx, by, bz, e, p;

    for (int i = 0; i < nDirection; i++) {
        rho = U[0][i];
        u = U[1][i] / rho;
        v = U[2][i] / rho;
        w = U[3][i] / rho;

        bx = U[4][i];
        by = U[5][i];
        bz = U[6][i];

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
    muscl.getLeftComponent(componentsCenter.rho, componentsLeft.rho);
    muscl.getLeftComponent(componentsCenter.u,   componentsLeft.u);
    muscl.getLeftComponent(componentsCenter.v,   componentsLeft.v);
    muscl.getLeftComponent(componentsCenter.w,   componentsLeft.w);
    muscl.getLeftComponent(componentsCenter.by,  componentsLeft.by);
    muscl.getLeftComponent(componentsCenter.bz,  componentsLeft.bz);
    muscl.getLeftComponent(componentsCenter.p,   componentsLeft.p); 

    for (int i = 0; i < nDirection; i++) {
        componentsLeft.bx[i] = componentsCenter.bx[i];
    }
}


void CalculateHalfComponents::calculateRightComponents()
{ 
    muscl.getRightComponent(componentsCenter.rho, componentsRight.rho);
    muscl.getRightComponent(componentsCenter.u,   componentsRight.u);
    muscl.getRightComponent(componentsCenter.v,   componentsRight.v);
    muscl.getRightComponent(componentsCenter.w,   componentsRight.w);
    muscl.getRightComponent(componentsCenter.by,  componentsRight.by);
    muscl.getRightComponent(componentsCenter.bz,  componentsRight.bz);
    muscl.getRightComponent(componentsCenter.p,   componentsRight.p);

    for (int i = 0; i < nDirection; i++) {
        componentsRight.bx[i] = componentsCenter.bx[i];
    }
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

