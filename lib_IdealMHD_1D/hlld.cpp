#include <vector>
#include "const.hpp"
#include "hlld.hpp"


HLLD::HLLD()
{
    CalculateHalfComponents calculateHalfComponents;

    componentsCenter.rho = std::vector<double>(nx, 0.0);
    componentsCenter.u = std::vector<double>(nx, 0.0);
    componentsCenter.v = std::vector<double>(nx, 0.0);
    componentsCenter.w = std::vector<double>(nx, 0.0);
    componentsCenter.bx = std::vector<double>(nx, 0.0);
    componentsCenter.by = std::vector<double>(nx, 0.0);
    componentsCenter.bz = std::vector<double>(nx, 0.0);
    componentsCenter.p = std::vector<double>(nx, 0.0);

    outerFanParameters.rho = std::vector<double>(nx, 0.0);
    outerFanParameters.u = std::vector<double>(nx, 0.0);
    outerFanParameters.v = std::vector<double>(nx, 0.0);
    outerFanParameters.w = std::vector<double>(nx, 0.0);
    outerFanParameters.bx = std::vector<double>(nx, 0.0);
    outerFanParameters.by = std::vector<double>(nx, 0.0);
    outerFanParameters.bz = std::vector<double>(nx, 0.0);
    outerFanParameters.e = std::vector<double>(nx, 0.0);

    innerFanParameters.rho = std::vector<double>(nx, 0.0);
    innerFanParameters.u = std::vector<double>(nx, 0.0);
    innerFanParameters.v = std::vector<double>(nx, 0.0);
    innerFanParameters.w = std::vector<double>(nx, 0.0);
    innerFanParameters.bx = std::vector<double>(nx, 0.0);
    innerFanParameters.by = std::vector<double>(nx, 0.0);
    innerFanParameters.bz = std::vector<double>(nx, 0.0);
    innerFanParameters.e = std::vector<double>(nx, 0.0);
}


void HLLD::calculateHLLDParametersForOuterFanParameters(
    const std::vector<std::vector<double>> U
)
{
    calculateHalfComponents.setPhysicalParameters(U);
    calculateHalfComponents.calculateLeftComponents();
    calculateHalfComponents.calculateRightComponents();

    Components leftComponents, rightComponents;
    leftComponents = calculateHalfComponents.getLeftComponents();
    rightComponents = calculateHalfComponents.getRightComponents();

    calculateHLLDSubParameters(
        leftComponents, 
        hLLDLeftParameters
    );
    calculateHLLDSubParameters(
        rightComponents, 
        hLLDRightParameters
    );

    
}

void HLLD::calculateHLLDSubParameters(
    const Components components, 
    HLLDParameters hLLDParameters
)
{
    double rho, u, v, w, bx, by, bz, p;
    double pT, e, cs, ca, va, cf;
    for (int i = 0; i < nx; i++) {
        rho = components.rho[i];
        u = components.u[i];
        v = components.v[i];
        w = components.w[i];
        bx = components.bx[i];
        by = components.by[i];
        bz = components.bz[i];
        p = components.p[i];

        pT = p + 0.5 * (bx * bx + by * by + bz * bz);
        e = p / (gamma_mhd - 1.0)
          + 0.5 * rho * (u * u + v * v + w * w)
          + 0.5 * (bx * bx + by * by + bz * bz);
        cs = sqrt(gamma_mhd * p / rho);
        ca = sqrt((bx * bx + by * by + bz * bz) / rho);
        va = sqrt(bx * bx / rho);
        cf = sqrt(0.5 * (cs * cs + ca * ca
           + sqrt((cs * cs + ca * ca) * (cs * cs + ca * ca)
           - 4.0 * cs * cs * va * va)));

        hLLDParameters.pT[i] = pT;
        hLLDParameters.e[i] = e;
        hLLDParameters.cs[i] = cs;
        hLLDParameters.ca[i] = ca;
        hLLDParameters.va[i] = va;
        hLLDParameters.cf[i] = cf;
    }
}


void HLLD::calculateHLLDParametersForInnerFanParameters()
{

}


Components HLLD::getOuterFanParameters()
{

}


Components HLLD::getInnerFanParameters()
{

}

