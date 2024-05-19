#include <vector>
#include "const.hpp"
#include "hlld.hpp"


inline double sign(double x)
{
    return (x > 0.0) - (x < 0.0);
}


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

    outerLeftFanParameters.rho = std::vector<double>(nx, 0.0);
    outerLeftFanParameters.u = std::vector<double>(nx, 0.0);
    outerLeftFanParameters.v = std::vector<double>(nx, 0.0);
    outerLeftFanParameters.w = std::vector<double>(nx, 0.0);
    outerLeftFanParameters.bx = std::vector<double>(nx, 0.0);
    outerLeftFanParameters.by = std::vector<double>(nx, 0.0);
    outerLeftFanParameters.bz = std::vector<double>(nx, 0.0);
    outerLeftFanParameters.e = std::vector<double>(nx, 0.0);

    outerRightFanParameters.rho = std::vector<double>(nx, 0.0);
    outerRightFanParameters.u = std::vector<double>(nx, 0.0);
    outerRightFanParameters.v = std::vector<double>(nx, 0.0);
    outerRightFanParameters.w = std::vector<double>(nx, 0.0);
    outerRightFanParameters.bx = std::vector<double>(nx, 0.0);
    outerRightFanParameters.by = std::vector<double>(nx, 0.0);
    outerRightFanParameters.bz = std::vector<double>(nx, 0.0);
    outerRightFanParameters.e = std::vector<double>(nx, 0.0);

    innerLeftFanParameters.rho = std::vector<double>(nx, 0.0);
    innerLeftFanParameters.u = std::vector<double>(nx, 0.0);
    innerLeftFanParameters.v = std::vector<double>(nx, 0.0);
    innerLeftFanParameters.w = std::vector<double>(nx, 0.0);
    innerLeftFanParameters.bx = std::vector<double>(nx, 0.0);
    innerLeftFanParameters.by = std::vector<double>(nx, 0.0);
    innerLeftFanParameters.bz = std::vector<double>(nx, 0.0);
    innerLeftFanParameters.e = std::vector<double>(nx, 0.0);

    innerRightFanParameters.rho = std::vector<double>(nx, 0.0);
    innerRightFanParameters.u = std::vector<double>(nx, 0.0);
    innerRightFanParameters.v = std::vector<double>(nx, 0.0);
    innerRightFanParameters.w = std::vector<double>(nx, 0.0);
    innerRightFanParameters.bx = std::vector<double>(nx, 0.0);
    innerRightFanParameters.by = std::vector<double>(nx, 0.0);
    innerRightFanParameters.bz = std::vector<double>(nx, 0.0);
    innerRightFanParameters.e = std::vector<double>(nx, 0.0);

    std::vector<std::vector<double>> flux(8, std::vector<double>(nx, 0.0));

}


void HLLD::calculateFlux(
    const std::vector<std::vector<double>> U
)
{
    setComponents(U);
    calculateHLLDParametersForMiddleFan();
    calculateHLLDParametersForInnerFan();
    //setFlux();

}



void HLLD::setComponents(
    const std::vector<std::vector<double>> U
)
{
    calculateHalfComponents.setPhysicalParameters(U);
    calculateHalfComponents.calculateLeftComponents();
    calculateHalfComponents.calculateRightComponents();

    componentsLeft = calculateHalfComponents.getLeftComponents();
    componentsRight = calculateHalfComponents.getRightComponents();
}


void HLLD::calculateHLLDParametersForOuterFan(
    const Components components
)
{
    double rho, u, v, w, bx, by, bz, p, e, pT;
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
    }
}


void HLLD::calculateHLLDParametersForMiddleFan()
{
    calculateHLLDSubParameters(
        componentsLeft, 
        hLLDLeftParameters
    );
    calculateHLLDSubParameters(
        componentsRight, 
        hLLDRightParameters
    );


    double SL, SR, SM, pT1, pTL, pTR;
    double rho1, u1, v1, w1, bx1, by1, bz1, e1, pT1L, pT1R;
    double rhoL, rhoR, uL, uR, cfL, cfR;
    for (int i = 0; i < nx; i++) {
        rhoL = componentsLeft.rho[i];
        rhoR = componentsRight.rho[i];
        uL = componentsLeft.u[i];
        uR = componentsRight.u[i];
        pTL = hLLDLeftParameters.pT[i];
        pTR = hLLDRightParameters.pT[i];
        cfL = hLLDLeftParameters.cf[i];
        cfR = hLLDRightParameters.cf[i];

        SL = std::min(uL, uR) - std::max(cfL, cfR);
        SR = std::max(uL, uR) + std::max(cfL, cfR);
        SL = std::min(SL, 0.0);
        SR = std::max(SR, 0.0);

        SM = ((SR - uR) * rhoR * uR - (SL - uL) * rhoL * uL - pTR + pTL)
           / ((SR - uR) * rhoR - (SL - uL) * rhoL + EPS);
        pT1 = ((SR - uR) * rhoR * pTL - (SL - uL) * rhoL * pTR
            + rhoL * rhoR * (SR - uR) * (SL - uL) * (uR - uL))
            / ((SR - uR) * rhoR - (SL - uL) * rhoL + EPS);
        pT1L = pTL;
        pT1R = pTR;

        hLLDLeftParameters.S[i] = SL;
        hLLDRightParameters.S[i] = SR;
        hLLDLeftParameters.SM[i] = SM;
        hLLDRightParameters.SM[i] = SM; //やらなくてもOK
        hLLDLeftParameters.pT1[i] = pT1L;
        hLLDRightParameters.pT1[i] = pT1R;
    }

    calculateHLLDParameters1(
        componentsLeft, hLLDLeftParameters, outerLeftFanParameters
    );
    calculateHLLDParameters1(
        componentsRight, hLLDRightParameters, outerRightFanParameters
    );
}


void HLLD::calculateHLLDParametersForInnerFan()
{
    double S1L, S1R, SM, rho1L, rho1R, bx1;
    for (int i = 0; i < nx; i++) {
        SM = hLLDLeftParameters.SM[i]; //RightでもOK
        rho1L = outerLeftFanParameters.rho[i];
        rho1R = outerRightFanParameters.rho[i];
        bx1 = outerLeftFanParameters.bx[i]; //RightでもOK

        S1L = SM - sqrt(bx1 * bx1 / rho1L);
        S1R = SM + sqrt(bx1 * bx1 / rho1R);

        hLLDLeftParameters.S1[i] = S1L;
        hLLDRightParameters.S1[i] = S1R;
    }

    calculateHLLDParameters2(
        outerLeftFanParameters, 
        outerRightFanParameters, 
        innerLeftFanParameters, 
        innerRightFanParameters
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


void HLLD::calculateHLLDParameters1(
    const Components components, 
    const HLLDParameters hlldParameters, 
    FanParameters outerFanParameters
)
{
    double rho, u, v, w, bx, by, bz, e, pT, pT1, S, SM;
    double rho1, u1, v1, w1, bx1, by1, bz1, e1;
    for (int i = 0; i < nx; i++) {
        rho = components.rho[i];
        u = components.u[i];
        v = components.v[i];
        w = components.w[i];
        bx = components.bx[i];
        by = components.by[i];
        bz = components.bz[i];

        e = hlldParameters.e[i];
        pT = hlldParameters.pT[i];
        pT1 = hlldParameters.pT1[i];
        S = hlldParameters.S[i];
        SM = hlldParameters.SM[i];

        rho1 = rho * (S - u) / (S - SM + EPS);
        u1 = SM;
        v1 = v - bx * by * (SM - u) / (rho * (S - u) * (S - SM) - bx * bx + EPS);
        w1 = w - bx * bz * (SM - u) / (rho * (S - u) * (S - SM) - bx * bx + EPS);
        bx1 = bx;
        by1 = by * (rho * (S - u) * (S - u) - bx * bx)
            / (rho * (S - u) * (S - SM) - bx * bx + EPS);
        bz1 = bz * (rho * (S - u) * (S - u) - bx * bx)
            / (rho * (S - u) * (S - SM) - bx * bx + EPS);
        e1 = ((S - u) * e - pT * u + pT1 * SM
           + bx * ((u * bx + v * by + w * bz) - (u1 * bx1 + v1 * by1 + w1 * bz1)))
           / (S - SM + EPS);
        
        outerFanParameters.rho[i] = rho1;
        outerFanParameters.u[i] = u1;
        outerFanParameters.v[i] = v1;
        outerFanParameters.w[i] = w1;
        outerFanParameters.bx[i] = bx1;
        outerFanParameters.by[i] = by1;
        outerFanParameters.bz[i] = bz1;
        outerFanParameters.e[i] = e1;
    }
}


void HLLD::calculateHLLDParameters2(
    const FanParameters outerLeftFanParameters,
    const FanParameters outerRightFanParameters, 
    FanParameters innerLeftFanParameters, 
    FanParameters innerRightFanParameters
)
{
    double rho1L, u1L, v1L, w1L, bx1L, by1L, bz1L, e1L;
    double rho1R, u1R, v1R, w1R, bx1R, by1R, bz1R, e1R;
    double SM;
    double rho2L, u2L, v2L, w2L, bx2L, by2L, bz2L, e2L;
    double rho2R, u2R, v2R, w2R, bx2R, by2R, bz2R, e2R;

    for (int i = 0; i < nx; i++) {
        rho1L = outerLeftFanParameters.rho[i];
        u1L = outerLeftFanParameters.u[i];
        v1L = outerLeftFanParameters.v[i];
        w1L = outerLeftFanParameters.w[i];
        bx1L = outerLeftFanParameters.bx[i];
        by1L = outerLeftFanParameters.by[i];
        bz1L = outerLeftFanParameters.bz[i];
        e1L = outerLeftFanParameters.e[i];
        
        rho1R = outerRightFanParameters.rho[i];
        u1R = outerRightFanParameters.u[i];
        v1R = outerRightFanParameters.v[i];
        w1R = outerRightFanParameters.w[i];
        bx1R = outerRightFanParameters.bx[i];
        by1R = outerRightFanParameters.by[i];
        bz1R = outerRightFanParameters.bz[i];
        e1R = outerRightFanParameters.e[i];

        rho2L = rho1L;
        rho2R = rho1R;
        u2L = SM;
        u2R = u2L;
        v2L = (sqrt(rho1L) * v1L + sqrt(rho1R) * v1R + (by1R - by1L) * sign(bx1L))
            / (sqrt(rho1L) + sqrt(rho1R) + EPS);
        v2R = v2L;
        w2L = (sqrt(rho1L) * w1L + sqrt(rho1R) * w1R + (bz1R - bz1L) * sign(bx1L))
            / (sqrt(rho1L) + sqrt(rho1R) + EPS);
        w2R = w2L;
        bx2L = bx1L;
        bx2R = bx1R;
        by2L = (sqrt(rho1L) * by1R + sqrt(rho1R) * by1L + sqrt(rho1L * rho1R) * (v1R - v1L) * sign(bx1L))
             / (sqrt(rho1L) + sqrt(rho1R) + EPS);
        by2R = by2L;
        bz2L = (sqrt(rho1L) * bz1R + sqrt(rho1R) * bz1L + sqrt(rho1L * rho1R) * (w1R - w1L) * sign(bx1L))
             / (sqrt(rho1L) + sqrt(rho1R) + EPS);
        bz2R = bz2L;
        e2L = e1L - sqrt(rho1L)
            * ((u1L * bx1L + v1L * by1L + w1L * bz1L) - (u2L * bx2L + v2L * by2L + w2L * bz2L))
            * sign(bx2L);
        e2R = e1R + sqrt(rho1R)
            * ((u1R * bx1R + v1R * by1R + w1R * bz1R) - (u2R * bx2R + v2R * by2R + w2R * bz2R))
            * sign(bx2R);
        

        innerLeftFanParameters.rho[i] = rho2L;
        innerLeftFanParameters.u[i] = u2L;
        innerLeftFanParameters.v[i] = v2L;
        innerLeftFanParameters.w[i] = w2L;
        innerLeftFanParameters.bx[i] = bx2L;
        innerLeftFanParameters.by[i] = by2L;
        innerLeftFanParameters.bz[i] = bz2L;
        innerLeftFanParameters.e[i] = e2L;

        innerRightFanParameters.rho[i] = rho2R;
        innerRightFanParameters.u[i] = u2R;
        innerRightFanParameters.v[i] = v2R;
        innerRightFanParameters.w[i] = w2R;
        innerRightFanParameters.bx[i] = bx2R;
        innerRightFanParameters.by[i] = by2R;
        innerRightFanParameters.bz[i] = bz2R;
        innerRightFanParameters.e[i] = e2R;
    }
}


void HLLD::setFlux()
{
    double rho, u, v, w, bx, by, bz, p, e;
    double pT;
    for (int i = 0; i < nx; i++) {
        rho = componentsLeft.rho[i];
        u = componentsLeft.u[i];
        v = componentsLeft.v[i];
        w = componentsLeft.w[i];
        bx = componentsLeft.bx[i];
        by = componentsLeft.by[i];
        bz = componentsLeft.bz[i];
        p = componentsLeft.p[i];
        e = p / (gamma_mhd - 1.0)
          + 0.5 * rho * (u * u + v * v + w * w)
          + 0.5 * (bx * bx + by * by + bz * bz);
    }
}
