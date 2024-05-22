#include <vector>
#include "const.hpp"
#include "calculate_half_components.hpp"


class CalculateHalfComponents;


struct FanParameters
{
    int nSize; 

    std::vector<double> rho;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> bx;
    std::vector<double> by;
    std::vector<double> bz;
    std::vector<double> e;
    std::vector<double> pT;

    FanParameters() : nSize(0) {}
    FanParameters(int nDirection);
};

struct HLLDParameters
{
    int nSize; 

    std::vector<double> pT;
    std::vector<double> pT1;
    std::vector<double> pT2;
    std::vector<double> e;
    std::vector<double> cs;
    std::vector<double> ca;
    std::vector<double> va;
    std::vector<double> cf;
    std::vector<double> S;
    std::vector<double> S1;
    std::vector<double> SM;

    HLLDParameters() : nSize(0) {}
    HLLDParameters(int nDirection);
};

struct Flux1D
{
    int nSize; 
    std::vector<std::vector<double>> flux;

    Flux1D() : nSize(0) {}
    Flux1D(int nDirection);
};




class HLLD
{
private:
    int nSize; 

    CalculateHalfComponents calculateHalfComponents;
    Components componentsCenter;
    Components componentsLeft;
    Components componentsRight;
    FanParameters outerLeftFanParameters;
    FanParameters outerRightFanParameters;
    FanParameters middleLeftFanParameters;
    FanParameters middleRightFanParameters;
    FanParameters innerLeftFanParameters;
    FanParameters innerRightFanParameters;
    HLLDParameters hLLDLeftParameters;
    HLLDParameters hLLDRightParameters;

    Flux1D flux;
    Flux1D fluxOuterLeft, fluxMiddleLeft, fluxInnerLeft;
    Flux1D fluxOuterRight, fluxMiddleRight, fluxInnerRight;

public:
    HLLD() : nSize(0) {}
    HLLD(int nDirection);

    void calculateFlux(
        const std::vector<std::vector<double>> U
    );

    Components getLeftComponents();
    Components getRightComponents();

    FanParameters getOuterLeftFanParameters();
    FanParameters getOuterRightFanParameters();
    FanParameters getMiddleLeftFanParameters();
    FanParameters getMiddleRightFanParameters();
    FanParameters getInnerLeftFanParameters();
    FanParameters getInnerRightFanParameters();

    HLLDParameters getHLLDLeftParameters();
    HLLDParameters getHLLDRightParameters();

    Flux1D getFlux();

private:
    double sign(double x);

    void setComponents(
        const std::vector<std::vector<double>> U
    );

    void calculateHLLDParametersForOuterFan();

    void calculateHLLDParametersForMiddleFan();

    void calculateHLLDParametersForInnerFan();

    void setFanParametersFromComponents(
        const Components components, 
        FanParameters& fanParameters
    );

    void calculateHLLDSubParametersForMiddleFan(
        const Components components, 
        const FanParameters outerFanParameters, 
        HLLDParameters& hLLDParameters
    );

    void calculateHLLDParameters1(
        const FanParameters outerFanParameters, 
        const HLLDParameters hLLDParameters, 
        FanParameters& middleFanParameters
    );

    void calculateHLLDParameters2();

    void setFlux(
        const FanParameters fanParameters, 
        Flux1D& flux
    );
};

