#include <vector>
#include "const.hpp"
#include "calculate_half_components.hpp"


class CalculateHalfComponents;

struct Components
{
    std::vector<double> rho;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> bx;
    std::vector<double> by;
    std::vector<double> bz;
    std::vector<double> p;
};

struct FanParameters
{
    std::vector<double> rho;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> bx;
    std::vector<double> by;
    std::vector<double> bz;
    std::vector<double> e;
    std::vector<double> pT;
};

struct HLLDParameters
{
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
};

struct Flux
{
    std::vector<std::vector<double>> flux;
};




class HLLD
{
private:
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

    Flux flux;
    Flux fluxLeft, fluxOuterLeft, fluxInnerLeft;
    Flux fluxRight, fluxOuterRight, fluxInnerRight;

public:
    HLLD();

    void calculateFlux(
        const std::vector<std::vector<double>> U
    );

private:
    double sign(double x);

    void setComponents(
        const std::vector<std::vector<double>> U
    );

    //void calculateHLLDParametersForOuterFan();

    void calculateHLLDParametersForMiddleFan();

    void calculateHLLDParametersForInnerFan();

    void calculateHLLDSubParameters(
        const Components components, 
        HLLDParameters hLLDParameters
    );

    void calculateHLLDParameters1(
        const Components components, 
        const HLLDParameters hLLDParameters, 
        FanParameters outerFanParameters
    );

    void calculateHLLDParameters2(
        const FanParameters outerLeftFanParameters,
        const FanParameters outerRightFanParameters, 
        FanParameters innerLeftFanParameters, 
        FanParameters innerRightFanParameters
    );

    void setFlux();
};

