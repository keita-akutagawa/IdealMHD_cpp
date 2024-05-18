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



class HLLD
{
private:
    CalculateHalfComponents calculateHalfComponents;
    Components componentsCenter;
    FanParameters outerFanParameters;
    FanParameters innerFanParameters;
    HLLDParameters hLLDLeftParameters;
    HLLDParameters hLLDRightParameters;

public:
    HLLD();

    void calculateHLLDParametersForOuterFanParameters(
        const std::vector<std::vector<double>> U
    );

    void calculateHLLDParametersForInnerFanParameters();

    Components getOuterFanParameters();

    Components getInnerFanParameters();

private:
    void calculateHLLDSubParameters(
        const Components components, 
        HLLDParameters hLLDParameters
    );
};

