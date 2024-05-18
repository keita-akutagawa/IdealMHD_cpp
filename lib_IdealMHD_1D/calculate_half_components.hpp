#include <vector>
#include "const.hpp"
#include "muscl.hpp"


class MUSCL;

class CalculateHalfComponents
{
private:
    MUSCL* muscl;

    std::vector<double> rho;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> bx;
    std::vector<double> by;
    std::vector<double> bz;
    std::vector<double> p;

    std::vector<double> rhoL;
    std::vector<double> uL;
    std::vector<double> vL;
    std::vector<double> wL;
    std::vector<double> bxL;
    std::vector<double> byL;
    std::vector<double> bzL;
    std::vector<double> pL;

    std::vector<double> rhoR;
    std::vector<double> uR;
    std::vector<double> vR;
    std::vector<double> wR;
    std::vector<double> bxR;
    std::vector<double> byR;
    std::vector<double> bzR;
    std::vector<double> pR;

    std::vector<double> tmpVector;

public:
    CalculateHalfComponents();

    void getPhysicalParameters(
        const std::vector<std::vector<double>> U
    );

    void calculateLeftComponents();

    void calculateRightComponents();
};

