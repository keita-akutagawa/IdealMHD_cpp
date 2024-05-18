#include <vector>
#include "const.hpp"
#include "muscl.hpp"


class MUSCL;

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


class CalculateHalfComponents
{
private:
    MUSCL* muscl;

    Components componentsCenter;
    Components componentsLeft;
    Components componentsRight;

    std::vector<double> tmpVector;

public:
    CalculateHalfComponents();

    void getPhysicalParameters(
        const std::vector<std::vector<double>> U
    );

    void calculateLeftComponents();

    void calculateRightComponents();

    Components getCenterComponents();

    Components getLeftComponents();

    Components getRightComponents();
};

