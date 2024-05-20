#include <vector>
#include <string>
#include "flux_solver.hpp"


class IdealMHD1D
{
private:
    FluxSolver fluxSolver;
    std::vector<std::vector<double>> U;

public:
    void setU(
        const std::vector<std::vector<double>> UInit
    ); 
    void oneStepRK2();
    void save(
        std::string directoryname, 
        std::string filenameWithoutStep, 
        int step
    );


private:
    void calculateDt();
};



