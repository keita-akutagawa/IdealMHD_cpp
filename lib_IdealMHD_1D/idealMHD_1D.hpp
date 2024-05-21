#include <vector>
#include <string>
#include "flux_solver.hpp"


class IdealMHD1D
{
private:
    FluxSolver fluxSolver;
    std::vector<std::vector<double>> U;
    std::vector<std::vector<double>> UBar;
    Flux fluxF;

public:
    IdealMHD1D();

    void initializeU(
        const std::vector<std::vector<double>> UInit
    ); 

    void oneStepRK2();

    void save(
        std::string directoryname, 
        std::string filenameWithoutStep, 
        int step
    );

    std::vector<std::vector<double>> getU();

private:
    void calculateDt();
};



