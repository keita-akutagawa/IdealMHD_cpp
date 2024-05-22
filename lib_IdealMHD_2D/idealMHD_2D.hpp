#include <vector>
#include <string>
#include "flux_solver.hpp"
#include "boundary.hpp"


class IdealMHD2D
{
private:
    FluxSolver fluxSolver;
    Flux2D fluxF, fluxG;
    std::vector<std::vector<std::vector<double>>> U;
    std::vector<std::vector<std::vector<double>>> UBar;
    Boundary boundary;

public:
    IdealMHD2D();

    void initializeU(
        const std::vector<std::vector<std::vector<double>>> UInit
    ); 

    void oneStepRK2();

    void save(
        std::string directoryname, 
        std::string filenameWithoutStep, 
        int step
    );

    std::vector<std::vector<std::vector<double>>> getU();

    void calculateDt();
};


