#include <vector>
#include <string>
#include "const.hpp"
#include "flux_solver.hpp"
#include "boundary.hpp"
#include "divB_cleaner.hpp"


class IdealMHD2D
{
private:
    FluxSolver fluxSolver;
    Flux2D flux2D;
    std::vector<std::vector<std::vector<double>>> U;
    std::vector<std::vector<std::vector<double>>> UBar;
    Boundary boundary;
    std::vector<std::vector<double>> bxOld;
    std::vector<std::vector<double>> byOld;
    DivBCleaner divBCleaner;

public:
    IdealMHD2D() :
        U(8, std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0))), 
        UBar(8, std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0)))
        {}

    void initializeU(
        const std::vector<std::vector<std::vector<double>>>& UInit
    ); 

    void oneStepRK2();

    void save(
        std::string directoryname, 
        std::string filenameWithoutStep, 
        int step
    );

    std::vector<std::vector<std::vector<double>>> getU();

    void calculateDt();

    bool checkCalculationIsCrashed();
};



