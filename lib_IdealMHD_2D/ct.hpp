#include <vector>
#include "flux_solver.hpp"


class CT
{
private:
    std::vector<std::vector<double>> EzVector;

public:
    void divBClean(
        const Flux2D& flux2D,  
        const std::vector<std::vector<double>>& bxOld,
        const std::vector<std::vector<double>>& byOld, 
        std::vector<std::vector<std::vector<double>>>& U
    );
};

