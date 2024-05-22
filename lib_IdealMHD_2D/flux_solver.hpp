#include "const.hpp"
#include "hlld.hpp"


struct Flux2D
{
    std::vector<std::vector<std::vector<double>>> fluxF;
    std::vector<std::vector<std::vector<double>>> fluxG;

    Flux2D();
};


class HLLD;

class FluxSolver
{
private:
    HLLD hLLDForF, hLLDForG;
    Flux1D flux1DForF, flux1DForG;
    Flux2D flux2D;
    std::vector<std::vector<double>> tmpVectorForF;
    std::vector<std::vector<double>> tmpVectorForG;
    std::vector<std::vector<std::vector<double>>> tmpFlux;

public:
    FluxSolver();

    Flux2D getFluxF(
        const std::vector<std::vector<std::vector<double>>> U
    );

    Flux2D getFluxG(
        const std::vector<std::vector<std::vector<double>>> U
    );

    void setFluxGToProperPosition();
};


