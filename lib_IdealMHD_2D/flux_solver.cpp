#include "flux_solver.hpp"


Flux2D::Flux2D()
{
    fluxF = std::vector(8, std::vector(nx, std::vector<double>(ny, 0.0)));
    fluxG = std::vector(8, std::vector(nx, std::vector<double>(ny, 0.0)));
}


FluxSolver::FluxSolver()
{
    nDirection = nx;
    hLLDForF = HLLD();
    flux1DForF = Flux1D();

    nDirection = ny;
    hLLDForF = HLLD();
    flux1DForG = Flux1D();

    //バグが起きたら壊れるように
    nDirection = -1;

    tmpVectorForF = std::vector(8, std::vector<double>(nx, 0.0));
    tmpVectorForG = std::vector(8, std::vector<double>(ny, 0.0));
}


Flux2D FluxSolver::getFluxF(
    const std::vector<std::vector<std::vector<double>>> U
)
{
    for (int j = 0; j < ny; j++) {
        for (int comp = 0; comp < 8; comp++) {
            for (int i = 0; i < nx; i++) {
                tmpVectorForF[comp][i] = U[comp][i][j];
            }
        }

        hLLDForF.calculateFlux(tmpVectorForF);
        flux1DForF = hLLDForF.getFlux();

        for (int comp = 0; comp < 8; comp++) {
            for (int i = 0; i < nx; i++) {
                flux2D.fluxF[comp][i][j] = flux1DForF.flux[comp][i];
            }
        }
    }

    return flux2D;
}


Flux2D FluxSolver::getFluxG(
    const std::vector<std::vector<std::vector<double>>> U
)
{
    for (int i = 0; i < nx; i++) {
        for (int comp = 0; comp < 8; comp++) {
            for (int j = 0; j < ny; j++) {
                tmpVectorForG[comp][j] = U[comp][i][j];
            }
        }

        hLLDForG.calculateFlux(tmpVectorForG);
        flux1DForG = hLLDForG.getFlux();

        for (int comp = 0; comp < 8; comp++) {
            for (int j = 0; j < ny; j++) {
                flux2D.fluxG[comp][i][j] = flux1DForG.flux[comp][j];
            }
        }
    }

    return flux2D;
}
