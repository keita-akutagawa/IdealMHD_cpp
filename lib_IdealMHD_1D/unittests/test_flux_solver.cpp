#include <gtest/gtest.h>
#include "../flux_solver.hpp"
#include "../const.hpp"
#include <algorithm>
#include <vector>


TEST(FluxSolverTest, ConstU)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    FluxSolver fluxSolver;
    Flux fluxF;

    fluxF = fluxSolver.getFluxF(U);

    for (int i = 0; i < nx; i++) {
        EXPECT_NE(fluxF.flux[0][i], 0.0); //flux.flux[4][i]は0.0
    }
}

