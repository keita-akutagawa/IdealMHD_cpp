#include <gtest/gtest.h>
#include "../hlld.hpp"
#include "../const.hpp"
#include <algorithm>
#include <vector>


TEST(HLLDTest, Constructor)
{
    HLLD hlld;
    Flux flux;

    flux = hlld.getFlux();

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            EXPECT_EQ(flux.flux[comp][i], 0.0);
        }
    }
}


TEST(HLLDTest, CalculateFluxSetComponents)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    HLLD hlld;
    Components components;

    hlld.calculateFlux(U);
    components = hlld.getComponents();

    for (int i = 0; i < nx; i++) {
        EXPECT_EQ(components.rho[i], 1.0); 
    }
}

TEST(HLLDTest, CalculateFluxSetFanParameters)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    HLLD hlld;
    FanParameters fanParameters;

    hlld.calculateFlux(U);
    fanParameters = hlld.getFanParameters();

    EXPECT_EQ(fanParameters.rho[10], 1.0); 
}

