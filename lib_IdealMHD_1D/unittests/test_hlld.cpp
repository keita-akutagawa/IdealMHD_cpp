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

TEST(HLLDTest, CalculateFluxConstructor)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 0.0));

    HLLD hlld;
    Flux flux;

    hlld.calculateFlux(U);
    flux = hlld.getFlux();

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            EXPECT_EQ(flux.flux[comp][i], 0.0);
        }
    }
}

TEST(HLLDTest, CalculateFluxConstU)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));


    HLLD hlld;
    Flux flux;

    hlld.calculateFlux(U);
    flux = hlld.getFlux();

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            EXPECT_EQ(flux.flux[comp][i], 0.0);
        }
    }
}

TEST(HLLDTest, CalculateFluxNotConstU)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            if ((i > 10) && (i < nx - 10)) {
                U[comp][i] = 2.0;
            }
        }
    }


    HLLD hlld;
    Flux flux;

    hlld.calculateFlux(U);
    flux = hlld.getFlux();

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            EXPECT_NE(flux.flux[comp][i], 0.0);
        }
    }
}

