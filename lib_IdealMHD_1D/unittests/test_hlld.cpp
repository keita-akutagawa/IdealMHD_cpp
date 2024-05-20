#include <gtest/gtest.h>
#include "../hlld.hpp"
#include "../const.hpp"
#include <algorithm>
#include <vector>


TEST(HLLDTest, CheckStructConstructor)
{
    HLLD hLLD;
    Components components;
    FanParameters fanParameters;
    HLLDParameters hLLDParameters;
    Flux flux;

    components = hLLD.getLeftComponents();
    fanParameters = hLLD.getOuterLeftFanParameters();
    hLLDParameters = hLLD.getHLLDLeftParameters();
    flux = hLLD.getFlux();

    for (int i = 0; i < nx; i++) {
        EXPECT_EQ(components.rho[i], 0.0);
        EXPECT_EQ(fanParameters.pT[i], 0.0);
        EXPECT_EQ(hLLDParameters.ca[i], 0.0);
    }

    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            EXPECT_EQ(flux.flux[comp][i], 0.0);
        }
    }
}

TEST(HLLDTest, SetComponents)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    HLLD hLLD;
    Components components;
    FanParameters fanParameters;
    HLLDParameters hLLDParameters;
    Flux flux;

    hLLD.calculateFlux(U);

    components = hLLD.getLeftComponents();

    for (int i = 0; i < nx; i++) {
        EXPECT_EQ(components.rho[i], 1.0); 
        EXPECT_EQ(components.bx[i], 1.0); 
    }
}

TEST(HLLDTest, CalculateFanParameters)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    HLLD hLLD;
    Components components;
    FanParameters fanParameters;
    HLLDParameters hLLDParameters;
    Flux flux;

    hLLD.calculateFlux(U);
    
    fanParameters = hLLD.getOuterLeftFanParameters();

    for (int i = 0; i < nx; i++) {
        EXPECT_EQ(fanParameters.rho[i], 1.0); 
        EXPECT_NE(fanParameters.pT[i], 0.0); 
    }
}

TEST(HLLDTest, CalculateHLLDParameters)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    HLLD hLLD;
    Components components;
    FanParameters fanParameters;
    HLLDParameters hLLDParameters;
    Flux flux;

    hLLD.calculateFlux(U);
    
    hLLDParameters = hLLD.getHLLDLeftParameters();

    for (int i = 0; i < nx; i++) {
        EXPECT_NE(hLLDParameters.pT[i], 0.0); 
        EXPECT_NE(hLLDParameters.S[i], 0.0); 
    }
}

TEST(HLLDTest, CalculateFlux)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    HLLD hLLD;
    Components components;
    FanParameters fanParameters;
    HLLDParameters hLLDParameters;
    Flux flux;

    hLLD.calculateFlux(U);
    
    flux = hLLD.getFlux();

    for (int i = 0; i < nx; i++) {
        EXPECT_NE(flux.flux[0][i], 0.0); //flux.flux[4][i]は0.0
    }
}



