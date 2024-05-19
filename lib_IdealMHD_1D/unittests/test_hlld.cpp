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

