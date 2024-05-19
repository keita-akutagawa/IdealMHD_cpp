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

    for (int i = 0; i < nx; i++) {
        EXPECT_EQ(flux.flux[0][i], 0.0);
    }
}

