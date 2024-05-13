#include <gtest/gtest.h>
#include "../minmod.cpp"
#include <algorithm>


TEST(MinmodTest, Zero)
{
    EXPECT_LE(std::abs(minmod(1.0, -0.1)), 1e-10);
    EXPECT_LE(std::abs(minmod(-1.0, 0.1)), 1e-10);
}

TEST(MinmodTest, Positive)
{
    EXPECT_LE(std::abs(minmod(1.0, 0.1) - 0.1), 1e-10);
    EXPECT_LE(std::abs(minmod(0.1, 1.0) - 0.1), 1e-10);
    EXPECT_LE(std::abs(minmod(-1.0, -0.1) - (-0.1)), 1e-10);
    EXPECT_LE(std::abs(minmod(-0.1, -1.0) - (-0.1)), 1e-10);
}


