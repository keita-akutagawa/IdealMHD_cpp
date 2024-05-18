#include <gtest/gtest.h>
#include "../calculate_half_components.cpp"
#include "../const.cpp"
#include <algorithm>
#include <vector>


TEST(CalculateHalfComponentsTest, Constructor)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    CalculateHalfComponents calculateHalfComponents;

    calculateHalfComponents.getPhysicalParameters(U);

    for (int i = 0; i < nx; i++) {
        EXPECT_EQ(U[0][i], 1.0);
    }
}



