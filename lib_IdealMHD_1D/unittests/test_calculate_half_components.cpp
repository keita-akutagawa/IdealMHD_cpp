#include <gtest/gtest.h>
#include "../calculate_half_components.hpp"
#include "../const.hpp"
#include <algorithm>
#include <vector>


TEST(CalculateHalfComponentsTest, Constructor)
{
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 1.0));

    CalculateHalfComponents calculateHalfComponents;

    calculateHalfComponents.getPhysicalParameters(U);

    Components components;
    components = calculateHalfComponents.getCenterComponents();

    for (int i = 0; i < nx; i++) {
        EXPECT_EQ(components.rho[i], 1.0);
    }
}



