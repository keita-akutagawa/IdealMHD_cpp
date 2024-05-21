#include <gtest/gtest.h>
#include "../idealMHD_1D.hpp"
#include "../const.hpp"
#include <algorithm>
#include <vector>


TEST(IdealMHD1D, Initialize)
{
    std::vector<std::vector<double>> UInit(8, std::vector<double>(nx, 1.0));
    std::vector<std::vector<double>> U(8, std::vector<double>(nx, 0.0));

    IdealMHD1D idealMHD1D;

    idealMHD1D.initializeU(UInit);

    U = idealMHD1D.getU();

    for (int i = 0; i < nx; i++) {
        EXPECT_EQ(U[0][i], 1.0);
    }
}


