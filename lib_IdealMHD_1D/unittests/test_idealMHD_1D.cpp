#include <gtest/gtest.h>
#include "../idealMHD_1D.hpp"
#include "../const.hpp"
#include <algorithm>
#include <vector>


TEST(IdealMHD1D, initializeU)
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

TEST(IdealMHD1D, calculateDt)
{
    double rho, u, v, w, bx, by, bz, p, e;
    rho = u = v = w = bx = by = bz = p = 1.0;
    e = p / (gamma_mhd - 1.0)
      + 0.5 * rho * (u * u + v * v + w * w)
      + 0.5 * (bx * bx + by * by + bz * bz);

    std::vector<std::vector<double>> UInit(8, std::vector<double>(nx, 0.0));
    for (int i = 0; i < nx; i++) {
        UInit[0][i] = rho;
        UInit[1][i] = u;
        UInit[2][i] = v;
        UInit[3][i] = w;
        UInit[4][i] = bx;
        UInit[5][i] = by;
        UInit[6][i] = bz;
        UInit[7][i] = e;
    }

    IdealMHD1D idealMHD1D;

    idealMHD1D.initializeU(UInit);

    idealMHD1D.calculateDt();

    std::cout << dt << std::endl;
}


