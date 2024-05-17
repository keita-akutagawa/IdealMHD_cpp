#include <gtest/gtest.h>
#include "../muscl.cpp"
#include <algorithm>
#include <vector>


TEST(MUSCLTest, CheckConstInputConstOutput)
{
    int n = 5;
    std::vector<double> q = std::vector<double>(5, 1.0);
    std::vector<double> qLeft = std::vector<double>(5, 0.0);
    std::vector<double> qRight = std::vector<double>(5, 0.0);

    MUSCL muscl;

    muscl.getLeftComponent(q, qLeft);
    muscl.getRightComponent(q, qRight);

    for (int i = 0; i < n; i++) {
        EXPECT_NEAR(qLeft[i], q[i], 1e-10);
        EXPECT_NEAR(qRight[i], q[i], 1e-10);
    }
}

TEST(MUSCLTest, CheckDeltaInputChangeOutput)
{
    int n = 5;
    std::vector<double> q = std::vector<double>(5, 0.0);
    std::vector<double> qLeft = std::vector<double>(5, 0.0);
    std::vector<double> qRight = std::vector<double>(5, 0.0);
    q[2] = 1.0;

    MUSCL muscl;

    muscl.getLeftComponent(q, qLeft);
    muscl.getRightComponent(q, qRight);

    for (int i = 0; i < n; i++) {
        EXPECT_NEAR(qLeft[i], q[i], 1e-10);
        EXPECT_NEAR(qRight[i], q[i], 1e-10);
    }
}



