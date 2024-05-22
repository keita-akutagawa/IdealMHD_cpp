#include "boundary.hpp"
#include "const.hpp"


void Boundary::periodicBoundary(
    std::vector<std::vector<std::vector<double>>>& U
)
{
    //他のクラスの関数は周期境界を想定して組んでいるので
    //何もしなくてよい
    return;
}


void Boundary::symmetricBoundary2ndX(
    std::vector<std::vector<std::vector<double>>>& U
)
{
    for (int comp = 0; comp < 8; comp++) {
        for (int j = 0; j < ny; j++) {
            U[comp][0][j] = U[comp][2][j];
            U[comp][1][j] = U[comp][2][j];
            U[comp][nx-1][j] = U[comp][nx-3][j];
            U[comp][nx-2][j] = U[comp][nx-3][j];
        }
    }
}


void Boundary::symmetricBoundary2ndY(
    std::vector<std::vector<std::vector<double>>>& U
)
{
    for (int comp = 0; comp < 8; comp++) {
        for (int i = 0; i < nx; i++) {
            U[comp][i][0] = U[comp][i][2];
            U[comp][i][1] = U[comp][i][2];
            U[comp][i][nx-1] = U[comp][i][nx-3];
            U[comp][i][nx-2] = U[comp][i][nx-3];
        }
    }
}

