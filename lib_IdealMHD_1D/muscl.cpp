#include <algorithm>
#include <vector>
#include "minmod.hpp"
#include "muscl.hpp"


std::vector<double> MUSCL::getLeftComponent(
    const std::vector<double> q, 
    std::vector<double> qLeft
)
{

    for (int i = 1; i < q.size() - 1; i++) {
        qLeft[i] = q[i] + 0.5 * minmod(q[i] - q[i-1], q[i+1] - q[i]);
    }

    //周期境界条件
    qLeft[0] = q[0] + 0.5 * minmod(
        q[0] - q[q.size()-1], q[1] - q[0]
        );
    qLeft[q.size()-1] = q[q.size()-1] + 0.5 * minmod(
        q[q.size()-1] - q[q.size()-2], q[0] - q[q.size()-1]
        );
    
    return qLeft;

}

std::vector<double> MUSCL::getRightComponent(
    const std::vector<double> q, 
    std::vector<double> qRight
)
{

    for (int i = 1; i < q.size() - 1; i++) {
        qRight[i] = q[i] - 0.5 * minmod(q[i+1] - q[i], q[i+2] - q[i+1]);
    }

    //周期境界条件
    qRight[0] = q[0] + 0.5 * minmod(
        q[1] - q[0], q[2] - q[1]
        );
    qRight[q.size()-1] = q[q.size()-1] + 0.5 * minmod(
        q[0] - q[q.size()-1], q[1] - q[0]
        );
    
    return qRight;

}

