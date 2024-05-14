#include <vector>

class MUSCL
{
private:

public:
    std::vector<double> getLeftComponent(
        const std::vector<double> q, 
        std::vector<double> qLeft
    );
    std::vector<double> getRightComponent(
        const std::vector<double> q, 
        std::vector<double> qRight
    );
};

