#include "const.hpp"
#include "ct.hpp"


class DivBCleaner
{
private:
    CT ct;

public:
    void divBClean(
        const Flux2D& flux2D,  
        const std::vector<std::vector<double>>& bxOld,
        const std::vector<std::vector<double>>& byOld, 
        std::vector<std::vector<std::vector<double>>>& U
    );
};

