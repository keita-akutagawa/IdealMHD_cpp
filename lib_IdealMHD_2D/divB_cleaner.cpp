#include "divB_cleaner.hpp"


void DivBCleaner::divBClean(
    const Flux2D& flux2D,  
    const std::vector<std::vector<double>>& bxOld,
    const std::vector<std::vector<double>>& byOld, 
    std::vector<std::vector<std::vector<double>>>& U
)
{
    ct.divBClean(flux2D, bxOld, byOld, U);
}


