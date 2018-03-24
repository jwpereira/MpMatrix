#pragma once 

#include "mpmatrix.hpp"

namespace momentmp {
    /**
     * @brief Cholesky decomposition for MpMatrix
     * 
     * @param[in] initial Address of source matrix
     * @param[out] lower Address of destination matrix
     * @return true (Exceptions aside) was able to perform decomposition
     * @return false (Exceptions aside) was unable to perform decomposition
     */
    bool cholesky(const MpMatrix &initial, MpMatrix &lower);
}
