#pragma once 

#include <functional>
#include "mpmatrix.hpp"

namespace momentmp {

    inline bool generate(MpMatrix &dest, std::function<fixedmpz(size_t, size_t)> generator) {
        return false;
    }

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
