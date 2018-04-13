#pragma once

#include "fixedmpz.hpp"
#include "mpmatrix.hpp"

namespace momentmp {

    inline void momentInitCol(MpArray &col) {
        const auto shift = col.getShift();
        const auto size  = col.getSize();

        const auto i = col.getId();

        for (size_t index = i; index < size; index++) {            
            col[index] = factorial((2 * (i + index)) + 1);
            col[index] <<= shift + 1;
            col[index].setShift(shift);
        }
    }

    inline void momentInit(MpMatrix &matrix) {
        std::for_each(matrix.begin(), matrix.end(), momentInitCol);
    }

    inline void cholesky_decompose(MpMatrix &matrix) {
        return;
    }

    inline void cholesky_update(const MpArray &leftCol, MpArray &currCol) {
        return;
    }
}