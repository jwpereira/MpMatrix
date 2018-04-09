#pragma once

#include "fixedmpz.hpp"
#include "mpmatrix.hpp"

namespace momentmp {

    inline void momentHalfGamma(fmp_t &dest, size_t i, size_t j) {
        const auto shift = dest.getShift();
        auto fac = factorial((2 * (j + i)) + 1);
        fac <<= (1 + shift);
        fac.setShift(shift);
        dest = fac;
    }

    inline void momentInitCol(MpColumn &prev, MpColumn &curr) {
        if (prev.getSize() != curr.getSize()) {
            throw std::runtime_error("Can't initialize column from another column of different size");
        }

        const auto limit = curr.getId() + 1;
        for (size_t i = 0; i < limit; i++) {
            curr[i] = prev[i + 1];
        }

        const auto i = curr.getId();
        for (size_t j = i; j < curr.getSize(); j++) {
            momentHalfGamma(curr[j], i, j);
        }
    }

    inline void momentInitCol0(MpColumn &col) {
        const auto shift = col.getShift();
        const auto size  = col.getSize();

        const auto i = col.getId();

        for (size_t j = 0; j < size; j++) {
            momentHalfGamma(col[j], i, j);
        }
    }

    inline void momentInit(MpMatrix &matrix) {
        std::for_each(matrix.begin(), matrix.end(), momentInitCol0);
    }

    inline void momentInit_reuse(MpMatrix &matrix) {
        // std::for_each(matrix.begin(), matrix.end(), momentInitCol0);
        const auto size = matrix.getCols();
        momentInitCol0(matrix[0]);
        
        if (size > 1) {
            for (size_t i = 1; i < size; i++) {
                momentInitCol(matrix[i-1], matrix[i]);
            }
        }
    }

    inline void cholesky_decompose(MpMatrix &matrix) {
        return;
    }

    inline void cholesky_update(const MpColumn &leftCol, MpColumn &currCol) {
        return;
    }
}