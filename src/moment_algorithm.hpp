#pragma once

#include "fixedmpz.hpp"
#include "mpmatrix.hpp"

namespace momentmp {

    inline void momentInitCol(MpArray &col) {
        const auto shift = col.getShift();
        const auto size  = col.size();

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
        auto dim = matrix.getRows();

        // procCol is the column currently being applied to every other column
        for (auto &procCol : matrix) {
            auto id = procCol.getId();
            auto start = id + 1;
            MpArray orig(procCol);

            // Replace procCol with all the values under diagonal with those values divided by
            // diagonal
            auto &diagonal = orig[id];
            for (size_t row = start; row < dim; row++) {
                procCol[row] /= diagonal;
            }

            // Apply procCol to all other columns to its right
            for (size_t col = start; col < dim; col++) {
                MpArray &destCol = matrix[col];

                // Going down the rows for each col, z' = z - yx
                for (size_t row = col; row < dim; row++) {
                    auto &z = destCol[row];
                    auto &y = orig[col];
                    auto &x = procCol[row];

                    z = z - (y * x);
                }
            }
        }
    }

    inline void cholesky_update(const MpArray &leftCol, MpArray &currCol) {
        return;
    }
}