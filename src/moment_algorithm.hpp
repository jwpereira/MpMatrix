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
        #pragma omp parallel for schedule(dynamic, 1)
        for (auto it = matrix.begin(); it < matrix.end(); it++) {
            momentInitCol(*it);
        }
    }

    inline void cholesky_decompose(MpMatrix &matrix) {
        auto dim = matrix.getDim();

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
                    const auto &y = orig[col];
                    const auto &x = procCol[row];

                    z = z - (y * x);
                }
            }
        }
    }

    inline void extract_diagonal(MpMatrix &src, MpArray &dest, bool replace=true) {
        if (dest.size() != src.getDim()) {
            throw std::runtime_error("Cannot extract to different size array");
        }

        auto shift = fmpshift(src.getShift());
        auto one = (1^shift);

        for (auto &proc : src) {
            auto id = proc.getId();
            dest[id] = proc[id];

            if (replace) {
                proc[id] = one;
            }
        }
    }

    inline void impose_diagonal(MpArray &diagonal, MpMatrix &dest) {
        if (diagonal.size() != dest.getDim()) {
            throw std::runtime_error("Cannot impose diagonal on different dimension matrix");
        }

        for (auto &proc : dest) {
            auto id = proc.getId();
            proc[id] = diagonal[id];
        }
    }

    inline void invert(MpMatrix &matrix) {
        auto dim = matrix.getDim();

        // procRow is the row currently being applied to all other rows
        for (auto &procRow : matrix) {
            auto id = procRow.getId();
            auto start = id + 1;

            for (size_t row = start; row < dim; row++) {
                auto &destRow = matrix[row];
                auto scale = destRow[id];

                for (size_t i = 0; i < dim; i++) {
                    if (i == id) {
                        destRow[i] = -destRow[i];
                    } else {
                        destRow[i] = destRow[i] - (procRow[i] * scale);
                    }
                }
            }
        }
    }

    inline void invert_diagonal(MpArray &diagonal) {
        auto shift = fmpshift(diagonal.getShift());
        auto one = 1^shift;

        // std::for_each(diagonal.begin(), diagonal.end(), [&](auto &element) {
        //     element = one / element;
        // });

        for (auto &elem : diagonal) {
            elem = one / elem;
        }
    }
}