#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include <vector>
#include "fixedmpz.hpp"

/**
 * @brief Namespace for the Multiple Precision Matrix Project
 * 
 * Classes and functions in this namespace are particularly geared towards the needs of the 
 * HankelHacker project. 
 */
namespace momentmp {
    using fmp_t = fixedmpz; ///< convenience alias to the underlying number type used for operations
    using fmp_shift_t = fmp_shift_t;   ///< alias to the scaling type really used for fixedmpz

    class MpColumn {
      private:
        std::vector<fmp_t> col;
        size_t rows, id;
        fmp_shift_t shift;
      public:
        MpColumn(size_t rows, size_t id, fmp_shift_t shift) noexcept 
                : rows(rows), id(id), shift(shift) {
            this->col = std::vector<fmp_t>(rows, fixedmpz(0, shift));
        }

        MpColumn(const MpColumn &other) = default;
        MpColumn(MpColumn &&other) = default;
        
        size_t getId() {
            return this->id;
        }

        size_t getId() const {
            return this->id;
        }

        size_t getSize() {
            return this->rows;
        }
        
        size_t getSize() const {
            return this->rows;
        }

        fmp_shift_t getShift() {
            return this->shift;
        }

        fmp_shift_t getShift() const {
            return this->shift;
        }

        void setId(size_t id) {
            this->id = id;
        }

        fmp_t &operator[](size_t row) {
            return this->col[row];
        }

        fmp_t &operator[](size_t row) const {
            return const_cast<fmp_t&>(this->col[row]);
        }
    };

    class MpMatrix {
      private:
        std::vector<MpColumn> matrix;
        size_t rows, cols;       ///< dimension * dimension = rows [we're working with square matricies]
        fmp_shift_t shift;  ///< for keeping track of the shift/precision factor across the matrix
      public:
        MpMatrix(size_t rows, size_t cols, fmp_shift_t shift) noexcept 
                : rows(rows), cols(cols), shift(shift) {
            this->matrix = std::vector<MpColumn>(cols, MpColumn(rows, 0, shift));
        }

        MpMatrix(const MpMatrix &other) = default;
        MpMatrix(MpMatrix &&other) = default;

        size_t getRows() {
            return this->rows;
        }

        size_t getRows() const {
            return this->rows;
        }

        size_t getCols() {
            return this->cols;
        }

        size_t getCols() const {
            return this->cols;
        }

        MpColumn &operator[](const size_t col) {
            return this->matrix[col];
        }

        MpColumn &operator[](const size_t col) const {
            return const_cast<MpColumn&>(this->matrix[col]);
        }

        friend std::ostream &operator<<(std::ostream &os, const MpMatrix &mp);
    };

    inline std::ostream &operator<<(std::ostream &os, const MpMatrix &mp) {
        std::ios::fmtflags initialFlags(os.flags());

        for (size_t i = 0; i < mp.getRows(); i++) {
            for (size_t j = 0; j < mp.getCols(); j++) {
                os << mp[i][j] << '\t';
            }
            os << '\n';
        }

        os.flags(initialFlags);

        return os;
    }
}
