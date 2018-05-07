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

    class MpArray {
      private:
        std::vector<fmp_t> col;
        size_t dim, id;
        fmp_shift_t shift;
      public:
        MpArray(size_t dim, fmp_shift_t shift, size_t id=-1) noexcept 
                : dim(dim), id(id), shift(shift) {
            this->col = std::vector<fmp_t>(dim, fmp_t(0, shift));
        }

        MpArray(const MpArray &other) = default;
        MpArray(MpArray &&other) = default;
        
        size_t getId() {
            return this->id;
        }

        size_t getId() const {
            return this->id;
        }

        size_t size() {
            return this->dim;
        }
        
        size_t size() const {
            return this->dim;
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

        /**
         * @brief Returns a begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) begin() {
            return this->col.begin();
        }

        /**
         * @brief Returns a begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) begin() const {
            return this->col.begin();
        }

        /**
         * @brief Returns a constant begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) cbegin() const {
            return this->col.cbegin();
        }

        /**
         * @brief Returns a begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) end() {
            return this->col.end();
        }

        /**
         * @brief Returns a begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) end() const {
            return this->col.end();
        }

        /**
         * @brief Returns a constant end() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) cend() const {
            return this->col.cend();
        }

        friend std::ostream &operator<<(std::ostream &os, MpArray &array);
    };

    inline std::ostream &operator<<(std::ostream &os, MpArray &array) {
        for (auto &e : array) {
            os << e << ' ';
        }

        return os << std::endl;
    }

    enum MpMatrixMode : bool { ROW_ORIENTED=true , COL_ORIENTED=false };

    class MpMatrix {
      private:
        std::vector<MpArray> matrix;
        size_t dim;         ///< dimension * dimension = rows [we're working with square matricies]
        fmp_shift_t shift;  ///< for keeping track of the shift/precision factor across the matrix
        MpMatrixMode mode;
      public:
        MpMatrix(size_t dim, fmp_shift_t shift, MpMatrixMode mode = COL_ORIENTED) noexcept 
                : dim(dim),  shift(shift) {
            this->matrix = std::vector<MpArray>(dim, MpArray(dim, shift, 0));
            for (size_t i = 0; i < dim; i++) {
                this->matrix[i].setId(i);
            }
            this->mode = mode;
        }

        MpMatrix(const MpMatrix &other) = default;
        MpMatrix(MpMatrix &&other) = default;

        size_t getDim() {
            return this->dim;
        }

        size_t getDim() const {
            return this->dim;
        }

        fmp_shift_t getShift() {
            return this->shift;
        }

        fmp_shift_t getShift() const {
            return this->shift;
        }

        decltype(auto) getMode() const {
            return this->mode;
        }

        /**
         * @brief Returns a begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) begin() {
            return this->matrix.begin();
        }

        /**
         * @brief Returns a begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) begin() const {
            return this->matrix.begin();
        }

        /**
         * @brief Returns a constant begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) cbegin() const {
            return this->matrix.cbegin();
        }

        /**
         * @brief Returns a begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) end() {
            return this->matrix.end();
        }

        /**
         * @brief Returns a begin() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) end() const {
            return this->matrix.end();
        }

        /**
         * @brief Returns a constant end() iterator from the internal vector class.
         * 
         * @return decltype(auto) 
         */
        decltype(auto) cend() const {
            return this->matrix.cend();
        }

        void setMode(MpMatrixMode mode) {
            this->mode = mode;
        }

        MpArray &operator[](const size_t col) {
            return this->matrix[col];
        }

        MpArray &operator[](const size_t col) const {
            return const_cast<MpArray&>(this->matrix[col]);
        }

        void clear() {
            auto zero = 0^fmpshift(this->getShift());
            for (auto &array : *this) {
                for (auto &elem : array) {
                    elem = zero;
                }
            }
        }

        void dumpVecDouble(std::vector<double> &dest) const {
            auto dim = this->getDim();

            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    auto &elem = matrix[i][j];
                    dest[i * dim + j] = elem.to_mpf().get_d();
                }
            }
        }

        friend std::ostream &operator<<(std::ostream &os, const MpMatrix &mp);
    };

    inline std::ostream &operator<<(std::ostream &os, const MpMatrix &mp) {
        std::ios::fmtflags initialFlags(os.flags());

        auto mode = mp.getMode();

        if (mode == COL_ORIENTED) {
            for (size_t i = 0; i < mp.getDim(); i++) {
                for (size_t j = 0; j < mp.getDim(); j++) {
                    os << mp[j][i] << '\t';
                }
                os << '\n';
            }
        } else if (mode == ROW_ORIENTED) {
            for (size_t i = 0; i < mp.getDim(); i++) {
                for (size_t j = 0; j < mp.getDim(); j++) {
                    os << mp[i][j] << '\t';
                }
                os << '\n';
            }
        }

        os.flags(initialFlags);

        return os;
    }

    inline void multiply(const MpMatrix &multiplicand, const MpMatrix &multiplier, MpMatrix &product) {
        if (multiplicand.getDim() != multiplier.getDim()) {
            throw std::runtime_error("Unable to multiply MpMatricies of different dimensions");
        }
        
        if (multiplicand.getDim() != product.getDim()) {
            throw std::runtime_error("Destination (Product) matrix must be of same dimension as factors");
        }

        auto dim = multiplicand.getDim();
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                for (size_t k = 0; k < dim; k++) {
                    product[i][j] += (multiplicand[i][k] * multiplier[k][j]);
                }
            }
        }
    }

    inline void transpose(MpMatrix &matrix) {
        auto dim = matrix.getDim();

        for (size_t n = 0; n < (dim - 1); n++) {
            for (size_t m = (n + 1); m < dim; m++) {
                auto temp = matrix[m][n];
                matrix[m][n] = matrix[n][m];
                matrix[n][m] = temp;
            }
        }
    }

    inline void reorient(MpMatrix &matrix) {
        transpose(matrix);

        auto mode = matrix.getMode();
        if (mode == COL_ORIENTED) {
            matrix.setMode(ROW_ORIENTED);
        } else if (mode == ROW_ORIENTED) {
            matrix.setMode(COL_ORIENTED);
        }
    }

    inline void reflect(MpMatrix &matrix) {
        auto dim = matrix.getDim();

        for (size_t n = 0; n < (dim - 1); n++) {
            for (size_t m = (n + 1); m < dim; m++) {
                auto temp = matrix[m][n];
                matrix[m][n] = matrix[n][m];
                // matrix[n][m] = temp;
            }
        }
    }
}
