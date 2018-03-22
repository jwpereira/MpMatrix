#pragma once

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

    /**
     * @brief Matrix class based class focused on containing and fmp_t elements.
     * 
     * MpMatrix is a square-matrix container-wrapper for fmp_t elements. It allows for individual
     * elements in the underlying container to be accessed using a reasonable syntax 
     * (i.e., <code>name_of_matrix(row, col)</code>).
     */
    class MpMatrix {
      private:
        std::vector<fmp_t> matrix;
        size_t dim;       ///< dimension * dimension = size [we're working with square matricies]
        fmp_shift_t shift;  ///< for keeping track of the shift/precision factor across the matrix
      public:
        MpMatrix(size_t dim, fmp_shift_t shift) : dim(dim), shift(shift) {
            this->matrix = std::vector<fmp_t>(dim * dim, fmp_t(0, shift));
        }

        template <typename Iterator>
        MpMatrix(size_t dim, fmp_shift_t shift, Iterator begin, Iterator end) : MpMatrix(dim, shift) {
            //std::copy(begin, end, this->matrix.begin());
            auto iter = this->matrix.begin();
            std::for_each(begin, end, [&](auto &val) {
                iter->setNumber(val);
                iter->setshift(val.getShift());
                iter++;
            });
        }

        fmp_t &operator()(const size_t row, const size_t col) {
            return this->matrix[(row * this->dim) + col];
        }

        fmp_t &operator()(const size_t row, const size_t col) const {
            return const_cast<fmp_t&>(this->matrix[(row * this->dim) + col]);
        }

        fmp_shift_t getShift() const {
            return this->shift;
        }

        size_t getDimension() const {
            return this->dim;
        }

        decltype(auto) begin() const {
            return this->matrix.begin();
        }

        decltype(auto) cbegin() const {
            return this->matrix.cbegin();
        }

        decltype(auto) end() const {
            return this->matrix.end();
        }

        decltype(auto) cend() const {
            return this->matrix.cend();
        }

        friend std::ostream &operator<<(std::ostream &os, const MpMatrix &matrix);
    };

    inline std::ostream &operator<<(std::ostream &os, const MpMatrix &matrix) {        
        //Capture the initial flags of the output stream
        std::ios::fmtflags initialFlags(os.flags());

        size_t counter = 0;
        for (auto &num : matrix) {
            counter++;
            os << num << '\t';
            if (counter % matrix.getDimension() == 0) {
                os << '\n';
            }
        }

        //Restore the initial flags of the output streams
        os.flags(initialFlags);

        return os << std::endl;
    }

}
