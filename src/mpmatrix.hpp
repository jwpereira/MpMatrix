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
     * elements in the underlying container to be accessed using a reasonable syntax (i.e.,
     * <code>name_of_matrix(row, col)</code>).
     */
    class MpMatrix {
      private:
        std::vector<fmp_t> matrix;
        size_t dim;       ///< dimension * dimension = size [we're working with square matricies]
        fmp_shift_t shift;  ///< for keeping track of the shift/precision factor across the matrix
      public:
        /**
         * @brief Construct a new MpMatrix object
         * 
         * @param dim Size of the matrix (n for an n*n matrix)
         * @param shift Amount each of the elements are. Per \link fixedmpz \endlink
         */
        MpMatrix(size_t dim, fmp_shift_t shift) : dim(dim), shift(shift) {
            this->matrix = std::vector<fmp_t>(dim * dim, fmp_t(0, shift));
        }

        /**
         * @brief Construct a new MpMatrix object
         *
         * Consider this like a copy constructor using Iterators (so as to allow loading data into
         * MpMatrix from other containers/means). For example: for testing purposes, one could make
         * a raw array of \link fixedmpz \endlink objects (such that the array is of an n*n square
         * matrix) and have an MpMatrix object be created from that. 
         *
         * @param dim Size of the matrix (n for an n*n matrix)
         * @param shift Amount each of the elements are. Per \link fixedmpz \endlink
         * @param begin Iterator from where copying starts
         * @param end Iterator from where copying ends
         */
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

        /**
         * @brief Returns a reference to the fmp_t element at position (row, col).
         * 
         * Operator overload for () allows for easy access to elements inside the matrix using 
         * an intuitive syntax that is close to MpMatrix[row][col].
         * 
         * @param row The row where the element is located.
         * @param col The column where the element is located. 
         * @return fmp_t& The address of the element at position (row, col).
         */
        fmp_t &operator()(const size_t row, const size_t col) {
            return this->matrix[(row * this->dim) + col];
        }
        
        /**
         * @brief This is a const version of <code>&operator()</code>
         */
        fmp_t &operator()(const size_t row, const size_t col) const {
            return const_cast<fmp_t&>(this->matrix[(row * this->dim) + col]);
        }

        /**
         * @brief Get the amount the fmp_t elements are shifted by
         * 
         * @return The shift amount.
         */
        fmp_shift_t getShift() const {
            return this->shift;
        }

        /**
         * @brief Get the n-dimension of the matrix
         * 
         * For an n*n square matrix (for which \link MpMatrix \endlink represents), this returns n.
         * 
         * @return size_t n-dimension
         */
        size_t getDimension() const {
            return this->dim;
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

        friend std::ostream &operator<<(std::ostream &os, const MpMatrix &matrix);
    };

    /**
     * @brief Basic overload for the << operator to allow simple formatted-output to std::ostream
     * 
     * Example usage for \link MpMatrix \endlink named <code>matrix</code>:
     * <code>std::cout << matrix;</code>
     * 
     * <strong>Note:</strong> As of now this <strong>does not</strong> set the precision of the 
     * outputted numbers! This means that rounded results may be printed instead!
     */
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
