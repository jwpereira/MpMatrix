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

        MpMatrix(const MpMatrix &other) = default;
        MpMatrix(MpMatrix &&other) = default;

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

        MpMatrix &operator=(const MpMatrix &other) = default;
        MpMatrix &operator=(MpMatrix &&other) = default;

        MpMatrix &operator+=(const MpMatrix &addend) {
            if (addend.dim != this->dim) {
                throw std::runtime_error("Unable to add together matricies of different dimensions");
            }
            
            std::transform(this->begin(), this->end(), addend.cbegin(), this->begin(), std::plus<fmp_t>());

            return *this;
        }

        MpMatrix &operator+=(const fmp_t &addend) {
            std::for_each(this->begin(), this->end(), [&](auto &fmp) {
                fmp += addend;
            });

            return *this;
        }

        MpMatrix &operator-=(const MpMatrix &subtrahend) {
            if (subtrahend.dim != this->dim) {
                throw std::runtime_error("Unable to subtract matricies of different dimensions");
            }
            
            std::transform(this->begin(), this->end(), subtrahend.cbegin(), this->begin(), std::minus<fmp_t>());

            return *this;
        }

        MpMatrix &operator-=(const fmp_t &subtrahend) {
            std::for_each(this->begin(), this->end(), [&](auto &fmp) {
                fmp -= subtrahend;
            });

            return *this;
        }

        MpMatrix &operator*=(const MpMatrix &multiplicand) {
            if (multiplicand.dim != this->dim) {
                throw std::runtime_error("MpMatrix multiplication requires same dimensions");
            }
            
            MpMatrix &multiplier = *this;
            MpMatrix ret(this->dim, this->shift);

            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    for (size_t k = 0; k < dim; k++) {
                        ret(i, j) += multiplier(i, k) * multiplicand(k, j);
                    }
                }
            }

            *this = ret;
            return *this;
        }

        MpMatrix &operator*=(const fmp_t &multiplicand) {
            std::for_each(this->begin(), this->end(), [&](auto &fmp) {
                fmp *= multiplicand;
            });

            return *this;
        }

        friend std::ostream &operator<<(std::ostream &os, const MpMatrix &matrix);
    };

    inline MpMatrix operator+(MpMatrix lhs, const MpMatrix &rhs) {
        lhs += rhs;
        return lhs;
    }

    inline MpMatrix operator+(MpMatrix lhs, const fmp_t &rhs) {
        lhs += rhs;
        return lhs;
    }

    inline MpMatrix operator-(MpMatrix lhs, const MpMatrix &rhs) {
        lhs -= rhs;
        return lhs;
    }

    inline MpMatrix operator-(MpMatrix lhs, const fmp_t &rhs) {
        lhs -= rhs;
        return lhs;
    }

    inline MpMatrix operator*(MpMatrix lhs, const MpMatrix &rhs) {
        lhs *= rhs;
        return lhs;
    }

    inline MpMatrix operator*(MpMatrix lhs, const fmp_t &rhs) {
        lhs *= rhs;
        return lhs;
    }


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
