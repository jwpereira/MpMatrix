/**
 * @brief For the high-precision matrix implementation
 *
 * MpMatrix uses the fixedmpz wrapper for GMP's mpz_class (itself a c++ wrapper for GMP's own mpz_t)
 * to store multiple-precision numbers in a matrix container format.
 *
 * @file mpmatrix.hpp
 * @author jwpereira
 */

#pragma once

#include <iostream>
#include <vector>

#include <gmpxx.h>
#include <omp.h>

#include "fixedmpz.hpp"

/**
 * @brief Namespace for the Multiple Precision Matrix Project
 *
 * Classes and functions in this namespace are particularly geared towards the needs of the
 * HankelHacker project.
 */
namespace momentmp {
    using fmp_t = fixedmpz; ///< convenience alias to the underlying number type used for operations
    using fmpz_shift_t = fmpz_shift_t;   ///< alias to the scaling type really used for fixedmpz

    /**
     * @brief MpArray is a wrapper for a vector of fixedmpz numbers
     */
    class MpArray {
      private:
        std::vector<fmp_t> col;
        size_t dim, id;
        fmpz_shift_t shift;
      public:
        MpArray(size_t dim, fmpz_shift_t shift, size_t id=-1) noexcept
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

        /**
         * @brief Return the size (or dim) of the MpArray
         */
        size_t size() {
            return this->dim;
        }

        /**
         * @brief Return the size (or dim) of the MpArray
         */
        size_t size() const {
            return this->dim;
        }

        /**
         * @brief Return the amount each element of the MpArray -should- be shifted by
         *
         * Note: MpArray does not enforce by how much any value in it is shifted.
         */
        fmpz_shift_t getShift() {
            return this->shift;
        }

        /**
         * @brief Return the amount each element of the MpArray -should- be shifted by
         *
         * Note: MpArray does not enforce by how much any value in it is shifted.
         */
        fmpz_shift_t getShift() const {
            return this->shift;
        }

        /**
         * @brief Manually set the id of the MpArray
         *
         * In a matrix, this would help each row/col (but not both row/col at the same time)
         * identify itself to another piece of code.
         */
        void setId(size_t id) {
            this->id = id;
        }

        /**
         * @brief Exposes the underlying vector's [] operator to access elements inside the MpArray
         */
        fmp_t &operator[](size_t row) {
            return this->col[row];
        }

        /**
         * @brief Exposes the underlying vector's [] operator to access elements inside the MpArray
         */
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

    /**
     * @brief Operator overload for having an MpArray be printed out using std::ostream
     *
     * Example usage : <code>std::cout << array;<\code>
     */
    inline std::ostream &operator<<(std::ostream &os, MpArray &array) {
        for (auto &e : array) {
            os << e << ' ';
        }

        return os << std::endl;
    }

    /**
     * @brief Signifies the two orientations an MpMatrix can take on
     *
     * As different matrix algorithms require different orientations of the data, this just serves
     * as a marker for which one an MpMatrix takes on. Especially hqelps with the printing of the
     * matrix.
     */
    enum MpMatrixMode : bool { ROW_ORIENTED=true , COL_ORIENTED=false };

    /**
     * @brief Matrix class based class focused on containing and accessing fmp_t elements.
     *
     * MpMatrix is a square-matrix container-wrapper for fmp_t elements. It allows for individual
     * elements in the underlying container to be accessed using a reasonable syntax (i.e.,
     * <code>name_of_matrix[row][col]</code> (if row-oriented) or
     * <code>name_of_matrix[col][row]</code>(if col-oriented)).
     */
    class MpMatrix {
      private:
        std::vector<MpArray> matrix;
        size_t dim;         ///< dimension * dimension = rows [we're working with square matricies]
        fmpz_shift_t shift;  ///< for keeping track of the shift/precision factor across the matrix
        MpMatrixMode mode;
      public:
        MpMatrix(size_t dim, fmpz_shift_t shift, MpMatrixMode mode = COL_ORIENTED) noexcept
                : dim(dim),  shift(shift) {
            this->matrix = std::vector<MpArray>(dim, MpArray(dim, shift, 0));
            for (size_t i = 0; i < dim; i++) {
                this->matrix[i].setId(i);
            }
            this->mode = mode;
        }

        MpMatrix(const MpMatrix &other) = default;
        MpMatrix(MpMatrix &&other) = default;

        /**
         * @brief Return the size of the MpArray.
         *
         * For an NxN matrix, this will return N.
         */
        size_t getDim() {
            return this->dim;
        }

        /**
         * @brief Return the size of the MpMatrix.
         *
         * For an NxN matrix, this will return N.
         */
        size_t getDim() const {
            return this->dim;
        }

        /**
         * @brief Return the amount each element of the MpMatrix -should- be shifted by
         *
         * Note: MpMatrix does not enforce by how much any value in it is shifted.
         */
        fmpz_shift_t getShift() {
            return this->shift;
        }

        /**
         * @brief Return the amount each element of the MpMatrix -should- be shifted by
         *
         * Note: MpMatrix does not enforce by how much any value in it is shifted.
         */
        fmpz_shift_t getShift() const {
            return this->shift;
        }

        /**
         * @brief Return the orientation mode currently set for the MpMatrix
         */
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

        /**
         * @brief Manually set the orientation of the MpMatrix
         *
         * Note: this does not do any transposing, use reorient() for that.
         */
        void setMode(MpMatrixMode mode) {
            this->mode = mode;
        }

        /**
         * @brief Exposes the underlying vector's [] operator to access elements inside the MpMatrix
         *
         * The return of this will be a reference to the desired MpArray inside the MpMatrix
         */
        MpArray &operator[](const size_t col) {
            return this->matrix[col];
        }

        /**
         * @brief Exposes the underlying vector's [] operator to access elements inside the MpMatrix
         *
         * The return of this will be a reference to the desired MpArray inside the MpMatrix
         */
        MpArray &operator[](const size_t col) const {
            return const_cast<MpArray&>(this->matrix[col]);
        }

        void clear() {
            auto zero = 0^fmpzshift(this->getShift());
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

    /**
     * @brief Operator overload for having an MpMatrix be printed out using std::ostream
     *
     * Example usage : <code>std::cout << matrix;<\code>
     */
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

    /**
     * @brief This function multiplies together two MpMatrix objects.
     *
     * Nothing particularly fancy, just schoolhouse multiplication.
     */
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

    /**
     * @brief Performs an in-place transpose of an MpMatrix
     *
     * per https://en.wikipedia.org/wiki/In-place_matrix_transposition#Square_matrices
     */
    inline void transpose(MpMatrix &matrix) {
        auto dim = matrix.getDim();

        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t n = 0; n < (dim - 1); n++) {
            for (size_t m = (n + 1); m < dim; m++) {
                auto temp = matrix[m][n];
                matrix[m][n] = matrix[n][m];
                matrix[n][m] = temp;
            }
        }
    }

    /**
     * @brief Transposes an MpMatrix and switches its orientation mode
     *
     * per https://en.wikipedia.org/wiki/In-place_matrix_transposition#Square_matrices
     */
    inline void reorient(MpMatrix &matrix) {
        transpose(matrix);

        auto mode = matrix.getMode();
        if (mode == COL_ORIENTED) {
            matrix.setMode(ROW_ORIENTED);
        } else if (mode == ROW_ORIENTED) {
            matrix.setMode(COL_ORIENTED);
        }
    }

    /**
     * @brief Given a lower triangular matrix, this reflects it across the diagonal to complete it
     *
     * Perhaps useful for symmetrical matrices
     */
    inline void reflect(MpMatrix &matrix) {
        auto dim = matrix.getDim();

        for (size_t n = 0; n < (dim - 1); n++) {
            for (size_t m = (n + 1); m < dim; m++) {
                auto temp = matrix[m][n];
                matrix[m][n] = matrix[n][m];
            }
        }
    }
}
