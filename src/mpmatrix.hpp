#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include <vector>

namespace mpmatrix {
    class MpMatrix {
      private:
        std::vector<mpz_class> matrix;
        size_t dim;       // dimension * dimension = size [we're working with square matricies]

      public:
        MpMatrix(size_t dim) : dim(dim) {
            this->matrix = std::vector<mpz_class>(dim * dim, mpz_class(0));
        }

        template <typename Iterator>
        MpMatrix(size_t dim, Iterator begin, Iterator end) : MpMatrix(dim) {
            std::copy(begin, end, this->matrix.begin());
        }

        mpz_class &operator()(size_t row, size_t col) {
            return this->matrix[(row * this->dim) + col];
        }

        mpz_class &operator()(size_t row, size_t col) const {
            return const_cast<mpz_class&>(this->matrix[(row * this->dim) + col]);
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

        // std::cout.precision(matrix.getPrecision());

        size_t counter = 0;
        for (auto &mp : matrix) {
            counter++;
            os << mp << '\t';
            if (counter % matrix.getDimension() == 0) {
                os << '\n';
            }
        }

        //Restore the initial flags of the output streams
        os.flags(initialFlags);

        return os << std::endl;
    }

}
