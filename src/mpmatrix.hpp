#pragma once

#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include <vector>

namespace mpmatrix {
    class MpMatrix {
      private:
        std::vector<mpf_class> matrix;
        size_t dim;       // dimension * dimension = size [we're working with square matricies]
        mp_bitcnt_t prec; // precision of each number

      public:
        MpMatrix(size_t dim, mp_bitcnt_t prec) : dim(dim), prec(prec) {
            this->matrix = std::vector<mpf_class>(dim * dim, mpf_class(0, prec));
        }

        template <typename Iterator>
        MpMatrix(size_t dim, mp_bitcnt_t prec, Iterator begin, Iterator end) : MpMatrix(dim, prec) {
            std::copy(begin, end, this->matrix.begin());
        }

        mpf_class &operator()(size_t row, size_t col) {
            return this->matrix[(row * this->dim) + col];
        }

        mpf_class &operator()(size_t row, size_t col) const {
            return const_cast<mpf_class&>(this->matrix[(row * this->dim) + col]);
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

        size_t getDim() const {
            return this->dim;
        }

        friend std::ostream &operator<<(std::ostream &os, const MpMatrix &matrix);
    };

    inline std::ostream &operator<<(std::ostream &os, const MpMatrix &matrix) {
        size_t counter = 0;
        for (auto &mp : matrix) {
            counter++;
            os << mp << '\t';
            if (counter % matrix.getDim() == 0) {
                os << '\n';
            }
        }

        return os << std::endl;
    }
}
