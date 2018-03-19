#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include <vector>
#include "fixedmpz.hpp"

namespace mpmatrix {
    using mp_t = fixedmpz;

    class MpMatrix {
      private:
        std::vector<mp_t> matrix;
        size_t dim;       // dimension * dimension = size [we're working with square matricies]
        mp_bitcnt_t scale;
      public:
        MpMatrix(size_t dim, mp_bitcnt_t scale) : dim(dim), scale(scale) {
            this->matrix = std::vector<mp_t>(dim * dim, mp_t(0, scale));
        }

        template <typename Iterator>
        MpMatrix(size_t dim, mp_bitcnt_t scale, Iterator begin, Iterator end) : MpMatrix(dim, scale) {
            //std::copy(begin, end, this->matrix.begin());
            auto iter = this->matrix.begin();
            std::for_each(begin, end, [&](auto &val) {
                iter->setNumber(val);
                iter->setScale(val.getScale());
                iter++;
            });
        }

        mp_t &operator()(size_t row, size_t col) {
            return this->matrix[(row * this->dim) + col];
        }

        mp_t &operator()(size_t row, size_t col) const {
            return const_cast<mp_t&>(this->matrix[(row * this->dim) + col]);
        }

        mp_bitcnt_t getScale() const {
            return this->scale;
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
        for (auto &mp_t : matrix) {
            counter++;
            os << mp_t << '\t';
            if (counter % matrix.getDimension() == 0) {
                os << '\n';
            }
        }

        //Restore the initial flags of the output streams
        os.flags(initialFlags);

        return os << std::endl;
    }

}
