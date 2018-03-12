#pragma once

#include <gmpxx.h>
#include <iostream>

namespace mpmatrix {
    class fixedmpz {
      private:
        mpz_class number;
        mp_bitcnt_t scale;

      public:
        fixedmpz(mpz_class number, mp_bitcnt_t scale) : number(number << scale), scale(scale) {}

        operator mpz_class() {
            return this->number;
        }

        friend std::ostream &operator<<(std::ostream &os, const fixedmpz fmp);
    };

    std::ostream &operator<<(std::ostream &os, const fixedmpz fmp) {
        return os << fmp.number;
    }
}