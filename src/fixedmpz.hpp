#pragma once

#include <gmpxx.h>
#include <iostream>

namespace mpmatrix {
    class fixedmpz {
      private:
        mpz_class number;
        mp_bitcnt_t dot;
        
      public:
        fixedmpz(mpz_class number, mp_bitcnt_t dot) : number(number), dot(dot) {}

        operator mpz_class() {
            return this->number;
        }

        friend std::ostream &operator<<(std::ostream &os, const fixedmpz fmp);
    };

    std::ostream &operator<<(std::ostream &os, const fixedmpz fmp) {
        return os << fmp.number;
    }
}