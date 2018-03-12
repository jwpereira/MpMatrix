#pragma once

#include <gmpxx.h>

namespace mpmatrix {
    class fixedmpz {
      private:
        mpz_class number;
        mp_bitcnt_t dot;
      public:
        fixedmpz(mpz_class number, mp_bitcnt_t dot) : number(number), dot(dot) {}

        operator mpz_class*() {
          return &(this->number);
        }
    };
}