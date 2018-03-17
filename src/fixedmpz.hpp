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

        mp_bitcnt_t getScale() const {
            return this->scale;
        }

        void setScale(mp_bitcnt_t scale) {
            this->scale = scale;
        }

        mpz_class &operator()() {
            return this->number;
        }

        operator mpz_class() {
            return this->number;
        }

        friend std::ostream &operator<<(std::ostream &os, const fixedmpz fmp);
    };

    inline std::ostream &operator<<(std::ostream &os, const fixedmpz fmp) {
        return os << fmp.number;
    }
}