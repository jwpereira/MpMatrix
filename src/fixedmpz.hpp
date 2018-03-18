#pragma once

#include <gmpxx.h>
#include <iostream>

namespace mpmatrix {
    class fixedmpz {
      private:
        mpz_class number;
        mp_bitcnt_t scale;

      public:
        fixedmpz(mpz_class number, mp_bitcnt_t scale) : number(number), scale(scale) {}
        fixedmpz(mpz_class number) : fixedmpz(number, 0) {}

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

        decltype(auto) get_mpz_t() {
            return this->number.get_mpz_t();
        }

        fixedmpz operator+=(const fixedmpz &addend) {
            this->number += addend.number;
            return *this;
        }

        fixedmpz operator-=(const fixedmpz &subtrahend) {
            this->number -= subtrahend.number;
            return *this;
        }

        fixedmpz operator*=(const fixedmpz &multiplicand) {
            this->number >>= this->scale;
            this->number *= (multiplicand.number >> multiplicand.scale);

            return *this;
        }

        fixedmpz operator/=(const fixedmpz &divisor) {
            this->number /= divisor.number;
            return *this;
        }

        fixedmpz operator>>=(const mp_bitcnt_t &amount) {
            this->number >>= amount;
            return *this;
        }

        fixedmpz operator<<=(const mp_bitcnt_t &amount) {
            this->number <<= amount;
            return *this;
        }

        friend std::ostream &operator<<(std::ostream &os, const fixedmpz fmp);
        friend void square(fixedmpz &rop, fixedmpz &op);
        friend decltype(auto) sqrt(fixedmpz &rop);
    };

    inline fixedmpz operator+(fixedmpz lhs, const fixedmpz &rhs) {
        lhs += rhs;
        return lhs;
    }

    inline fixedmpz operator-(fixedmpz lhs, const fixedmpz &rhs) {
        lhs -= rhs;
        return lhs;
    }

    inline fixedmpz operator*(fixedmpz lhs, const fixedmpz &rhs) {
        lhs *= rhs;
        return lhs;
    }

    inline fixedmpz operator/(fixedmpz lhs, const fixedmpz &rhs) {
        lhs /= rhs;
        return lhs;
    }

    inline fixedmpz operator>>(fixedmpz lhs, const mp_bitcnt_t &rhs) {
        lhs >>= rhs;
        return lhs;
    }

    inline fixedmpz operator<<(fixedmpz lhs, const mp_bitcnt_t &rhs) {
        lhs <<= rhs;
        return lhs;
    }

    inline std::ostream &operator<<(std::ostream &os, const fixedmpz fmp) {
        return os << fmp.number << ">>" << fmp.scale;
    }

    inline void square(fixedmpz &rop, fixedmpz &op) {
        op >>= op.scale;
        mpz_pow_ui(rop.get_mpz_t(), op.get_mpz_t(), 2);
    }

    inline decltype(auto) sqrt(fixedmpz &rop) {
        rop >>= rop.scale;
        return sqrt(rop.number);
    }

}
