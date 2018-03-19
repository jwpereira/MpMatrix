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

        decltype(auto) get_mpz_t() {
            return this->number.get_mpz_t();
        }

        decltype(auto) get_mpz_t() const {
            return this->number.get_mpz_t();
        }

        void setScale(mp_bitcnt_t scale) {
            this->scale = scale;
        }

        void setNumber(mpz_class number) {
            this->number = number;
        }

        mpz_class &operator()() {
            return this->number;
        }

        mpz_class &operator()() const {
            return const_cast<mpz_class&>(this->number);
        }

        operator mpz_class() {
            return this->number;
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
            this->number >>= (this->scale / 2);
            this->number *= (multiplicand.number >> multiplicand.scale / 2);

            return *this;
        }

        fixedmpz operator/=(const fixedmpz &divisor) {
            this->number <<= this->scale;
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
        // friend fixedmpz sqrt(const fixedmpz &rop);
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
        mpf_class scaled(0, fmp.scale);
        mpf_set_z(scaled.get_mpf_t(), fmp.number.get_mpz_t());
        scaled >>= fmp.scale;
        
        // return os << fmp.number << '(' << std::hex << fmp.number << std::dec << ')'
        //     << ">>" << fmp.scale << '=' << scaled
        //     << '(' << std::hex << scaled << ')' << std::dec;
        return os << scaled;
    }

    inline fixedmpz sq(const fixedmpz &op) {
        return op * op;
    }

    inline fixedmpz sqrt(const fixedmpz &rop) {
        fixedmpz ret = rop << rop.getScale();
        ret() = sqrt(ret());
        return ret;
    }
}
