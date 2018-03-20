#pragma once

#include <gmpxx.h>
#include <iostream>
#include <string>

namespace mpmatrix {
    using fmpz_scale = mp_bitcnt_t;

    class fixedmpz {
      private:
        mpz_class number;
        fmpz_scale scale;

      public:
        fixedmpz(mpz_class number, fmpz_scale scale) : number(number), scale(scale) {}
        fixedmpz(mpz_class number) : fixedmpz(number, 0) {}

        fmpz_scale getScale() const {
            return this->scale;
        }

        decltype(auto) get_mpz_t() {
            return this->number.get_mpz_t();
        }

        decltype(auto) get_mpz_t() const {
            return this->number.get_mpz_t();
        }

        void setScale(fmpz_scale scale) {
            this->scale = scale;
        }

        void setNumber(mpz_class number) {
            this->number = number;
        }

        mpf_class to_mpf() const {
            mpf_class scaled(0, this->scale);
            mpf_set_z(scaled.get_mpf_t(), this->get_mpz_t());
            scaled >>= this->scale;

            return scaled;
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

        fixedmpz operator>>=(const fmpz_scale &amount) {
            this->number >>= amount;
            return *this;
        }

        fixedmpz operator<<=(const fmpz_scale &amount) {
            this->number <<= amount;
            return *this;
        }

        friend std::ostream &operator<<(std::ostream &os, const fixedmpz fmp);
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

    inline fixedmpz operator>>(fixedmpz lhs, const fmpz_scale &rhs) {
        lhs >>= rhs;
        return lhs;
    }

    inline fixedmpz operator<<(fixedmpz lhs, const fmpz_scale &rhs) {
        lhs <<= rhs;
        return lhs;
    }

    inline std::ostream &operator<<(std::ostream &os, const fixedmpz fmp) {
        return os << fmp.to_mpf();
    }

    inline fixedmpz sq(const fixedmpz &op) {
        return op * op;
    }

    inline fixedmpz sqrt(const fixedmpz &rop) {
        fixedmpz ret = rop << rop.getScale();
        ret() = sqrt(ret());
        return ret;
    }

    class fmpz_scale_t {
      public:
        const fmpz_scale scale;
        constexpr fmpz_scale_t(unsigned long long scale) : scale(scale) {}
    };

    constexpr inline fmpz_scale_t operator""_fmpz(unsigned long long literal) {
        return fmpz_scale_t(literal);
    }

    inline fixedmpz operator^(unsigned long literal, fmpz_scale_t scale) {
        mpz_class number(literal);
        number <<= scale.scale;
        return fixedmpz(number, scale.scale);
    }
}
