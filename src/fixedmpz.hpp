#pragma once

#include <gmpxx.h>
#include <iostream>
#include <string>

namespace momentmp {
    // forward declaration for the class below; allows the aliases immediately following to work.
    class fixedmpz;

    using fmp_t = fixedmpz; ///< convenience alias to the underlying number type used for operations
    using fmp_scale_t = mp_bitcnt_t;   ///< alias to the scaling type really used for fixedmpz

    /// mpz_class based class purposed for fixed-precision arithmetic.
    /** 
     * Use this to have mpz's that automatically get shifted into place.
     * Has overloaded operators for convenience.
     */
    class fixedmpz {
      private:
        mpz_class number;
        fmp_scale_t scale;

      public:
        fixedmpz(mpz_class number, fmp_scale_t scale) : number(number), scale(scale) {}
        fixedmpz(mpz_class number) : fixedmpz(number, 0) {}

        fmp_scale_t getScale() const {
            return this->scale;
        }

        decltype(auto) get_mpz_t() {
            return this->number.get_mpz_t();
        }

        decltype(auto) get_mpz_t() const {
            return this->number.get_mpz_t();
        }

        void setScale(fmp_scale_t scale) {
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

        fixedmpz operator>>=(const fmp_scale_t &amount) {
            this->number >>= amount;
            return *this;
        }

        fixedmpz operator<<=(const fmp_scale_t &amount) {
            this->number <<= amount;
            return *this;
        }

        friend std::ostream &operator<<(std::ostream &os, const fixedmpz &fmp);
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

    inline fixedmpz operator>>(fixedmpz lhs, const fmp_scale_t &rhs) {
        lhs >>= rhs;
        return lhs;
    }

    inline fixedmpz operator<<(fixedmpz lhs, const fmp_scale_t &rhs) {
        lhs <<= rhs;
        return lhs;
    }

    inline std::ostream &operator<<(std::ostream &os, const fixedmpz &fmp) {
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

    class fmpz_adapter {
      public:
        const fmp_scale_t scale;
        constexpr fmpz_adapter(unsigned long long scale) : scale(scale) {}
    };

    constexpr inline fmpz_adapter operator""_fmpz(unsigned long long literal) {
        return fmpz_adapter(literal);
    }

    inline fixedmpz operator^(unsigned long literal, fmpz_adapter scale) {
        mpz_class number(literal);
        number <<= scale.scale;
        return fixedmpz(number, scale.scale);
    }
}
