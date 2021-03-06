/**
 * @brief Wrapper class for GMP's mpz_class geared towards fixed-point arithmetic
 *
 * @file fixedmpz.hpp
 * @author jwpereira
 */

#pragma once

#include <gmpxx.h>
#include <iostream>
#include <string>

namespace momentmp {
    // forward declaration for the class below; allows the aliases immediately following to work.
    class fixedmpz;

    using fmp_t = fixedmpz; ///< convenience alias to the underlying number type used for operations

    /** Alias for fixedmpz's underlying type representing the where the "dot" goes (in base 2) */
    using fmpz_shift_t = mp_bitcnt_t;

    /**
     * @brief mpz_class based class purposed for fixed-precision arithmetic.
     *
     * Use this to have mpz's that automatically get shifted into place.
     * Has overloaded operators (such as +, -, *, /, <<, >>) for convenience. This class essentially wraps GMP's mpz_class,
     * allowing it to be used for fixed point arithmetic.
     */
    class fixedmpz {
      private:
        mpz_class number;
        fmpz_shift_t shift;

      public:
        fixedmpz(mpz_class number, fmpz_shift_t shift) noexcept : number(number), shift(shift) {}
        fixedmpz(mpz_class number) noexcept : fixedmpz(number, 0) {}
        fixedmpz(const fixedmpz &other) = default;
        fixedmpz(fixedmpz &&other) = default;

        /// Returns how much the the underlying mpz_class is shifted by
        fmpz_shift_t getShift() const {
            return this->shift;
        }

        /**
         * @brief Internal mpz_class's get_mpz_t()
         */
        decltype(auto) get_mpz_t() {
            return this->number.get_mpz_t();
        }

        /**
         * @brief Internal mpz_class's get_mpz_t()
         */
        decltype(auto) get_mpz_t() const {
            return this->number.get_mpz_t();
        }

        void setShift(fmpz_shift_t shift) {
            this->shift = shift;
        }

        void setNumber(mpz_class number) {
            this->number = number;
        }

        /**
         * @brief Returns floating point representation as mpf_class
         */
        mpf_class to_mpf() const {
            mpf_class shiftd(0, this->shift);
            mpf_set_z(shiftd.get_mpf_t(), this->get_mpz_t());
            shiftd >>= this->shift;

            return shiftd;
        }

        /**
         * @brief Return underlying mpz_class
         */
        mpz_class &operator()() {
            return this->number;
        }

        /**
         * @brief Return underlying mpz_class
         */
        mpz_class &operator()() const {
            return const_cast<mpz_class&>(this->number);
        }

        /**
         * @brief Allows this object to be cast to mpz_class
         *
         * Returns internal mpz_class
         */
        operator mpz_class() {
            return this->number;
        }

        fixedmpz &operator=(const fixedmpz &other) = default;
        fixedmpz &operator=(fixedmpz &&other) = default;

        fixedmpz &operator+=(const fixedmpz &addend) {
            this->number += addend.number;
            return *this;
        }

        fixedmpz &operator-=(const fixedmpz &subtrahend) {
            this->number -= subtrahend.number;
            return *this;
        }

        fixedmpz &operator*=(const fixedmpz &multiplier) {
            this->number *= multiplier.number;
            this->number >>= this->shift;
            return *this;
        }

        fixedmpz &operator/=(const fixedmpz &divisor) {
            this->number <<= this->shift;
            this->number /= divisor.number;
            return *this;
        }

        fixedmpz &operator>>=(const fmpz_shift_t &amount) {
            this->number >>= amount;
            return *this;
        }

        fixedmpz &operator<<=(const fmpz_shift_t &amount) {
            this->number <<= amount;
            return *this;
        }

        fixedmpz operator-() const {
            fixedmpz ret(*this);
            mpz_neg(ret.get_mpz_t(), this->get_mpz_t());
            return ret;
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

    inline fixedmpz operator>>(fixedmpz lhs, const fmpz_shift_t &rhs) {
        lhs >>= rhs;
        return lhs;
    }

    inline fixedmpz operator<<(fixedmpz lhs, const fmpz_shift_t &rhs) {
        lhs <<= rhs;
        return lhs;
    }

    inline std::ostream &operator<<(std::ostream &os, const fixedmpz &fmp) {
        return os << fmp.to_mpf();
    }

    /**
     * @brief Returns the square of a fixedmpz number in fixedmpz format
     *
     * Performs the squaring by multiplying input * input.
     *
     * @param[in] op fixedmpz number to be squared.
     */
    inline fixedmpz sq(const fixedmpz &op) {
        return op * op;
    }

    /**
     * @brief Returns the square root of a fixedmpz number in fixedmpz format
     *
     * @param[in] op fixedmpz number to take the square root of.
     */
    inline fixedmpz sqrt(const fixedmpz &op) {
        fixedmpz ret = op << op.getShift();
        ret() = sqrt(ret());
        return ret;
    }

    /**
     * @brief Returns the factorial result of the number (unsigned int) passed in
     *
     * Simply forwards the calculation to mpz's mpz_fac_ui function.
     *
     * @param base number to factorialize
     * @return fixedmpz factorial of base
     */
    inline fixedmpz factorial(unsigned long base) {
        fixedmpz ret(0);
        mpz_fac_ui(ret.get_mpz_t(), base);
        return ret;
    }

    /**
     * @brief Utility class designed purely to make the _fmpz literal operator happen
     *
     * This class really should not be used for anything as it is just made to allow
     * <code>_fmpz</code> to work. What it is actually doing is serving as a wrapping type for \link
     * fmpz_shift_t \endlink. The <code>_fmpz</code> suffix turns an unsigned long
     */
    class fmpz_adapter {
      public:
        const fmpz_shift_t shift;
        constexpr fmpz_adapter(unsigned long long shift) : shift(shift) {}
    };

    constexpr inline fmpz_adapter operator""_fmpz(unsigned long long literal) {
        return fmpz_adapter(literal);
    }

    /**
     * @brief Use _fmpz as an operator for literal values to be instantly made into fixedmpz values.
     *
     * Example usage: 1^256_fmpz will essentially make a fixedmpz of value 1, shifted 256 places.
     * Because it uses fixedmpz, it will still be legible as 1.
     */
    inline fixedmpz operator^(unsigned long literal, fmpz_adapter shift) {
        mpz_class number(literal);
        number <<= shift.shift;
        return fixedmpz(number, shift.shift);
    }

    inline fmpz_adapter fmpzshift(unsigned long shift) {
        return fmpz_adapter(shift);
    }
}
