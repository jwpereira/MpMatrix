#include <gmpxx.h>
#include "mpmatrix.hpp"
#include "mpm_algorithm.hpp"

using namespace mpmatrix;

bool mpmatrix::cholesky(const MpMatrix &initial, MpMatrix &lower) {
    if (initial.getDimension() != lower.getDimension()) {
        return false;
    }

    size_t dim = initial.getDimension();

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j <= i; j++) {
            mpz_class sum(0_mpz);
            if (j == i) {
                for (size_t k = 0; k < j; k++) {
                    mpz_class addend(0_mpz);
                    mpz_pow_ui(addend.get_mpz_t(), lower(j, k).get_mpz_t(), 2);
                    sum += addend;
                }
                lower(j, j) = sqrt(initial(j, j) - sum);
            } else {
                for (size_t k = 0; k < j; k++) {
                    sum += (lower(i, k) * lower(j, k));
                }
                lower(i, j) = (initial(i, j) - sum) / lower(j, j);
            }
        }
    }

    return true;
}