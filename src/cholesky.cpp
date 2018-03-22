#include <gmpxx.h>
#include "mpmatrix.hpp"
#include "mpmatrix_algorithm.hpp"

using namespace momentmp;
using momentmp::MpMatrix;
/**
 * This function takes in an initial matrix and outputs the lower triangular matrix that is the
 * result of a Cholesky decomposition on it. 
 */
bool momentmp::cholesky(const MpMatrix &initial, MpMatrix &lower) {
    if (initial.getDimension() != lower.getDimension()) {
        return false;
    }

    size_t dim = initial.getDimension();

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j <= i; j++) {
            fmp_t sum(0, initial.getShift());
            if (j == i) {
                for (size_t k = 0; k < j; k++) {
                    sum += sq(lower(j, k));
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