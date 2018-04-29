#include <iostream>
#include <gmpxx.h>
#include "fixedmpz.hpp"
#include "mpmatrix.hpp"
#include "moment_algorithm.hpp"

using namespace momentmp;

bool sample_gen_ij(fmp_t &dest, size_t i, size_t j) {
    auto shift = fmpshift(dest.getShift());
    dest = ((i+1)^shift) / ((j+1)^shift);
    return true;
}

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    fmp_shift_t m_shift = 256;

    // size_t dim = 4;
    // MpMatrix m(dim, dim, m_shift);
    // fixedmpz m_raw[] = {
    //     18^256_fmpz, 22^256_fmpz,  54^256_fmpz,  42^256_fmpz,  
    //     22^256_fmpz, 70^256_fmpz,  86^256_fmpz,  62^256_fmpz,
    //     54^256_fmpz, 86^256_fmpz, 174^256_fmpz, 134^256_fmpz, 
    //     42^256_fmpz, 62^256_fmpz, 134^256_fmpz, 106^256_fmpz
    // };


    auto dim = 3;
    MpMatrix m(dim, dim, m_shift);

    fixedmpz m_raw[] = {
        4^256_fmpz,     0^256_fmpz,     0^256_fmpz, 
        12^256_fmpz,    37^256_fmpz,    0^256_fmpz,
        -(16^256_fmpz), -(43^256_fmpz), 98^256_fmpz 
    };

    for (size_t col = 0; col < dim; col++) {
        for (size_t row = 0; row < dim; row++) {
            m[col][row] = m_raw[row * dim + col];
        }
    }

    std::cout << "before:\n";
    std::cout << m;

    cholesky_decompose(m);
    std::cout << "after:\n";
    std::cout << m;

    transpose(m);
    std::cout << "after transpose:\n";
    std::cout << m;

    return 0;
}
