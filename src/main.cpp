#include <iostream>
#include <gmpxx.h>
#include "fixedmpz.hpp"
#include "mpmatrix.hpp"

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

    size_t m_dim = 4;
    fmp_shift_t m_shift = 256;
    fixedmpz m_raw[] = {
        18^256_fmpz, 22^256_fmpz,  54^256_fmpz,  42^256_fmpz,  
        22^256_fmpz, 70^256_fmpz,  86^256_fmpz,  62^256_fmpz,
        54^256_fmpz, 86^256_fmpz, 174^256_fmpz, 134^256_fmpz, 
        42^256_fmpz, 62^256_fmpz, 134^256_fmpz, 106^256_fmpz
    };

    MpMatrix m(2, 2, 256);
    m[0][0] = 100^256_fmpz;
    std::cout << m[0][0] << std::endl;

    return 0;
}
