#include <iostream>
#include <gmpxx.h>
#include "fixedmpz.hpp"
#include "mpmatrix_algorithm.hpp"
#include "mpmatrix.hpp"
#include "prettymp.hpp"

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

    MpMatrix m(m_dim, m_shift, m_raw, (m_raw + (m_dim * m_dim)));

    std::cout << "Before:\n";
    std::cout << m << std::endl;

    MpMatrix m_cholesky(m_dim, m_shift);

    if (cholesky(m, m_cholesky)) {
        std::cout << "After:\n";
        std::cout << m_cholesky << std::endl;
    } else {
        std::cout << "An error occured:\n";
    }    

    std::cout << "Debug:\n";
    std::cout << PrecPrint(m_cholesky, (64 >> 2));

    MpMatrix m_i_over_j(m_dim, m_shift);
    apply(m_i_over_j, sample_gen_ij);

    std::cout << PrecPrint(m_i_over_j, 12);   

    MpMatrix m_i_plus_j(m_dim, m_shift);
    apply(m_i_plus_j, [](auto &fmp, auto i, auto j) {
        auto shift = fmpshift(fmp.getShift());
        fmp = (i^shift) + (j^shift);
        return true;
    }); 

    std::cout << PrecPrint(m_i_plus_j, 12); 

    return 0;
}
