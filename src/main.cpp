#include <iostream>
#include <gmpxx.h>
#include "fixedmpz.hpp"
#include "mpm_algorithm.hpp"
#include "mpmatrix.hpp"
#include "prettymp.hpp"

using namespace momentmp;

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    size_t m_dim = 4;
    fmpz_scale m_scale = 256;
    fixedmpz m_raw[] = {
        18^256_fmpz, 22^256_fmpz,  54^256_fmpz,  42^256_fmpz,  
        22^256_fmpz, 70^256_fmpz,  86^256_fmpz,  62^256_fmpz,
        54^256_fmpz, 86^256_fmpz, 174^256_fmpz, 134^256_fmpz, 
        42^256_fmpz, 62^256_fmpz, 134^256_fmpz, 106^256_fmpz
    };

    MpMatrix m(m_dim, m_scale, m_raw, (m_raw + (m_dim * m_dim)));

    std::cout << "Before:\n";
    std::cout << m << std::endl;

    MpMatrix m_cholesky(m_dim, m_scale);

    if (cholesky(m, m_cholesky)) {
        std::cout << "After:\n";
        std::cout << m_cholesky << std::endl;
    } else {
        std::cout << "An error occured:\n";
    }    

    std::cout << "Debug:\n";
    std::cout << PrecPrint(m_cholesky, (64 >> 2));

    return 0;
}
