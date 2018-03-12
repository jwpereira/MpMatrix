#include <iostream>
#include "fixedmpz.hpp"
#include <gmpxx.h>
#include "mpmatrix.hpp"
#include "mpm_algorithm.hpp"
#include "prettymp.hpp"

using namespace mpmatrix;

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    fixedmpz numero(1_mpz, 256);

    std::cout << numero << std::endl;

    /*
    size_t m_dim = 4;
    mp_bitcnt_t m_prec = 8;
    mpz_class m_raw[] = {
        mpz_class(18_mpz), mpz_class(22_mpz), mpz_class(54_mpz) , mpz_class(42_mpz),  
        mpz_class(22_mpz), mpz_class(70_mpz), mpz_class(86_mpz) , mpz_class(62_mpz),
        mpz_class(54_mpz), mpz_class(86_mpz), mpz_class(174_mpz), mpz_class(134_mpz), 
        mpz_class(42_mpz), mpz_class(62_mpz), mpz_class(134_mpz), mpz_class(106_mpz)
    };

    MpMatrix m(m_dim, m_raw, (m_raw + (m_dim * m_dim)));

    std::cout << "Before:\n";
    std::cout << m << std::endl;

    mp_bitcnt_t m_cholesky_prec = 256;
    MpMatrix m_cholesky(m_dim);

    if (cholesky(m, m_cholesky)) {
        std::cout << "After:\n";
        std::cout << m_cholesky << std::endl;
    } else {
        std::cout << "An error occured:\n";
    }    

    std::cout << "Debug:\n";
    std::cout << DebugPrint(m_cholesky);
    */

    return 0;
}
