#include <iostream>
#include <gmpxx.h>
#include "fixedmpz.hpp"
#include "mpm_algorithm.hpp"
#include "mpmatrix.hpp"
#include "prettymp.hpp"

using namespace mpmatrix;

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    std::cout << "THE SQUARE ROOT OF 174 IS" << sqrt(fixedmpz(174)) << std::endl;

    fixedmpz a(4096_mpz, 3);
    fixedmpz b(64_mpz, 3);

    fixedmpz c = a * b;

    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "c: " << c << std::endl;

    fixedmpz base(256, 2);
    fixedmpz square(0, 2);
    fixedmpz square_root(0, 2);

    std::cout << "base: " << base << std::endl;
    std::cout << "square: " << sq(base) << std::endl;
    std::cout << "sqrt: " << mpmatrix::sqrt(base) << std::endl;

    size_t m_dim = 4;

    fixedmpz testo(18_mpz << 32, 32);

    std::cout << testo << std::endl;

   // exit(0);

    mp_bitcnt_t m_scale = 64;
    fixedmpz m_raw[] = {
        fixedmpz(18_mpz << m_scale, m_scale), fixedmpz(22_mpz << m_scale, m_scale), fixedmpz( 54_mpz << m_scale, m_scale), fixedmpz( 42_mpz << m_scale, m_scale),  
        fixedmpz(22_mpz << m_scale, m_scale), fixedmpz(70_mpz << m_scale, m_scale), fixedmpz( 86_mpz << m_scale, m_scale), fixedmpz( 62_mpz << m_scale, m_scale),
        fixedmpz(54_mpz << m_scale, m_scale), fixedmpz(86_mpz << m_scale, m_scale), fixedmpz(174_mpz << m_scale, m_scale), fixedmpz(134_mpz << m_scale, m_scale), 
        fixedmpz(42_mpz << m_scale, m_scale), fixedmpz(62_mpz << m_scale, m_scale), fixedmpz(134_mpz << m_scale, m_scale), fixedmpz(106_mpz << m_scale, m_scale)
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
    std::cout << DebugPrint(m_cholesky);

    return 0;
}
