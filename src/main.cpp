#include <iostream>
#include <gmpxx.h>
#include "mpmatrix.hpp"

using mpmatrix::MpMatrix;

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    size_t m_dim = 4;
    mp_bitcnt_t m_prec = 8;
    mpf_class m_raw[] = {
        mpf_class(18_mpf, m_prec), mpf_class(22_mpf, m_prec), mpf_class(54_mpf, m_prec), mpf_class(42_mpf, m_prec),  
        mpf_class(22_mpf, m_prec), mpf_class(70_mpf, m_prec), mpf_class(86_mpf, m_prec), mpf_class(62_mpf, m_prec),
        mpf_class(54_mpf, m_prec), mpf_class(86_mpf, m_prec), mpf_class(174_mpf,m_prec), mpf_class(134_mpf,m_prec), 
        mpf_class(42_mpf, m_prec), mpf_class(62_mpf, m_prec), mpf_class(134_mpf,m_prec), mpf_class(106_mpf,m_prec)
    };

    MpMatrix m(m_dim, 8, m_raw, (m_raw + (m_dim * m_dim)));

    std::cout << "Before:\n";
    std::cout << m << std::endl;

<<<<<<< HEAD
    mp_bitcnt_t m_cholesky_prec = 256;
    MpMatrix m_cholesky(m_dim, m_cholesky_prec);

    if (mpmatrix::cholesky(m, m_cholesky)) {
        std::cout << "After:\n";
        std::cout << m_cholesky << std::endl;
    } else {
        std::cout << "An error occured:\n";
    }    
=======
    mp_bitcnt_t m_cholesky_prec = mpmatrix::bitsFromDigits(1000);
    MpMatrix m_cholesky(m_dim, m_cholesky_prec);
    mpmatrix::cholesky(m, m_cholesky);

    std::cout << "After:\n";
    std::cout << m_cholesky << std::endl;
>>>>>>> 601bc4f... updated main demo

    return 0;
}
