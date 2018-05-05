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

void demo(MpMatrix &m, size_t dim, fmp_shift_t m_shift) {
    std::cout << "original:\n";
    std::cout << m;

    cholesky_decompose(m);
    std::cout << "cholesky:\n";
    std::cout << m;

    MpMatrix lt(m);
    transpose(lt);
    std::cout << "after transpose:\n";
    std::cout << lt;

    reorient(m);
    std::cout << "reoriented:\n";
    std::cout << m;

    MpArray diagonals(dim, m_shift);
    extract_diagonal(m, diagonals);
    std::cout << "diagonals:\n";
    std::cout << diagonals;
    std::cout << "lower:\n";
    std::cout << m;

    // reorient(m);
    invert(m);
    std::cout << "inverted:\n";
    std::cout << m;

    invert_diagonal(diagonals);
    std::cout << "diagonals inverted:\n";
    std::cout << diagonals;

    invert(m);
    std::cout << "inverted inverted:\n";
    std::cout << m;

    invert_diagonal(diagonals);
    std::cout << "diagonals inverted inverted:\n";
    std::cout << diagonals;
}

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    fmp_shift_t m_shift = 256;

    // size_t dim = 4;
    // MpMatrix m(dim, m_shift);
    // fixedmpz m_raw[] = {
    //     18^256_fmpz, 22^256_fmpz,  54^256_fmpz,  42^256_fmpz,  
    //     22^256_fmpz, 70^256_fmpz,  86^256_fmpz,  62^256_fmpz,
    //     54^256_fmpz, 86^256_fmpz, 174^256_fmpz, 134^256_fmpz, 
    //     42^256_fmpz, 62^256_fmpz, 134^256_fmpz, 106^256_fmpz
    // };

    auto dim = 3;
    MpMatrix m(dim, m_shift);
    fixedmpz m_raw[] = {
        4^256_fmpz,     0^256_fmpz,     0^256_fmpz, 
        12^256_fmpz,    37^256_fmpz,    0^256_fmpz,
        -(16^256_fmpz), -(43^256_fmpz), 98^256_fmpz 
    };

    // size_t dim = 4;
    // MpMatrix m(dim, m_shift);
    // fixedmpz m_raw[] = {
    //     1^256_fmpz, 0^256_fmpz, 0^256_fmpz, 0^256_fmpz,  
    //     4^256_fmpz, 1^256_fmpz, 0^256_fmpz, 0^256_fmpz,
    //     6^256_fmpz, 2^256_fmpz, 1^256_fmpz, 0^256_fmpz, 
    //     3^256_fmpz, 7^256_fmpz, 4^256_fmpz, 1^256_fmpz
    // };

    for (size_t col = 0; col < dim; col++) {
        for (size_t row = 0; row < dim; row++) {
            m[col][row] = m_raw[row * dim + col];
        }
    }

    // demo(m, dim, m_shift);

    MpMatrix l(m);
    cholesky_decompose(l);
    std::cout << "cholesky ld:\n" << l << std::endl;

    MpArray diagonal(dim, m_shift);
    extract_diagonal(l, diagonal);

    reorient(l);

    MpMatrix lt(l);
    transpose(lt);

    std::cout << "m:\n" << m << std::endl;
    std::cout << "l:\n" << l << std::endl;
    std::cout << "lt:\n" << lt << std::endl;
    std::cout << "diagonal: " << diagonal << std::endl;

    MpMatrix d(dim, m_shift, ROW_ORIENTED);
    impose_diagonal(diagonal, d);
    std::cout << "ld:\n" << l << std::endl;

    MpMatrix ld(dim, m_shift, ROW_ORIENTED);
    multiply(l, d, ld);

    MpMatrix ldlt(dim, m_shift, ROW_ORIENTED);
    multiply(ld, lt, ldlt);
    
    std::cout << "ldlt:\n" << ldlt << std::endl;

    return 0;
}
