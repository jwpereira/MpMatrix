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

void demo_multiply(MpMatrix &m, size_t dim, fmp_shift_t m_shift) {
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
}

void the_project(MpMatrix &m) {
    auto dim = m.getDim();
    auto shift = m.getShift();

    // Initialize the matrix with the seeding function
    momentInit(m);
    std::cout << "Original Matrix:\n" << m << '\n';

    // Perform cholesky decomposition on the matrix
    cholesky_decompose(m);
    std::cout << "After Cholesky:\n" << m << '\n';

    // All this really does is help me keep the math straight lol
    auto &l = m;

    // Extract diagonals and impose them onto new matrix. Since we're inverting everything we'll
    // invert the diagonal here itself.
    MpArray diagonal(dim, shift);
    extract_diagonal(m, diagonal);
    MpMatrix d(dim, shift, ROW_ORIENTED);
    impose_diagonal(diagonal, d);

    // just for printing purposes:
    std::cout << "L:\n" << l;
    std::cout << "D:\n" << d;
    l.setMode(ROW_ORIENTED);
    std::cout << "Lt:\n" << l;
    l.setMode(COL_ORIENTED);

    // We'll first take the inverse of L to get L'
    reorient(l);    // first get L into row-oriented form
    invert(l);
    
    MpMatrix lt_inverse(l);
    transpose(lt_inverse);

    // Here we're actually inverting the diagonal
    invert_diagonal(diagonal);
    impose_diagonal(diagonal, d);

    std::cout << "\nInverted:\n";
    std::cout << "L':\n" << l;
    std::cout << "D':\n" << d;
    std::cout << "Lt':\n" << lt_inverse;

    // for max clarity, for me
    auto &l_inverse = l;
    auto &d_inverse = d;

    // Now multiply together to get inverse of m
    MpMatrix ltd_inverse(dim, shift, ROW_ORIENTED);
    multiply(lt_inverse, d_inverse, ltd_inverse);

    MpMatrix m_inverse(dim, shift, ROW_ORIENTED); // m' = (lt)'d'l'
    multiply(ltd_inverse, l_inverse, m_inverse);

    std::cout << "\nInverse of m:\n" << m_inverse;
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

    // auto dim = 3;
    // MpMatrix m(dim, m_shift);
    // fixedmpz m_raw[] = {
    //     4^256_fmpz,     0^256_fmpz,     0^256_fmpz, 
    //     12^256_fmpz,    37^256_fmpz,    0^256_fmpz,
    //     -(16^256_fmpz), -(43^256_fmpz), 98^256_fmpz 
    // };

    // size_t dim = 4;
    // MpMatrix m(dim, m_shift);
    // fixedmpz m_raw[] = {
    //     1^256_fmpz, 0^256_fmpz, 0^256_fmpz, 0^256_fmpz,  
    //     4^256_fmpz, 1^256_fmpz, 0^256_fmpz, 0^256_fmpz,
    //     6^256_fmpz, 2^256_fmpz, 1^256_fmpz, 0^256_fmpz, 
    //     3^256_fmpz, 7^256_fmpz, 4^256_fmpz, 1^256_fmpz
    // };

    // for (size_t col = 0; col < dim; col++) {
    //     for (size_t row = 0; row < dim; row++) {
    //         m[col][row] = m_raw[row * dim + col];
    //     }
    // }

    // demo(m, dim, m_shift);
    // demo_multiply(m, dim, m_shift);

    auto dim = 3;
    MpMatrix m(dim, m_shift);

    the_project(m);

    return 0;
}
