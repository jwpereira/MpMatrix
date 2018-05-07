#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "demo.hpp"
#include "fixedmpz.hpp"
#include "mpmatrix.hpp"
#include "moment_algorithm.hpp"

using namespace momentmp;

void print_matrix(std::vector<double> &matrix, size_t dim) {
    std::cout << "\nAs doubles:\n";
    for (size_t i = 0; i < (dim * dim); i++) {
        if (i > 1 && i % dim == 0) {
            std::cout << '\n';
        } 
        std::cout << matrix[i] << '\t';
    }
    std::cout << '\n';
}

void inversion(MpMatrix &m, MpMatrix &m_inverse) {
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
    auto &l_inverse = l;    // for max clarity, for me
    
    MpMatrix lt_inverse(l);
    transpose(lt_inverse);

    // Here we're actually inverting the diagonal
    invert_diagonal(diagonal);
    impose_diagonal(diagonal, d);
    auto &d_inverse = d;    // for max clarity, for me

    std::cout << "\nInverted:\n";
    std::cout << "L':\n" << l_inverse;
    std::cout << "D':\n" << d_inverse;
    std::cout << "Lt':\n" << lt_inverse;

    // Now multiply together to get inverse of m
    MpMatrix ltd_inverse(dim, shift, ROW_ORIENTED);
    multiply(lt_inverse, d_inverse, ltd_inverse);
    multiply(ltd_inverse, l_inverse, m_inverse);    // m' = (lt)'d'l'

    std::cout << "\nInverse of m:\n" << m_inverse;
}

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    fmp_shift_t m_shift = 256;

    auto dim = 3;
    MpMatrix m(dim, m_shift, COL_ORIENTED), m_inverse(dim, m_shift, ROW_ORIENTED);

    inversion(m, m_inverse);

    std::vector<double> m_inverse_as_doubles(dim * dim);
    m_inverse.dump_vec_double(m_inverse_as_doubles);

    print_matrix(m_inverse_as_doubles, dim);

    return 0;
}
