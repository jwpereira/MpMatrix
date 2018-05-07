#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "demo.hpp"
#include "eigen.hpp"
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
    MpMatrix source(dim, m_shift, COL_ORIENTED);

    // Initialize the matrix with the seeding function
    momentInit(source);
    MpMatrix m(source), m_inverse(dim, m_shift, ROW_ORIENTED); // copy original before it gets changed
    std::cout << "Original Matrix:\n" << source << '\n';

    inversion(m, m_inverse);

    std::vector<double> m_inverse_as_doubles(dim * dim);
    m_inverse.dumpVecDouble(m_inverse_as_doubles);

    print_matrix(m_inverse_as_doubles, dim);

    double smallest_eigenvalue = get_eigenvalue(m_inverse, SMALLEST);
    double inverse_of_largest_eigenvalue = 1.0 / get_eigenvalue(source, LARGEST);

    std::cout << "Smallest eigenvalue: " << smallest_eigenvalue << '\n';
    std::cout << "Inverse of largest: " << inverse_of_largest_eigenvalue << '\n';

    return 0;
}
