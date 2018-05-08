#include <iostream>
#include <vector>

#include <gmpxx.h>

#include "demo.hpp"
#include "eigen.hpp"
#include "fixedmpz.hpp"
#include "mpmatrix.hpp"
#include "moment_algorithm.hpp"

using namespace momentmp;

const size_t INV_DIM = 10;

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
    // std::cout << "After Cholesky:\n" << m << '\n';

    // All this really does is help me keep the math straight lol
    auto &l = m;

    // Extract diagonals and impose them onto new matrix. Since we're inverting everything we'll
    // invert the diagonal here itself.
    MpArray diagonal(dim, shift);
    extract_diagonal(m, diagonal);

    std::cout << "last_diagonal: " << diagonal[diagonal.size() - 1] << std::endl;

    // We'll first take the inverse of L to get L'
    reorient(l);    // first get L into row-oriented form
    invert(l);
    auto &l_inverse = l;    // for max clarity, for me
    
    MpMatrix lt_inverse(l);
    transpose(lt_inverse);

    auto zero = 0^fmpshift(shift);
    for (size_t i = 0; (i < INV_DIM && i < dim); i++) {
        for (size_t j = 0; (j < INV_DIM && j < dim); j++) {
            auto sum = zero;
            for (size_t k = 0; k < dim; k++) {
                sum += lt_inverse[i][k] * l_inverse[k][j] / diagonal[k];
            }
            m_inverse[i][j] = sum;
        }
    }
}

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    fmp_shift_t m_shift = 4096;

    auto dim = 100;
    MpMatrix source(dim, m_shift, COL_ORIENTED);

    // Initialize the matrix with the seeding function
    momentInit(source);
    MpMatrix m(source), m_inverse(INV_DIM, m_shift, ROW_ORIENTED); // copy original before it gets changed
    // std::cout << "Original Matrix:\n" << source << '\n';

    inversion(m, m_inverse);

    std::vector<double> m_inverse_section_as_doubles(INV_DIM * INV_DIM);
    m_inverse.dumpVecDouble(m_inverse_section_as_doubles);

    print_matrix(m_inverse_section_as_doubles, INV_DIM);

    // double smallest_eigenvalue = get_eigenvalue(source, SMALLEST);
    double inverse_of_largest_eigenvalue = 1.0 / get_eigenvalue(m_inverse, LARGEST);

    // std::cout << "Smallest eigenvalue: " << smallest_eigenvalue << '\n';
    std::cout << "Inverse of largest: " << inverse_of_largest_eigenvalue << '\n';

    return 0;
}
