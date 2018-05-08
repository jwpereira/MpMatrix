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
const bool DEBUG = true;

void print_matrix(std::vector<double> &matrix, size_t dim) {
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

    if (DEBUG) std::cerr << "Cholesky-decompose input matrix... ";
    // Perform cholesky decomposition on the matrix
    cholesky_decompose(m);
    if (DEBUG) std::cerr << "done!\n";

    // All this really does is help me keep the math straight lol
    auto &l = m;

    // Extract diagonals and impose them onto new matrix. Since we're inverting everything we'll
    // invert the diagonal here itself.
    if (DEBUG) std::cerr << "Extracting diagonals... ";
    MpArray diagonal(dim, shift);
    extract_diagonal(m, diagonal);
    if (DEBUG) std::cerr << "done!\n";

    std::cout << "last diagonal: " << diagonal[diagonal.size() - 1] << std::endl;

    // We'll first take the inverse of L to get L'
    if (DEBUG) std::cerr << "Transposing L into row-oriented form... ";
    reorient(l);    if (DEBUG) std::cout << "done!\n";    // first get L into row-oriented form    
    if (DEBUG) std::cerr << "Inverting L to get L'... ";
    invert(l);      if (DEBUG) std::cout << "done!\n";
    auto &l_inverse = l;    // for max clarity, for me
    
    if (DEBUG) std::cerr << "Transposing L' to get (Lt)'... ";
    MpMatrix lt_inverse(l);
    transpose(lt_inverse);
    if (DEBUG) std::cerr << "done!\n";

    if (DEBUG) std::cerr << "Creating first" << INV_DIM << "x" << INV_DIM << "of inverse of M... ";
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
    if (DEBUG) std::cerr << "done!\n";
}

int main(int argc, char *argv[]) {
    // Not using printf, therefore no need to have cout sync with stdio ->
    // better performance
    std::ios_base::sync_with_stdio(false);

    fmp_shift_t m_shift = 4096;

    auto dim = 500;
    MpMatrix source(dim, m_shift, COL_ORIENTED);

    // Initialize the matrix with the seeding function
    if (DEBUG) std::cerr << "Generating source matrix... ";
    momentInit(source);
    MpMatrix m(source), m_inverse(INV_DIM, m_shift, ROW_ORIENTED); // copy original before it gets changed
    if (DEBUG) std::cerr << "done!\n";

    if (DEBUG) std::cerr << "Finding inverse of source matrix:\n";
    inversion(m, m_inverse);

    if (DEBUG) std::cerr << "Extracting largest eigenvalue... ";
    double inverse_of_largest_eigenvalue = 1.0 / get_eigenvalue(m_inverse, LARGEST);
    if (DEBUG) std::cerr << "done!\n";

    // std::cout << "Smallest eigenvalue: " << smallest_eigenvalue << '\n';
    std::cout << "Inverse of largest: " << inverse_of_largest_eigenvalue << '\n';

    return 0;
}
