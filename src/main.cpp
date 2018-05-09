/**
 * @brief Main driver file for the entire HankelHacker project
 * 
 * This file puts it all together.
 * 
 * @file main.cpp
 * @author jwpereira
 */

#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>

#include <gmpxx.h>

#include "demo.hpp"
#include "eigen.hpp"
#include "fixedmpz.hpp"
#include "mpmatrix.hpp"
#include "moment_algorithm.hpp"

using namespace momentmp;

/**
 * @brief This is the size of the matrix to be handed over to the GSL eigensolver
 */
const size_t INV_DIM = 10;

/**
 * @brief Boolean for whether to print log statements or not to stderr
 */
const bool DEBUG = true;

/**
 * @brief Convenience function for printing out a vector-based matrix
 * 
 * This is purely for matrices that are actually std::vector<double>. MpMatrix has its own
 * << overload that enables it to be printed out!
 */
void print_matrix(std::vector<double> &matrix, size_t dim) {
    for (size_t i = 0; i < (dim * dim); i++) {
        if (i > 1 && i % dim == 0) {
            std::cout << '\n';
        }
        std::cout << matrix[i] << '\t';
    }
    std::cout << '\n';
}

/**
 * @brief Main routine for inverting the source matrix
 *
 * This function has the matrix decomposed into essentially LDLT form (via cholesky decomposition).
 * The decomposition turns the input matrix into L with D superimposed on it. D is extracted out of
 * L (leaving 1s in L's diagonal). D is not an actual MpMatrix but rather an MpArray, simply because
 * it has all zeros except for the diagonal itself. By default, the input matrix is in
 * column-oriented form such that the cholesky decomposition works. However, because the rest of the
 * operations being performed will be more optimal in row-oriented form, L is reoriented to
 * row-oriented form. From there, L is inverted to get L'. L' is then copied and then transposed to
 * get (Lt)'.
 *
 * Typically by glove's rule, we could get the original matrix's inverse by three matrix
 * multiplications: M'=(Lt)'D'L'. However, since we will really only be interested in the largest
 * eigenvalue, we can focus on solely the first 10x10 (or the value of INV_DIM) upper left section
 * of the inverted matrix. This is also why D is not inverted, as 1/D simply means dividing by D.
 *
 * Input matrix m is changed by this procedure, m_inverse will hold the 10x10 (or otherwised defined
 * size) of the inverse of m. m_inverse should be passed onto the eigensolver.
 */
void inversion(MpMatrix &m, MpMatrix &m_inverse) {
    auto dim = m.getDim();
    auto shift = m.getShift();

    // Perform cholesky decomposition on the matrix
        if (DEBUG) std::cerr << "Cholesky-decompose input matrix... ";
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

    // We'll take the inverse of L to get L'
        if (DEBUG) std::cerr << "Transposing L into row-oriented form... ";
    reorient(l);                                            // first get L into row-oriented form
        if (DEBUG) std::cerr << "done!\n";    
        if (DEBUG) std::cerr << "Inverting L to get L'... ";
    invert(l);      
        if (DEBUG) std::cerr << "done!\n";
    auto &l_inverse = l;    // for max clarity, for me

    // Then we'll take the transpose of that to get (Lt)'
        if (DEBUG) std::cerr << "Transposing L' to get (Lt)'... ";
    MpMatrix lt_inverse(l);
    transpose(lt_inverse);
        if (DEBUG) std::cerr << "done!\n";

    // Calculate out the first 10x10 (or INV_DIMxINV_DIM) of the inverse of the source matrix
        if (DEBUG) std::cerr << "Creating first " << INV_DIM << "x" << INV_DIM << " of inverse of M... ";
    auto zero = 0^fmpzshift(shift);
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

    // HankelHacker must be launched with 2 extra arguments <dim> and <shift>
    if (argc < 3) {
        std::cerr << "Error: Missing arguments.\n";
        std::cerr << "Usage: hankelhacker <dimension of source> <shift amount>\n";
        return -1;
    }

    // Take command line arguments and store them
    auto dim = strtoul(argv[1], NULL, 10);
    fmpz_shift_t m_shift = strtoul(argv[2], NULL, 10);

    std::cout << "Size of matrix: " << dim << " by " << dim << "\n";
    std::cout << "Shift: " << m_shift << "\n";

    // Since time is of interest, note the start time
    auto start_time = std::chrono::high_resolution_clock::now();

    // Initialize the matrix with the seeding function
        if (DEBUG) std::cerr << "Generating source matrix... ";
    MpMatrix m(dim, m_shift, COL_ORIENTED);
    momentInit(m);
    MpMatrix m_inverse(INV_DIM, m_shift, ROW_ORIENTED);
        if (DEBUG) std::cerr << "done!\n";

    // Invert the source matrix
        if (DEBUG) std::cerr << "Finding inverse of source matrix:\n";
    inversion(m, m_inverse);

    // Extract the largest eigenvalue
        if (DEBUG) std::cerr << "Extracting largest eigenvalue... ";
    double inverse_of_largest_eigenvalue = 1.0 / get_eigenvalue(m_inverse, LARGEST);
        if (DEBUG) std::cerr << "done!\n";

    std::cout << "Inverse of largest: " << std::setprecision(15) <<  std::scientific
              << inverse_of_largest_eigenvalue << '\n';

    // Print out the time it took to do all these operations by taking the difference in times
    auto finish_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = finish_time - start_time;
    std::cout << "Completed in " << elapsed_time.count() << " seconds\n";

    return 0;
}
