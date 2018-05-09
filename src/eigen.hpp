/**
 * @brief Functions for extracting the eigenvalues out of an MpMatrix
 *
 * @file eigen.hpp
 * @author jwpereira
 */

#pragma once

#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "mpmatrix.hpp"

namespace momentmp {

    /**
     * @brief Take in an MpMatrix and return an std::vector of its eigenvalues
     *
     * First the MpMatrix is brought down to just double-precision, and then it is handed off to
     * the GSL eigensolver for it to do its magic.
     */
    inline void eigen_solve(const MpMatrix &matrix, std::vector<double> &eigenvalues) {
        auto dim = matrix.getDim();
        std::vector<double> m(dim * dim);
        matrix.dumpVecDouble(m);

        gsl_vector *eval = gsl_vector_alloc(dim);
        gsl_eigen_symm_workspace *workspace = gsl_eigen_symm_alloc(dim);
        gsl_matrix_view gsl_m;

        gsl_m = gsl_matrix_view_array(&m[0], dim, dim);
        gsl_eigen_symm(&gsl_m.matrix, eval, workspace);
        gsl_eigen_symm_free(workspace);

        for (size_t index = 0; index < dim; index++) {
            eigenvalues[index] = gsl_vector_get(eval, index);
        }
    }

    /**
     * @brief For this project, we would only be interested in the largest/smallest eigenvalue
     */
    enum EigenMode : bool { SMALLEST=true, LARGEST=false };

    /**
     * @brief Calls eigen_solve, but then returns either the smallest/largest eigenvalue it finds
     */
    inline double get_eigenvalue(const MpMatrix &matrix, const EigenMode mode) {
        auto dim = matrix.getDim();
        std::vector<double> eigenvalues(dim);

        eigen_solve(matrix, eigenvalues);

        if (mode == SMALLEST) {
            return *(std::min_element(eigenvalues.begin(), eigenvalues.end()));
        } else if (mode == LARGEST) {
            return *(std::max_element(eigenvalues.begin(), eigenvalues.end()));
        }

        // Hopefully this never happens
        throw std::runtime_error("Desired output neither SMALLEST nor LARGEST");
    }
}
