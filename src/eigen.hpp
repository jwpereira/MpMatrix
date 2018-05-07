#pragma once

#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "mpmatrix.hpp"

namespace momentmp {
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

    enum EigenMode : bool { SMALLEST=true, LARGEST=false };

    inline double get_eigenvalue(const MpMatrix &matrix, const EigenMode mode) {
        auto dim = matrix.getDim();
        std::vector<double> eigenvalues(dim);

        eigen_solve(matrix, eigenvalues);
        double est = eigenvalues[0];

        if (mode == SMALLEST) {
            for (size_t index = 1; index < dim; index++) {
                if (eigenvalues[index] < est) {
                    est = eigenvalues[index];
                }
            }
        } else if (mode == LARGEST) {
            for (size_t index = 1; index < dim; index++) {
                if (eigenvalues[index] > est) {
                    est = eigenvalues[index];
                }
            }
        }

        return est;
    }
}
