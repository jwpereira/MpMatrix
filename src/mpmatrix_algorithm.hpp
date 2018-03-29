#pragma once 

#include <functional>
#include "mpmatrix.hpp"

namespace momentmp {

    /**
     * @brief Conjugate transpose a matrix
     * 
     * TODO: Look into more efficient methods of transposing matricies
     * 
     * @param[in] initial source matrix
     * @param[out] lower destination matrix
     */
    inline bool transpose(const MpMatrix &src, MpMatrix &dest) {
        if (src.getDimension() != dest.getDimension()) {
            throw std::runtime_error("Unable to transpose square matrix to different size matrix");
        }
        
        auto dim = src.getDimension();
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                dest(j, i) = src(i, j);
            }
        }

        return true;
    }

    inline MpMatrix transpose(const MpMatrix &src) {
        MpMatrix dest(src.getDimension(), src.getShift());
        transpose(src, dest);
        return dest;
    }

    /**
     * @brief Function wrapper alias for the \link apply \endlink function
     * 
     * functions that can be assigned to applicator_f have signatures similar to:
     * <code>function(\link fmp_t \endlink&, size_t row, size_t col)</code>
     */
    using applicator_f = std::function<bool(fmp_t&, size_t, size_t)>;

    /**
     * @brief Applies a function (with specific arguments) at every row,col in an mpmatrix.
     *
     * The function to be applied must have a signature similar to : function(\link fmp_t \endlink&,
     * size_t row, size_t col). This means that it takes in an element of the destination matrix by
     * reference 
     *
     */
    inline bool apply(MpMatrix &dest, const applicator_f generator) {
        const size_t dim = dest.getDimension();
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                bool result = generator(dest(i, j), i, j);
            }
        }

        return true;
    }

    /**
     * @brief Applicator functor that turns an \link MpMatrix \endlink into an identity matrix
     */
    inline applicator_f identity = [](auto &fmp, auto row, auto col) {
        if (row == col) {
            fmp = (1^fmpshift(fmp.getShift()));
        } else {
            fmp = (0^fmpshift(fmp.getShift()));
        }
        return true;
    };

    /**
     * @brief Cholesky decomposition for MpMatrix
     * 
     * @param[in] initial source matrix
     * @param[out] lower destination matrix
     * @return true (Exceptions aside) was able to perform decomposition
     * @return false (Exceptions aside) was unable to perform decomposition
     */
    bool cholesky(const MpMatrix &initial, MpMatrix &lower);
}
