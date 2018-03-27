#pragma once 

#include <functional>
#include "mpmatrix.hpp"

namespace momentmp {

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
     * @param[in] initial Address of source matrix
     * @param[out] lower Address of destination matrix
     * @return true (Exceptions aside) was able to perform decomposition
     * @return false (Exceptions aside) was unable to perform decomposition
     */
    bool cholesky(const MpMatrix &initial, MpMatrix &lower);
}
