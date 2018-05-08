/**
 * @brief This file is just a demo for random functionality.
 *
 * Maybe it'll have its use some day, but now is not its time.
 *
 * @file demo.hpp
 * @author jwpereira
 */

#pragma once

#include <iostream>
#include <gmpxx.h>
#include "fixedmpz.hpp"
#include "mpmatrix.hpp"
#include "moment_algorithm.hpp"

using namespace momentmp;

void demo(MpMatrix &m, size_t dim, fmpz_shift_t m_shift) {
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

void demo_multiply(MpMatrix &m, size_t dim, fmpz_shift_t m_shift) {
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