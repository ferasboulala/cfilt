/**
 * Copyright 2020 Feras Boulala <ferasboulala@gmail.com>
 *
 * This file is part of cfilt.
 *
 * cfilt is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cfilt is distributed in the hope it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cfilt. If not, see <https://www.gnu.org/licenses/>.
 */

#include "cfilt/util.h"
#include "utest.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

int
test_cfilt_matrix_invert(void)
{
    gsl_matrix* src = gsl_matrix_alloc(2, 2);
    gsl_matrix* dst = gsl_matrix_alloc(2, 2);
    gsl_matrix* sol = gsl_matrix_alloc(2, 2);
    gsl_permutation* perm = gsl_permutation_alloc(2);

    gsl_matrix_set(src, 0, 0, 1.0);
    gsl_matrix_set(src, 0, 1, 2.0);
    gsl_matrix_set(src, 1, 0, 3.0);
    gsl_matrix_set(src, 1, 1, 4.0);

    gsl_matrix_set(sol, 0, 0, 4.0);
    gsl_matrix_set(sol, 0, 1, -2.0);
    gsl_matrix_set(sol, 1, 0, -3.0);
    gsl_matrix_set(sol, 1, 1, 1.0);
    gsl_matrix_scale(sol, -0.5);

    UTEST_EXEC_ASSERT(cfilt_matrix_invert, src, dst, perm);
    UTEST_EXEC_ASSERT(cfilt_matrix_cmp_tol, dst, sol, 0.1);

    gsl_matrix_free(src);
    gsl_matrix_free(dst);
    gsl_matrix_free(sol);
    gsl_permutation_free(perm);

    src = gsl_matrix_alloc(3, 3);
    dst = gsl_matrix_alloc(3, 3);
    sol = gsl_matrix_alloc(3, 3);
    perm = gsl_permutation_alloc(3);

    gsl_matrix_set_identity(src);
    gsl_matrix_set_identity(sol);

    UTEST_EXEC_ASSERT(cfilt_matrix_invert, src, dst, perm);
    UTEST_EXEC_ASSERT(cfilt_matrix_cmp_tol, dst, sol, 0.1);

    gsl_matrix_free(src);
    gsl_matrix_free(dst);
    gsl_matrix_free(sol);
    gsl_permutation_free(perm);

    return GSL_SUCCESS;
}

int
test_cfilt_matrix_tri_zero(void)
{
    gsl_matrix* src = gsl_matrix_alloc(3, 3);
    gsl_matrix* sol = gsl_matrix_alloc(3, 3);

    gsl_matrix_set_identity(src);
    gsl_matrix_set_identity(sol);

    gsl_matrix_set(src, 0, 1, 1.0);
    gsl_matrix_set(src, 0, 2, 1.0);
    gsl_matrix_set(src, 1, 2, 1.0);

    UTEST_EXEC_ASSERT(cfilt_matrix_tri_zero, src, 1);
    UTEST_EXEC_ASSERT(cfilt_matrix_cmp, src, sol);

    gsl_matrix_set_identity(src);
    gsl_matrix_set(src, 1, 0, 1.0);
    gsl_matrix_set(src, 2, 0, 1.0);
    gsl_matrix_set(src, 2, 1, 1.0);

    UTEST_EXEC_ASSERT(cfilt_matrix_tri_zero, src, 0);
    UTEST_EXEC_ASSERT(cfilt_matrix_cmp, src, sol);

    gsl_matrix_free(src);
    gsl_matrix_free(sol);

    return GSL_SUCCESS;
}

int
test_cfilt_matrix_cmp(void)
{
    gsl_matrix* a = gsl_matrix_alloc(3, 3);
    gsl_matrix* b = gsl_matrix_alloc(3, 3);

    gsl_matrix_set_zero(a);
    gsl_matrix_set_identity(b);

    UTEST_EXEC_ASSERT_(cfilt_matrix_cmp, a, b);

    gsl_matrix_set_zero(b);

    UTEST_EXEC_ASSERT(cfilt_matrix_cmp, a, b);

    gsl_matrix_free(b);
    b = gsl_matrix_alloc(3, 2);

    UTEST_EXEC_ASSERT_(cfilt_matrix_cmp, a, b);

    gsl_matrix_free(a);
    gsl_matrix_free(b);

    return GSL_SUCCESS;
}

int
test_cfilt_matrix_cmp_tol(void)
{
    gsl_matrix* a = gsl_matrix_alloc(3, 3);
    gsl_matrix* b = gsl_matrix_alloc(3, 3);

    gsl_matrix_set_zero(a);
    gsl_matrix_set_all(b, 0.01);

    UTEST_EXEC_ASSERT(cfilt_matrix_cmp_tol, a, b, 0.1);

    gsl_matrix_set_all(b, 1.0);

    UTEST_EXEC_ASSERT_(cfilt_matrix_cmp_tol, a, b, 0.1);

    gsl_matrix_free(b);
    b = gsl_matrix_alloc(3, 2);

    UTEST_EXEC_ASSERT_(cfilt_matrix_cmp_tol, a, b, 0.1);

    gsl_matrix_free(a);
    gsl_matrix_free(b);

    return GSL_SUCCESS;
}

int
test_cfilt_vector_cmp(void)
{
    gsl_vector* a = gsl_vector_alloc(3);
    gsl_vector* b = gsl_vector_alloc(3);

    gsl_vector_set_zero(a);
    gsl_vector_set_all(b, 1.0);

    UTEST_EXEC_ASSERT_(cfilt_vector_cmp, a, b);

    gsl_vector_set_zero(b);

    UTEST_EXEC_ASSERT(cfilt_vector_cmp, a, b);

    gsl_vector_free(b);
    b = gsl_vector_alloc(2);

    UTEST_EXEC_ASSERT_(cfilt_vector_cmp, a, b);

    gsl_vector_free(a);
    gsl_vector_free(b);

    return GSL_SUCCESS;
}

int
test_cfilt_vector_cmp_tol(void)
{
    gsl_vector* a = gsl_vector_alloc(3);
    gsl_vector* b = gsl_vector_alloc(3);

    gsl_vector_set_zero(a);
    gsl_vector_set_all(b, 0.05);

    UTEST_EXEC_ASSERT(cfilt_vector_cmp_tol, a, b, 0.1);

    gsl_vector_set_zero(b);

    UTEST_EXEC_ASSERT(cfilt_vector_cmp_tol, a, b, 0.1);

    gsl_vector_set_all(b, 1.0);

    UTEST_EXEC_ASSERT_(cfilt_vector_cmp_tol, a, b, 0.1);

    gsl_vector_free(b);
    b = gsl_vector_alloc(2);

    UTEST_EXEC_ASSERT_(cfilt_vector_cmp_tol, a, b, 0.1);

    gsl_vector_free(a);
    gsl_vector_free(b);

    return GSL_SUCCESS;
}

int
main(void)
{
    RUN_TEST(test_cfilt_matrix_invert);
    RUN_TEST(test_cfilt_matrix_tri_zero);
    RUN_TEST(test_cfilt_matrix_cmp);
    RUN_TEST(test_cfilt_matrix_cmp_tol);
    RUN_TEST(test_cfilt_vector_cmp);
    RUN_TEST(test_cfilt_vector_cmp_tol);

    return GSL_SUCCESS;
}
