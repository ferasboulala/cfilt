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

#include <string.h>

#include <sys/types.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

/**
 * Some of the functions in this module could be optimized. The choice of
 * implementation
 * is not entirely arbitrary. I opted for simplicity and readability.
 *
 * These functions are not meant to be used outside of cfilt modules and most
 * of them are used in unit testing anyway (comparison functions for instance).
 */

int
cfilt_matrix_invert(gsl_matrix* src, gsl_matrix* dst, gsl_permutation* perm)
{
    int signum;
    EXEC_ASSERT(gsl_linalg_LU_decomp, src, perm, &signum);
    EXEC_ASSERT(gsl_linalg_LU_invert, src, perm, dst);

    return GSL_SUCCESS;
}

int
cfilt_matrix_tri_zero(gsl_matrix* src, int upper)
{
    if (!upper)
    {
        for (size_t j = 0; j < src->size2 - 1; ++j)
        {
            for (size_t i = j + 1; i < src->size1; ++i)
            {
                gsl_matrix_set(src, i, j, 0.0);
            }
        }
    }
    else
    {
        for (size_t i = 0; i < src->size1 - 1; ++i)
        {
            for (size_t j = i + 1; j < src->size2; ++j)
            {
                gsl_matrix_set(src, i, j, 0.0);
            }
        }
    }

    return GSL_SUCCESS;
}

int
cfilt_matrix_cmp(gsl_matrix* a, gsl_matrix* b)
{
    return cfilt_matrix_cmp_tol(a, b, 0.0);
}

int
cfilt_matrix_cmp_tol(const gsl_matrix* a, const gsl_matrix* b, const double tol)
{
    if (a->size1 != b->size2 || a->size2 != b->size2 || a->tda != b->tda)
    {
        return GSL_EBADLEN;
    }

    for (size_t i = 0; i < a->size1; ++i)
    {
        for (size_t j = 0; j < a->size2; ++j)
        {
            const double x = gsl_matrix_get(a, i, j);
            const double y = gsl_matrix_get(b, i, j);

            if (!IS_EQ_TOL(x, y, tol))
            {
                return GSL_EFAILED;
            }
        }
    }

    return GSL_SUCCESS;
}

int
cfilt_vector_cmp(const gsl_vector* a, const gsl_vector* b)
{
    return cfilt_vector_cmp_tol(a, b, 0.0);
}

int
cfilt_vector_cmp_tol(const gsl_vector* a, const gsl_vector* b, const double tol)
{
    if (a->size != b->size || a->stride != b->stride)
    {
        return GSL_EBADLEN;
    }

    for (size_t i = 0; i < a->size; ++i)
    {
        const double x = gsl_vector_get(a, i);
        const double y = gsl_vector_get(b, i);

        if (!IS_EQ_TOL(x, y, tol))
        {
            return GSL_EFAILED;
        }
    }

    return GSL_SUCCESS;
}

void
cfilt_fprintf_matrix_rows(FILE* file, const gsl_matrix *mat)
{
    if (mat->size1 * mat->size2 == 0)
    {
        return;
    }

    for (size_t i = 0; i < mat->size1; ++i)
    {
        gsl_vector_const_view view = gsl_matrix_const_row(mat, i);
        cfilt_fprintf_vector_row(file, &view.vector);
    }
}

void
cfilt_fprintf_vector_row(FILE* file, const gsl_vector *vec)
{
    if (vec->size == 0)
    {
        return;
    }

    for (size_t i = 0; i < vec->size - 1; ++i)
    {
        fprintf(file, "%f,", gsl_vector_get(vec, i));
    }

    fprintf(file, "%f\n", gsl_vector_get(vec, vec->size-1));
}

void
cfilt_printf_matrix_rows(const gsl_matrix *mat)
{
    cfilt_fprintf_matrix_rows(stdout, mat);
}

void
cfilt_printf_vector_row(const gsl_vector *vec)
{
    cfilt_fprintf_vector_row(stdout, vec);
}
