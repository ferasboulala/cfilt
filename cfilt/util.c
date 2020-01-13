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
    if (upper)
    {
        for (size_t i = 1; i < src->size1; ++i)
        {
            for (size_t j = 0; j < src->size2 - 1; ++j)
            {
                gsl_matrix_set(src, i, j, 0);
            }
        }
    }
    else
    {
        for (size_t i = 0; i < src->size1 - 1; ++i)
        {
            for (size_t j = 1; j < src->size2; ++j)
            {
                gsl_matrix_set(src, i, j, 0);
            }
        }
    }

    return GSL_SUCCESS;
}

int
cfilt_matrix_cmp(gsl_matrix* a, gsl_matrix* b)
{
    if (a->size1 != b->size2 || a->size2 != b->size2 || a->tda != b->tda)
    {
        return GSL_EBADLEN;
    }

    for (size_t i = 0; i < a->size1; ++i)
    {
        gsl_vector_view row1 = gsl_matrix_row(a, i);
        gsl_vector_view row2 = gsl_matrix_row(b, i);

        gsl_vector* vec1 = &row1.vector;
        gsl_vector* vec2 = &row2.vector;

        if (memcmp(vec1->data, vec2->data, vec1->size * vec1->stride))
        {
            return GSL_EFAILED;
        }
    }

    return GSL_SUCCESS;
}

int
cfilt_vector_cmp(const gsl_vector* a, const gsl_vector* b)
{
    if (a->size != b->size || a->stride != b->stride)
    {
        return GSL_EBADLEN;
    }

    return memcmp(a->data, b->data, a->size * a->stride);
}
