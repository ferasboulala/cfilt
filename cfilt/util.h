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

#ifndef CFILT_UTIL_H_
#define CFILT_UTIL_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

#define FREE_IF_NOT_NULL(p, func)                                              \
    if (p)                                                                     \
        func(p);
#define V_FREE_IF_NOT_NULL(v) FREE_IF_NOT_NULL(v, gsl_vector_free)
#define M_FREE_IF_NOT_NULL(m) FREE_IF_NOT_NULL(m, gsl_matrix_free)

#define M_ALLOC_ASSERT(p, n, m, func, ...)                                     \
    do                                                                         \
    {                                                                          \
        p = gsl_matrix_alloc((n), (m));                                        \
        if (p == NULL)                                                         \
        {                                                                      \
            func(__VA_ARGS__);                                                 \
            return GSL_ENOMEM;                                                 \
        }                                                                      \
    } while (0);

#define V_ALLOC_ASSERT(p, n, func, ...)                                        \
    do                                                                         \
    {                                                                          \
        p = gsl_vector_alloc(n);                                               \
        if (p == NULL)                                                         \
        {                                                                      \
            func(__VA_ARGS__);                                                 \
            return GSL_ENOMEM;                                                 \
        }                                                                      \
    } while (0);

#define EXEC_ASSERT(func, ...)                                                 \
    do                                                                         \
    {                                                                          \
        const int status_ = func(__VA_ARGS__);                                 \
        if (status_ != GSL_SUCCESS)                                            \
        {                                                                      \
            return status_;                                                    \
        }                                                                      \
    } while (0);

int cfilt_matrix_invert(gsl_matrix* src, gsl_matrix* dst,
                        gsl_permutation* perm);

int cfilt_matrix_tri_zero(gsl_matrix* src, int upper);

int cfilt_matrix_cmp(const gsl_matrix* a, const gsl_matrix* b);

int cfilt_vector_cmp(const gsl_vector* a, const gsl_vector* b);

#endif // CFILT_UTIL_H_
