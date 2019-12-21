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

#include "cfilt/kalman.h"
#include "cfilt/util.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

#include <string.h>

#define V_ALLOC_ASSERT_(p, n)                                                  \
    V_ALLOC_ASSERT(p, n, cfilt_kalman_filter_free, filt)
#define M_ALLOC_ASSERT_(p, n, m)                                               \
    M_ALLOC_ASSERT(p, n, m, cfilt_kalman_filter_free, filt)

int
cfilt_kalman_filter_alloc(cfilt_kalman_filter* filt, const size_t n,
                          const size_t m, const size_t k)
{
    memset(filt, 0, sizeof(cfilt_kalman_filter));

    M_ALLOC_ASSERT_(filt->F, n, n);
    M_ALLOC_ASSERT_(filt->P_, n, n);
    M_ALLOC_ASSERT_(filt->B, n, m);
    M_ALLOC_ASSERT_(filt->Q, n, n);
    M_ALLOC_ASSERT_(filt->P, n, n);
    M_ALLOC_ASSERT_(filt->H, k, n);
    M_ALLOC_ASSERT_(filt->R, k, k);
    M_ALLOC_ASSERT_(filt->K, n, k);

    V_ALLOC_ASSERT_(filt->x, n);
    V_ALLOC_ASSERT_(filt->x_, n);
    V_ALLOC_ASSERT_(filt->z, k);
    V_ALLOC_ASSERT_(filt->u, m);
    V_ALLOC_ASSERT_(filt->y, k);

    M_ALLOC_ASSERT_(filt->_FP, n, n);
    M_ALLOC_ASSERT_(filt->_PH_T, n, k);
    M_ALLOC_ASSERT_(filt->_PH_T_R, k, k);
    M_ALLOC_ASSERT_(filt->_inv, k, k);
    M_ALLOC_ASSERT_(filt->_I, n, n);

    filt->_perm = gsl_permutation_alloc(k);
    if (filt->_perm == NULL)
    {
        return GSL_ENOMEM;
    }

    return GSL_SUCCESS;
}

void
cfilt_kalman_filter_free(cfilt_kalman_filter* filt)
{
    M_FREE_IF_NOT_NULL(filt->F);
    M_FREE_IF_NOT_NULL(filt->B);
    M_FREE_IF_NOT_NULL(filt->Q);
    M_FREE_IF_NOT_NULL(filt->P);
    M_FREE_IF_NOT_NULL(filt->H);
    M_FREE_IF_NOT_NULL(filt->R);
    M_FREE_IF_NOT_NULL(filt->K);

    V_FREE_IF_NOT_NULL(filt->x);
    V_FREE_IF_NOT_NULL(filt->x_);
    V_FREE_IF_NOT_NULL(filt->z);
    V_FREE_IF_NOT_NULL(filt->u);
    V_FREE_IF_NOT_NULL(filt->y);

    M_FREE_IF_NOT_NULL(filt->_FP);
    M_FREE_IF_NOT_NULL(filt->_PH_T);
    M_FREE_IF_NOT_NULL(filt->_PH_T_R);
    M_FREE_IF_NOT_NULL(filt->_inv);
    M_FREE_IF_NOT_NULL(filt->_I);

    if (filt->_perm)
    {
        gsl_permutation_free(filt->_perm);
    }
}

int
cfilt_kalman_filter_predict(cfilt_kalman_filter* filt)
{
    // x_ = Fx + Bu
    EXEC_ASSERT(gsl_blas_dgemv, CblasNoTrans, 1.0, filt->F, filt->x, 0.0,
                filt->x_);
    EXEC_ASSERT(gsl_blas_dgemv, CblasNoTrans, 1.0, filt->B, filt->u, 1.0,
                filt->x_);

    // P_ = FPF^T + Q
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasNoTrans, 1.0, filt->F,
                filt->P, 0.0, filt->_FP);
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasTrans, 1.0, filt->_FP,
                filt->F, 0.0, filt->P_);
    EXEC_ASSERT(gsl_matrix_add, filt->P_, filt->Q);

    return GSL_SUCCESS;
}

int
cfilt_kalman_filter_update(cfilt_kalman_filter* filt)
{
    // K = P_H^T(HP_H^T + R)^-1
    // _PH_T_R is used to avoid changing R
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasTrans, 1.0, filt->P_,
                filt->H, 0.0, filt->_PH_T);
    EXEC_ASSERT(gsl_matrix_memcpy, filt->_PH_T_R, filt->R);
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasNoTrans, 1.0, filt->H,
                filt->_PH_T, 1.0, filt->_PH_T_R);
    EXEC_ASSERT(cfilt_matrix_invert, filt->_PH_T_R, filt->_inv, filt->_perm);
    EXEC_ASSERT(gsl_matrix_memcpy, filt->_PH_T_R, filt->_inv);
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasNoTrans, 1.0, filt->_PH_T,
                filt->_PH_T_R, 0.0, filt->K);

    // y = z - Hx_
    // y is used to avoid changing z
    EXEC_ASSERT(gsl_vector_memcpy, filt->y, filt->z);
    EXEC_ASSERT(gsl_blas_dgemv, CblasNoTrans, -1.0, filt->H, filt->x_, 1.0,
                filt->y);

    // x = x_ + Ky
    // x_ is copied over to x to avoid changing x_
    EXEC_ASSERT(gsl_vector_memcpy, filt->x, filt->x_);
    EXEC_ASSERT(gsl_blas_dgemv, CblasNoTrans, 1.0, filt->K, filt->y, 1.0,
                filt->x);

    // P = (I - KH)P_
    gsl_matrix_set_identity(filt->_I);
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasNoTrans, -1.0, filt->K,
                filt->H, 1.0, filt->_I);
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasNoTrans, 1.0, filt->_I,
                filt->P_, 0.0, filt->P);

    return GSL_SUCCESS;
}
