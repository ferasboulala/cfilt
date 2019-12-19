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

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

#include <string.h>

#define FREE_IF_NOT_NULLV(v)                                                                                           \
    if (v)                                                                                                             \
        gsl_vector_free(v);
#define FREE_IF_NOT_NULLM(m)                                                                                           \
    if (m)                                                                                                             \
        gsl_matrix_free(m);

int
cfilt_kalman_alloc(cfilt_kalman_filter* filt, const size_t n, const size_t m, const size_t k)
{
    memset(filt, 0, sizeof(cfilt_kalman_filter));

    filt->F = gsl_matrix_calloc(n, n);
    filt->P_ = gsl_matrix_calloc(n, n);
    filt->B = gsl_matrix_calloc(n, m);
    filt->Q = gsl_matrix_calloc(n, n);
    filt->P = gsl_matrix_calloc(n, n);
    filt->H = gsl_matrix_calloc(k, n);
    filt->R = gsl_matrix_calloc(k, k);
    filt->K = gsl_matrix_calloc(n, k);

    filt->x = gsl_vector_calloc(n);
    filt->x_ = gsl_vector_calloc(n);
    filt->z = gsl_vector_calloc(k);
    filt->u = gsl_vector_calloc(m);
    filt->y = gsl_vector_calloc(k);

    filt->_FP = gsl_matrix_alloc(n, n);
    filt->_PH_T = gsl_matrix_alloc(n, k);
    filt->_PH_T_R = gsl_matrix_alloc(k, k);
    filt->_inv = gsl_matrix_alloc(k, k);
    filt->_perm = gsl_permutation_alloc(k);
    filt->_I = gsl_matrix_alloc(n, n);

    if (!filt->F || !filt->B || !filt->Q || !filt->P || !filt->H || !filt->R || !filt->x || !filt->x_ || !filt->z ||
        !filt->u || !filt->y || !filt->_FP || !filt->_PH_T || !filt->_I || !filt->_perm || !filt->_inv)
    {
        GSL_ERROR("failed to allocate space for kalman filter matrices", GSL_ENOMEM);
    }

    return GSL_SUCCESS;
}

void
cfilt_kalman_free(cfilt_kalman_filter* filt)
{
    FREE_IF_NOT_NULLM(filt->F);
    FREE_IF_NOT_NULLM(filt->B);
    FREE_IF_NOT_NULLM(filt->Q);
    FREE_IF_NOT_NULLM(filt->P);
    FREE_IF_NOT_NULLM(filt->H);
    FREE_IF_NOT_NULLM(filt->R);
    FREE_IF_NOT_NULLM(filt->K);

    FREE_IF_NOT_NULLV(filt->x);
    FREE_IF_NOT_NULLV(filt->x_);
    FREE_IF_NOT_NULLV(filt->z);
    FREE_IF_NOT_NULLV(filt->u);
    FREE_IF_NOT_NULLV(filt->y);

    FREE_IF_NOT_NULLM(filt->_FP);
    FREE_IF_NOT_NULLM(filt->_PH_T);
    FREE_IF_NOT_NULLM(filt->_PH_T_R);
    FREE_IF_NOT_NULLM(filt->_inv);
    FREE_IF_NOT_NULLM(filt->_I);

    if (filt->_perm)
    {
        gsl_permutation_free(filt->_perm);
    }
}

int
cfilt_kalman_predict(cfilt_kalman_filter* filt)
{
    // x_ = Fx + Bu
    if (gsl_blas_dgemv(CblasNoTrans, 1.0, filt->F, filt->x, 0.0, filt->x_))
    {
        GSL_ERROR("failed to compute Fx + [0] => x_", GSL_EFAILED);
    }

    if (gsl_blas_dgemv(CblasNoTrans, 1.0, filt->B, filt->u, 1.0, filt->x_))
    {
        GSL_ERROR("failed to compute Bu + x_ => x_", GSL_EFAILED);
    }

    // P_ = FPF^T + Q
    if (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, filt->F, filt->P, 0.0, filt->_FP))
    {
        GSL_ERROR("failed to compute FP + [0] => _FP", GSL_EFAILED);
    }

    if (gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, filt->_FP, filt->F, 0.0, filt->P_))
    {
        GSL_ERROR("failed to compute _FPF^T + [0] => P_", GSL_EFAILED);
    }

    if (gsl_matrix_add(filt->P_, filt->Q))
    {
        GSL_ERROR("failed to add Q into P_", GSL_EFAILED);
    }

    return GSL_SUCCESS;
}

#include <stdio.h>
int
cfilt_kalman_update(cfilt_kalman_filter* filt)
{
    // K = P_H^T(HP_H^T + R)^-1
    // _PH_T_R is used to avoid changing R
    if (gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, filt->P_, filt->H, 0.0, filt->_PH_T))
    {
        GSL_ERROR("failed to compute P_H^T + [0] => _PH_T", GSL_EFAILED);
    }

    if (gsl_matrix_memcpy(filt->_PH_T_R, filt->R))
    {
        GSL_ERROR("failed to copy R into _PH_T_R", GSL_EFAILED);
    }

    if (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, filt->H, filt->_PH_T, 1.0, filt->_PH_T_R))
    {
        GSL_ERROR("failed to compute H _PH_T + _PH_T_R => _PH_T_R", GSL_EFAILED);
    }

    // TODO : Make the matrix inversion a function if it is reused somewhere else
    int signum;
    if (gsl_linalg_LU_decomp(filt->_PH_T_R, filt->_perm, &signum))
    {
        GSL_ERROR("failed to compute the LU decomposition of _PH_T_R", GSL_EFAILED);
    }

    if (gsl_linalg_LU_invert(filt->_PH_T_R, filt->_perm, filt->_inv))
    {
        GSL_ERROR("failed to compute the inverse of R", GSL_EFAILED);
    }

    if (gsl_matrix_memcpy(filt->_PH_T_R, filt->_inv))
    {
        GSL_ERROR("failed to copy inv into _PH_T_R", GSL_EFAILED);
    }

    if (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, filt->_PH_T, filt->_PH_T_R, 0.0, filt->K))
    {
        GSL_ERROR("failed to compute P_H^T _PH_T_R", GSL_EFAILED);
    }

    // y = z - Hx_
    // y is used to avoid changing z
    if (gsl_vector_memcpy(filt->y, filt->z))
    {
        GSL_ERROR("failed to copy z into y", GSL_EFAILED);
    }

    if (gsl_blas_dgemv(CblasNoTrans, -1.0, filt->H, filt->x_, 1.0, filt->y))
    {
        GSL_ERROR("failed to compute -Hx_ + y => y", GSL_EFAILED);
    }

    // x = x_ + Ky
    // x_ is copied over to x to avoid changing x_
    if (gsl_vector_memcpy(filt->x, filt->x_))
    {
        GSL_ERROR("failed to copy x_ into x", GSL_EFAILED);
    }

    if (gsl_blas_dgemv(CblasNoTrans, 1.0, filt->K, filt->y, 1.0, filt->x))
    {
        GSL_ERROR("failed to compute Ky + x => x", GSL_EFAILED);
    }

    // P = (I - KH)P_
    gsl_matrix_set_identity(filt->_I);

    if (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, filt->K, filt->H, 1.0, filt->_I))
    {
        GSL_ERROR("failed to compute -KH + I => I", GSL_EFAILED);
    }

    if (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, filt->_I, filt->P_, 0.0, filt->P))
    {
        GSL_ERROR("failed to compute IP_ => P", GSL_EFAILED);
    }

    return GSL_SUCCESS;
}
