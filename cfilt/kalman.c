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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <string.h>

#define FREE_IF_NO_NULLV(v) if (v) gsl_vector_free(v);
#define FREE_IF_NO_NULLM(m) if (m) gsl_matrix_free(m);


int cfilt_kalman_alloc(cfilt_kalman_filter *filt, const size_t n, const size_t m, const size_t k)
{
    memset(filt, 0, sizeof(cfilt_kalman_filter));

    filt->F = gsl_matrix_calloc(n, n);
    filt->B = gsl_matrix_calloc(n, m);
    filt->Q = gsl_matrix_calloc(n, n);
    filt->P = gsl_matrix_calloc(n, n);
    filt->H = gsl_matrix_calloc(k, n);
    filt->R = gsl_matrix_calloc(k, k);
    
    filt->x = gsl_vector_calloc(n);
    filt->x_ = gsl_vector_calloc(n);
    filt->z = gsl_vector_calloc(k);
    filt->u = gsl_vector_calloc(m);

    filt->_FP = gsl_matrix_alloc(n, n);
    filt->_P_H_T = gsl_matrix_alloc(n, k);

    if (!filt->F || !filt->B || !filt->Q || !filt->P || !filt->H || !filt->R || !filt->x || !filt->x_ || !filt->z || !filt->_FP || !filt->_P_H_T)
    {
        GSL_ERROR("failed to allocate space for kalman filter matrices", GSL_ENOMEM);
    }

    return GSL_SUCCESS;
}

void cfilt_kalman_free(cfilt_kalman_filter *filt)
{
    FREE_IF_NO_NULLM(filt->F);
    FREE_IF_NO_NULLM(filt->B);
    FREE_IF_NO_NULLM(filt->Q);
    FREE_IF_NO_NULLM(filt->P);
    FREE_IF_NO_NULLM(filt->H);
    FREE_IF_NO_NULLM(filt->R);

    FREE_IF_NO_NULLV(filt->x);
    FREE_IF_NO_NULLV(filt->x_);
    FREE_IF_NO_NULLV(filt->z);
    FREE_IF_NO_NULLV(filt->u);

    FREE_IF_NO_NULLM(filt->_FP);
    FREE_IF_NO_NULLM(filt->_P_H_T);
}

int cfilt_kalman_predict(cfilt_kalman_filter *filt)
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

    if (gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, filt->_FP, filt->F, 1.0, filt->Q))
    {
        GSL_ERROR("failed to compute _FPF^T + Q => Q", GSL_EFAILED);
    }

    if (gsl_matrix_memcpy(filt->P_, filt->Q))
    {
        GSL_ERROR("failed to copy Q into P", GSL_EFAILED);
    }

    return GSL_SUCCESS;
}

int cfilt_kalman_update(__attribute__((unused)) cfilt_kalman_filter *filt)
{
    // K = P_H^T(HP_H^T + R)^-1
    if (gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, filt->P_, filt->H, 0.0, filt->_P_H_T))
    {
        GSL_ERROR("failed to compute P_H^T + [0] => _P_H_T", GSL_EFAILED);
    }

    if (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, filt->H, filt->_P_H_T, 1.0, filt->R))
    {
        GSL_ERROR("failed to compute H _P_H_T + R => R", GSL_EFAILED);
    }

    if (gsl_blas_dgemm(

    // y = z - Hx_
    if (gsl_blas_dgemv(CblasNoTrans, -1.0, filt->H, filt->x_, 1.0, filt->z))
    {
        GSL_ERROR("failed to compute -Hx_ + z => z", GSL_EFAILED);
    }
    // x = x_ + Ky
    // P = (I - KH)P_
    
    return GSL_SUCCESS;
}
