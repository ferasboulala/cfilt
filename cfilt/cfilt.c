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

#include "cfilt/cfilt.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

#include <string.h>

#define FREE_IF_NOT_NULL(m)                                                                                            \
    if (m)                                                                                                             \
        gsl_matrix_free(m);

int
cfilt_process_cov_discrete_white_noise(gsl_matrix* tau, const double sigma, gsl_matrix* Q)
{
    if (gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tau, tau, sigma, Q))
    {
        GSL_ERROR("failed to compute Tau Tau^T => Q", GSL_EFAILED);
    }

    return GSL_SUCCESS;
}

static void
cfilt_mahalanobis_free(gsl_matrix* x_copy, gsl_matrix* mu_copy, gsl_matrix* cov_inv, gsl_matrix* mahalanobis,
                       gsl_permutation* perm)
{
    FREE_IF_NOT_NULL(x_copy);
    FREE_IF_NOT_NULL(mu_copy);
    FREE_IF_NOT_NULL(cov_inv);
    FREE_IF_NOT_NULL(mahalanobis);
    if (perm)
    {
        gsl_permutation_free(perm);
    }
}

int
cfilt_mahalanobis(gsl_vector* x, gsl_vector* mu, gsl_matrix* cov, double* res)
{
    const int N = x->size;
    gsl_matrix* x_copy = gsl_matrix_alloc(N, 1);
    gsl_matrix* mu_copy = gsl_matrix_alloc(N, 1);
    gsl_matrix* cov_inv = gsl_matrix_alloc(N, N);
    gsl_matrix* mahalanobis = gsl_matrix_alloc(1, 1);
    gsl_permutation* perm = gsl_permutation_alloc(N);

    if (!x_copy || !mu_copy || !cov_inv || !mahalanobis || !perm)
    {
        cfilt_mahalanobis_free(x_copy, mu_copy, cov_inv, mahalanobis, perm);
        GSL_ERROR("failed to allocate memory for mahalanobis computation", GSL_ENOMEM);
    }

    memcpy(x_copy->data, x->data, N * sizeof(double));
    memcpy(mu_copy->data, mu->data, N * sizeof(double));

    // x - mu
    if (gsl_matrix_sub(x_copy, mu_copy))
    {
        GSL_ERROR("failed to compute x_copy - mu_copy => x_copy", GSL_EFAILED);
    }

    // Q^(-1)
    int signum;
    if (gsl_linalg_LU_decomp(cov, perm, &signum))
    {
        GSL_ERROR("failed to compute the LU decomposition of the covariance matrix", GSL_EFAILED);
    }

    if (gsl_linalg_LU_invert(cov, perm, cov_inv))
    {
        GSL_ERROR("failed to compute the inverse of the covariance matrix", GSL_EFAILED);
    }

    // (x - mu)Q^T^(-1)
    if (gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, x_copy, cov_inv, 0, mu_copy))
    {
        GSL_ERROR("failed to compute (x - mu)^T Q^-1 => mu_copy", GSL_EFAILED);
    }

    // (x - mu)Q^-1(x - mu)
    if (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mu_copy, x_copy, 0.0, mahalanobis))
    {
        GSL_ERROR("failed to compute the final step to mahalanobis distance", GSL_EFAILED);
    }

    *res = *(double*)mahalanobis->data;

    cfilt_mahalanobis_free(x_copy, mu_copy, cov_inv, mahalanobis, perm);

    return GSL_SUCCESS;
}
