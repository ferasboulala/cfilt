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
#include "cfilt/util.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

#include <math.h>
#include <string.h>

int
cfilt_discrete_white_noise(gsl_matrix* tau, const double sigma, gsl_matrix* Q)
{
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasTrans, 1.0, tau, tau, 0, Q);
    gsl_matrix_scale(Q, sigma);

    return GSL_SUCCESS;
}

static void
cfilt_mahalanobis_free(gsl_matrix* x_copy, gsl_matrix* mu_copy,
                       gsl_matrix* xmuq, gsl_matrix* cov_inv,
                       gsl_matrix* mahalanobis, gsl_permutation* perm)
{
    M_FREE_IF_NOT_NULL(x_copy);
    M_FREE_IF_NOT_NULL(mu_copy);
    M_FREE_IF_NOT_NULL(cov_inv);
    M_FREE_IF_NOT_NULL(xmuq);
    M_FREE_IF_NOT_NULL(mahalanobis);

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
    gsl_matrix* xmuq = gsl_matrix_alloc(1, N);
    gsl_matrix* cov_inv = gsl_matrix_alloc(N, N);
    gsl_matrix* mahalanobis = gsl_matrix_alloc(1, 1);
    gsl_permutation* perm = gsl_permutation_alloc(N);

    if (!x_copy || !mu_copy || !cov_inv || !mahalanobis || !perm)
    {
        cfilt_mahalanobis_free(x_copy, mu_copy, xmuq, cov_inv, mahalanobis,
                               perm);

        return GSL_ENOMEM;
    }

    // TODO : Make this a util function
    gsl_vector_view x_copy_view = gsl_matrix_column(x_copy, 0);
    gsl_vector *x_copy_vector = &x_copy_view.vector;
    EXEC_ASSERT(gsl_vector_memcpy, x_copy_vector, x);

    gsl_vector_view mu_copy_view = gsl_matrix_column(mu_copy, 0);
    gsl_vector *mu_copy_vector = &mu_copy_view.vector;
    EXEC_ASSERT(gsl_vector_memcpy, mu_copy_vector, mu);

    // x - mu
    EXEC_ASSERT(gsl_matrix_sub, x_copy, mu_copy);

    // Q^(-1)
    EXEC_ASSERT(cfilt_matrix_invert, cov, cov_inv, perm);

    // (x - mu)^TQ^(-1)
    EXEC_ASSERT(gsl_blas_dgemm, CblasTrans, CblasNoTrans, 1.0, x_copy, cov_inv,
                0.0, xmuq);

    // (x - mu)Q^-1(x - mu)
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasNoTrans, 1.0, xmuq, x_copy,
                0.0, mahalanobis);

    *res = sqrt(*(double*)mahalanobis->data);

    cfilt_mahalanobis_free(x_copy, mu_copy, xmuq, cov_inv, mahalanobis, perm);

    return GSL_SUCCESS;
}

int
cfilt_norm_estimated_error_squared(gsl_vector* x_, gsl_matrix* cov, double* res)
{
    // Here, x_ = x - x_estimation
    // We do not compute the difference here because it would require an
    // intermediary result which requires memory allocation.
    // Unlike the previous function, this one is simple enough not to do it.
    gsl_vector* zero = gsl_vector_alloc(x_->size);
    if (zero == NULL)
    {
        return GSL_ENOMEM;
    }

    gsl_vector_set_zero(zero);

    EXEC_ASSERT(cfilt_mahalanobis, x_, zero, cov, res);

    *res *= *res;

    gsl_vector_free(zero);

    return GSL_SUCCESS;
}
