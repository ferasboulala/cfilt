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

#include "cfilt/sigma.h"
#include "cfilt/util.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#define M_ALLOC_ASSERT_VDM(p, n, m)                                            \
    M_ALLOC_ASSERT(p, n, m, cfilt_sigma_generator_van_der_merwe_free, *gen)
#define V_ALLOC_ASSERT_VDM(v, n)                                               \
    V_ALLOC_ASSERT(v, n, cfilt_sigma_generator_van_der_merwe_free, *gen)
#define VDM(n) (2 * (n) + 1)

typedef struct
{
    cfilt_sigma_generator_common_ _common;
    double alpha;
    double beta;
    double kappa;
    double lambda;

    gsl_matrix* _chol;
} cfilt_sigma_generator_van_der_merwe;

static void
cfilt_sigma_generator_van_der_merwe_free(
  cfilt_sigma_generator_van_der_merwe* gen)
{
    M_FREE_IF_NOT_NULL(gen->_chol);
}

static int
cfilt_sigma_generator_van_der_merwe_alloc(
  cfilt_sigma_generator_van_der_merwe** gen, const size_t n, const double alpha,
  const double beta, const double kappa)
{
    *gen = malloc(sizeof(cfilt_sigma_generator_van_der_merwe));
    if (*gen == NULL)
    {
        return GSL_ENOMEM;
    }

    cfilt_sigma_generator_van_der_merwe vdm;
    vdm._common.type = CFILT_SIGMA_VAN_DER_MERWE;
    vdm._common.n = n;

    vdm.alpha = alpha;
    vdm.beta = beta;
    vdm.kappa = kappa;
    vdm.lambda = pow(alpha, 2) * (n + kappa) - n;

    M_ALLOC_ASSERT_VDM(vdm._common.points, VDM(n), n);
    M_ALLOC_ASSERT_VDM(vdm._chol, n, n);

    V_ALLOC_ASSERT_VDM(vdm._common.mu_weights, VDM(n));
    V_ALLOC_ASSERT_VDM(vdm._common.sigma_weights, VDM(n));

    const double weight = 1.0 / (2.0 * (n + vdm.lambda));
    gsl_vector_set_all(vdm._common.mu_weights, weight);
    gsl_vector_set_all(vdm._common.sigma_weights, weight);

    gsl_vector_set(vdm._common.mu_weights, 0, vdm.lambda / (vdm.lambda + n));
    gsl_vector_set(vdm._common.sigma_weights, 0,
                   gsl_vector_get(vdm._common.mu_weights, 0) + 1 -
                     pow(alpha, 2) + beta);

    memcpy(*gen, &vdm, sizeof(cfilt_sigma_generator_van_der_merwe));

    return GSL_SUCCESS;
}

static int
cfilt_sigma_generator_van_der_merwe_generate(
  cfilt_sigma_generator_van_der_merwe* gen, const gsl_vector* mu,
  const gsl_matrix* cov)
{
    // X_0
    gsl_vector_view first_row = gsl_matrix_row(gen->_common.points, 0);
    gsl_vector* first_point = &first_row.vector;
    gsl_vector_memcpy(first_point, mu);

    // sqrt( (n + lambda) * cov )
    EXEC_ASSERT(gsl_matrix_memcpy, gen->_chol, cov);
    gsl_matrix_scale(gen->_chol, gen->_common.n + gen->lambda);
    EXEC_ASSERT(gsl_linalg_cholesky_decomp1, gen->_chol);

    // gsl will return a lower triangular matrix and stores the original
    // matrix in the upper matrix. We must zero out that upper matrix
    // but whether or not the result is upper or lower is irrelevant for
    // as long as we pick either columns or rows, respectively.
    EXEC_ASSERT(cfilt_matrix_tri_zero, gen->_chol, 1);

    // mu +/- variance
    for (size_t i = 0; i < gen->_common.n; ++i)
    {
        // The following operations work because matrices are row major
        // And we are working with rows because the cholesky decomposition
        // returned a lower triangular matrix
        gsl_vector_view row = gsl_matrix_row(gen->_chol, i);
        gsl_vector* src = &row.vector;

        gsl_vector_view row1 = gsl_matrix_row(gen->_common.points, i + 1);
        gsl_vector* dst1 = &row1.vector;

        // +
        gsl_vector_memcpy(dst1, src);
        gsl_vector_add(dst1, mu);

        gsl_vector_view row2 =
          gsl_matrix_row(gen->_common.points, i + 1 + gen->_common.n);
        gsl_vector* dst2 = &row2.vector;

        // -
        gsl_vector_memcpy(dst2, src);
        gsl_vector_sub(dst2, mu);
        gsl_vector_scale(dst2, -1);
    }

    return GSL_SUCCESS;
}

int
cfilt_sigma_generator_alloc(const cfilt_sigma_generator_type type,
                            cfilt_sigma_generator** gen, const size_t n, ...)
{
    if (n == 0)
    {
        GSL_ERROR("cannot initialize a generator of dimension 0", GSL_EINVAL);
    }

    va_list valist;
    switch (type)
    {
        case CFILT_SIGMA_VAN_DER_MERWE:
            va_start(valist, n);

            const double alpha = va_arg(valist, double);
            const double beta = va_arg(valist, double);
            const double kappa = va_arg(valist, double);

            va_end(valist);

            EXEC_ASSERT(cfilt_sigma_generator_van_der_merwe_alloc,
                        (cfilt_sigma_generator_van_der_merwe**)gen, n, alpha,
                        beta, kappa);
            break;
        default:
            GSL_ERROR("Invalid sigma generator type", GSL_EINVAL);
    }

    return GSL_SUCCESS;
}

void
cfilt_sigma_generator_free(cfilt_sigma_generator* gen)
{
    M_FREE_IF_NOT_NULL(gen->points);

    V_FREE_IF_NOT_NULL(gen->mu_weights);
    V_FREE_IF_NOT_NULL(gen->sigma_weights);

    switch (gen->type)
    {
        case CFILT_SIGMA_VAN_DER_MERWE:
            cfilt_sigma_generator_van_der_merwe_free(
              (cfilt_sigma_generator_van_der_merwe*)gen);
            break;
    }

    free(gen);
}

int
cfilt_sigma_generator_generate(cfilt_sigma_generator* gen, const gsl_vector* mu,
                               const gsl_matrix* cov)
{
    switch (gen->type)
    {
        case CFILT_SIGMA_VAN_DER_MERWE:
            EXEC_ASSERT(cfilt_sigma_generator_van_der_merwe_generate,
                        (cfilt_sigma_generator_van_der_merwe*)gen, mu, cov);
            break;
        default:
            GSL_ERROR("Could not recognize sigma generator type", GSL_EINVAL);
    }

    return GSL_SUCCESS;
}
