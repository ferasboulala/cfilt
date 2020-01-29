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
#include "utest.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_vector.h>

int
test_cfilt_sigma_generator_alloc_van_der_merwe(void)
{
    cfilt_sigma_generator* gen;
    UTEST_EXEC_ASSERT(cfilt_sigma_generator_alloc, CFILT_SIGMA_VAN_DER_MERWE,
                      &gen, 2, 0.5, 2.0, 1.0);
    UTEST_ASSERT(
      gen->type == CFILT_SIGMA_VAN_DER_MERWE,
      "The wrong type of generator was created. Asked for van der merwe");

    cfilt_sigma_generator_free(gen);

    return GSL_SUCCESS;
}

int
test_cfilt_sigma_generator_alloc(void)
{
    RUN_TEST(test_cfilt_sigma_generator_alloc_van_der_merwe);

    return GSL_SUCCESS;
}

int
test_cfilt_sigma_generator_generate_van_der_merwe(void)
{
    // Van Der MerweCA
    // Reasonable parameters are 0 < alpha < 1, beta = 2 and kappa = 3 - n
    const int n = 2;
    const double alpha = 0.5;
    const double beta = 2.0;
    const double kappa = 3 - n;
    const double lambda = alpha * alpha * (n + kappa) - n;

    gsl_matrix* cov = gsl_matrix_alloc(n, n);
    gsl_vector* mu = gsl_vector_alloc(n);

    gsl_matrix_set(cov, 0, 0, 1.0);
    gsl_matrix_set(cov, 1, 1, 2.0);
    gsl_vector_set(mu, 0, 0.0);
    gsl_vector_set(mu, 1, 1.0);

    // Expected results
    gsl_matrix* expected_points = gsl_matrix_alloc(n * n + 1, n);
    gsl_vector* expected_mu_weights = gsl_vector_alloc(n * n + 1);
    gsl_vector* expected_sigma_weights = gsl_vector_alloc(n * n + 1);

    gsl_vector_set_all(expected_sigma_weights, 1.0 / 2.0 / (n + lambda));
    gsl_vector_memcpy(expected_mu_weights, expected_sigma_weights);
    gsl_vector_set(expected_mu_weights, 0, lambda / (n + lambda));
    gsl_vector_set(expected_sigma_weights, 0,
                   lambda / (n + lambda) + 1 + beta - alpha * alpha);

    // First point is just mu
    gsl_vector_view row_view = gsl_matrix_row(expected_points, 0);
    gsl_vector* row = &row_view.vector;
    gsl_vector_memcpy(row, mu);
    // The rest are picked from the Cholesky decomposition matrix's rows

    cfilt_sigma_generator* gen;
    UTEST_EXEC_ASSERT(cfilt_sigma_generator_alloc, CFILT_SIGMA_VAN_DER_MERWE,
                      &gen, n, alpha, beta, kappa);
    UTEST_EXEC_ASSERT(cfilt_sigma_generator_generate, gen, mu, cov);

    // The first weight must be the largest
    double first_weight = gsl_vector_get(gen->mu_weights, 0);
    UTEST_ASSERT(first_weight == lambda / (n + lambda),
                 "First weight is not the largest for mu weights");
    first_weight = gsl_vector_get(gen->sigma_weights, 0);
    UTEST_ASSERT(first_weight ==
                   lambda / (n + lambda) + 1 + beta - pow(alpha, 2),
                 "First weight is not the largest for sigma weights");

    // We sort the weights to compare them directly as vectors
    gsl_sort_vector(expected_mu_weights);
    gsl_sort_vector(expected_sigma_weights);
    gsl_sort_vector(gen->mu_weights);
    gsl_sort_vector(gen->sigma_weights);

    UTEST_EXEC_ASSERT(cfilt_vector_cmp_tol, gen->mu_weights,
                      expected_mu_weights, 0.1);
    UTEST_EXEC_ASSERT(cfilt_vector_cmp_tol, gen->sigma_weights,
                      expected_sigma_weights, 0.1);
    row_view = gsl_matrix_row(gen->points, 0);
    UTEST_EXEC_ASSERT(cfilt_vector_cmp_tol, &row_view.vector, mu, 0.1);

    cfilt_sigma_generator_free(gen);
    gsl_matrix_free(cov);
    gsl_vector_free(mu);
    gsl_matrix_free(expected_points);
    gsl_vector_free(expected_mu_weights);
    gsl_vector_free(expected_sigma_weights);

    return GSL_SUCCESS;
}

int
test_cfilt_sigma_generator_generate(void)
{
    RUN_TEST(test_cfilt_sigma_generator_generate_van_der_merwe);

    return GSL_SUCCESS;
}

int
main(void)
{
    RUN_TEST(test_cfilt_sigma_generator_alloc);
    RUN_TEST(test_cfilt_sigma_generator_generate);

    return GSL_SUCCESS;
}
