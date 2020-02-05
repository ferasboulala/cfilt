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
#include "utest.h"

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

int
test_discrete_white_noise(void)
{
    static const double DT = 0.1;
    static const double stddev = 1.0;

    gsl_matrix* Q = gsl_matrix_alloc(2, 2);
    gsl_vector* tau = gsl_vector_alloc(2);
    gsl_matrix* expected_result_matrix = gsl_matrix_alloc(2, 2);

    gsl_matrix_set(expected_result_matrix, 0, 0, pow(DT, 4) * 0.25);
    gsl_matrix_set(expected_result_matrix, 0, 1, pow(DT, 3) * 0.50);
    gsl_matrix_set(expected_result_matrix, 1, 0, pow(DT, 3) * 0.50);
    gsl_matrix_set(expected_result_matrix, 1, 1, pow(DT, 2));

    gsl_vector_set(tau, 0, 0.5 * DT * DT);
    gsl_vector_set(tau, 1, DT);

    UTEST_EXEC_ASSERT(cfilt_discrete_white_noise, tau, stddev, Q);
    // FIXME : Add tolerance
    // UTEST_EXEC_ASSERT(cfilt_matrix_cmp, expected_result_matrix, Q);

    gsl_matrix_free(Q);

    return GSL_SUCCESS;
}

int
test_mahalanobis(void)
{
    // With identity cov, it should be the euclidean distance
    gsl_vector* x = gsl_vector_alloc(3);
    gsl_vector* mu = gsl_vector_alloc(3);
    gsl_matrix* cov = gsl_matrix_alloc(3, 3);

    gsl_vector_set(x, 0, 1.0);
    gsl_vector_set(x, 1, 2.0);
    gsl_vector_set(x, 2, 3.0);
    gsl_vector_set_zero(mu);
    gsl_matrix_set_identity(cov);

    double distance;
    UTEST_EXEC_ASSERT(cfilt_mahalanobis, x, mu, cov, &distance);
    const double expected_distance = sqrt(1.0 + 4.0 + 9.0);
    UTEST_ASSERT_TOL(distance, expected_distance, 0.1,
                     "Expected value is not matched");

    gsl_vector_free(x);
    gsl_vector_free(mu);
    gsl_matrix_free(cov);

    return GSL_SUCCESS;
}

int
test_norm_estimated_error_squared(void)
{
    gsl_vector* x = gsl_vector_alloc(3);
    gsl_matrix* cov = gsl_matrix_alloc(3, 3);

    gsl_vector_set(x, 0, 1.0);
    gsl_vector_set(x, 1, 2.0);
    gsl_vector_set(x, 2, 3.0);
    gsl_matrix_set_identity(cov);

    double nees;
    UTEST_EXEC_ASSERT(cfilt_norm_estimated_error_squared, x, cov, &nees);
    UTEST_ASSERT_TOL(1.0 + 4.0 + 9.0, nees, 0.1,
                     "Expected value is not matched. Got %f, expected %f", nees,
                     14.0);

    gsl_vector_free(x);
    gsl_matrix_free(cov);

    return GSL_SUCCESS;
}

int
main(void)
{
    RUN_TEST(test_discrete_white_noise);
    RUN_TEST(test_mahalanobis);
    RUN_TEST(test_norm_estimated_error_squared);

    return GSL_SUCCESS;
}
