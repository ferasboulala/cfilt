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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <stdio.h>

#define DT 0.1

int
main(void)
{
    // Testing the discrete white noise covariance matrix generator
    
    gsl_matrix *tau = gsl_matrix_alloc(2, 1);
    gsl_matrix_set(tau, 0, 0, 1.0 / 2 * DT * DT);
    gsl_matrix_set(tau, 1, 0, DT);

    const double sigma = 5;

    gsl_matrix *Q = gsl_matrix_alloc(2, 2);

    if (cfilt_process_cov_discrete_white_noise(tau, sigma, Q))
    {
        fprintf(stderr, "Failed to compute the covariance matrix\n");
        goto cleanup;
    }

    if (gsl_matrix_fprintf(stderr, Q, "%f"))
    {
        fprintf(stderr, "Failed to print out the covariance matrix\n");
        goto cleanup;
    }

    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, 1);
    gsl_vector_set(x, 1, 2);

    gsl_vector *mu = gsl_vector_alloc(2);
    gsl_vector_set_zero(mu);

    gsl_matrix *cov = gsl_matrix_alloc(2, 2);
    gsl_matrix_set_identity(cov);

    double res;

    if (cfilt_mahalanobis(x, mu, cov, &res))
    {
        fprintf(stderr, "Failed to compute the mahalanobis distance\n");
    }
    else
    {
        fprintf(stderr, "mahalanobis distance: %f\n", res);
    }

    gsl_vector_free(x);
    gsl_vector_free(mu);
    gsl_matrix_free(cov);

cleanup:
    gsl_matrix_free(tau);
    gsl_matrix_free(Q);

    return 0;
}
