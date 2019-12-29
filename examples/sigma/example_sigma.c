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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <stdio.h>

#define N 2
#define ALPHA 0.5
#define BETA 2.0
#define KAPPA 3.0 - N

/**
 * This test emulates an entity moving in a straight line, in 2D. Its sensors
 * yield
 * position and velocity.
 */
int
main(int argc, char** argv)
{
    cfilt_sigma_generator* generator;
    if (cfilt_sigma_generator_alloc(CFILT_SIGMA_VAN_DER_MERWE, &generator, N,
                                    ALPHA, BETA, KAPPA))
    {
        fprintf(
          stderr,
          "Could not allocate memory for the Van Der Merwe sigma generator\n");
        return -1;
    }

    gsl_matrix* cov = gsl_matrix_alloc(N, N);
    gsl_matrix_set_identity(cov);
    gsl_vector* mu = gsl_vector_alloc(N);
    gsl_vector_set_all(mu, 1);

    if (cfilt_sigma_generator_generate(generator, mu, cov))
    {
        fprintf(stderr, "Could not generate sigma points\n");
        goto cleanup;
    }

    fprintf(stderr, "mu_weights=\n");
    gsl_vector_fprintf(stderr, generator->mu_weights, "%f");

    fprintf(stderr, "\nsigma_weights=\n");
    gsl_vector_fprintf(stderr, generator->sigma_weights, "%f");

    fprintf(stderr, "\npoints=\n");
    gsl_matrix_fprintf(stderr, generator->points, "%f");

cleanup:
    cfilt_sigma_generator_free(generator);
    gsl_matrix_free(cov);
    gsl_vector_free(mu);

    return 0;
}
