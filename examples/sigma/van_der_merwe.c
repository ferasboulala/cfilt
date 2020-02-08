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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <stdio.h>
#include <sys/types.h>

void
print_usage(char* prog_name)
{
    fprintf(stderr, "Usage : %s alpha beta kappa\n", prog_name);
}

int
main(int argc, char** argv)
{
    if (argc != 4)
    {
        print_usage(argv[0]);
        return -1;
    }

    const double alpha = atof(argv[1]);
    const double beta = atof(argv[2]);
    const double kappa = atof(argv[3]);
    const size_t N = 2;

    cfilt_sigma_generator* generator;
    if (cfilt_sigma_generator_alloc(CFILT_SIGMA_VAN_DER_MERWE, &generator, N,
                                    alpha, beta, kappa))
    {
        fprintf(
          stderr,
          "Could not allocate memory for the Van Der Merwe sigma generator\n");
        return -1;
    }

    gsl_matrix* cov = gsl_matrix_alloc(N, N);
    gsl_matrix_set_identity(cov);
    gsl_matrix_scale(cov, 3.0);
    gsl_matrix_set(cov, 0, 1, 1.0);
    gsl_matrix_set(cov, 1, 0, 1.0);

    gsl_vector* mu = gsl_vector_alloc(N);
    gsl_vector_set_zero(mu);

    if (cfilt_sigma_generator_generate(generator, mu, cov))
    {
        fprintf(stderr, "Could not generate sigma points\n");
        return -1;
    }

    printf("points\n");
    cfilt_printf_matrix_rows(generator->points);

    printf("mu_weights\n");
    cfilt_printf_vector_row(generator->mu_weights);

    printf("sigma_weights\n");
    cfilt_printf_vector_row(generator->sigma_weights);

    printf("mu\n");
    cfilt_printf_vector_row(mu);

    printf("cov\n");
    cfilt_printf_matrix_rows(cov);

    cfilt_sigma_generator_free(generator);
    gsl_matrix_free(cov);
    gsl_vector_free(mu);

    return 0;
}
