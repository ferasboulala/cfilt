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

int
main(void)
{
    double x[] = { 1.0, 2.0, 3.0 };
    const int N = sizeof(x) / sizeof(double);

    gsl_vector_view x_view = gsl_vector_view_array(x, N);
    gsl_vector* x_vec = &x_view.vector;

    double mu[] = { 3.0, 2.0, 1.0 };

    gsl_vector_view mu_view = gsl_vector_view_array(mu, N);
    gsl_vector* mu_vec = &mu_view.vector;

    gsl_matrix* cov = gsl_matrix_alloc(N, N);
    gsl_matrix_set_identity(cov);
    gsl_matrix_scale(cov, 2.0);

    double distance;

    if (cfilt_mahalanobis(x_vec, mu_vec, cov, &distance))
    {
        fprintf(stderr,
                "An error occured while computing the mahalanobis distance\n");
        return -1;
    }

    printf("The result is %f\n", distance);

    gsl_matrix_free(cov);

    return 0;
}
