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

#include <math.h>
#include <stdio.h>

#define DT 0.1

int
main(void)
{
    double tau[] = { 0.5 * pow(DT, 2), DT };
    const int N = sizeof(tau) / sizeof(double);

    gsl_vector_view tau_view = gsl_vector_view_array(tau, N);
    gsl_vector* tau_vec = &tau_view.vector;

    gsl_matrix* Q = gsl_matrix_alloc(N, N);
    double variance = 1.0;

    if (cfilt_discrete_white_noise(tau_vec, variance, Q) != GSL_SUCCESS)
    {
        fprintf(stderr, "An error occured while computing the discrete white "
                        "noise covariance matrix\n");
        return -1;
    }

    gsl_matrix_fprintf(stdout, Q, "%f");

    gsl_matrix_free(Q);

    return 0;
}
