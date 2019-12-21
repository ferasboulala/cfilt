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

#ifndef CFILT_H_
#define CFILT_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    double mean;
    double var;
} cfilt_gauss;

int cfilt_discrete_white_noise(gsl_matrix* tau, const double sigma,
                               gsl_matrix* Q);

int cfilt_mahalanobis(gsl_vector* x, gsl_vector* mu, gsl_matrix* cov,
                      double* res);

int cfilt_norm_estimated_error_squared(gsl_vector* x_, gsl_matrix* cov,
                                       double* res);

#ifdef __cplusplus
}
#endif

#endif // CFILT_H_
