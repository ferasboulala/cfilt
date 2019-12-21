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

#ifndef KALMAN_H_
#define KALMAN_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * n            : Number of variables tracked
 * m            : Number of control inputs
 * k            : Number of measurement variables
 *
 * F (n x n)    : State transition matrix (or state transition function)
 * B (n x m)    : Control matrix
 * Q (n x n)    : Process covariance matrix (noise)
 * P (n x n)    : State covariance matrix
 * P_(n x n)    : State estimate covariance matrix
 * K (n x k)    : Kalman gain matrix
 *
 * H (k x n)    : Measurement matrix
 * R (k x k)    : Measurement covariance matrix (noise)
 *
 * x (n x 1)    : State vector
 * x_(n x 1)    : State estimate vector
 * z (k x 1)    : Measurement vector
 * u (m x 1)    : Control input vector
 * y (k x 1)    : Residual vector
 */

typedef struct
{
    gsl_vector* x;
    gsl_vector* x_;
    gsl_vector* z;
    gsl_vector* u;
    gsl_vector* y;

    gsl_matrix* F;
    gsl_matrix* B;
    gsl_matrix* Q;
    gsl_matrix* P;
    gsl_matrix* P_;
    gsl_matrix* H;
    gsl_matrix* R;
    gsl_matrix* K;

    // Intermediary results
    gsl_matrix* _FP;
    gsl_matrix* _PH_T;
    gsl_matrix* _PH_T_R;
    gsl_matrix* _inv;
    gsl_permutation* _perm;
    gsl_matrix* _I;

} cfilt_kalman_filter;

int cfilt_kalman_filter_alloc(cfilt_kalman_filter* filt, const size_t n, const size_t m, const size_t k);

void cfilt_kalman_filter_free(cfilt_kalman_filter* filt);

int cfilt_kalman_filter_predict(cfilt_kalman_filter* filt);

int cfilt_kalman_filter_update(cfilt_kalman_filter* filt);

#ifdef __cplusplus
}
#endif

#endif // KALMAN_H_
