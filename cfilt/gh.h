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

#ifndef GH_H_
#define GH_H_

/**
 * A simple implementation of the gh filter (also known as the alpha beta
 * filter).
 * It tracks the values of an arbitrary amount of state variables that are
 * related as such:
 *      x_(i+1) = d(x_i)/dt
 *
 * There are assumptions made with this filters. First, when tracking an set of
 * state variables,
 * all intermediary variables will be tracked. For example, if we wish to track
 * position and
 * acceleration, velocity will be tracked too. This is enforced by requiring a
 * dimension for
 * the filter (dim := difference between the highest and lowest order variables
 * + 1).
 *
 * When predicting the value of a state variable, we use all the higher order
 * variables to do so.
 * That is, x_pred_i = x_i + dt * x_(i+1) + dt^2 * x_(i+2) * 1/2 + ... +
 * dt^(n-i-1) * x_n * 1/(n-i-1)!.
 * For x_n, it is assumed that the value is constant only during the prediction.
 *
 * When updating the value of a state variable, the value of the residual is
 * computed as:
 *      1) The variable received data: Residual = x_pred_i - _zi
 *      2) If not:
 *          a) The immediate lower order variable received data: Residual =
 * x_pred_i - (x_pred_(i-1) - _z(i - 1)) / dt
 *          b) If not: Residual = 0.
 * In all cases, the update equation is x_i = x_i - gh_i * residual
 *
 * Runtime of functions with n being the filter dimensions:
 * init : Theta(n)
 * free : Theta(n)
 * write : O(1)
 * predict : Theta(n^2)
 * update : Theta(n)
 */

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    double* gh;
    double* x;
    double* x_pred;
    size_t dim;

    char* _upd;
    double* _z;
    void* _ptr;

} cfilt_gh_filter;

int cfilt_gh_alloc(cfilt_gh_filter* filt, const size_t dim);

void cfilt_gh_free(cfilt_gh_filter* filt);

void cfilt_gh_write(cfilt_gh_filter* filt, const double val, const size_t ord);

void cfilt_gh_predict(cfilt_gh_filter* filt, const double dt);

void cfilt_gh_update(cfilt_gh_filter* filt, const double dt);

#ifdef __cplusplus
}
#endif

#endif // GH_H_
