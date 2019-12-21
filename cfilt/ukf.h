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

#ifndef UKF_H_
#define UKF_H_

#include "cfilt/sigma.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    gsl_vector* x_;
    gsl_vector* x;

    gsl_matrix* P_;
    gsl_matrix* P;
    gsl_matrix* R;
    gsl_matrix* K;

    int (*F)(void*, void*);
    int (*H)(void*, void*);

    cfilt_sigma_generator* gen;

} cfilt_ukf;

int cfilt_ukf_alloc(cfilt_ukf* filt, const size_t n, const size_t m,
                    const size_t k, int (*F)(void*, void*),
                    int (*H)(void*, void*), cfilt_sigma_generator* gen);

void cfilt_ukf_free(cfilt_ukf* filt);

int cfilt_ukf_predict(cfilt_ukf* filt);

int cfilt_ukf_update(cfilt_ukf* filt);

#ifdef __cplusplus
}
#endif

#endif // UKF_
