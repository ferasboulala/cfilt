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

#ifndef CFILT_SIGMA_H_
#define CFILT_SIGMA_H_

#include <sys/types.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { CFILT_SIGMA_VAN_DER_MERWE = 0 } cfilt_sigma_generator_type;

typedef struct
{
    cfilt_sigma_generator_type type;
    gsl_matrix* points;
    gsl_vector* mu_weights;
    gsl_vector* sigma_weights;
    size_t n;
} cfilt_sigma_generator_common_;

typedef cfilt_sigma_generator_common_ cfilt_sigma_generator;

int cfilt_sigma_generator_alloc(cfilt_sigma_generator_type type,
                                cfilt_sigma_generator** gen, const size_t n,
                                ...);

void cfilt_sigma_generator_free(cfilt_sigma_generator* gen);

int cfilt_sigma_generator_generate(cfilt_sigma_generator* gen,
                                   const gsl_vector* mu, const gsl_matrix* cov);

#ifdef __cplusplus
}
#endif

#endif // CFILT_SIGMA_H_
