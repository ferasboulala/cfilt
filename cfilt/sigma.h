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

#include <gsl/gsl_matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { CFILT_VAN_DER_MERWE = 0 } cfilt_sigma_generator_type;

typedef struct
{
    cfilt_sigma_generator_type type;
    gsl_matrix* points;
} cfilt_sigma_generator_common_;

typedef cfilt_sigma_generator_common_ cfilt_sigma_generator;

typedef struct
{
    cfilt_sigma_generator_common_ common_;
    double alpha;
    double beta;
    double kappa;
} cfilt_sigma_generator_van_der_merwe;

int cfilt_sigma_alloc(cfilt_sigma_generator_type, cfilt_sigma_generator** gen);

void cfilt_sigma_free(cfilt_sigma_generator* gen);

int cfilt_sigma_generate(cfilt_sigma_generator* gen, gsl_matrix** points);

#ifdef __cplusplus
}
#endif

#endif // CFILT_SIGMA_H_
