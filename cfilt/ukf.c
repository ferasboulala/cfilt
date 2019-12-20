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

#include "cfilt/ukf.h"
#include "cfilt/common.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <string.h>
#include <sys/types.h>

#define V_ALLOC_ASSERT_(p, n) V_ALLOC_ASSERT(p, n, cfilt_ukf_free, filt)
#define M_ALLOC_ASSERT_(p, n, m) M_ALLOC_ASSERT(p, n, m, cfilt_ukf_free, filt)

int
cfilt_ukf_alloc(cfilt_ukf* filt, const size_t n, const size_t m, const size_t k, int (*F)(void*, void*),
                int (*H)(void*, void*), cfilt_sigma_generator* gen)
{
    memset(filt, 0, sizeof(cfilt_ukf));

    V_ALLOC_ASSERT_(filt->x_, n);
    V_ALLOC_ASSERT_(filt->x, n);

    M_ALLOC_ASSERT_(filt->P_, n, n);
    M_ALLOC_ASSERT_(filt->P, n, n);
    M_ALLOC_ASSERT_(filt->R, k, k);
    M_ALLOC_ASSERT_(filt->K, n, k);

    filt->F = F;
    filt->H = H;
    filt->gen = gen;

    return GSL_SUCCESS;
}

void
cfilt_ukf_free(cfilt_ukf* filt)
{
    V_FREE_IF_NOT_NULL(filt->x_);
    V_FREE_IF_NOT_NULL(filt->x);

    M_FREE_IF_NOT_NULL(filt->P_);
    M_FREE_IF_NOT_NULL(filt->P);
    M_FREE_IF_NOT_NULL(filt->R);
    M_FREE_IF_NOT_NULL(filt->K);
}

int
cfilt_ukf_predict(cfilt_ukf* filt)
{
    return GSL_SUCCESS;
}

int
cfilt_ukf_update(cfilt_ukf* filt)
{
    return GSL_SUCCESS;
}
