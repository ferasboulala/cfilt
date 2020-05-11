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
#include "cfilt/sigma.h"
#include "cfilt/util.h"

#include <gsl/gsl_blas.h>

#include <string.h>

#define V_ALLOC_ASSERT_(p, n) V_ALLOC_ASSERT(p, n, cfilt_ukf_free, filt)
#define M_ALLOC_ASSERT_(p, n, m) M_ALLOC_ASSERT(p, n, m, cfilt_ukf_free, filt)

int
cfilt_ukf_sanity_check(const size_t n, const size_t m, const size_t k, cfilt_sigma_generator *gen)
{
    if (n * m * k == 0 || n == 1)
    {
        GSL_ERROR(
          "n, m and k must be positive integers. n must be greater than 1",
          GSL_EINVAL);
    }

    if (F == NULL || H == NULL)
    {
        GSL_ERROR("F and H must be non null pointers to a function",
                  GSL_EINVAL);
    }

    if (gen == NULL || gen->points->size2 != n)
    {
        GSL_ERROR(
          "Sigma generator dimensionality does not match filter's dimensions or NULL was given",
          GSL_EINVAL);
    }

    return GSL_SUCCESS;
}

int
cfilt_ukf_alloc(cfilt_ukf* filt, const size_t n, const size_t m, const size_t k,
                int (*F)(cfilt_ukf*, void*), int (*H)(cfilt_ukf*, void*), int (*X_MEAN)(cfilt_ukf*, void*),
                int (*Z_MEAN)(cfilt_ukf*, void*), int (*X_DIFF)(cfilt_ukf*, void*),
                int (*Y_DIFF)(cfilt_ukf*, void*), int (*Z_DIFF)(cfilt_ukf*, void*),
                cfilt_sigma_generator* gen)
{
    memset(filt, 0, sizeof(cfilt_ukf));

    EXEC_ASSERT(cfilt_ukf_realloc, filt, n, m, k, gen, 0);

    filt->F = F;
    filt->H = H;
    filt->X_MEAN = X_MEAN;
    filt->Z_MEAN = Z_MEAN;
    filt->X_DIFF = X_DIFF;

    return GSL_SUCCESS;
}

int
cfilt_ukf_matrix_realloc(cfilt_ukf* filt, gsl_matrix** a, const size_t n, const size_t m, const int keep_values)
{
    if (!filt->allocated_once)
    {
        M_ALLOC_ASSERT_(*a, n, m);
        return GSL_SUCCESS;
    }

    if (cfilt_matrix_realloc(a, n, m, keep_values))
    {
        cfilt_ukf_free(filt);
        return GSL_ENOMEM;
    }

    return GSL_SUCCESS;
}

int
cfilt_ukf_vector_realloc(cfilt_ukf* filt, gsl_vector** v, const size_t n, const int keep_values)
{
    if (!filt->allocated_once)
    {
        V_ALLOC_ASSERT_(*v, n);
        return GSL_SUCCESS;
    }

    if (cfilt_vector_realloc(v, n, keep_values))
    {
        cfilt_ukf_free(filt);
        return GSL_ENOMEM;
    }

    return GSL_SUCCESS;
}

int
cfilt_ukf_realloc(cfilt_ukf* filt, const size_t n, const size_t m, const size_t k, cfilt_sigma_generator* gen, const int keep_values)
{
    EXEC_ASSERT(cfilt_ukf_sanity_check, n, m, t, gen);

    filt->gen = gen;

    // predict step
    EXEC_ASSERT(cfilt_ukf_vector_realloc, filt, &filt->x, n, keep_values);

    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, &filt->P_, n, n, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, &filt->Y, gen->points->size1, n, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, &filt->Q_, n, n, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, &filt->_Y_x, 1, n, keep_values);

    // update step
    EXEC_ASSERT(cfilt_ukf_vector_realloc, filt, &filt->z, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_vector_realloc, filt, &filt->x, n, keep_values);
    EXEC_ASSERT(cfilt_ukf_vector_realloc, filt, &filt->u_z, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_vector_realloc, filt, &filt->y, k, keep_values);

    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->P, n, n, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->R, k, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->P_z, k, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->K, n, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->Z, gen->points->size1, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->_P_z_inv, n, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->_Z_u, 1, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->_Y_x_Z_u, n, k, keep_values);
    EXEC_ASSERT(cfilt_ukf_matrix_realloc, filt, filt->_K_P_z, n, k, keep_values);

    if (!filt->allocated_once)
    {
        filt->_perm = gsl_permutation_alloc(k);
        if (filt->_perm == NULL)
        {
            cfilt_ukf_free(filt);
            return GSL_ENOMEM;
        }
    }
    else
    {
        if (cfilt_permutation_realloc(&filt->_perm, k))
        {
            cfilt_ukf_free(filt);
            return GSL_ENOMEM;
        }
    }

    return GSL_SUCCESS;
}

void
cfilt_ukf_free(cfilt_ukf* filt)
{
    V_FREE_IF_NOT_NULL(filt->x_);
    V_FREE_IF_NOT_NULL(filt->x);
    V_FREE_IF_NOT_NULL(filt->u_z);
    V_FREE_IF_NOT_NULL(filt->y);
    V_FREE_IF_NOT_NULL(filt->z);

    M_FREE_IF_NOT_NULL(filt->P_);
    M_FREE_IF_NOT_NULL(filt->Y);
    M_FREE_IF_NOT_NULL(filt->Q);
    M_FREE_IF_NOT_NULL(filt->_Y_x);
    M_FREE_IF_NOT_NULL(filt->P);
    M_FREE_IF_NOT_NULL(filt->R);
    M_FREE_IF_NOT_NULL(filt->K);
    M_FREE_IF_NOT_NULL(filt->Z);
    M_FREE_IF_NOT_NULL(filt->P_z);
    M_FREE_IF_NOT_NULL(filt->_P_z_inv);
    M_FREE_IF_NOT_NULL(filt->_Z_u);
    M_FREE_IF_NOT_NULL(filt->_K_P_z);
    M_FREE_IF_NOT_NULL(filt->_Y_x_Z_u);

    if (filt->_perm != NULL)
    {
        gsl_permutation_free(filt->_perm);
    }
}

int
cfilt_ukf_predict(cfilt_ukf* filt, void* ptr)
{
    // Y = f(X)
    EXEC_ASSERT(cfilt_sigma_generator_generate, filt->gen, filt->x, filt->P);
    EXEC_ASSERT(filt->F, filt, ptr);

    // x_ = sum_i [mu_weight * Y_i]
    if (filt->X_MEAN == NULL)
    {
        gsl_vector_set_zero(filt->x_);
        for (size_t i = 0; i < filt->gen->points->size1; ++i)
        {
            const double weight = gsl_vector_get(filt->gen->mu_weights, i);
            gsl_vector_view row = gsl_matrix_row(filt->Y, i);
            gsl_vector* point = &row.vector;

            EXEC_ASSERT(gsl_vector_axpby, weight, point, 1.0, filt->x_);
        }
    }
    else
    {
        EXEC_ASSERT(filt->X_MEAN, filt->Y, filt->gen, filt->x_, ptr);
    }

    // P_ = sum_i [(Y - x_)(Y - x_)^T] + Q
    gsl_matrix_set_zero(filt->P_);
    for (size_t i = 0; i < filt->gen->points->size1; ++i)
    {
        const double weight = gsl_vector_get(filt->gen->sigma_weights, i);
        gsl_vector_view row = gsl_matrix_row(filt->Y, i);
        gsl_vector* point = &row.vector;

        gsl_vector_view dst_row = gsl_matrix_row(filt->_Y_x, 0);
        gsl_vector* dst = &dst_row.vector;

        EXEC_ASSERT(gsl_vector_memcpy, dst, filt->x_);

        if (filt->X_SUB != NULL)
        {
            EXEC_ASSERT(filt->X_SUB, filt->dst, filt->point, ptr);
        }
        else
        {
            EXEC_ASSERT(gsl_vector_sub, dst, point);
            EXEC_ASSERT(gsl_vector_scale, dst, -1);
        }

        EXEC_ASSERT(gsl_blas_dgemm, CblasTrans, CblasNoTrans, weight,
                    filt->_Y_x, filt->_Y_x, 1, filt->P_);
    }

    return GSL_SUCCESS;
}

int
cfilt_ukf_update(cfilt_ukf* filt, void* ptr)
{
    // Z = h(Y)
    EXEC_ASSERT(filt->H, filt, ptr);

    // u_z = sum_i [mu_weight * Z_i]
    if (filt->Z_MEAN == NULL)
    {
        gsl_vector_set_zero(filt->u_z);
        for (size_t i = 0; i < filt->gen->points->size1; ++i)
        {
            const double weight = gsl_vector_get(filt->gen->mu_weights, i);
            gsl_vector_view row = gsl_matrix_row(filt->Z, i);
            gsl_vector* point = &row.vector;

            EXEC_ASSERT(gsl_vector_axpby, weight, point, 1.0, filt->u_z);
        }
    }
    else
    {
        EXEC_ASSERT(filt->Z_MEAN, filt->Z, filt->gen, filt->u_z, ptr);
    }

    // y = z - u_z
    if (filt->Y_DIFF == NULL)
    {
        EXEC_ASSERT(gsl_vector_memcpy, filt->y, filt->z);
        EXEC_ASSERT(gsl_vector_sub, filt->y, filt->u_z);
    }
    else
    {
        EXEC_ASSERT(filt->Y_DIFF, filt->z, filt->u_z, filt->y, ptr);
    }

    // P_z = sum_i [sigma_weight * (z_i - u)(z_i - u)^T] + R
    gsl_matrix_set_zero(filt->P_z);
    for (size_t i = 0; i < filt->gen->points->size1; ++i)
    {
        const double weight = gsl_vector_get(filt->gen->sigma_weights, i);
        gsl_vector_view row = gsl_matrix_row(filt->Z, i);
        gsl_vector* point = &row.vector;

        gsl_vector_view dst_row = gsl_matrix_row(filt->_Z_u, 0);
        gsl_vector* dst = &dst_row.vector;

        EXEC_ASSERT(gsl_vector_memcpy, dst, filt->u_z);

        if (filt->Z_DIFF != NULL)
        {
            EXEC_ASSERT(filt->Z_DIFF, dst, point, ptr);
        }
        else
        {
            EXEC_ASSERT(gsl_vector_sub, dst, point);
            EXEC_ASSERT(gsl_vector_scale, dst, -1);
        }

        EXEC_ASSERT(gsl_blas_dgemm, CblasTrans, CblasNoTrans, weight,
                    filt->_Z_u, filt->_Z_u, 1.0, filt->P_z);
    }

    // K = sum_i [sigma_weight * (Y_i - x_)(Z_i - u_z)^T]P_z^(-1)
    gsl_matrix_set_zero(filt->_Y_x_Z_u);
    for (size_t i = 0; i < filt->gen->points->size1; ++i)
    {
        const double weight = gsl_vector_get(filt->gen->sigma_weights, i);

        gsl_vector_view row_z = gsl_matrix_row(filt->Z, i);
        gsl_vector* point_z = &row_z.vector;
        gsl_vector_view row_y = gsl_matrix_row(filt->Y, i);
        gsl_vector* point_y = &row_y.vector;

        gsl_vector_view dst_row_z = gsl_matrix_row(filt->_Z_u, 0);
        gsl_vector* dst_z = &dst_row_z.vector;
        gsl_vector_view dst_row_y = gsl_matrix_row(filt->_Y_x, 0);
        gsl_vector* dst_y = &dst_row_y.vector;

        EXEC_ASSERT(gsl_vector_memcpy, dst_z, filt->u_z);
        EXEC_ASSERT(gsl_vector_memcpy, dst_y, filt->x_);

        if (X_DIFF != NULL)
        {
            EXEC_ASSERT(filt->X_DIFF, dst_y, point_y, ptr);
        }
        else
        {
            EXEC_ASSERT(gsl_vector_sub, dst_y, point_y);
            EXEC_ASSERT(gsl_vector_scale, dst_y, -1);
        }

        if (Z_DIFF != NULL)
        {
            EXEC_ASSERT(filt->Z_DIFF, dst_z, point_z, ptr);
        }
        else
        {
            EXEC_ASSERT(gsl_vector_sub, dst_z, point_z);
            EXEC_ASSERT(gsl_vector_scale, dst_z, -1);
        }

        EXEC_ASSERT(gsl_blas_dgemm, CblasTrans, CblasNoTrans, weight,
                    filt->_Y_x, filt->_Z_u, 1.0, filt->_Y_x_Z_u);
    }

    EXEC_ASSERT(cfilt_matrix_invert, filt->P_z, filt->_P_z_inv, filt->_perm);
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasNoTrans, 1.0, filt->_Y_x_Z_u,
                filt->_P_z_inv, 0.0, filt->K);

    // x = x_ + Ky
    if (filt->X_UPDT == NULL)
    {
        EXEC_ASSERT(gsl_vector_memcpy, filt->x, filt->x_);
        EXEC_ASSERT(gsl_blas_dgemv, CblasNoTrans, 1.0, filt->K, filt->y, 1.0,
                    filt->x);
    }
    else
    {
        EXEC_ASSERT(filt->X_UPDT, filt->x_, filt->y, filt->x, ptr);
    }

    // P = P_ - KP_zK^T
    EXEC_ASSERT(gsl_matrix_memcpy, filt->P, filt->P_);
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasNoTrans, 1.0, filt->K,
                filt->P_z, 0.0, filt->_K_P_z);
    EXEC_ASSERT(gsl_blas_dgemm, CblasNoTrans, CblasTrans, -1.0, filt->_K_P_z,
                filt->K, 1.0, filt->P);

    return GSL_SUCCESS;
}
