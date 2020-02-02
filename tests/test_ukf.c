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

#include "cfilt/sigma.h"
#include "cfilt/ukf.h"
#include "utest.h"

#include <gsl/gsl_errno.h>

#define __unused__ __attribute__((unused))

int
no_op_func(__unused__ cfilt_ukf* filt, __unused__ void* ptr)
{
    return GSL_SUCCESS;
}

int
test_cfilt_ukf_alloc(void)
{
    cfilt_ukf filt;
    cfilt_sigma_generator* gen;
    UTEST_EXEC_ASSERT(cfilt_sigma_generator_alloc, CFILT_SIGMA_VAN_DER_MERWE,
                      &gen, 3, 0.5, 2.0, 1.0);

    gsl_error_handler_t* hdl = gsl_set_error_handler_off();
    UTEST_EXEC_ASSERT_(cfilt_ukf_alloc, &filt, 0, 3, 3, no_op_func, no_op_func,
                       gen);
    UTEST_EXEC_ASSERT_(cfilt_ukf_alloc, &filt, 3, 3, 3, NULL, no_op_func, gen);
    UTEST_EXEC_ASSERT_(cfilt_ukf_alloc, &filt, 2, 3, 3, no_op_func, no_op_func,
                       gen);
    gsl_set_error_handler(hdl);

    UTEST_EXEC_ASSERT(cfilt_ukf_alloc, &filt, 3, 3, 3, no_op_func, no_op_func,
                      gen);
    cfilt_ukf_free(&filt);
    cfilt_sigma_generator_free(gen);

    return GSL_SUCCESS;
}

int
test_cfilt_ukf_predict(void)
{
    cfilt_ukf filt;
    cfilt_sigma_generator* gen;
    UTEST_EXEC_ASSERT(cfilt_sigma_generator_alloc, CFILT_SIGMA_VAN_DER_MERWE,
                      &gen, 3, 0.5, 2.0, 1.0);
    UTEST_EXEC_ASSERT(cfilt_ukf_alloc, &filt, 3, 3, 3, no_op_func, no_op_func,
                      gen);

    // The only matrix that must be set for the sigma generator callback
    gsl_matrix_set_identity(filt.P);

    UTEST_EXEC_ASSERT(cfilt_ukf_predict, &filt, NULL);
    UTEST_EXEC_ASSERT(cfilt_ukf_predict, &filt, NULL);
    UTEST_EXEC_ASSERT(cfilt_ukf_predict, &filt, NULL);

    cfilt_ukf_free(&filt);
    cfilt_sigma_generator_free(gen);

    return GSL_SUCCESS;
}

int
test_cfilt_ukf_update(void)
{
    cfilt_ukf filt;
    cfilt_sigma_generator* gen;
    UTEST_EXEC_ASSERT(cfilt_sigma_generator_alloc, CFILT_SIGMA_VAN_DER_MERWE,
                      &gen, 3, 0.5, 2.0, 1.0);
    UTEST_EXEC_ASSERT(cfilt_ukf_alloc, &filt, 3, 3, 3, no_op_func, no_op_func,
                      gen);

    gsl_matrix_set_identity(filt.P);
    UTEST_EXEC_ASSERT(cfilt_ukf_predict, &filt, NULL);
    UTEST_EXEC_ASSERT(cfilt_ukf_update, &filt, NULL);
    gsl_matrix_set_identity(filt.P);
    UTEST_EXEC_ASSERT(cfilt_ukf_predict, &filt, NULL);
    gsl_matrix_set_identity(filt.P);
    UTEST_EXEC_ASSERT(cfilt_ukf_predict, &filt, NULL);
    gsl_matrix_set_identity(filt.P);
    UTEST_EXEC_ASSERT(cfilt_ukf_predict, &filt, NULL);
    gsl_matrix_set_identity(filt.P);
    UTEST_EXEC_ASSERT(cfilt_ukf_update, &filt, NULL);

    cfilt_ukf_free(&filt);
    cfilt_sigma_generator_free(gen);

    return GSL_SUCCESS;
}

int
main(void)
{
    RUN_TEST(test_cfilt_ukf_alloc);
    RUN_TEST(test_cfilt_ukf_predict);
    RUN_TEST(test_cfilt_ukf_update);

    return GSL_SUCCESS;
}
