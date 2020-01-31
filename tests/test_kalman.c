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

#include "cfilt/kalman.h"
#include "utest.h"

#include <gsl/gsl_errno.h>

int
test_kalman_filter_alloc(void)
{
    cfilt_kalman_filter filt;

    gsl_error_handler_t* hdl = gsl_set_error_handler_off();
    UTEST_EXEC_ASSERT_(cfilt_kalman_filter_alloc, &filt, 0, 3, 3);
    gsl_set_error_handler(hdl);

    cfilt_kalman_filter_alloc(&filt, 3, 3, 3);
    cfilt_kalman_filter_free(&filt);

    return GSL_SUCCESS;
}

int
test_kalman_filter_predict(void)
{
    cfilt_kalman_filter filt;
    UTEST_EXEC_ASSERT(cfilt_kalman_filter_alloc, &filt, 3, 3, 3);

    // No matrix inversions in the prediction phase which means the values
    // in the matrices don't need to be set

    UTEST_EXEC_ASSERT(cfilt_kalman_filter_predict, &filt);
    UTEST_EXEC_ASSERT(cfilt_kalman_filter_predict, &filt);

    cfilt_kalman_filter_free(&filt);

    return GSL_SUCCESS;
}

int
test_kalman_filter_update(void)
{
    cfilt_kalman_filter filt;
    UTEST_EXEC_ASSERT(cfilt_kalman_filter_alloc, &filt, 3, 3, 3);

    // The actual values are not tested. We leave that for examples as
    // unit tests would take way too long to write and would serve
    // no purpose considering all those filters do is call gsl matrix
    // operations functions.
    gsl_matrix_set_identity(filt.H);
    gsl_matrix_set_identity(filt.P_);
    gsl_matrix_set_identity(filt.R);

    UTEST_EXEC_ASSERT(cfilt_kalman_filter_update, &filt);
    UTEST_EXEC_ASSERT(cfilt_kalman_filter_update, &filt);
    UTEST_EXEC_ASSERT(cfilt_kalman_filter_update, &filt);
    UTEST_EXEC_ASSERT(cfilt_kalman_filter_update, &filt);
    UTEST_EXEC_ASSERT(cfilt_kalman_filter_update, &filt);

    cfilt_kalman_filter_free(&filt);

    return GSL_SUCCESS;
}

int
main(void)
{
    RUN_TEST(test_kalman_filter_alloc);
    RUN_TEST(test_kalman_filter_predict);
    RUN_TEST(test_kalman_filter_update);

    return GSL_SUCCESS;
}
