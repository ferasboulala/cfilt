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
some_function(__unused__ cfilt_ukf* filt, __unused__ void* ptr)
{
    // Dummy function for the allocation unit test
    return GSL_SUCCESS;
}

int
test_ukf_alloc(void)
{
    cfilt_ukf filt;
    cfilt_sigma_generator* gen;
    UTEST_EXEC_ASSERT(cfilt_sigma_generator_alloc, CFILT_SIGMA_VAN_DER_MERWE,
                      &gen, 3, 0.5, 2.0, 1.0);

    gsl_error_handler_t* hdl = gsl_set_error_handler_off();
    UTEST_EXEC_ASSERT_(cfilt_ukf_alloc, &filt, 0, 3, 3, some_function,
                       some_function, gen);
    UTEST_EXEC_ASSERT_(cfilt_ukf_alloc, &filt, 3, 3, 3, NULL, some_function,
                       gen);
    UTEST_EXEC_ASSERT_(cfilt_ukf_alloc, &filt, 2, 3, 3, some_function,
                       some_function, gen);
    gsl_set_error_handler(hdl);

    UTEST_EXEC_ASSERT(cfilt_ukf_alloc, &filt, 3, 3, 3, some_function,
                      some_function, gen);
    cfilt_ukf_free(&filt);
    cfilt_sigma_generator_free(gen);

    return GSL_SUCCESS;
}

int
main(void)
{
    RUN_TEST(test_ukf_alloc);

    return GSL_SUCCESS;
}
