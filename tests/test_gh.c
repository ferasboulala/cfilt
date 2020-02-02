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

#include "cfilt/gh.h"
#include "utest.h"

#include <gsl/gsl_errno.h>

int
test_cfilt_gh_alloc(void)
{
    // Nothing much to do here as it does only memory allocation.
    cfilt_gh_filter filt;
    gsl_error_handler_t* hdl = gsl_set_error_handler_off();
    UTEST_EXEC_ASSERT_(cfilt_gh_alloc, &filt, 0);
    gsl_set_error_handler(hdl);

    UTEST_EXEC_ASSERT(cfilt_gh_alloc, &filt, 3);

    cfilt_gh_free(&filt);

    return GSL_SUCCESS;
}

int
test_cfilt_gh_predict(void)
{
    // As for all filters, we do not test values. They are way
    // too hard to track so we leave that to examples. Correctness
    // will be asserted with the code just not crashing because
    // of segfaults and whatnot. This comment will not be repeated
    // and so I assume whoever cares enough to read my tests, they
    // read them in order of filter complexity (gh is the simplest).
    cfilt_gh_filter filt;

    cfilt_gh_alloc(&filt, 3);
    cfilt_gh_predict(&filt, 0.1);
    cfilt_gh_predict(&filt, 0.0);
    cfilt_gh_free(&filt);

    cfilt_gh_alloc(&filt, 1);
    cfilt_gh_predict(&filt, 0.1);
    cfilt_gh_free(&filt);

    return GSL_SUCCESS;
}

int
test_cfilt_gh_update(void)
{
    cfilt_gh_filter filt;

    cfilt_gh_alloc(&filt, 3);
    cfilt_gh_write(&filt, 0.0, 0.1);
    // With and without updates
    cfilt_gh_update(&filt, 0.1);
    cfilt_gh_update(&filt, 0.0);

    cfilt_gh_free(&filt);

    return GSL_SUCCESS;
}

int
main(void)
{
    RUN_TEST(test_cfilt_gh_alloc);
    RUN_TEST(test_cfilt_gh_predict);
    RUN_TEST(test_cfilt_gh_update);

    return GSL_SUCCESS;
}
