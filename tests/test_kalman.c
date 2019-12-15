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

#include <gsl/gsl_randist.h>

#include <stdio.h>
#include <time.h>

#define N_STEPS 1000
#define DT 0.1
#define V_X 1.0
#define V_Y 2.0
#define V_X_NOISE 0.1 * V_X
#define V_Y_NOISE 0.1 * V_Y
#define Q_VAR 0.1

/**
 * This test emulates an entity moving in a straight line, in 2D. Its sensors yield
 * position and velocity.
 */
int
main(void)
{
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));

    cfilt_kalman_filter filt;
    if (cfilt_kalman_alloc(&filt, 4, 1, 2))
    {
        fprintf(stderr, "Could not allocate kalman filter memory\n");
        goto cleanup;
    }

    // Setting up the state transition matrix F
    // [1 0 0 0
    //  0 1 0 0
    //  0 0 1 0
    //  0 0 0 1]
    gsl_matrix_set_identity(filt.F);

    // Setting up the control matrix B
    // [0
    //  V_X
    //  0
    //  V_Y]
    gsl_matrix_set_zero(filt.B);
    gsl_matrix_set(filt.B, 1, 0, V_X);
    gsl_matrix_set(filt.B, 3, 0, V_Y);

    // Setting up the control input vector u
    // [dt]
    gsl_vector_set(filt.u, 0, DT);

    // Setting up the process covariance matrix Q
    gsl_matrix_set_identity(filt.Q);
    gsl_matrix_scale(filt.Q, Q_VAR);

    // Setting up the measurement matrix H
    // [0 1 0 0
    //  0 0 0 1]
    gsl_matrix_set_zero(filt.H);
    gsl_matrix_set(filt.H, 0, 1, 1);
    gsl_matrix_set(filt.H, 1, 3, 1);

    // Setting up the measurement covariance matrix R
    // [V_X_NOISE 0
    //     0  V_Y_NOISE]
    gsl_matrix_set_zero(filt.R);
    gsl_matrix_set(filt.R, 0, 0, V_X_NOISE);
    gsl_matrix_set(filt.R, 1, 1, V_Y_NOISE);

    // Initializing the state vector x
    // [0 V_X 0 V_Y]^T
    gsl_vector_set_zero(filt.x);
    gsl_vector_set(filt.x, 1, V_X);
    gsl_vector_set(filt.x, 3, V_Y);

    // Initializing the covariance matrix P
    gsl_matrix_set_identity(filt.P);
    gsl_matrix_scale(filt.P, 0.5);

    printf("x_,dx_,y_,dy_,P_[0,0].P_[0,1],P_[1,0],P_[1,1]\n");

    for (int i = 0; i < N_STEPS; ++i)
    {
        if (cfilt_kalman_predict(&filt))
        {
            fprintf(stderr, "An error occured with the prediction step\n");
            break;
        }

        printf("%f,%f,%f,%f,", gsl_vector_get(filt.x_, 0), gsl_vector_get(filt.x_, 1), gsl_vector_get(filt.x_, 2),
               gsl_vector_get(filt.x_, 3));
        printf("%f,%f,%f,%f\n", gsl_matrix_get(filt.P_, 0, 0), gsl_matrix_get(filt.P_, 0, 1),
               gsl_matrix_get(filt.P_, 1, 0), gsl_matrix_get(filt.P_, 1, 1));

        // Because we are only predicting, x_ = x and P_ = P
        gsl_vector_memcpy(filt.x, filt.x_);
        gsl_matrix_memcpy(filt.P, filt.P_);
    }

cleanup:
    cfilt_kalman_free(&filt);
    gsl_rng_free(rng);

    return 0;
}
