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

#define X0 100
#define Y0 50

/**
 * This test emulates an entity moving in a straight line, in 2D. Its sensors yield
 * position and velocity.
 */
int
main(int argc, char** argv)
{
    if (argc != 10)
    {
        fprintf(stderr, "Usage : test_kalman N_STEPS DT V_X V_Y X_NOISE Y_NOISE A_X A_Y Q_VAR\n");
        return -1;
    }

    const int N_STEPS = atoi(argv[1]);
    const double DT = atof(argv[2]);
    const double V_X = atof(argv[3]);
    const double V_Y = atof(argv[4]);
    const double X_NOISE = atof(argv[5]);
    const double Y_NOISE = atof(argv[6]);
    const double A_X = atof(argv[7]);
    const double A_Y = atof(argv[8]);
    const double Q_VAR = atof(argv[9]);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));

    cfilt_kalman_filter filt;
    if (cfilt_kalman_alloc(&filt, 4, 2, 2))
    {
        fprintf(stderr, "Could not allocate kalman filter memory\n");
        goto cleanup;
    }

    // Setting up the state transition matrix F
    // [1 DT 0 0
    //  0  1 0 0
    //  0  0 1 DT
    //  0  0 0 1]
    gsl_matrix_set_identity(filt.F);
    gsl_matrix_set(filt.F, 0, 1, DT);
    gsl_matrix_set(filt.F, 2, 3, DT);

    // Setting up the control matrix B
    // [0 0
    //  1 0
    //  0 0
    //  0 1]
    gsl_matrix_set_zero(filt.B);
    gsl_matrix_set(filt.B, 1, 0, 1);
    gsl_matrix_set(filt.B, 3, 1, 1);

    // Setting up the control input vector u (changes)
    // [V_X
    //  V_Y]
    gsl_vector_set(filt.u, 0, V_X);
    gsl_vector_set(filt.u, 1, V_Y);

    // Setting up the process covariance matrix Q
    gsl_matrix_set_identity(filt.Q);
    gsl_matrix_scale(filt.Q, Q_VAR);

    // Setting up the measurement matrix H
    // [1 0 0 0
    //  0 0 1 0]
    gsl_matrix_set_zero(filt.H);
    gsl_matrix_set(filt.H, 0, 0, 1);
    gsl_matrix_set(filt.H, 1, 2, 1);

    // Setting up the measurement covariance matrix R
    // [X_NOISE 0
    //     0  Y_NOISE]
    gsl_matrix_set_zero(filt.R);
    gsl_matrix_set(filt.R, 0, 0, X_NOISE);
    gsl_matrix_set(filt.R, 1, 1, Y_NOISE);

    // Initializing the state vector x
    // [0 V_X 0 V_Y]^T
    gsl_vector_set_zero(filt.x);
    gsl_vector_set(filt.x, 1, V_X);
    gsl_vector_set(filt.x, 3, V_Y);

    // Initializing the covariance matrix P
    gsl_matrix_set_zero(filt.P);
    gsl_matrix_set(filt.P, 0, 0, X_NOISE);
    gsl_matrix_set(filt.P, 2, 2, Y_NOISE);

    printf("x_,dx_,y_,dy_,x,x_var,dx,y,y_var,dy,x_real,dx_real,y_real,dy_real\n");

    double x = X0;
    double y = Y0;
    double v_x = V_X;
    double v_y = V_Y;

    for (int i = 0; i < N_STEPS; ++i)
    {
        if (cfilt_kalman_predict(&filt))
        {
            fprintf(stderr, "An error occured with the prediction step\n");
            break;
        }

        printf("%f,%f,%f,%f,", gsl_vector_get(filt.x_, 0), gsl_vector_get(filt.x_, 1), gsl_vector_get(filt.x_, 2),
               gsl_vector_get(filt.x_, 3));

        const double x_noise = gsl_ran_gaussian(rng, X_NOISE);
        const double y_noise = gsl_ran_gaussian(rng, Y_NOISE);

        x += DT * v_x;
        y += DT * v_y;
        v_x += DT * A_X;
        v_y += DT * A_Y;

        gsl_vector_set(filt.z, 0, x + x_noise);
        gsl_vector_set(filt.z, 1, y + y_noise);

        if (cfilt_kalman_update(&filt))
        {
            fprintf(stderr, "An error occured with the update step\n");
            break;
        }

        gsl_vector_set(filt.u, 0, DT * A_X);
        gsl_vector_set(filt.u, 1, DT * A_Y);

        printf("%f,%f,%f,%f,%f,%f,", gsl_vector_get(filt.x, 0), gsl_matrix_get(filt.P, 0, 0), gsl_vector_get(filt.x, 1), gsl_vector_get(filt.x, 2), gsl_matrix_get(filt.P, 2, 2),
               gsl_vector_get(filt.x, 3));

        printf("%f,%f,%f,%f\n", x, v_x, y, v_y);
    }

cleanup:
    cfilt_kalman_free(&filt);
    gsl_rng_free(rng);

    return 0;
}
