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

#include "cfilt/kalman1d.h"

#include <gsl/gsl_randist.h>

#include <stdio.h>
#include <time.h>

#define N_STEPS 1000
#define DT 0.1
#define X0 0.0
#define V0 10.0
#define X_NOISE 2.0
#define V_NOISE 2.0

/**
 * This test emulates an entity moving in a straight line. Its sensors yield
 * position and velocity. Both state variables will be tracked independently.
 */
int
main(void)
{
    double real_position = X0;
    double real_velocity = V0;

    struct cfilt_gauss position = {.mean = X0, .var = X_NOISE * X_NOISE };
    struct cfilt_gauss velocity = {.mean = V0, .var = V_NOISE * V_NOISE };
    struct cfilt_gauss acceleration = {.mean = 0, .var = V_NOISE * V_NOISE };

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));

    printf("x_pred,v_pred,x,v,x_real,v_real,e_x,e_v\n");

    for (int i = 0; i < N_STEPS; ++i)
    {
        struct cfilt_gauss position_pred;
        struct cfilt_gauss velocity_pred;

        const struct cfilt_gauss dx = {.mean = velocity.mean * DT, .var = velocity.var * DT * DT };
        const struct cfilt_gauss dv = {.mean = acceleration.mean, .var = acceleration.var * DT * DT };

        cfilt_kalman1d_predict(&position_pred, position, dx);
        cfilt_kalman1d_predict(&velocity_pred, velocity, dv);

        real_position += real_velocity * DT;
        real_velocity += acceleration.mean * DT;

        const double x_noise = gsl_ran_gaussian(rng, X_NOISE);
        const double v_noise = gsl_ran_gaussian(rng, V_NOISE);

        struct cfilt_gauss position_z = {.mean = real_position + x_noise, X_NOISE * X_NOISE };
        struct cfilt_gauss velocity_z = {.mean = real_velocity + v_noise, V_NOISE * V_NOISE };

        cfilt_kalman1d_update(&position, position_pred, position_z);
        cfilt_kalman1d_update(&velocity, velocity_pred, velocity_z);

        const double position_error = real_position - position.mean;
        const double velocity_error = real_velocity - velocity.mean;

        printf("%f/%f,%f/%f,%f/%f,%f/%f,%f,%f,%f,%f\n", position_pred.mean, position_pred.var, velocity_pred.mean,
               velocity_pred.var, position.mean, position.var, velocity.mean, velocity.var, real_position,
               real_velocity, position_error, velocity_error);
    }

    gsl_rng_free(rng);

    return 0;
}
