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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define WHEEL_BASE 0.5
#define N_LANDMARKS 3

#define STATE_DIM 3
#define CONTROL_DIM 2
#define SENSOR_DIM N_LANDMARKS * 2 // distance and angle

#define ALPHA 0.00001
#define BETA 2.0
#define KAPPA 0.0

#define DT 0.01

int
F(cfilt_ukf* filt, void* ptr)
{
    const double dt = *((double*)ptr);

    gsl_vector* state = filt->x;
    const double heading = gsl_vector_get(state, 2);

    gsl_vector* control = filt->u_z;
    const double vel = gsl_vector_get(control, 0);
    const double steering_angle = gsl_vector_get(control, 1);

    const double dist = vel * dt;

    double estimate_state[] = { 0.0, 0.0, 0.0 };
    gsl_vector_view view = gsl_vector_view_array(estimate_state, 3);
    gsl_vector* estimate_state_vec = &view.vector;

    if (abs(steering_angle) > 0.0001)
    {
        const double tan_steering_angle = tan(steering_angle);
        const double beta = (dist / WHEEL_BASE) * tan_steering_angle;
        const double r = WHEEL_BASE / tan_steering_angle;
        const double sinh = sin(heading);
        const double sinh_ = sin(heading + beta);
        const double cosh = cos(heading);
        const double cosh_ = cos(heading + beta);

        estimate_state[0] = r * sinh_ - r * sinh;
        estimate_state[1] = r * cosh - r * cosh_;
        estimate_state[2] = beta;
    }
    else
    {
        estimate_state[0] = dist * cos(heading);
        estimate_state[1] = dist * sin(heading);
        estimate_state[2] = 0;
    }

    return gsl_vector_memcpy(filt->x_, filt->x) ||
           gsl_vector_add(filt->x_, estimate_state_vec);
}

double
normalize_angle(const double angle)
{
    double angle_ = fmod(angle, M_PI * 2);
    if (angle_ > M_PI)
    {
        angle_ -= 2 * M_PI;
    }

    return angle_;
}

int
H(cfilt_ukf* filt, void* ptr)
{
    gsl_matrix* landmarks = (gsl_matrix*)ptr;

    for (size_t i = 0; i < landmarks->size1; ++i)
    {
        gsl_vector_view landmark_view = gsl_matrix_row(landmarks, i);
        gsl_vector* landmark = &landmark_view.vector;
        const double p_x = gsl_vector_get(landmark, 0);
        const double p_y = gsl_vector_get(landmark, 1);

        const double x = gsl_vector_get(filt->x, 0);
        const double y = gsl_vector_get(filt->x, 1);
        const double theta = gsl_vector_get(filt->x, 2);

        const double dist = sqrt(pow(p_x - x, 2) + pow(p_y - y, 2));
        const double angle = atan2(p_y - y, p_x - x);
        const double normalized_angle = normalize_angle(angle - theta);

        double h_i[] = { dist, normalized_angle };
    }

    return GSL_SUCCESS;
}

int
main(int argc, char** argv)
{
    cfilt_sigma_generator* gen;
    if (cfilt_sigma_generator_alloc(CFILT_SIGMA_VAN_DER_MERWE, &gen, STATE_DIM,
                                    ALPHA, BETA, KAPPA) != GSL_SUCCESS)
    {
        fprintf(stderr,
                "An error occured creating the sigma point generator\n");
        return -1;
    }

    cfilt_ukf filt;
    if (cfilt_ukf_alloc(&filt, STATE_DIM, CONTROL_DIM, SENSOR_DIM, F, H, gen) !=
        GSL_SUCCESS)
    {
        fprintf(stderr,
                "An error occured creating the uscented kalman filter\n");
        return -1;
    }

    const double dt = DT;

    cfilt_ukf_free(&filt);
    cfilt_sigma_generator_free(gen);

    return 0;
}
