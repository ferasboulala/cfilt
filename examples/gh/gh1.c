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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double
rand_minus_1_to_minus_1(void)
{
    const int x = rand();
    return 2.0 * ((double)x / (RAND_MAX - 1) - 0.5);
}

double
scaled_rand(const double scale)
{
    return scale * rand_minus_1_to_minus_1();
}

double
uniform_noise(const double scale)
{
    return scaled_rand(scale);
}

/**
 * This test emulates an entity moving in a straight line. Its sensors yield
 * position and velocity but we wish to track acceleration too.
 */
int
main(int argc, char** argv)
{
    if (argc != 11)
    {
        fprintf(stderr, "Usage : %s N_STEPS DT X_NOISE_MAG V_NOISE_MAG X0 V0 "
                        "A0 GH0 GH1 GH2\n",
                argv[0]);
        return -1;
    }

    int arg_i = 1;

    const int N_STEPS = atoi(argv[arg_i++]);
    const double DT = atof(argv[arg_i++]);

    const double X_NOISE_MAG = atof(argv[arg_i++]);
    const double V_NOISE_MAG = atof(argv[arg_i++]);

    const double X0 = atof(argv[arg_i++]);
    const double V0 = atof(argv[arg_i++]);
    const double A0 = atof(argv[arg_i++]);
    const double GH0 = atof(argv[arg_i++]);
    const double GH1 = atof(argv[arg_i++]);
    const double GH2 = atof(argv[arg_i++]);

    cfilt_gh_filter filter;
    cfilt_gh_alloc(&filter, 3);

    srand(time(NULL));

    double position = X0;
    double velocity = V0;
    const double acceleration = A0; // Constant acceleration

    filter.x[0] = 0;
    filter.x[1] = 0;
    filter.x[2] = 0;

    filter.gh[0] = GH0;
    filter.gh[1] = GH1;
    filter.gh[2] = GH2;

    printf("x,v,a,x_pred,v_pred,a_pred,z_x,z_v,x_,v_,a_\n");
    for (int i = 0; i < N_STEPS; ++i)
    {
        cfilt_gh_predict(&filter, DT);
        const double x_pred = filter.x_pred[0];
        const double v_pred = filter.x_pred[1];
        const double a_pred = filter.x_pred[2];

        const double z_x = position + uniform_noise(X_NOISE_MAG);
        const double z_v = velocity + uniform_noise(V_NOISE_MAG);

        cfilt_gh_write(&filter, z_x, 0);
        cfilt_gh_write(&filter, z_v, 1);

        cfilt_gh_update(&filter, DT);
        const double x_ = filter.x[0];
        const double v_ = filter.x[1];
        const double a_ = filter.x[2];

        position += velocity * DT;
        velocity += acceleration * DT;

        printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", position, velocity,
               acceleration, x_pred, v_pred, a_pred, z_x, z_v, x_, v_, a_);
    }

    cfilt_gh_free(&filter);

    return 0;
}
