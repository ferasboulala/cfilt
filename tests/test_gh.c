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

#define RAND() ((double)rand() / (RAND_MAX - 1) - 0.5)
#define SCRAND(scale) ((scale)*RAND())
#define UNOISE(scale) SCRAND((scale))

#define X0 0.0

/**
 * This test emulates an entity moving in a straight line. Its sensors yield
 * position and velocity but we wish to track acceleration aswell.
 */
int
main(int argc, char** argv)
{
    if (argc != 10)
    {
        fprintf(stderr, "Usage : test_gh N_STEPS DT X_NOISE V_NOISE A0 V0 GH0 GH1 GH2\n");
        return -1;
    }

    const int N_STEPS = atoi(argv[1]);
    const double DT = atof(argv[2]);
    const double X_NOISE = atof(argv[3]);
    const double V_NOISE = atof(argv[4]);
    const double A0 = atof(argv[5]);
    const double V0 = atof(argv[6]);
    const double GH0 = atof(argv[7]);
    const double GH1 = atof(argv[8]);
    const double GH2 = atof(argv[9]);

    cfilt_gh_filter filter;
    cfilt_gh_alloc(&filter, 3);

    srand(time(NULL));

    double position = X0;
    double velocity = V0;
    double acceleration = A0;

    filter.x[0] = position;
    filter.x[1] = velocity;
    filter.x[2] = acceleration;

    filter.gh[0] = GH0;
    filter.gh[1] = GH1;
    filter.gh[2] = GH2;

    printf("x,v,a,x_pred,v_pred,a_pred,z_x,z_v,x_,v_,a_,e_x,e_v,e_a\n");

    for (int i = 0; i < N_STEPS; ++i)
    {
        cfilt_gh_predict(&filter, DT);
        const double x_pred = filter.x_pred[0];
        const double v_pred = filter.x_pred[1];
        const double a_pred = filter.x_pred[2];

        const double z_x = position + UNOISE(X_NOISE);
        const double z_v = velocity + UNOISE(V_NOISE);

        cfilt_gh_write(&filter, z_x, 0);
        cfilt_gh_write(&filter, z_v, 1);

        cfilt_gh_update(&filter, DT);
        const double x = filter.x[0];
        const double v = filter.x[1];
        const double a = filter.x[2];

        position += velocity * DT;
        velocity += acceleration * DT;

        const double e_x = position - filter.x[0];
        const double e_v = velocity - filter.x[1];
        const double e_a = acceleration - filter.x[2];

        printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", position, velocity, acceleration, x_pred, v_pred, a_pred,
               z_x, z_v, x, v, a, e_x, e_v, e_a);
    }

    cfilt_gh_free(&filter);

    return 0;
}
