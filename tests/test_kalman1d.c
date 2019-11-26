#include "cfilt/kalman1d.h"
#include "common.h"

#include <string.h>

#define N_STEPS 1000
#define DT 0.1
#define X0 0.0
#define V0 1.0
#define VEL 10.0
#define V_MAX 100.0
#define X_NOISE 1.0
#define V_NOISE 1.0

/**
 * This test emulates an entity moving in a straight line. Its sensors yield
 * position and velocity. Both state variables will be tracked independently.
 */

// TODO : Use gsl's functions to get gaussian noise and add it to common.h

int
main(void)
{
    srand(time(NULL));

    double real_position = X0;
    double real_velocity = V0;

    struct cfilt_gauss position = { .mean = X0, .var = X_NOISE * X_NOISE };
    struct cfilt_gauss velocity = { .mean = V0, .var = V_NOISE * V_NOISE };
    struct cfilt_gauss dummy_acceleration;
    bzero(&dummy_acceleration, sizeof(struct cfilt_gauss));

    printf("x_pred,v_pred,x,v,x_real,v_real,e_x,e_v\n");

    for (int i = 0; i < N_STEPS; ++i)
    {
        struct cfilt_gauss position_pred;
        struct cfilt_gauss velocity_pred;

        cfilt_kalman1d_predict(&position_pred, position, velocity);
        cfilt_kalman1d_predict(&velocity_pred, velocity, dummy_acceleration);

        real_position += real_velocity * DT;
        real_velocity += 0;

        struct cfilt_gauss position_z = { .mean = real_position + GNOISE(X_NOISE), X_NOISE * X_NOISE };
        struct cfilt_gauss velocity_z = { .mean = real_velocity + GNOISE(V_NOISE), V_NOISE * V_NOISE };

        cfilt_kalman1d_update(&position, position_pred, position_z);
        cfilt_kalman1d_update(&velocity, velocity_pred, velocity_z);

        const double position_error = real_position - position.mean;
        const double velocity_error = real_velocity - velocity.mean;

        printf("%f/%f,%f/%f,%f/%f,%f/%f,%f,%f,%f,%f\n", position_pred.mean, position_pred.var, velocity_pred.mean, velocity_pred.var, position.mean, position.var, velocity.mean, velocity.var, real_position, real_velocity, position_error, velocity_error);
    }

    return 0;
}
