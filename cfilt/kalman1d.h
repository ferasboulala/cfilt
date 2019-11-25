#ifndef KALMAN1D_H
#define KALMAN1D_H

#include "common.h"

/**
 * Simple implementation of a 1D Kalman Filter.
 * All gaussians are represented by the structure cfilt_gauss (a pair of mean and variance).
 */

// Predicts the state variable x according to the rate of change dx.
void cfilt_kalman1d_predict(struct cfilt_gauss* x_pred, struct cfilt_gauss x, struct cfilt_gauss dx);

// Updates the state variable x according to the sensor data z.
void cfilt_kalman1d_update(struct cfilt_gauss* x, struct cfilt_gauss x_pred, struct cfilt_gauss z);

#endif // KALMAN1D_H
