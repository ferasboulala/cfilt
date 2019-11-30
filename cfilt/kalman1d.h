#ifndef KALMAN1D_H
#define KALMAN1D_H

#include "cfilt/cfilt.h"

/**
 * Simple implementation of a 1D Kalman Filter.
 * All gaussians are represented by the structure cfilt_gauss (a pair of mean and variance).
 */

#ifdef __cplusplus
extern "C" {
#endif

void cfilt_kalman1d_predict(struct cfilt_gauss* x_pred, struct cfilt_gauss x, struct cfilt_gauss dx);

void cfilt_kalman1d_update(struct cfilt_gauss* x, struct cfilt_gauss x_pred, struct cfilt_gauss z);

#ifdef __cplusplus
}
#endif

#endif // KALMAN1D_H
