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
 * You should have received a copy of the GNU GEneral Public License
 * along with cfilt. If not, see <https://www.gnu.org/licenses/>.
 */

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
