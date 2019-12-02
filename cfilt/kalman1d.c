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

void
cfilt_kalman1d_predict(struct cfilt_gauss* x_pred, struct cfilt_gauss x, struct cfilt_gauss dx)
{
    x_pred->mean = x.mean + dx.mean;
    x_pred->var = x.var + dx.var;
}

void
cfilt_kalman1d_update(struct cfilt_gauss* x, struct cfilt_gauss x_pred, struct cfilt_gauss z)
{
    const double residual = z.mean - x_pred.mean;
    const double kalman_gain = x_pred.var / (x_pred.var + z.var);
    x->mean = kalman_gain * residual + x_pred.mean;
    x->var = (1 - kalman_gain) * x_pred.var;
}
