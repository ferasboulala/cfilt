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
    x->var = x_pred.var * z.var / (x_pred.var + z.var);
}
