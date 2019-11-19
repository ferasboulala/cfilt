#include "gh.h"

#include <stdlib.h>
#include <string.h>

void gh_alloc(struct gh_filter *filt, size_t dim)
{
    filt->dim = dim;
    filt->gh = malloc(dim * sizeof(double));
    filt->x = malloc(dim * sizeof(double));
    filt->x_pred = malloc(dim * sizeof(double));
    filt->z_ = malloc(dim * sizeof(double));
    filt->upd_ = calloc(dim, sizeof(char));
}

void gh_free(struct gh_filter *filt)
{
    free(filt->gh);
    free(filt->x);
    free(filt->x_pred);
    free(filt->z_);
    free(filt->upd_);
    memset(filt, 0, sizeof(struct gh_filter));
}

void gh_write(struct gh_filter *filt, double val, size_t order)
{
    filt->z_[order] = val;
    filt->upd_[order] = 1;
}

void gh_predict(struct gh_filter *filt, double dt)
{
    for (size_t i = 0; i < filt->dim - 1; ++i)
    {
        filt->x_pred[i] = filt->x[i];
        unsigned int denum = 1;
        double dt_ = dt;
        for (size_t j = i + 1; j < filt->dim; ++j)
        {
            filt->x_pred[i] += filt->x[j] * dt_ / denum;
            dt_ *= dt;
            denum *= (denum + 1);
        }
    }
    filt->x_pred[filt->dim - 1] = filt->x[filt->dim - 1];
}

void gh_update(struct gh_filter *filt, double dt)
{
    filt->x[0] += filt->upd_[0] ? (filt->gh[0] * (filt->z_[0] - filt->x_pred[0])) : filt->x_pred[0];
    filt->upd_[0] = 0;
    for (size_t i = filt->dim - 1; i > 0; --i)
    {
        filt->x[i] = filt->x_pred[i];
        double residual = 0.0;
        if (filt->upd_[i])
        {
            residual = filt->z_[i] - filt->x_pred[i];
            filt->upd_[i] = 0;
        }
        else if (filt->upd_[i - 1])
        {
            residual = (filt->z_[i - 1] - filt->x[i - 1]) / dt - filt->x_pred[i];
        }

        filt->x[i] += filt->gh[i] * residual;
    }
}
