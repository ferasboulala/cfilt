#include "cfilt/gh.h"

#include <gsl/gsl_errno.h>

#include <stdlib.h>
#include <string.h>

int
gh_alloc(struct gh_filter* filt, size_t dim)
{
    filt->dim = dim;
    filt->ptr_ = calloc(1, 4 * dim * sizeof(double) + dim * sizeof(char));
    if (filt->ptr_ == NULL)
    {
        GSL_ERROR("failed to allocate space for gh filter", GSL_ENOMEM);
    }

    filt->gh = filt->ptr_;
    filt->x = (void*)filt->gh + dim * sizeof(double);
    filt->x_pred = (void*)filt->x + dim * sizeof(double);
    filt->z_ = (void*)filt->x_pred + dim * sizeof(double);
    filt->upd_ = (void*)filt->z_ + dim * sizeof(double);

    return GSL_SUCCESS;
}

void
gh_free(struct gh_filter* filt)
{
    free(filt->ptr_);
    memset(filt, 0, sizeof(struct gh_filter));
}

void
gh_write(struct gh_filter* filt, double val, size_t order)
{
    filt->z_[order] = val;
    filt->upd_[order] = 1;
}

void
gh_predict(struct gh_filter* filt, double dt)
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

void
gh_update(struct gh_filter* filt, double dt)
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
