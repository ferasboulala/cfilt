#ifndef _GH_H
#define _GH_H

/**
 * A simple implementation of the gh filter (also known as the alpha beta filter).
 * It tracks the values of an arbitrary amount of state variables that are related as such:
 *      x_(i+1) = d(x_i)/dt
 *
 * There are assumptions made with this filters. First, when tracking an set of state variables,
 * all intermediary variables will be tracked. For example, if we wish to track position and
 * acceleration, velocity will be tracked too. This is enforced by requiring a dimension for
 * the filter (dim := difference between the highest and lowest order variables + 1).
 *
 * When predicting the value of a state variable, we use all the higher order variables to do so.
 * That is, x_pred_i = x_i + dt * x_(i+1) + dt^2 * x_(i+2) * 1/2 + ... + dt^(n-i-1) * x_n * 1/(n-i-1)!.
 * For x_n, it is assumed that the value is constant only during the prediction.
 *
 * When updating the value of a state variable, the value of the residual is computed as:
 *      1) The variable received data: Residual = x_pred_i - z_i
 *      2) If not:
 *          a) The immediate lower order variable received data: Residual = x_pred_i - (x_pred_(i-1) - z_(i - 1)) / dt
 *          b) If not: Residual = 0.
 * In all cases, the update equation is x_i = x_i - gh_i * residual
 *
 * Runtime of functions with n being the filter dimensions:
 * init : Theta(n)
 * free : Theta(n)
 * write : O(1)
 * predict : Theta(n^2)
 * update : Theta(n)
 */

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

struct gh_filter
{
    double* gh;
    double* x;
    double* x_pred;
    size_t dim;

    char* upd_;
    double* z_;
    void* ptr_;
};

int gh_alloc(struct gh_filter* filt, size_t dim);

void gh_free(struct gh_filter* filt);

void gh_write(struct gh_filter* filt, double val, size_t ord);

void gh_predict(struct gh_filter* filt, double dt);

void gh_update(struct gh_filter* filt, double dt);

#ifdef __cplusplus
}
#endif

#endif // _GH_H
