import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd

import io
import os
import pathlib
import subprocess
import sys

BIN_DIR = os.path.join(pathlib.Path(__file__).resolve().parents[1], 'bin')

def get_program_output(binary, *args):
    """ Runs a command with arguments and returns the output as a list of lines. """
    bin_path = os.path.join(BIN_DIR, binary)
    proc = subprocess.Popen([bin_path, *args], stdout=subprocess.PIPE)
    
    return list(map(lambda l: l.strip(), io.TextIOWrapper(proc.stdout, encoding='utf-8')))

def get_program_output_df(binary, *args):
    """ Runs a command with arguments and returns the output as a dataframe. """
    bin_path = os.path.join(BIN_DIR, binary)
    proc = subprocess.Popen([bin_path, *args], stdout=subprocess.PIPE)

    return pd.read_csv(proc.stdout)

def csv_line_to_list(line):
    """ Converts comma separated floating values to a list of floats. """
    return np.array(list(map(lambda v: float(v), line.split(','))))

def csv_line_to_np_array(line):
    """ Converts comma separated floating values to a numpy array. """
    return np.array(csv_line_to_list(line))

def generate_ellipse_2D(ax, mean, cov, n_std=1, **kwargs):
    """ Returns an ellipse object given a mean, a covariance and the number of stddev. """
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse(
            (0, 0),
            width=ell_radius_x * 2,
            height=ell_radius_y * 2,
            alpha=0.2,
            facecolor='none',
            **kwargs,
    )

    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std

    mean_x, mean_y = mean[0], mean[1]

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)

    return ellipse

def stringify(param):
    """ Converts a sequence of objects into their string representation. """
    return list(map(lambda p: str(p), param))

def extract_results_van_der_merwe(generator):
    """ Parses the output of the Van Der Merwe example program.  """
    assert next(generator) == 'points'
    next_line = next(generator)
    X, Y = [], []
    while next_line != 'mu_weights':
        data = csv_line_to_list(next_line)
        X.append(data[0])
        Y.append(data[1])
        next_line = next(generator)

    m = len(X)

    assert next_line == 'mu_weights'
    mu_weights = csv_line_to_np_array(next(generator))
    assert m == len(mu_weights)

    assert next(generator) == 'sigma_weights'
    sigma_weights = csv_line_to_np_array(next(generator))
    assert m == len(sigma_weights)

    assert next(generator) == "mu"
    mu = csv_line_to_np_array(next(generator))

    assert next(generator) == "cov"
    cov = []
    for _ in range(2):
        row = csv_line_to_list(next(generator))
        cov.append(row)
    cov = np.array(cov)

    return X, Y, mu_weights, sigma_weights, mu, cov

def plot_result_van_der_merwe(ax, X, Y, mu_weights, sigma_weights, mu, cov, alpha, beta, kappa):
    """ Plots the results of the Van Der Merwe example program. """
    ellipse = generate_ellipse_2D(ax, mu, cov, 1, edgecolor='firebrick')
    ax.add_patch(ellipse)
    ellipse = generate_ellipse_2D(ax, mu, cov, 2, edgecolor='fuchsia')
    ax.add_patch(ellipse)
    ellipse = generate_ellipse_2D(ax, mu, cov, 3, edgecolor='blue')
    ax.add_patch(ellipse)

    scale = np.abs(mu_weights) / np.max(np.abs(mu_weights)) * 500
    ax.scatter(X, Y, s=scale)

    title = 'Sigma points alpha={}, beta={}, kappa={}'.format(alpha, beta, kappa)
    title = title.format(title)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')

def plot_results_van_der_merwe(results, parameters):
    assert len(results) == len(parameters)

    fig, axes = plt.subplots(1, len(results))
    for ax, result, params in zip(axes, results, parameters):
        plot_result_van_der_merwe(ax, *result, *params)

    plt.show()

def van_der_merwe():
    """ Van Der Merwe Sigma Point Generation Example Program main function. """
    parameters = [
            (0.3, 2, 3),
            (1.0, 2, 3),
            (0.5, 2, -1),
    ]
    results = []
    for param in parameters:
        str_params = stringify(param)
        result = get_program_output('van_der_merwe', *str_params)
        generator = iter(result)
        result = extract_results_van_der_merwe(generator)
        results.append(result)

    plot_results_van_der_merwe(results, parameters)

def plot_result_constant_accel(ax, df, *params):
    x = df['x']
    v = df['v']
    a = df['a']

    x_pred = df['x_pred']
    v_pred = df['v_pred']
    a_pred = df['a_pred']

    z_x = df['z_x']
    z_v = df['z_v']

    x_ = df['x_']
    v_ = df['v_']
    a_ = df['a_']

    x_ax, v_ax, a_ax = ax

    x_ax.plot(x)
    x_ax.plot(x_pred)
    x_ax.plot(z_x)
    x_ax.plot(x_)
    x_ax.legend(['x', 'x_pred', 'z_x', 'x_'])

    v_ax.plot(v)
    v_ax.plot(v_pred)
    v_ax.plot(z_v)
    v_ax.plot(v_)
    v_ax.legend(['v', 'v_pred', 'z_v', 'v_'])

    a_ax.plot(a)
    a_ax.plot(a_pred)
    a_ax.plot(a_)
    a_ax.legend(['a', 'a_pred', 'a_'])

def plot_results_constant_accel(results, parameters):
    assert len(results) == len(parameters)

    fig, axes = plt.subplots(len(results), 3)
    for ax, result, params in zip(axes, results, parameters):
        plot_result_constant_accel(ax, result, *params)

    plt.show()

def constant_accel():
    """ GH filter Example Program main function. """
    'N_STEPS DT X_NOISE_MAG V_NOISE_MAG X0 V0 A0 GH0 GH1 GH2'
    parameters = [
            (100, 0.1, 1, 1, 0, 0, 1, 0.5, 0.5, 0.5),
            (100, 0.1, 1, 1, 500, 100, 1, 0.5, 0.5, 0.5),
            (100, 0.1, 5, 5, 0, 0, 1, 0.5, 0.5, 0.5),
    ]

    results = []
    for param in parameters:
        str_params = stringify(param)
        result = get_program_output_df('constant_accel', *str_params)
        results.append(result)

    plot_results_constant_accel(results, parameters)

def main():
    if len(sys.argv) != 2:
        print("Usage: python {} NAME_OF_THE_EXAMPLE_PROGRAM".format(sys.argv[0]))
        sys.exit(1)

    program = sys.argv[1]
    print(program)
    if program == "van_der_merwe":
        van_der_merwe()
    elif program == "constant_accel":
        constant_accel()
    else:
        print("Unknown program")
        sys.exit(1)

if __name__ == '__main__':
    main()
