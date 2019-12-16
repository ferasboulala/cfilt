import matplotlib.pyplot as plt
import pandas as pd
import time
import os

df = pd.read_csv('test_kalman1d.csv')

plt, axes = plt.subplots(2, 1)
x_plot, v_plot = axes

steps = [i for i in range(len(df))]

x_plot.errorbar(steps, df['x_pred'], yerr=df['x_pred_var'])
x_plot.errorbar(steps, df['z_x'], yerr=df['z_x_var'])
x_plot.errorbar(steps, df['x'], yerr=df['x_var'])
x_plot.plot(df['x_real'])
x_plot.legend(['$\widetilde{x}$', '$z_x$', '$\hat{x}$', 'x'])
x_plot.set_xlabel('step')
x_plot.set_ylabel('position')
x_plot.set_title('Position tracking')

v_plot.errorbar(steps, df['v_pred'], yerr=df['v_pred_var'])
v_plot.errorbar(steps, df['z_v'], yerr=df['z_v_var'])
v_plot.errorbar(steps, df['v'], yerr=df['v_var'])
v_plot.plot(df['v_real'])
v_plot.legend(['$\widetilde{v}$', '$z_v$', '$\hat{v}$', 'v'])
v_plot.set_xlabel('step')
v_plot.set_ylabel('velocity')
v_plot.set_title('Velocity tracking')

plt.show()
