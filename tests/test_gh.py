import matplotlib.pyplot as plt
import pandas as pd
import time
import os

df = pd.read_csv('test_gh.csv')

plt, axes = plt.subplots(3, 1)
x_plot, v_plot, a_plot = axes

x_plot.plot(df['x'])
x_plot.plot(df['z_x'])
x_plot.plot(df['x_'])
x_plot.legend(['x', '$z_x$', '$\hat{x}$'])
x_plot.set_xlabel('step')
x_plot.set_ylabel('position')
x_plot.set_title('Position tracking')

v_plot.plot(df['v'])
v_plot.plot(df['z_v'])
v_plot.plot(df['v_'])
v_plot.legend(['v', '$z_v$', '$\hat{v}$'])
v_plot.set_xlabel('step')
v_plot.set_ylabel('velocity')
v_plot.set_title('Velocity tracking')

a_plot.plot(df['a'])
a_plot.plot(df['a_'])
a_plot.legend(['a', '$\hat{a}$'])
a_plot.set_xlabel('step')
a_plot.set_ylabel('acceleration')
a_plot.set_title('Acceleration tracking')

plt.show()
