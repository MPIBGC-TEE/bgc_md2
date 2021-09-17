#https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
import matplotlib.pyplot as plt
import numpy as np

xs = np.linspace(0, 2 * np.pi, 100)
ys1 = np.sin(xs)
ys2 = np.sin(xs + np.pi / 16)

fig, ax = plt.subplots(figsize=(9.2, 5))
ax.set_title("blue chasing red")
ax.plot(xs,ys1,color='red', label='sin(x)')
ax.plot(xs,ys2,color='blue', label='sin(x+pi/16)')
ax.legend()
fig.savefig('myfigure.pdf')
