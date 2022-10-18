# small example to set ticks
import matplotlib as mpl
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure(tight_layout=True)  
gs = fig.add_gridspec(
    1,
    3,
    width_ratios=(1, 1, 1),
)
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(
    gs[0,0],
    projection='3d',
    box_aspect=[.5,.5,1]
)
theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
x = r * np.sin(theta)
y = r * np.cos(theta)
ax.plot(x, y, z, label='parametric curve')
ax.plot(x, y, z*2, label='test')
ax.legend()

ax.set_xlabel('$X$', fontsize=20)
ax.set_ylabel('$Y$')
ax.yaxis._axinfo['label']['space_factor'] = 3.0
# set z ticks and labels
ax.set_zticks([-2, 0, 2])
ax.set_zticklabels(["a", 0, "b"])
ax.set_xticks([min(x), 0, max(x)])
ax.set_xticklabels(["x_min", 0, "x_max"])
ax.set_yticks([min(y), 0, max(y)])
ax.set_yticklabels(["y_min", 0, "y_max"])
# change fontsize
for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(10)
# disable auto rotation
ax.zaxis.set_rotate_label(False) 
ax.set_zlabel('$\gamma$', fontsize=30, rotation = 0)
plt.show()

#enter image description here
