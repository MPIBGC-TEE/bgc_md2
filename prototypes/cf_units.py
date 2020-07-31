from cf_units import Unit
import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol


c = Unit("deg_c")
k = Unit("deg_k")
c.convert(0, k)

# Now we test the actual use case
# for a plot.

second = Unit("second")
minute = Unit("minute")
meter = Unit("meter")
xs = [i * second for i in range(1, 10)]
ys = [meter / x for x in xs]


def plot_with_units(ax, xt, yt):
    xs, x_unit = xt
    ys, y_unit = yt
    xnums = [x.convert(1, x_unit) for x in xs]
    ynums = [y.convert(1, y_unit) for y in ys]
    print(xnums, ynums)
    ax.set_xlabel(str(x_unit))
    ax.set_ylabel(str(y_unit))
    ax.plot(xnums, ynums)


def auto_plot_with_units(ax, xt, yt):
    raise Exception(
        """
        it is not possible to derive the standard unit from
        from the values. See documentation of cf_unit
        """
    )


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot_with_units(ax, (xs, minute), (ys, meter / minute))

plt.show()
