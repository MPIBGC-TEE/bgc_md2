from bokeh.plotting import figure, show
import numpy as np
from functools import reduce
from helpers import with_arrow_head


# create a new plot with a title and axis labels
p = figure(
    title="Simple line example", 
    x_axis_label='x', 
    y_axis_label='y'
)
# prepare some data
xs_1 = [ 3, 4, 5, 6, 7]
ys_1 = [ 1, 1, 1, 1, 1]
xs_2 = [ 3, 4, 5, 6, 7]
ys_2 = [ 2, 2, 2, 2, 2]
xss_1,yss_1=with_arrow_head(xs_1,ys_1,0.2,0.5)
xss_2,yss_2=with_arrow_head(xs_2,ys_2)
p.multi_line(
    [xss_1,xss_2],
    [yss_1,yss_2],
    line_color=["blue","red"], 
    line_width=4,
)
show(p)
