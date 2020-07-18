import matplotlib.pyplot as plt
from sympy import Symbol,Function,latex,prod,sin,pi
from sympy.physics.units import (
        convert_to,
        Quantity,
        length,
        time,
        day,
        second,
        minute,
        kelvin,
        meter,
        kilometer,
        kilogram)
from sympy.physics.units.systems import SI
from sympy.physics.units.systems.si import dimsys_SI
# only fixed qunatities (numbers times unit like in a parameter set)

from bgc_md2.described_quantities import describedQuantity
from bgc_md2.resolve.mvars import ParameterDict



a=Quantity("a")
SI.set_quantity_dimension(a,length)

t=Quantity("t")
SI.set_quantity_dimension(t,time)

res=a**2/t
# we can now determine the physical dimension of res
print(SI.get_dimensional_expr(res))



# In parameter dicts we can describe values along with units
a_val=Quantity("a_val")
SI.set_quantity_dimension(a_val,length)
SI.set_quantity_scale_factor(a_val,5*meter)

parameter_dict = {a:5*meter,t:4*second}

res_val=res.subs(parameter_dict)
# we can now determine the physical dimension of res_val
# and check it against the expected
print(SI.get_dimensional_expr(res_val))
#dimsys_SI.equivalent_dims(res,res_val)



# Now we test the actual use case
# for a plot.
xs=[i*second for i in range(1,10)]
ys=[meter/x for x in xs]

def plot_with_units(ax,xt,yt):
    xs,x_unit = xt
    ys,y_unit = yt
    xnums=[x/x_unit for x in xs]
    ynums=[y/y_unit for y in ys]
    ax.set_xlabel(str(x_unit))
    ax.set_ylabel(str(y_unit))
    ax.plot(xnums,ynums)

def auto_plot_with_units(ax,xs,ys):
    x_unit = prod((1.0*xs[0].n()).args[1:])
    y_unit = prod((1.0*ys[0].n()).args[1:])
    plot_with_units(ax,(xs,x_unit),(ys,y_unit)) 


fig = plt.figure()
#ax = fig.add_axes((0,0,1,1))
ax = fig.add_subplot(1,1,1)
#plot_with_units(ax,(xs,second),(ys,meter/second))
#auto_plot_with_units(ax,xs,ys)
#plt.show()



v=Function('v')
m=Quantity('m')
# create an expression containing a function
E=m/2*v(t)**2



# a function with uits
def v_unit(t):
    assert(
        dimsys_SI.equivalent_dims(
            SI.get_dimensional_expr(t),
            time
        )
    )
    omega=2*pi/second
    V_0=20*meter/second
    V_range=5*meter/second
    return V_0+V_range*sin(omega*t)

def v_num(t):
    omega=2*pi
    V_0=20
    V_range=5
    return V_0+V_range*sin(omega*t)


#print(v_unit(3*day))

#create a funcdict
funcDictUnit={v(t):v_unit}
funcDictNum={v(t):v_num}

parDictUnit={m:5*kilogram}
parDictNum ={m:5}

E_unit=E.subs(funcDictUnit).subs(parDictUnit)
E_num=E.subs(funcDictNum).subs(parDictNum)

E_num 

