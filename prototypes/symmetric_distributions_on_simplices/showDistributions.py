import numpy as np
import matplotlib.pyplot as plt
fig=plt.figure()
g=np.random.default_rng()

def dont_draw_x(ind):
    # the index argument is actually ignored and is just there to be able to use map
    y,z=g.random(2)
    return np.array([1-y-z, y, z]) if y+z <= 1 else dont_draw_x(ind)

def dont_draw_y(ind):
    # the index argument is actually ignored and is just there to be able to use map
    x,z=g.random(2)
    return np.array([x, 1-x-z, z]) if x+z <= 1 else dont_draw_y(ind)

def dont_draw_z(ind):
    # the index argument is actually ignored and is just there to be able to use map
    x,y=g.random(2)
    return np.array([x, y, 1-x-y]) if x+y <= 1 else dont_draw_z(ind)


mean,sigma=0.5,0.05
def dont_draw_x_normal(ind):
    # the index argument is actually ignored and is just there to be able to use map
    y,z=g.normal(loc=[mean,mean],scale=[sigma,sigma])
    return np.array([1-y-z, y, z]) if y>0 and z >0 and y+z <= 1 else dont_draw_x(ind)

def dont_draw_y_normal(ind):
    # the index argument is actually ignored and is just there to be able to use map
    x,z=np.abs(g.normal(loc=[mean,mean],scale=[sigma,sigma]))
    return np.array([x, 1-x-z, z]) if x+z <= 1 else dont_draw_y(ind)

def dont_draw_z_normal(ind):
    # the index argument is actually ignored and is just there to be able to use map
    x,y=np.abs(g.normal(loc=[mean,mean],scale=[sigma,sigma]))
    return np.array([x, y, 1-x-y]) if x+y <= 1 else dont_draw_z(ind)

def sampler_maker(draw):
    def sampler(n):
        return np.stack(list(map(draw,range(n))),axis=0)

    return sampler

n = 50000
n_bins = 20
uf_samplers = list(map(sampler_maker, [dont_draw_x, dont_draw_y, dont_draw_z]))

res_wox,res_woy,res_woz = map(lambda f:f(n),uf_samplers)

n_samplers = list(map(sampler_maker, [dont_draw_x_normal, dont_draw_y_normal, dont_draw_z_normal]))

res_wox_n,res_woy_n,res_woz_n = map(lambda f:f(n),n_samplers)
res_dirichlet =np.random.dirichlet(alpha=[1,1,1],size=(n))             

fig=plt.figure(figsize=(15,10))
axs=fig.subplots(3,4)

axs[0,0].set_title('draw y,z , x=1-y-z')
axs[0,1].set_title('draw_x,z , y=1-z-x')
axs[0,2].set_title('draw_x,y , z=1-x-y')
axs[0,3].set_title('dirichlet')
x_labels = ('x','y','z')
for r in range(3):
    pwox=axs[r,0]
    pwox.hist(res_wox[:,r],bins=n_bins)
    pwox.set_xlabel(x_labels[r])
    pwox.set_ylabel('frequency')

    pwoy=axs[r,1]
    pwoy.hist(res_woy[:,r],bins=n_bins)
    pwox.set_xlabel(x_labels[r])
    
    pwoz=axs[r,2]
    pwoz.hist(res_woz[:,r],bins=n_bins)
    pwox.set_xlabel(x_labels[r])

    pdirichlet=axs[r,3]
    pdirichlet.hist(res_dirichlet[:,r],bins=n_bins)
    pdirichlet.set_xlabel(x_labels[r])
fig.suptitle('Uniform distributions on simplex')   
fig.savefig('uniform_histograms.pdf')

fig=plt.figure(figsize=(15,10))
axs=fig.subplots(3,4)

axs[0,0].set_title('draw y,z, compute x=1-y-z')
axs[0,1].set_title('draw x,z, compute y=1-z-x')
axs[0,2].set_title('draw x,y, compute z=1-x-y')
axs[0,3].set_title('dirichlet')
x_labels = ('x','y','z')
for r in range(3):
    pwox=axs[r,0]
    pwox.hist(res_wox_n[:,r],bins=n_bins)
    pwox.set_xlabel(x_labels[r])
    pwox.set_ylabel('frequency')

    pwoy=axs[r,1]
    pwoy.hist(res_woy_n[:,r],bins=n_bins)
    pwox.set_xlabel(x_labels[r])
    
    pwoz=axs[r,2]
    pwoz.hist(res_woz_n[:,r],bins=n_bins)
    pwox.set_xlabel(x_labels[r])

    pdirichlet=axs[r,3]
    pdirichlet.hist(res_dirichlet[:,r],bins=n_bins)
    pdirichlet.set_xlabel(x_labels[r])
fig.suptitle('Normal distributions on simplex')   
fig.savefig('normal_histograms.pdf')
