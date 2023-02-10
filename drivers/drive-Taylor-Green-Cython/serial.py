""" Taylor-Green Vortex
"""
from phi.torch.flow import *
import matplotlib.pyplot as plt
assert backend.default_backend().set_default_device('GPU')

def taylor_green_pressure(x):
    return math.sum(math.cos(2 * x * VORTEX_COUNT), 'vector') / 4

def taylor_green_velocity(x):
    sin = math.sin(VORTEX_COUNT * x)
    cos = math.cos(VORTEX_COUNT * x)
    return stack({
        'x': sin.vector['x'] * cos.vector['y'] * cos.vector['z'],
        'y': -cos.vector['x'] * sin.vector['y'] * cos.vector['z'],
        'z': 0*sin.vector['z']},
        dim=channel('vector'))


@jit_compile
def step(v, p, dt):
    v = advect.mac_cormack(v, v, dt)
    v, p = fluid.make_incompressible(v, solve=Solve('auto', 1e-5, 0, x0=p))
    vort = field.curl(v.at_centers().z[n//4].vector[:2])

    return v, p, vort

n = 32
DOMAIN = dict(x=n, y=n, z=n, bounds=Box(x=2*PI, y=2*PI, z=2*PI), extrapolation=extrapolation.PERIODIC)
VORTEX_COUNT = 1
dt = 0.05
t = 0.

v0 = StaggeredGrid(taylor_green_velocity, **DOMAIN)  # also works with CenteredGrid
v0 = tensor(2*[v0.values], batch("init"), channel(vector='x,y,z'))
v0 = StaggeredGrid(v0, **DOMAIN)

noise_both = 1e-2 * Noise().at(v0)
noise_offset = 1e-5 * StaggeredGrid(math.random_normal(v0), **DOMAIN)

# v0 += noise_both    # adds the same noise to both, in order to break up symmettry in initial condition
v0 += noise_offset  # adds different noise to each field in batch

sim_velocity, sim_pressure = fluid.make_incompressible(v0)
vort = field.curl(sim_velocity.at_centers().z[n//4].vector[:2])

diff = []
cmap="plasma"

fig = plt.figure()
gs = fig.add_gridspec(2,2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])

im1 = ax1.imshow(vort.init[0].values.numpy('y,x'), cmap=cmap,  vmin=-4, vmax=4)

for _ in range(300):
    sim_velocity, sim_pressure, vort = step(sim_velocity, sim_pressure, dt)
    diff.append(field.l2_loss(sim_velocity.init[0] - sim_velocity.init[1])/math.sqrt(n**3*2))

    vort_diff = vort.init[0].values.numpy('y,x') - vort.init[1].values.numpy('y,x')
    im1.set_data(vort.init[1].values.numpy('y,x'))
    ax2.imshow(vort_diff)
    ax3.semilogy(diff)
    plt.pause(0.01)
    ax2.clear()
    ax3.clear()
   
plt.show()
    