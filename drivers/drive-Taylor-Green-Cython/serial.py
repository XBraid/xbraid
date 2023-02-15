""" Taylor-Green Vortex
"""
from phi.torch.flow import *
from progressbar import progressbar
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
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

def norm(v):
    math.sqrt(field.l2_loss(v)).numpy('x')[0]

@jit_compile
def step(v, p, dt):
    # diffusion with strang splitting
    v = advect.mac_cormack(v, v, dt/2)
    v = diffuse.explicit(v, 1/RE, dt, substeps=1)
    v = advect.mac_cormack(v, v, dt/2)
    # No diffusion
    # v = advect.mac_cormack(v, v, dt)
    v, p = fluid.make_incompressible(v, solve=Solve('auto', 1e-5, 0, x0=p))

    return v, p

mindt = 0.1
def richardson_step(v, p, dt, errtol=1e-4, ord=1):
    errtol *= field.l2_loss(v).sum # make errtol relative
    v2h, p2h = step(v, p, dt)
    vh, ph = step(v, p, dt/2)
    vh, ph = step(vh, ph, dt/2)
    e = (vh - v2h)/(2**ord - 1)
    diff = math.sqrt(field.l2_loss(e).sum)
    if diff <= errtol or dt < mindt:
        return vh + e, p
    else:
        print(f"failed errtol({errtol}) with h={dt}")
        # refine time-step size and substep:
        vh, ph = richardson_step(v, p, dt/2)
        vh, ph = richardson_step(vh, ph, dt/2)
        return vh, ph

n = 64
RE = 1600
DOMAIN = dict(x=n, y=n, z=n, bounds=Box(x=2*PI, y=2*PI, z=2*PI), extrapolation=extrapolation.PERIODIC)
VORTEX_COUNT = 1
dt = .2
t = 0.

v0 = StaggeredGrid(tensor(np.load(f"./checkpoint_{n}.npy"), spatial("x,y,z"), channel(vector="x,y,z")), **DOMAIN)
# v0 = StaggeredGrid(taylor_green_velocity, **DOMAIN)  # also works with CenteredGrid
v0 = tensor(2*[v0.values], batch("init"), channel(vector='x,y,z'))
v0 = StaggeredGrid(v0, **DOMAIN)

noise_both = 1e-2 * Noise().at(v0)
noise_offset = 1e-5 * StaggeredGrid(math.random_normal(v0), **DOMAIN)

# v0 += noise_both    # adds the same noise to both, in order to break up symmetry in initial condition
v0 += noise_offset  # adds different noise to each field in batch

sim_velocity, sim_pressure = fluid.make_incompressible(v0)
vort = field.curl(sim_velocity.at_centers().z[n//4].vector[:2])
vort_diff = vort.init[0].values.numpy('y,x') - vort.init[1].values.numpy('y,x')

diff = [math.sqrt(field.l2_loss(sim_velocity.init[0] - sim_velocity.init[1])/(n**3*2)).numpy('x')[0]]
time = [0]
cmap="plasma"

fig = plt.figure()
gs = fig.add_gridspec(2,2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])
ax3.set_xlabel("time")
ax3.set_ylabel("norm difference")
ax1.set_title("vorticity")
ax2.set_title("difference")
fig.tight_layout()

im1 = ax1.imshow(vort.init[0].values.numpy('y,x'), cmap=cmap,  vmin=-4, vmax=4)
im2 = ax2.imshow(vort_diff, cmap=cmap)
line, = ax3.semilogy(time, diff)

def update(frame):
    global sim_velocity, sim_pressure
    print(f"frame: {frame}")
    sim_velocity, sim_pressure = richardson_step(sim_velocity, sim_pressure, dt, ord=2)
    # sim_velocity, sim_pressure = step(sim_velocity, sim_pressure, dt)
    diff.append(math.sqrt(field.l2_loss(sim_velocity.init[0] - sim_velocity.init[1])/(n**3*2)).numpy('x')[0])
    time.append(time[-1] + dt)
    
    vort = field.curl(sim_velocity.at_centers().z[n//4].vector[:2])
    vort_diff = vort.init[0].values.numpy('y,x') - vort.init[1].values.numpy('y,x')
    vort = vort.init[0].values.numpy('y,x')
    im1.set_data(vort)
    im1.set_clim(vort.min(), vort.max())
    im2.set_data(vort_diff)
    im2.set_clim(vort_diff.min(), vort_diff.max())
    line.set_data(time, diff)
    ax3.set_ylim(min(diff), max(diff))
    ax3.set_xlim(0., time[-1])
    # line, = ax3.semilogy(diff)
    return im1, im2, line

# for frame in progressbar(range(200), redirect_stdout=True):
#     update(frame)
#     plt.pause(0.01)

frames = range(200)
frames = progressbar(frames, redirect_stdout=True)
ani = FuncAnimation(fig, update, frames=frames, blit=True)
ani.save(f"./turbulence_{n}.gif")
# ani.save(f"./turbulence_{n}_check1.gif")

# checkpoint:
np.save(f"./checkpoint_{n}.npy", sim_velocity.init[0].data.numpy('x,y,z,vector'))
