
""" Taylor-Green Vortex
"""
from phi.torch.flow import *
assert backend.default_backend().set_default_device('GPU')

def taylor_green_pressure(x):
    return math.sum(math.cos(2 * x * VORTEX_COUNT), 'vector') / 4 * math.exp(-4 * VORTEX_COUNT ** 2 * t / RE)

def taylor_green_velocity(x):
    sin = math.sin(VORTEX_COUNT * x)
    cos = math.cos(VORTEX_COUNT * x)
    return math.exp(-2 * VORTEX_COUNT ** 2 * t / RE) * stack({
        'x': sin.vector['x'] * cos.vector['y'] * cos.vector['z'],
        'y': -cos.vector['x'] * sin.vector['y'] * cos.vector['z'],
        'z': 0*sin.vector['z']},
        dim=channel('vector'))


@jit_compile
def step(v, p, dt):
    v = advect.mac_cormack(v, v, dt)
    v, p = fluid.make_incompressible(v, solve=Solve('auto', 1e-5, 0, x0=p))
    return v, p

n = 128
DOMAIN = dict(x=n, y=n, z=n, bounds=Box(x=2*PI, y=2*PI, z=2*PI), extrapolation=extrapolation.PERIODIC)
VORTEX_COUNT = 1
RE = 1600.  # Reynolds number for analytic function
dt = vis.control(0.01)
t = 0.

sim_pressure = CenteredGrid(taylor_green_pressure, **DOMAIN)
sim_velocity = StaggeredGrid(taylor_green_velocity, **DOMAIN)  # also works with CenteredGrid
lyap = StaggeredGrid(Noise(), **DOMAIN)

sim_velocity, sim_pressure = fluid.make_incompressible(sim_velocity)
vort = field.curl(sim_velocity.at_centers().z[n//4].vector[:2])

viewer = view(vort, namespace=globals())
for _ in viewer.range():
    t += dt
    sim_velocity, sim_pressure = step(sim_velocity, sim_pressure, dt)
    vort = field.curl(sim_velocity.at_centers().z[n//4].vector[:2])
