
"""Karman Vortex Street
Air flow around a static cylinder.
Vortices start appearing after a couple of hundred steps.
"""
# from phi.flow import *  # minimal dependencies
from phi.torch.flow import *
import matplotlib.pyplot as plt
# from phi.tf.flow import *
# from phi.jax.flow import *
assert backend.default_backend().set_default_device('GPU')

SPEED = 20.
OFFSET = 1e-6
v0 = tensor([(SPEED, 0), (SPEED + OFFSET, 0)], batch("speed"), channel(vector='x,y'))

nx, ny = (256, 128)
DOMAIN = {"extrapolation": extrapolation.BOUNDARY, 'x': nx, 'y': ny, "bounds": Box(x=128, y=64)}
velocity = StaggeredGrid(v0, **DOMAIN)
noise = StaggeredGrid(Noise(), **DOMAIN)
velocity += noise
CYLINDER = Obstacle(geom.infinite_cylinder(x=15, y=32, radius=5, inf_dim=None))
BOUNDARY_MASK = StaggeredGrid(Box(x=(-INF, 0.5), y=None), velocity.extrapolation, velocity.bounds, velocity.resolution)
vort = field.curl(velocity)
pressure = None

@jit_compile  # Only for PyTorch, TensorFlow and Jax
def step(v, p, dt=.1):
    v = advect.semi_lagrangian(v, v, dt)
    v = v * (1 - BOUNDARY_MASK) + BOUNDARY_MASK * (SPEED, 0)
    return fluid.make_incompressible(v, [CYLINDER], Solve('auto', 1e-5, 0, x0=p))

diff = []
cmap="plasma"
fig, (ax1, ax2, ax3) = plt.subplots(3,1)
im1 = ax1.imshow(vort.speed[0].values.numpy('y,x'), cmap=cmap, vmin=-2*SPEED, vmax=2*SPEED)
line = ax3.semilogy([])

for _ in range(100):
    for _ in range(10):
        velocity, pressure = step(velocity, pressure)
    diff.append(field.l2_loss(velocity.speed[0] - velocity.speed[1])/math.sqrt(nx*ny*2))
    vort = field.curl(velocity)
    vort_diff = vort.speed[0].values.numpy('y,x') - vort.speed[1].values.numpy('y,x')
    im1.set_data(vort.speed[0].values.numpy('y,x'))
    ax2.imshow(vort_diff)
    ax3.semilogy(diff)
    plt.pause(0.01)
    ax2.clear()
    ax3.clear()
    
plt.figure()
plt.semilogy(diff)
plt.show()