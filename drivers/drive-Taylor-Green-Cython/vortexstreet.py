
"""Karman Vortex Street
Air flow around a static cylinder.
Vortices start appearing after a couple of hundred steps.
"""
# from phi.flow import *  # minimal dependencies
from phi.torch.flow import *
import matplotlib.pyplot as plt
from progressbar import progressbar
# from phi.tf.flow import *
# from phi.jax.flow import *
# assert backend.default_backend().set_default_device('GPU')

SPEED = 60.
OFFSET = 1e-6
v0 = tensor([(SPEED, 0), (SPEED, 0)], batch("init"), channel(vector='x,y'))

nx, ny = (128, 64)
DOMAIN = {"extrapolation": extrapolation.BOUNDARY, 'x': nx, 'y': ny, "bounds": Box(x=40, y=20)}
velocity = StaggeredGrid(v0, **DOMAIN)
noise_both = Noise().at(velocity)
noise_offset = OFFSET * StaggeredGrid(math.random_normal(velocity), **DOMAIN)
velocity += noise_both
CYLINDER = Obstacle(geom.infinite_cylinder(x=10, y=10, radius=3, inf_dim=None))
BOX = Obstacle(geom.Cuboid(x=10, y=10))
BOUNDARY_MASK = StaggeredGrid(Box(x=(-INF, 0.5), y=None), velocity.extrapolation, velocity.bounds, velocity.resolution)
vort = field.curl(velocity.with_extrapolation(extrapolation.BOUNDARY))
pressure = None


@jit_compile  # Only for PyTorch, TensorFlow and Jax
def step(v, p, dt=.05):
    v = advect.semi_lagrangian(v, v, dt)
    v = v * (1 - BOUNDARY_MASK) + BOUNDARY_MASK * (SPEED, 0)
    v, p = fluid.make_incompressible(v, [CYLINDER], Solve('auto', 1e-5, 0, x0=p))
    # v, p = fluid.make_incompressible(v, [BOX], Solve('auto', 1e-5, 0, x0=p))
    return v, p

# reach mean flow
print("stepping to mean flow")
for _ in progressbar(range(100)):
    velocity, pressure = step(velocity, pressure)
    
velocity += noise_offset
diff = []
cmap="plasma"
fig, (ax1, ax2, ax3) = plt.subplots(3,1)
im1 = ax1.imshow(vort.init[0].values.numpy('y,x'), cmap=cmap, vmin=-2*SPEED, vmax=2*SPEED)
line = ax3.semilogy([])

for _ in range(200):
    for _ in range(1):
        velocity, pressure = step(velocity, pressure)
    diff.append(field.l2_loss(velocity.init[0] - velocity.init[1])/math.sqrt(nx*ny*2))
    vort = field.curl(velocity.with_extrapolation(extrapolation.BOUNDARY))
    vort_diff = vort.init[0].values.numpy('y,x') - vort.init[1].values.numpy('y,x')
    im1.set_data(vort.init[0].values.numpy('y,x'))
    ax2.imshow(vort_diff)
    ax3.semilogy(diff)
    plt.pause(0.01)
    ax2.clear()
    ax3.clear()
    
plt.figure()
plt.semilogy(diff)
plt.show()