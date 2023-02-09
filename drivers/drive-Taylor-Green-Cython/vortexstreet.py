
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

SPEED = vis.control(20.)
velocity = StaggeredGrid((SPEED, 0), extrapolation.BOUNDARY, x=256, y=128, bounds=Box(x=128, y=64))
CYLINDER = Obstacle(geom.infinite_cylinder(x=15, y=32, radius=5, inf_dim=None))
BOUNDARY_MASK = StaggeredGrid(Box(x=(-INF, 0.5), y=None), velocity.extrapolation, velocity.bounds, velocity.resolution)
pressure = None


@jit_compile  # Only for PyTorch, TensorFlow and Jax
def step(v, p, dt=1.):
    v = advect.semi_lagrangian(v, v, dt)
    v = v * (1 - BOUNDARY_MASK) + BOUNDARY_MASK * (SPEED, 0)
    return fluid.make_incompressible(v, [CYLINDER], Solve('auto', 1e-5, 0, x0=p))


for _ in range(1000):
    velocity, pressure = step(velocity, pressure)
    vorticity = field.curl(velocity)
    plt.imshow(vorticity.values)
    plt.draw()
    plt.pause(0.01)
    plt.clf()
    