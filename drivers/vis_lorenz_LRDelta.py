import numpy as np
import matplotlib.pyplot as plt

try:
    plt.style.use('https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pacoty.mplstyle')
except:
    pass

def read_csv(fname):
    u = []
    with open(fname, 'r') as f:
        for line in f:
            u.append(np.array(line.split(','), dtype=np.double))
    return np.stack(u)

if __name__=="__main__":
    u = read_csv("lorenz-LRDelta.out")
    tan = read_csv("lorenz-LRDelta-lv.out")

    # reshape the tangent vectors
    dim = tan.shape[-1] // 3
    tan = tan.reshape((tan.shape[0], dim, 3))

    ax = plt.axes(projection='3d')
    ax.plot3D(*u.T)
    for d in range(dim):
        ax.quiver(*u.T, *tan[:, d, :].T, color=f"C{1+d}")
    plt.show()