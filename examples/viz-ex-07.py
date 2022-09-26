import numpy as np
import matplotlib.pyplot as plt

def read_data(fname):
    u = []
    with open(fname, 'r') as f:
        for line in f:
            u.append(np.array(line.split(' ')[1:], dtype=np.double))
    return np.stack(u)

if __name__=="__main__":
    u = read_data("ex-07.out")
    tan = read_data("ex-07-lv.out")

    # reshape the tangent vectors
    dim = tan.shape[-1] // 3
    tan = tan.reshape((tan.shape[0], dim, 3))

    # plot the trajectory and tangent vectors in 3D
    ax = plt.axes(projection='3d')
    ax.plot3D(*u.T)
    for d in range(dim):
        ax.quiver(*u[::16].T, *tan[::16, d, :].T, color=f"C{1+d}")
    plt.show()