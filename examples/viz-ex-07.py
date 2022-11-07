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
    dim = tan.shape[-1] // 4
    tan = tan.reshape((tan.shape[0], dim, 4))

    # separate out the exponents
    exps = tan[:, :, 0]
    tan  = tan[:, :, 1:]

    # find the time-average of the exponents
    dt = 4*np.log(10.)/0.9/len(exps)
    exps /= dt
    mean_exp = exps.mean(axis=0)
    # mean_exp = exps.sum(axis=0)/(4*np.log(10.)/0.9) # assuming tstop = 4 T_\lambda.
    print(f"average Lyapunov exponents: {mean_exp}")

    # plot the trajectory and tangent vectors in 3D
    plt.figure(figsize=(8,8))
    ax = plt.axes(projection='3d')
    ax.plot3D(*u.T)
    for d in range(dim):
        ax.quiver(*u[::16].T, *tan[::16, d, :].T, color=f"C{1+d}", length=1.7, label=fr"$\lambda_{d} = {mean_exp[d]}$")
        # ax.quiver(*u.T, *tan[:, d, :].T, color=f"C{1+d}", length=1.1, label=fr"$\lambda_{d} = {exps[d]}$")
    plt.legend()
    plt.tight_layout()
    plt.show()