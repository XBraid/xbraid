from math import ceil
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from vis_lorenz_LRDelta import read_csv

def baseline(u):
    v = np.copy(u)
    for i in range(len(v)):
        v[i] -= np.mean(v[i])
    return v

def int_v(vec):
    out = np.zeros_like(vec)
    for i, v in enumerate(vec):
        out[i] = out[i-1] + vec[i]

    return out

def integrate(u):
    out = np.zeros_like(u)
    for i, vec in enumerate(u):
        out[i] = int_v(vec)
        
    return baseline(out)

if __name__=="__main__":
    u = read_csv("drive-ks.out")

    nt, nx = u.shape
    cf = ceil(len(u)/nt)
    l = nx
    T_lyap = np.log(10)/0.1
    Tf = 4*T_lyap

    # plot trajectories
    ratio = (5, 5)
    cmap = "plasma"
    # cmap = "cool_r"
    stride = max(nt//( nx*ratio[0]//ratio[1] )//cf, 1)
    extent = (0, l, Tf/T_lyap, 0)
    fig, axs = plt.subplots(1, 1, figsize=(ratio[0], ratio[1]))
    interp = "bilinear"
    axs.imshow(integrate(u[::cf*stride]), cmap=cmap, extent=extent, aspect="auto", interpolation=interp)
    axs.set_ylabel("$t_\lambda$")
    fig.tight_layout()
    axs.set_xlabel("u", fontsize=32)
    axs.grid()
    plt.show()

    # error norm + estimate lyapunov exponent
    # e_th = np.linalg.norm(u[cf::cf]-v[1:], axis=1)
    # # e = np.linalg.norm(u[::cf]-v, axis=1)
    # tc = np.linspace(0, Tf, num=len(e_th))
    # res = linregress(tc, np.log(e_th))
    # print(f"lambda_max = {res.slope}")
    # print(f"t_lambda = {np.log(10)/res.slope}")
    # t_lyap = tc*res.slope
    # plt.figure(figsize=(10, 8))
    # plt.semilogy(t_lyap/np.log(10), e_th, label=r"$\theta$ method")
    # plt.semilogy(t_lyap/np.log(10), np.exp(t_lyap + res.intercept), linestyle="dashed", label=f"$\lambda = {res.slope}$")
    # # plt.semilogy(t, e)
    # # plt.semilogy(t, np.exp(res.slope*t + res.intercept), 'r--')
    # plt.legend()
    # plt.show()