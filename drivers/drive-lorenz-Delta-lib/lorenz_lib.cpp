#include "lorenz_lib.hpp"

VEC f_lorenz(VEC u)
{
    double x, y, z;
    x = u[0];
    y = u[1];
    z = u[2];

    VEC out;
    out << sigma * (y - x),
        x * (rho - z) - y,
        x * y - beta * z;

    return out;
}

MAT f_lorenz_du(VEC u)
{
    double x, y, z;
    x = u[0];
    y = u[1];
    z = u[2];

    MAT out;
    out << -sigma, sigma, 0,
        rho - z, -1, -x,
        y, x, -beta;
    return out;
}

VEC euler(VEC u, double dt)
{
    return u + dt * f_lorenz(u);
}

MAT euler_du(VEC u, double dt)
{
    return Matrix3d::Identity() + dt * f_lorenz_du(u);
}

VEC theta1(const VEC u, const VEC guess, double dt, double theta, MAT *P_tan, int newton_iters, double tol)
{
    if (theta == 1.)
    {
        if (P_tan)
        {
            *P_tan = euler_du(u, dt);
        } // compute linearization
        return euler(u, dt);
    }
    VEC k1 = f_lorenz(u);
    VEC ustop = guess;
    VEC rhs;
    MAT A;
    bool max = false;

    for (int i = 0; i < newton_iters; i++)
    {
        rhs = ustop - u - dt * (theta * k1 + (1 - theta) * f_lorenz(ustop));

        // Newton iter
        if (rhs.norm() <= tol)
        {
            break;
        } // tolerance checking
        A = MAT::Identity() - dt * (1 - theta) * f_lorenz_du(ustop);
        ustop -= A.partialPivLu().solve(rhs);
        max = (i == newton_iters - 1);
    }
    // if (max)
    // {
    //     std::cout << "WARNING: max iters (" << newton_iters << ") reached in theta1!\n";
    // }
    if (P_tan)
    {
        A = MAT::Identity() - dt * (1 - theta) * f_lorenz_du(ustop);
        *P_tan = A.partialPivLu().solve(MAT::Identity() + dt * theta * f_lorenz_du(u));
    }
    return ustop;
}

MAT theta1_du(VEC u, VEC ustop, double dt, double theta)
{
    if (theta == 1.)
    {
        return euler_du(u, dt);
    }
    MAT rhs;
    MAT lhs;
    rhs = MAT::Identity() + dt * theta * f_lorenz_du(u);
    lhs = MAT::Identity() - dt * (1 - theta) * f_lorenz_du(ustop);
    return lhs.partialPivLu().solve(rhs);
}

VEC rk4(VEC u, double dt, MAT *P_tan)
{
    VEC k1, k2, k3, k4;
    VEC u1, u2, u3;

    k1 = dt * f_lorenz(u);
    u1 = u + 1. / 2. * k1;
    k2 = dt * f_lorenz(u1);
    u2 = u + 1. / 2. * k2;
    k3 = dt * f_lorenz(u2);
    u3 = u + k3;
    k4 = dt * f_lorenz(u3);

    if (P_tan)
    {
        // store K1 in passed pointer
        *P_tan = dt * f_lorenz_du(u);

        MAT &K1 = *P_tan;
        MAT K2, K3, K4;
        K2 = dt * f_lorenz_du(u1) * (MAT::Identity() + K1 / 2);
        K3 = dt * f_lorenz_du(u2) * (MAT::Identity() + K2 / 2);
        K4 = dt * f_lorenz_du(u3) * (MAT::Identity() + K3);
        *P_tan = MAT::Identity() + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
    }
    return u + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

VEC theta4(VEC u, VEC ustop, double dt, double theta, MAT *P_tan, int newton_iters, double tol)
{
    VEC k1, k2, k3, k4;
    VEC u1, u2, u3;

    k1 = dt * f_lorenz(u);
    // the matrix du1/du is also computed here, to be used to find K2 later
    u1 = theta1(u, u, dt / 2, (1. + theta) / 2, P_tan, newton_iters, tol);
    k2 = dt * f_lorenz(u1);
    u2 = u + (1. - theta) / 4 * k1 + (1. + theta) / 4 * k2;
    k3 = dt * f_lorenz(u2);
    u3 = u + k3;
    k4 = dt * f_lorenz(u3);

    if (P_tan)
    {
        // P_tan already stores part of K2 (i.e. du1/du)
        MAT &K2 = *P_tan;
        MAT K1, K3, K4;
        // MAT &tmp = K4; // just use this to store temporary values

        K1 = dt * f_lorenz_du(u);
        // tmp = dt * f_lorenz_du(u1);
        // K2 = tmp * (MAT::Identity() + (1 + theta) / 4 * K1 + (1 - theta) / 4 * tmp * K2);
        K2 = dt * f_lorenz_du(u1) * K2;
        K3 = dt * f_lorenz_du(u2) * (MAT::Identity() + (1 - theta) / 4 * K1 + (1 + theta) / 4 * K2);
        K4 = dt * f_lorenz_du(u3) * (MAT::Identity() + K3);
        *P_tan = MAT::Identity() + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
    }
    return u + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

void pack_array(std::ofstream &f, VEC u)
{
    f << u[0] << ',' << u[1] << ',' << u[2] << '\n';
}

int intpow(int base, int exp)
{
    int out = 1;
    for (int i = 0; i < exp; i++)
    {
        out *= base;
    }
    return out;
}

// int main()
// {
//     // test time-stepping
//     int nt = 1024;
//     int Tf_lyap = 2;
//     double Tf = Tf_lyap * T_lyap;
//     double dt = Tf / (nt - 1);

//     // point near the attractor
//     VEC u0{-1.8430428, -0.07036326, 23.15614636};

//     // evolve coordinates and write to file
//     std::ofstream f;
//     f.open("lorenz.csv");
//     pack_array(f, u0);
//     for (int i = 1; i < nt; i++)
//     {
//         u0 = euler(u0, dt);
//         pack_array(f, u0);
//     }
//     f.close();
//     return 0;
// }
