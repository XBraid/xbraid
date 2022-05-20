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
    if (max)
    {
        std::cout << "WARNING: max iters (" << newton_iters << ") reached in theta1!\n";
    }
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
