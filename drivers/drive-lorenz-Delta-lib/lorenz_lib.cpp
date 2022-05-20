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
    out <<  -sigma, sigma,     0,
           rho - z,    -1,    -x,
                 y,     x, -beta;
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

VEC theta1(VEC u, VEC ustop, double dt, double theta=0.5, int newton_iters=10, double tol=1e-10)
{
    VEC k1 = f_lorenz(u);
    VEC rhs;
    MAT A;
    for (int i=0; i < newton_iters; i++)
    {
        rhs = ustop - u - dt*(theta*k1 + (1-theta)*f_lorenz(ustop));
        if (rhs.norm() <= tol)
        {
            return ustop;
        }
        // Newton iter
        A = MAT::Identity() - dt*(1-theta)*f_lorenz_du(ustop);
        ustop -= A.partialPivLu().solve(rhs);
    }
    return ustop;
}

MAT theta1_du(VEC u, VEC ustop, double dt, double theta=0.5)
{
    MAT rhs;
    MAT lhs;
    rhs = MAT::Identity() + dt*theta*f_lorenz_du(u);
    lhs = MAT::Identity() - dt*(1-theta)*f_lorenz_du(ustop);
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
