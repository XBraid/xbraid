#include "lorenz_lib.hpp"

VEC f_lorenz(VEC u)
{
    double x, y, z;
    x = u[0];
    y = u[1];
    z = u[2];

    VEC out{0., 0., 0.};
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
