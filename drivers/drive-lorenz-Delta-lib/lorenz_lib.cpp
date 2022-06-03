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
        }
        return euler(u, dt);
    }
    VEC k1 = f_lorenz(u);
    VEC ustop = guess;
    VEC rhs;
    MAT A;
    // bool max = false;

    for (int i = 0; i < newton_iters; i++)
    {
        rhs = ustop - u - dt * (theta * k1 + (1 - theta) * f_lorenz(ustop));

        // Newton iter
        if (rhs.norm() <= tol)
        {
            // std::cout << "newton iters: " << i << '\n';
            break;
        } // tolerance checking
        A = MAT::Identity() - dt * (1 - theta) * f_lorenz_du(ustop);
        ustop -= A.partialPivLu().solve(rhs);
        // max = (i == newton_iters - 1);
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

VEC crank_nicolson(VEC u, VEC ustop, double dt, MAT *P_tan, int newton_iters, double tol)
{
    return theta1(u, ustop, dt, 0.5, P_tan, newton_iters, tol);
}

VEC theta2(VEC u, VEC ustop, double dt, double th_A, double th_B, double th_C, MAT *P_tan, int newton_iters, double tol)
{
    double th_Cs, a11, a12, a21, a22;
    Eigen::Vector<double, 2 * DIM> k;
    VEC u1, u2;
    Eigen::Matrix<double, 2 * DIM, 2 * DIM> A;

    th_Cs = 1. - th_A - th_B - th_C;
    a11 = (th_B + th_C) / 2.;
    a12 = -(th_C) / 2.;
    a21 = (th_A + th_B + th_C) / 2. + th_Cs;
    a22 = (th_A + th_C) / 2.;

    // TODO: find a better initial guess
    k.setZero(); // safest initial guess

    for (int i = 0; i < newton_iters; i++)
    {
        Eigen::Vector<double, 2 * DIM> rhs;
        u1 = u + a11 * k.head(DIM) + a12 * k.tail(DIM);
        u2 = u + a21 * k.head(DIM) + a22 * k.tail(DIM);

        rhs << k.head(DIM) - dt * f_lorenz(u1),
            k.tail(DIM) - dt * f_lorenz(u2);

        if (rhs.norm() <= tol)
        {
            break;
        }

        // block matrix
        MAT B1 = f_lorenz_du(u1);
        MAT B2 = f_lorenz_du(u2);
        A << MAT::Identity() - dt * a11 * B1, -dt * a12 * B1,
            -dt * a21 * B2, MAT::Identity() - dt * a22 * B2;

        k -= A.partialPivLu().solve(rhs);
    }
    if (P_tan)
    {
        Eigen::Matrix<double, 2 * DIM, DIM> RHS;
        MAT K1 = dt * f_lorenz_du(u1);
        MAT K2 = dt * f_lorenz_du(u2);

        RHS << K1, K2;
        A << MAT::Identity() - a11 * K1, -a12 * K1,
            -a21 * K2, MAT::Identity() - a22 * K2;
        RHS = A.partialPivLu().solve(RHS);
        K1 = RHS.topRows(DIM);
        K2 = RHS.bottomRows(DIM);

        *P_tan = MAT::Identity() + (K1 + K2) / 2.;
    }
    return u + (k.head(DIM) + k.tail(DIM)) / 2.;
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

VEC theta4(VEC u, VEC &uguess, double dt, double theta, MAT *P_tan, int newton_iters, double tol)
{
    VEC k1, k2, k3, k4;
    VEC u1, u2, u3;

    k1 = dt * f_lorenz(u);

    // get an initial guess for u1
    // k2 = uguess - u; // this is already 2nd order accurate!
    // u1 = u + (1. + theta) / 4 * k1 + (1. - theta) / 4 * k2;
    // or, use the previous computed value as the guess
    if (uguess.isZero())
    {
        uguess = u; // this is more stable when no good guess is available
    }
    u1 = uguess;

    // the matrix du1/du is also computed here, to be used to find K2 later
    u1 = theta1(u, u1, dt / 2, (1. + theta) / 2, P_tan, newton_iters, tol);
    k2 = dt * f_lorenz(u1);
    u2 = u + (1. - theta) / 4 * k1 + (1. + theta) / 4 * k2;
    k3 = dt * f_lorenz(u2);
    u3 = u + k3;
    k4 = dt * f_lorenz(u3);

    // std::cout << "guess accuracy: " << (u1 - uguess).norm() << '\n';
    uguess = u1; // save the value of u1 for next time

    if (P_tan)
    {
        // P_tan already stores part of K2 (i.e. du1/du)
        MAT &K2 = *P_tan;
        MAT K1, K3, K4;

        K1 = dt * f_lorenz_du(u);
        K2 = dt * f_lorenz_du(u1) * K2;
        K3 = dt * f_lorenz_du(u2) * (MAT::Identity() + (1 - theta) / 4 * K1 + (1 + theta) / 4 * K2);
        K4 = dt * f_lorenz_du(u3) * (MAT::Identity() + K3);
        *P_tan = MAT::Identity() + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
    }
    return u + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

Eigen::Vector<double, 3 * DIM> lobatto3rhs(VEC u, Eigen::Vector<double, 3 * DIM> k, double dt, Eigen::Matrix<double, 3, 3> A, VEC &u1, VEC &u2, VEC &u3)
{
    Eigen::Vector<double, 3 * DIM> rhs;
    u1 = u + A(0, 0) * k.head(DIM) + A(0, 1) * k.segment(DIM, DIM) + A(0, 2) * k.tail(DIM);
    u2 = u + A(1, 0) * k.head(DIM) + A(1, 1) * k.segment(DIM, DIM) + A(1, 2) * k.tail(DIM);
    u3 = u + A(2, 0) * k.head(DIM) + A(2, 1) * k.segment(DIM, DIM) + A(2, 2) * k.tail(DIM);

    rhs << k.head(DIM) - dt * f_lorenz(u1),
        k.segment(DIM, DIM) - dt * f_lorenz(u2),
        k.tail(DIM) - dt * f_lorenz(u3);
    return rhs;
}

double line_search(double &phi_0, Eigen::Vector<double, 3 * DIM> p, Eigen::Vector<double, 3 * DIM> &k, Eigen::Vector<double, 3 * DIM> &rhs, VEC u, double dt, Eigen::Matrix<double, 3, 3> A, VEC &u1, VEC &u2, VEC &u3)
{
    const double c1 = 0.2;
    double gam1 = 0.1;
    double gam2 = 0.5;
    double alpha = 1.;
    double phi_1, phi_prime;
    Eigen::Vector<double, 3 * DIM> k_new, rhs_new;

    k_new = k + alpha * p;
    rhs_new = lobatto3rhs(u, k_new, dt, A, u1, u2, u3);
    phi_1 = rhs_new.norm();
    phi_prime = rhs.dot(p);

    int maxiters = 4;
    int i = 0;

    // Wolf-Armijo condition
    while ((i < maxiters) && (phi_1 * phi_1 > phi_0 * phi_0 + 2 * c1 * alpha * phi_prime))
    {
        // approximately optimize alpha using quadratic approximation to phi
        double a_new = -alpha * phi_prime / (phi_1 * phi_1 - phi_0 * phi_0 - 2 * phi_prime * alpha);
        if (a_new < gam1 * alpha)
        {
            alpha = gam1 * alpha;
        }
        else
        {
            alpha = a_new;
        }

        k_new = k + alpha * p;
        rhs_new = lobatto3rhs(u, k_new, dt, A, u1, u2, u3);
        phi_1 = rhs_new.norm();
        i++;
    }

    k = k_new;
    rhs = rhs_new;
    phi_0 = phi_1;
    return alpha;
}

VEC thetaLobatto3(VEC u, VEC ustop, double dt, double th_A, double th_B, double th_C, MAT *P_tan, int newton_iters, double tol, bool linesearch)
{
    const int stages = 3;
    Eigen::Vector<double, stages * DIM> k, p;
    VEC u1, u2, u3;
    Eigen::Matrix<double, stages * DIM, stages * DIM> LHS;
    Eigen::Matrix<double, stages, stages> A;
    double norm;

    // set up coefficient matrix
    {
        double a11 = (th_B + th_C) / 6;
        double a12 = -th_B / 6 - th_C / 3;
        double a13 = th_C / 6;

        double a21 = 1 / 4. - th_A / 24 - (th_B + th_C) / 12;
        double a22 = 1 / 4. + (th_A + th_B) / 12 + th_C / 6;
        double a23 = -th_A / 24 - th_C / 12;

        double a31 = (th_A + th_B + th_C) / 6;
        double a32 = 1 - (th_A + th_C) / 3 - th_B / 6;
        double a33 = (th_A + th_C) / 6;
        A << a11, a12, a13,
            a21, a22, a23,
            a31, a32, a33;
    }

    k.setZero(); // safest initial guess?
    p.setZero(); // search direction for line search

    bool max = false;
    for (int i = 0; i < newton_iters; i++)
    {
        Eigen::Vector<double, stages * DIM> rhs;

        // backtracking linesearch algorithm
        if (linesearch && i > 0)
        {
            // this function finds and returns the near-optimal step size,
            // and also overwrites k, rhs, u1, u2, u3 and norm with the latest values
            double alph = line_search(norm, p, k, rhs, u, dt, A, u1, u2, u3);
            // std::cout << "line search step size: " << alph << '\n';
        }
        else
        {
            k += p;
            rhs = lobatto3rhs(u, k, dt, A, u1, u2, u3);
            norm = rhs.norm();
        }

        if (norm <= tol)
        {
            break;
        }
        max = (i == newton_iters - 1);

        // block matrix
        MAT B1 = dt * f_lorenz_du(u1);
        MAT B2 = dt * f_lorenz_du(u2);
        MAT B3 = dt * f_lorenz_du(u3);
        LHS << MAT::Identity() - A(0, 0) * B1, -A(0, 1) * B1, -A(0, 2) * B1,
            -A(1, 0) * B2, MAT::Identity() - A(1, 1) * B2, -A(1, 2) * B3,
            -A(2, 0) * B3, -A(2, 1) * B3, MAT::Identity() - A(2, 2) * B3;

        // this is the search direction given by Newton's method
        p = -LHS.partialPivLu().solve(rhs);
    }
    if (P_tan)
    {
        if (max)
        {
            u1 = u + A(0, 0) * k.head(DIM) + A(0, 1) * k.segment(DIM, DIM) + A(0, 2) * k.tail(DIM);
            u2 = u + A(1, 0) * k.head(DIM) + A(1, 1) * k.segment(DIM, DIM) + A(1, 2) * k.tail(DIM);
            u3 = u + A(2, 0) * k.head(DIM) + A(2, 1) * k.segment(DIM, DIM) + A(2, 2) * k.tail(DIM);
        }
        Eigen::Matrix<double, stages * DIM, DIM> RHS;
        MAT K1 = dt * f_lorenz_du(u1);
        MAT K2 = dt * f_lorenz_du(u2);
        MAT K3 = dt * f_lorenz_du(u3);

        RHS << K1, K2, K3;
        LHS << MAT::Identity() - A(0, 0) * K1, -A(0, 1) * K1, -A(0, 2) * K1,
            -A(1, 0) * K2, MAT::Identity() - A(1, 1) * K2, -A(1, 2) * K3,
            -A(2, 0) * K3, -A(2, 1) * K3, MAT::Identity() - A(2, 2) * K3;
        RHS = LHS.partialPivLu().solve(RHS);
        K1 = RHS.topRows(DIM);
        K2 = RHS.middleRows(DIM, DIM);
        K3 = RHS.bottomRows(DIM);

        *P_tan = MAT::Identity() + (K1 + 4 * K2 + K3) / 6.;
    }
    return u + (k.head(DIM) + 4 * k.segment(DIM, DIM) + k.tail(DIM)) / 6.;
}
// VEC thetaLobatto3(VEC u, VEC ustop, double dt, double th_A, double th_B, double th_C, MAT *P_tan, int newton_iters, double tol)
// {
//     const int stages = 3;
//     double a11, a12, a13, a21, a22, a23, a31, a32, a33;
//     Eigen::Vector<double, stages * DIM> k;
//     VEC u1, u2, u3;
//     Eigen::Matrix<double, stages * DIM, stages * DIM> A;

//     a11 = (th_B + th_C) / 6;
//     a12 = -th_B / 6 - th_C / 3;
//     a13 = th_C / 6;

//     a21 = 1 / 4. - th_A / 24 - (th_B + th_C) / 12;
//     a22 = 1 / 4. + (th_A + th_B) / 12 + th_C / 6;
//     a23 = -th_A / 24 - th_C / 12;

//     a31 = (th_A + th_B + th_C) / 6;
//     a32 = 1 - (th_A + th_C) / 3 - th_B / 6;
//     a33 = (th_A + th_C) / 6;

//     k.setZero(); // safest initial guess?

//     bool max = false;
//     for (int i = 0; i < newton_iters; i++)
//     {
//         Eigen::Vector<double, stages * DIM> rhs;
//         u1 = u + a11 * k.head(DIM) + a12 * k.segment(DIM, DIM) + a13 * k.tail(DIM);
//         u2 = u + a21 * k.head(DIM) + a22 * k.segment(DIM, DIM) + a23 * k.tail(DIM);
//         u3 = u + a31 * k.head(DIM) + a32 * k.segment(DIM, DIM) + a33 * k.tail(DIM);

//         rhs << k.head(DIM) - dt * f_lorenz(u1),
//                k.segment(DIM, DIM) - dt * f_lorenz(u2),
//                k.tail(DIM) - dt * f_lorenz(u3);

//         if (rhs.norm() <= tol)
//         {
//             break;
//         }

//         // block matrix
//         MAT B1 = dt * f_lorenz_du(u1);
//         MAT B2 = dt * f_lorenz_du(u2);
//         MAT B3 = dt * f_lorenz_du(u3);
//         A << MAT::Identity() - a11 * B1, -a12 * B1, -a13 * B1,
//             -a21 * B2, MAT::Identity() - a22 * B2, -a23 * B3,
//             -a31 * B3, -a32 * B3, MAT::Identity() - a33 * B3;

//         k -= A.partialPivLu().solve(rhs);
//         max = (i == newton_iters - 1);
//     }
//     if (P_tan)
//     {
//         if (max)
//         {
//             u1 = u + a11 * k.head(DIM) + a12 * k.segment(DIM, DIM) + a13 * k.tail(DIM);
//             u2 = u + a21 * k.head(DIM) + a22 * k.segment(DIM, DIM) + a23 * k.tail(DIM);
//             u3 = u + a31 * k.head(DIM) + a32 * k.segment(DIM, DIM) + a33 * k.tail(DIM);
//         }
//         Eigen::Matrix<double, stages * DIM, DIM> RHS;
//         MAT K1 = dt * f_lorenz_du(u1);
//         MAT K2 = dt * f_lorenz_du(u2);
//         MAT K3 = dt * f_lorenz_du(u3);

//         RHS << K1, K2, K3;
//         A << MAT::Identity() - a11 * K1, -a12 * K1, -a13 * K1,
//             -a21 * K2, MAT::Identity() - a22 * K2, -a23 * K3,
//             -a31 * K3, -a32 * K3, MAT::Identity() - a33 * K3;
//         RHS = A.partialPivLu().solve(RHS);
//         K1 = RHS.topRows(DIM);
//         K2 = RHS.middleRows(DIM, DIM);
//         K3 = RHS.bottomRows(DIM);

//         *P_tan = MAT::Identity() + (K1 + 4*K2 + K3) / 6.;
//     }
//     return u + (k.head(DIM) + 4*k.segment(DIM, DIM) + k.tail(DIM)) / 6.;
// }

void GramSchmidt(COL_MAT &A)
{ // orthonormalize A in place using gram-schmidt
    // normalize the first column
    A.col(0).normalize();

    for (Eigen::Index i = 1; i < A.cols(); i++)
    {
        // subtract the orthogonal projection of the columns to the left
        for (Eigen::Index j = 0; j < i; j++)
        {
            A.col(i) -= A.col(i).dot(A.col(j)) * A.col(j);
        }

        // normalize this column
        A.col(i).normalize();
    }
}

void pack_array(std::ofstream &f, VEC u)
{
    f << u[0] << ',' << u[1] << ',' << u[2] << '\n';
}

void pack_darray(std::ofstream &f, COL_MAT u)
{
    f << u(0, 0) << ',' << u(1, 0) << ',' << u(2, 0);
    for (Eigen::Index i = 1; i < u.cols(); i++)
    {
        f << ',' << u(0, i) << ',' << u(1, i) << ',' << u(2, i);
    }
    f << '\n';
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

void bf_pack_help(double *buf, const COL_MAT u, const size_t obj_size, size_t &bf_size)
{
   size_t i_bf = bf_size;
   for (size_t i = 0; i < obj_size; i++)
   {
      buf[i_bf] = u(i);
      i_bf++;
   }
   bf_size += obj_size;
}

void bf_pack_help(double *buf, const VEC u, const size_t obj_size, size_t &bf_size)
{
   size_t i_bf = bf_size;
   for (size_t i = 0; i < obj_size; i++)
   {
      buf[i_bf] = u(i);
      i_bf++;
   }
   bf_size += obj_size;
}

void bf_pack_help(double *buf, const MAT u, const size_t obj_size, size_t &bf_size)
{
   size_t i_bf = bf_size;
   for (size_t i = 0; i < obj_size; i++)
   {
      buf[i_bf] = u(i);
      i_bf++;
   }
   bf_size += obj_size;
}

void bf_unpack_help(double *buf, COL_MAT &u, const size_t obj_size, size_t &bf_size)
{
   size_t i_bf = bf_size;
   for (size_t i = 0; i < obj_size; i++)
   {
      u(i) = buf[i_bf];
      i_bf++;
   }
   bf_size += obj_size;
}

void bf_unpack_help(double *buf, VEC &u, const size_t obj_size, size_t &bf_size)
{
   size_t i_bf = bf_size;
   for (size_t i = 0; i < obj_size; i++)
   {
      u(i) = buf[i_bf];
      i_bf++;
   }
   bf_size += obj_size;
}

void bf_unpack_help(double *buf, MAT &u, const size_t obj_size, size_t &bf_size)
{
   size_t i_bf = bf_size;
   for (size_t i = 0; i < obj_size; i++)
   {
      u(i) = buf[i_bf];
      i_bf++;
   }
   bf_size += obj_size;
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
