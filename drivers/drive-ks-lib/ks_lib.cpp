#include "ks_lib.hpp"


// internal functions

void _setdiag(SPMAT &A, const VEC &u)
{
    std::vector<T> coefficients;
    coefficients.reserve(u.size());
    for (Index i = 0; i < u.size(); i++)
    {
        coefficients.push_back(T(i, i, u[i]));
    }
    A.setFromTriplets(coefficients.begin(), coefficients.end());
}

void _setup_bmat(SPMAT &out, const SPMAT &A, const SPMAT &B, const SPMAT &C, const SPMAT &D, const int nx)
{
    int nnz = A.nonZeros() + B.nonZeros() + C.nonZeros() + D.nonZeros();
    out.reserve(nnz);

    std::vector<T> nz;
    nz.reserve(nnz);

    // get nonzeros from blocks
    for (Index j = 0; j < nx; j++) // column major
    {
        for (SPMAT::InnerIterator it(A, j); it; ++it)
        {
            nz.push_back(T(it.row(), it.col(), it.value()));
        }

        for (SPMAT::InnerIterator it(B, j); it; ++it)
        {
            nz.push_back(T(it.row(), it.col() + nx, it.value()));
        }

        for (SPMAT::InnerIterator it(C, j); it; ++it)
        {
            nz.push_back(T(it.row() + nx, it.col(), it.value()));
        }

        for (SPMAT::InnerIterator it(D, j); it; ++it)
        {
            nz.push_back(T(it.row() + nx, it.col() + nx, it.value()));
        }
    }

    out.setFromTriplets(nz.begin(), nz.end());
    out.makeCompressed();
}

void _setup_col_bmat(SPMAT &out, const SPMAT &A1, const SPMAT &A2, const int nx)
{
    int nnz = A1.nonZeros() + A2.nonZeros();
    out.reserve(nnz);

    std::vector<T> nz;
    nz.reserve(nnz);

    // get nonzeros from blocks
    for (Index j = 0; j < nx; j++) // column major
    {
        for (SPMAT::InnerIterator it(A1, j); it; ++it)
        {
            nz.push_back(T(it.row(), it.col(), it.value()));
        }

        for (SPMAT::InnerIterator it(A2, j); it; ++it)
        {
            nz.push_back(T(it.row() + nx, it.col(), it.value()));
        }
    }

    out.setFromTriplets(nz.begin(), nz.end());
    out.makeCompressed();
}

Index _clip(Index i, int n)
{
    if (i < 0)
    {
        return i + n;
    }
    return i % n;
}

// exported functions
SPMAT circulant_from_stencil(Stencil stencil, int n)
{
    // e.g.
    // SPMAT Dxx = circulant_from_stencil({1, -2, 1}, 16)
    // yields a circulant matrix implementing the FD stencil {1, -2, 1}
    SPMAT A(n, n);
    int sten_width = stencil.size();
    int ord = sten_width / 2;

    // make triplets (int row, int col, double entry)
    std::vector<T> coefficients;
    coefficients.reserve(n * sten_width);
    for (Index i = 0; i < n; i++)
    {
        for (Index j_sten = 0; j_sten < sten_width; j_sten++)
        {
            Index j = _clip(i + j_sten - ord, n);
            coefficients.push_back(T(i, j, stencil[j_sten]));
        }
    }

    A.setFromTriplets(coefficients.begin(), coefficients.end());
    return A;
}

KSDiscretization::KSDiscretization(int nx_, double length, Stencil d1, Stencil d2, Stencil d4)
{
    nx = nx_;
    dx = length/nx;
    Dx = 1/dx * circulant_from_stencil(d1, nx);
    L  = 1/(dx*dx) * circulant_from_stencil(d2, nx);
    L += 1/(dx*dx*dx*dx) * circulant_from_stencil(d4, nx);
}

VEC KSDiscretization::f_ks(const VEC& u) const
{
    return - L * u - u.cwiseProduct(Dx * u);
}

SPMAT KSDiscretization::f_ks_du(const VEC& u) const
{
    SPMAT U(nx, nx);
    SPMAT N(nx, nx);
    _setdiag(U, u);
    _setdiag(N, Dx * u);
    return -(L + N + U * Dx).pruned();
}

VEC theta2(const VEC &u, const KSDiscretization& disc, double dt, double th_A, double th_B, double th_C, MAT *P_tan, int newton_iters, double tol)
{
    const int stages = 2;
    int nx = u.size();
    double th_Cs, a11, a12, a21, a22;
    th_Cs = 1. - th_A - th_B - th_C;
    a11 = (th_B + th_C) / 2;
    a12 = -(th_C) / 2;
    a21 = (th_A + th_B + th_C) / 2 + th_Cs;
    a22 = (th_A + th_C) / 2;

    VEC k(VEC::Zero(stages * nx)); // safest initial guess
    VEC u1(u), u2(u);

    SPMAT A(stages * nx, stages * nx);
    Eigen::UmfPackLU<SPMAT> solver; // best option I have tried so far!
    SPMAT eye(nx, nx);
    eye.setIdentity();
    VEC rhs(VEC::Zero(stages * nx));
    for (size_t i = 0; i < newton_iters; i++)
    {
        u1 = u + a11 * k.head(nx) + a12 * k.tail(nx);
        u2 = u + a21 * k.head(nx) + a22 * k.tail(nx);

        rhs << k.head(nx) - dt * disc.f_ks(u1),
               k.tail(nx) - dt * disc.f_ks(u2);

        if (rhs.norm() <= tol)
        {
            break;
        }

        // solve block system
        SPMAT B1 = dt * disc.f_ks_du(u1);
        SPMAT B2 = dt * disc.f_ks_du(u2);
        _setup_bmat(A, eye - a11 * B1, -a12 * B1, -a21 * B2, eye - a22 * B2, nx);

        solver.compute(A);
        k -= solver.solve(rhs);
    }
    if (P_tan)
    {
        VEC tmp(stages * nx); // just used to store solution vectors
        SPMAT RHS(stages * nx, nx);
        SPMAT K1 = dt * disc.f_ks_du(u1);
        SPMAT K2 = dt * disc.f_ks_du(u2);

        _setup_col_bmat(RHS, K1, K2, nx);
        _setup_bmat(A, eye - a11 * K1, -a12 * K1, -a21 * K2, eye - a22 * K2, nx);
        solver.compute(A); // does the factorization

        // solve for P_tan one column at a time (reusing factorization)
        for (Index j = 0; j < nx; j++)
        {
            tmp = solver.solve(RHS.col(j));
            P_tan->col(j) = eye.col(j) + (tmp.head(nx) + tmp.tail(nx)) / 2.;
        }
    }
    return u + (k.head(nx) + k.tail(nx)) / 2;
}

void pack_array(std::ofstream &f, const VEC &u)
{
    for (Index i = 0; i < u.size() - 1; i++)
    {
        f << u(i) << ',';
    }
    f << u(Eigen::last) << '\n';
}

int intpow(const int base, const int exp)
{
    int out = 1;
    for (int i = 0; i < exp; i++)
    {
        out *= base;
    }
    return out;
}

void bf_pack_help(double *buf, const MAT& u, const size_t obj_size, size_t &bf_size)
{
   size_t i_bf = bf_size;
   for (size_t i = 0; i < obj_size; i++)
   {
      buf[i_bf] = u(i);
      i_bf++;
   }
   bf_size += obj_size;
}

void bf_pack_help(double *buf, const VEC& u, const size_t obj_size, size_t &bf_size)
{
   size_t i_bf = bf_size;
   for (size_t i = 0; i < obj_size; i++)
   {
      buf[i_bf] = u(i);
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