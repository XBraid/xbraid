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

void _setup_nxnbmat(SPMAT &out, const std::vector<std::vector<SPMAT>> blocks, const int nx)
{
    // e.g.
    // SPMAT A(Identity(nx, nx)), B(Zero(nx, nx)), C(Zero(nx, nx)), D(Identity(nx, nx));
    // _setup_nxnbmat(A, {{A, B}, {C, D}}, nx);

    int nnz = 0;
    std::vector<T> nz;
    for (auto &&blockrow : blocks)
    {
        for (auto &&block : blockrow)
        {
            nnz += block.nonZeros();
        }
    }
    out.reserve(nnz);
    nz.reserve(nnz);

    Index brows = blocks.size();
    Index bcols = blocks[0].size();
    for (Index block_i = 0; block_i < brows; block_i++)
    {
        for (Index block_j = 0; block_j < bcols; block_j++)
        {
            for (Index j = 0; j < nx; j++)
            {
                for (SPMAT::InnerIterator it(blocks[block_i][block_j], j); it; ++it)
                {
                    nz.push_back(T(it.row() + block_i * nx, it.col() + block_j * nx, it.value()));
                }
            }
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

ArrayXd _pascals_row(int n)
{
    if (n == 1)
    {
        return ArrayXd::Ones(1);
    }
    ArrayXd out = ArrayXd::Zero(n);
    ArrayXd up = _pascals_row(n - 1);
    out.head(n - 1) = up;
    out.tail(n - 1) += up;
    // std::cout << out.transpose() << '\n';
    return out;
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

SPMAT circulant_from_stencil(ArrayXd stencil, int n)
{
    return circulant_from_stencil(Stencil(stencil.data(), stencil.data() + stencil.size()), n);
}

VEC FourierMode(int wavenum, int nx, double len)
{
    if (wavenum == 0)
    {
        return VEC::Ones(nx) / (double)nx;
    }
    VEC x = VEC::LinSpaced(nx, 0., len - len / nx);
    x *= 2 * wavenum * M_PI / len;
    return Eigen::sin(x.array());
}

void setFourierMatrix(MAT &A, const int nx, const double len)
{
    for (Index j = 0; j < A.cols(); j++)
    {
        A.col(j) = FourierMode(j + 1, nx, len);
        A.col(j).normalize();
    }
}

VEC smoothed_noise(int nx, int width)
{
    VEC out(nx);
    out.setRandom();

    ArrayXd stencil = _pascals_row(2 * width - 1);
    stencil /= std::pow(2, (2 * width - 1));

    SPMAT G = circulant_from_stencil(stencil, nx);
    return G * out;
}

double inf_norm(const VEC &u)
{
    return u.lpNorm<Eigen::Infinity>();
}

KSDiscretization::KSDiscretization(int nx_, double length, Stencil d1, Stencil d2, Stencil d4)
{
    nx = nx_;
    len = length;
    dx = len / nx;
    Dx = 1 / dx * circulant_from_stencil(d1, nx);
    L = 1 / (dx * dx) * circulant_from_stencil(d2, nx);
    L += 1 / (dx * dx * dx * dx) * circulant_from_stencil(d4, nx);
}

VEC KSDiscretization::f_ks(const VEC &u) const
{
    return -L * u - u.cwiseProduct(Dx * u);
}

SPMAT KSDiscretization::f_ks_du(const VEC &u) const
{
    SPMAT U(nx, nx);
    SPMAT N(nx, nx);
    _setdiag(U, u);
    _setdiag(N, Dx * u);
    return -(L + N + U * Dx).pruned();
}

void getGuessTheta2(VEC &guess, const VEC &u, const VEC &ustop, const KSDiscretization &disc, double dt)
{
    const int stages = 2;
    int nx = u.size();
    assert(guess.size() == stages * nx);

    guess.tail(nx) = dt * disc.f_ks(ustop);
    guess.head(nx) = 2 * (ustop - u) - guess.tail(nx);
}

VEC theta2(const VEC &u, VEC &guess, const KSDiscretization &disc, double dt, double th_A, double th_B, double th_C, MAT *P_tan, int newton_iters, double tol, double *err_est)
{
    const int stages = 2;
    int nx = u.size();
    assert(guess.size() == stages * nx);

    double th_Cs, a11, a12, a21, a22;
    th_Cs = 1. - th_A - th_B - th_C;
    a11 = (th_B + th_C) / 2;
    a12 = -(th_C) / 2;
    a21 = (th_A + th_B + th_C) / 2 + th_Cs;
    a22 = (th_A + th_C) / 2;

    VEC k(guess);
    VEC p(guess);
    VEC u1(u), u2(u);

    SPMAT A(stages * nx, stages * nx);
    Eigen::UmfPackLU<SPMAT> solver; // best option I have tried so far!
    SPMAT eye(nx, nx);
    eye.setIdentity();
    VEC rhs(stages * nx);
    for (int i = 0; i < newton_iters; i++)
    {
        u1 = u + a11 * k.head(nx) + a12 * k.tail(nx);
        u2 = u + a21 * k.head(nx) + a22 * k.tail(nx);

        rhs << k.head(nx) - dt * disc.f_ks(u1),
            k.tail(nx) - dt * disc.f_ks(u2);

        if (i > 0 && inf_norm(p) <= tol)
        {
            break;
        }

        // solve block system
        SPMAT B1 = dt * disc.f_ks_du(u1);
        SPMAT B2 = dt * disc.f_ks_du(u2);
        _setup_nxnbmat(A,
                       {{eye - a11 * B1, -a12 * B1},
                        {-a21 * B2, eye - a22 * B2}},
                       nx);

        solver.compute(A);
        p = solver.solve(rhs);
        k -= p;
    }
    if (P_tan)
    {
        // here we assume that the columns of P_tan contain tangent vectors to propagate
        // and we propagate them implicitly without forming the full lin. of Phi.
        // VEC tmp(stages * nx); // just used to store solution vectors
        // VEC rhs(stages * nx);
        MAT tmp(stages * nx, P_tan->cols());
        MAT rhs(stages * nx, P_tan->cols());
        SPMAT K1 = dt * disc.f_ks_du(u1);
        SPMAT K2 = dt * disc.f_ks_du(u2);

        _setup_nxnbmat(A, {{eye - a11 * K1, -a12 * K1}, {-a21 * K2, eye - a22 * K2}}, nx);
        solver.compute(A); // does the factorization

        // solve for P_tan one column at a time
        // for (Index j = 0; j < P_tan->cols(); j++)
        // {
        //     rhs << K1 * P_tan->col(j), K2 * P_tan->col(j);
        //     tmp = solver.solve(rhs);
        //     P_tan->col(j) += (tmp.head(nx) + tmp.tail(nx)) / 2.;
        // }
        rhs << K1 * (*P_tan), K2 * (*P_tan);
        tmp = solver.solve(rhs);
        *P_tan += (tmp.topRows(nx) + tmp.bottomRows(nx)) / 2.;
    }

    // std::cout << "guess accuracy: " << (k - guess).norm() << '\n';
    guess = k;
    u2 = u + (k.head(nx) + k.tail(nx)) / 2;
    if (err_est)
    {
        // use embedded 1st order method to estimate (3rd order) local discretization error
        u1 = u + k.head(nx);
        *err_est = (u1 - u2).lpNorm<Eigen::Infinity>();
    }
    return u2;
}

void getGuessTheta4(VEC &guess, const VEC &u, const VEC &ustop, const KSDiscretization &disc, double dt, double th_A, double th_B, double th_C)
{
    const int stages = 3;
    int nx = u.size();
    assert(guess.size() == stages * nx);

    // naive guess
    guess.head(nx) = dt * disc.f_ks(u);
    guess.segment(nx, nx) = ustop - u;
    guess.tail(nx) = dt * disc.f_ks(ustop);
}

VEC theta4(const VEC &u, VEC &guess, const KSDiscretization &disc, double dt, double th_A, double th_B, double th_C, MAT *P_tan, int newton_iters, double tol)
{
    const int stages = 3;
    int nx = u.size();
    assert(guess.size() == stages * nx);
    Eigen::Matrix<double, stages, stages> A;
    double norm;

    VEC k(guess);
    VEC u1(u), u2(u), u3(u);

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

    SPMAT LHS(stages * nx, stages * nx);
    Eigen::UmfPackLU<SPMAT> solver;
    SPMAT eye(nx, nx);
    eye.setIdentity();
    VEC rhs(stages * nx);
    for (int i = 0; i < newton_iters; i++)
    {
        u1 = u + A(0, 0) * k.head(nx) + A(0, 1) * k.segment(nx, nx) + A(0, 2) * k.tail(nx);
        u2 = u + A(1, 0) * k.head(nx) + A(1, 1) * k.segment(nx, nx) + A(1, 2) * k.tail(nx);
        u3 = u + A(2, 0) * k.head(nx) + A(2, 1) * k.segment(nx, nx) + A(2, 2) * k.tail(nx);
        rhs << k.head(nx) - dt * disc.f_ks(u1),
            k.segment(nx, nx) - dt * disc.f_ks(u2),
            k.tail(nx) - dt * disc.f_ks(u3);
        norm = rhs.norm();

        if (norm <= tol)
        {
            // std::cout << i << " newton iters\n";
            break;
        }

        // block matrix
        SPMAT B1 = dt * disc.f_ks_du(u1);
        SPMAT B2 = dt * disc.f_ks_du(u2);
        SPMAT B3 = dt * disc.f_ks_du(u3);
        _setup_nxnbmat(LHS,
                       {{eye - A(0, 0) * B1, -A(0, 1) * B1, -A(0, 2) * B1},
                        {-A(1, 0) * B2, eye - A(1, 1) * B2, -A(1, 2) * B3},
                        {-A(2, 0) * B3, -A(2, 1) * B3, eye - A(2, 2) * B3}},
                       nx);

        solver.compute(LHS);
        k -= solver.solve(rhs);
    }
    if (P_tan)
    {
        VEC tmp(stages * nx); // stores solution vectors
        SPMAT K1 = dt * disc.f_ks_du(u1);
        SPMAT K2 = dt * disc.f_ks_du(u2);
        SPMAT K3 = dt * disc.f_ks_du(u3);

        _setup_nxnbmat(LHS,
                       {{eye - A(0, 0) * K1, -A(0, 1) * K1, -A(0, 2) * K1},
                        {-A(1, 0) * K2, eye - A(1, 1) * K2, -A(1, 2) * K3},
                        {-A(2, 0) * K3, -A(2, 1) * K3, eye - A(2, 2) * K3}},
                       nx);
        solver.compute(LHS); // does the factorization

        // solve for P_tan one column at a time
        for (Index j = 0; j < P_tan->cols(); j++)
        {
            rhs << K1 * P_tan->col(j), K2 * P_tan->col(j), K3 * P_tan->col(j);
            tmp = solver.solve(rhs); // does the triangular solves
            P_tan->col(j) += (tmp.head(nx) + 4 * tmp.segment(nx, nx) + tmp.tail(nx)) / 6.;
        }
    }
    std::cout << "guess accuracy: " << (k - guess).norm() << '\n';
    guess = k;
    return u + (k.head(nx) + 4 * k.segment(nx, nx) + k.tail(nx)) / 6.;
}

void GramSchmidt(MAT &A)
{ // orthonormalize A in place using modified gram-schmidt
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

void pack_array(std::ofstream &f, const VEC &u)
{
    for (Index i = 0; i < u.size() - 1; i++)
    {
        f << u(i) << ',';
    }
    f << u(Eigen::last) << '\n';
}

void pack_darray(std::ofstream &f, const MAT &u)
{
    using Eigen::last;
    for (Index j = 0; j < u.cols() - 1; j++)
    {
        for (Index i = 0; i < u.rows(); i++)
        {
            f << u(i, j) << ',';
        }
        // f << u(j, last) << ';';
    }
    pack_array(f, u.col(u.cols() - 1));
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

void bf_pack_help(double *buf, const MAT &u, const size_t obj_size, size_t &bf_size)
{
    size_t i_bf = bf_size;
    for (size_t i = 0; i < obj_size; i++)
    {
        buf[i_bf] = u(i);
        i_bf++;
    }
    bf_size += obj_size;
}

void bf_pack_help(double *buf, const VEC &u, const size_t obj_size, size_t &bf_size)
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