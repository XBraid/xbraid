#ifndef KS_LIB
#define KS_LIB

// std lib
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>

// linear algebra includes
#include "../Eigen/Dense"
#include "../Eigen/Sparse"
#include "../Eigen/UmfPackSupport" // need to install SuiteSparse

// want the flexibility of dynamic arrays to change problem size at run-time
using Eigen::ArrayXd;
using Eigen::Index;
typedef Eigen::VectorXd VEC;
typedef Eigen::MatrixXd MAT;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SPMAT;
typedef Eigen::Triplet<double> T;
typedef std::vector<double> Stencil;

// Lyapunov time
const double T_lyap{log(10) / 0.1};

// exported functions
class KSDiscretization
{
public:
    SPMAT L;
    SPMAT Dx;
    int nx;
    double len;
    double dx;

    // constructor
    KSDiscretization() : L(SPMAT()), Dx(SPMAT()), nx(0), len(0.), dx(0.) {} // Never use this.
    KSDiscretization(int nx_, double length, Stencil d1, Stencil d2, Stencil d4);

    // implements discretization:
    VEC f_ks(const VEC& u) const;
    SPMAT f_ks_du(const VEC& u) const;
};

// KS helpers
SPMAT circulant_from_stencil(std::vector<double> stencil, int n);
SPMAT circulant_from_stencil(ArrayXd stencil, int n);
VEC smoothed_noise(int nx, int width);
VEC FourierMode(int wavenum, int nx);
void setFourierMatrix(MAT &A);
void setFourierMatrix(MAT &A, int rows, int cols);
void GramSchmidt(MAT& A);

// time steppers
VEC theta2(const VEC &u, const VEC &ustop, const KSDiscretization& disc, double dt, double th_A = 0., double th_B = 0., double th_C = 0., MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-9);
VEC theta4(const VEC &u, const VEC &ustop, const KSDiscretization& disc, double dt, double th_A = 0., double th_B = 0., double th_C = 0., MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-9);

// braid helpers
int intpow(int base, int exp);
void pack_array(std::ofstream &f, const VEC &u);
void pack_darray(std::ofstream &f, const MAT &u);
void bf_pack_help(double *buf, const MAT& u, const size_t obj_size, size_t &bf_size);
void bf_pack_help(double *buf, const VEC& u, const size_t obj_size, size_t &bf_size);
void bf_unpack_help(double *buf, MAT &u, const size_t obj_size, size_t &bf_size);
void bf_unpack_help(double *buf, VEC &u, const size_t obj_size, size_t &bf_size);

#endif