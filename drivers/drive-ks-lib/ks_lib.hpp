#ifndef KS_LIB
#define KS_LIB

// std lib
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>

// linear algebra solver includes
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/UmfPackSupport" // need to install SuiteSparse

// want the flexibility of dynamic arrays to dynamically change problem size
using Eigen::Index;
typedef Eigen::VectorXd VEC;
typedef Eigen::MatrixXd MAT;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SPMAT;
typedef Eigen::Triplet<double> T;
typedef std::vector<double> Stencil;


// exported functions
class KSDiscretization
{
public:
    SPMAT L;
    SPMAT Dx;
    int nx;
    double dx;

    // constructor
    KSDiscretization() : L(SPMAT()), Dx(SPMAT()), nx(0), dx(0.) {} // default. Never use this.
    KSDiscretization(int nx_, double length, Stencil d1, Stencil d2, Stencil d4);

    // implements discretization:
    VEC f_ks(const VEC& u) const;
    SPMAT f_ks_du(const VEC& u) const;
};

SPMAT circulant_from_stencil(std::vector<double> stencil, int n);
VEC theta2(const VEC &u, const KSDiscretization& disc, double dt, double th_A = 0., double th_B = 0., double th_C = 1., MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-10);
void pack_array(std::ofstream &f, const VEC &u);
int intpow(int base, int exp);

void bf_pack_help(double *buf, const MAT& u, const size_t obj_size, size_t &bf_size);
void bf_pack_help(double *buf, const VEC& u, const size_t obj_size, size_t &bf_size);
void bf_unpack_help(double *buf, MAT &u, const size_t obj_size, size_t &bf_size);
void bf_unpack_help(double *buf, VEC &u, const size_t obj_size, size_t &bf_size);

#endif