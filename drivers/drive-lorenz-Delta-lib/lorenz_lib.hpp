#ifndef LORENZ_LIB
#define LORENZ_LIB
#include <iostream>
#include <fstream>
#include <cmath>
#include "../Eigen/Dense" // expects Eigen library in xbraid/drivers/

#define DIM 3
typedef Eigen::Vector3d VEC;
typedef Eigen::Matrix3d MAT;
typedef Eigen::Matrix3Xd COL_MAT; // this is a 3xN matrix with 3d columns
using Eigen::Matrix3d;

// standard chaotic choices for Lorenz parameters
const double rho{28.};
const double sigma{10.};
const double beta{8. / 3.};

// Lyapunov time
const double T_lyap{log(10) / 0.9};

// functions
VEC f_lorenz(const VEC &u);
MAT f_lorenz_du(const VEC &u);

// first order
VEC euler(VEC u, double dt);
MAT euler_du(VEC u, double dt);
VEC theta1(const VEC u, const VEC guess, double dt, double theta = 0.5, MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-10);

// second order
VEC crank_nicolson(VEC u, VEC ustop, double dt, MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-10);
VEC theta2(VEC u, VEC ustop, double dt, double th_A = 0., double th_B = 0., double th_C = 0., MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-10);

// fourth order
VEC rk4(VEC u, double dt, MAT *P_tan = nullptr);
VEC theta4(VEC u, VEC &uguess, double dt, double theta = .2, MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-10);
VEC thetaLobatto3(VEC u, VEC ustop, double dt, double th_A, double th_B, double th_C, MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-10, bool linesearch=true);

void GramSchmidt(COL_MAT &A);
void pack_array(std::ofstream &f, VEC u);
void pack_darray(std::ofstream &f, COL_MAT u);
int intpow(int base, int exp);

// braid helpers
void bf_pack_help(double *buf, const COL_MAT u, const size_t obj_size, size_t &bf_size);
void bf_pack_help(double *buf, const VEC u, const size_t obj_size, size_t &bf_size);
void bf_pack_help(double *buf, const MAT u, const size_t obj_size, size_t &bf_size);
void bf_unpack_help(double *buf, COL_MAT &u, const size_t obj_size, size_t &bf_size);
void bf_unpack_help(double *buf, VEC &u, const size_t obj_size, size_t &bf_size);
void bf_unpack_help(double *buf, MAT &u, const size_t obj_size, size_t &bf_size);

#endif