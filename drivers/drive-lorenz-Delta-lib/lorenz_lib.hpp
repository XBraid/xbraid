#ifndef LORENZ_LIB
#define LORENZ_LIB
#include <iostream>
#include <fstream>
#include <cmath>
#include "Eigen/Dense" // expects Eigen library in this directory

#define VECSIZE 3
typedef Eigen::Vector3d VEC;
typedef Eigen::Matrix3d MAT;
typedef Eigen::Matrix3Xd COL_MAT; // this is a 3xN matrix with 3d columns
using Eigen::Matrix3d;

// standard chaotic choices for Lorenz parameters
const double rho { 28. };
const double sigma { 10. };
const double beta { 8. / 3. };

// Lyapunov time
const double T_lyap { log(10) / 0.9 };

// functions
VEC f_lorenz(VEC u);
MAT f_lorenz_du(VEC u);

// first order
VEC euler(VEC u, double dt);
MAT euler_du(VEC u, double dt);
VEC theta1(const VEC u, const VEC guess, double dt, double theta=0.5, MAT *P_tan=nullptr, int newton_iters=10, double tol=1e-10);
MAT theta1_du(VEC u, VEC ustop, double dt, double theta=0.5);

// fourth order
VEC rk4(VEC u, double dt, MAT *P_tan = nullptr);
VEC theta4(VEC u, VEC ustop, double dt, double theta = .2, MAT *P_tan = nullptr, int newton_iters = 10, double tol = 1e-10);

void GramSchmidt(COL_MAT &A);
void pack_array(std::ofstream &f, VEC u);
void pack_darray(std::ofstream &f, COL_MAT u);
int intpow(int base, int exp);

#endif