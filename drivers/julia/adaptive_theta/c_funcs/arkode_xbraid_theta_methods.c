/* theta_esdirk2 */

#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

void _theta_esdirk2_lhs(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * rhs[0] + th[0];
}

void _theta_esdirk2_lhs_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(1.0);
}

void _theta_esdirk2_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = RCONST(0.0);
  A[1] = RCONST(0.0);
  A[2] = RCONST(1.0) + RCONST(-1.0) * th[0];
  A[3] = th[0];
}

void _theta_esdirk2_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[0];
  b[1] = th[0];
}

void _theta_esdirk2_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = RCONST(0.0);
  c[1] = RCONST(1.0);
}

/* theta_sdirk2 */

void _theta_sdirk2_lhs(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * SUNRpowerI(th[0], 2) + RCONST(-1.0) * rhs[0] + RCONST(2.0) * th[0];
}

void _theta_sdirk2_lhs_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(2.0) + RCONST(-2.0) * th[0];
}

void _theta_sdirk2_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = th[0];
  A[1] = RCONST(0.0);
  A[2] = RCONST(1.0) + RCONST(-1.0) * th[0];
  A[3] = th[0];
}

void _theta_sdirk2_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[0];
  b[1] = th[0];
}

void _theta_sdirk2_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = th[0];
  c[1] = RCONST(1.0);
}

/* theta_esdirk3 */

void _theta_esdirk3_lhs(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * rhs[0] + RCONST(-1.0) * th[0] * th[1] + RCONST(-1.0) * th[1] * th[2] + th[0] + th[1];
  phi[1] = RCONST(-1.0) * rhs[1] + SUNRpowerI(th[1], 2) + RCONST(-1.0) * SUNRpowerI(th[1], 2) * th[0] + RCONST(-1.0) * SUNRpowerI(th[1], 2) * th[2] + th[0];
  phi[2] = RCONST(-1.0) * rhs[2] + SUNRpowerI(th[0], 2) + RCONST(-2.0) * SUNRpowerI(th[0], 2) * th[1] + RCONST(2.0) * th[0] * th[1] + RCONST(-2.0) * th[0] * th[1] * th[2];
}

void _theta_esdirk3_lhs_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(1.0) + RCONST(-1.0) * th[1];
  phi_J[1] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  phi_J[2] = RCONST(-1.0) * th[1];
  phi_J[3] = RCONST(1.0) + RCONST(-1.0) * SUNRpowerI(th[1], 2);
  phi_J[4] = RCONST(2.0) * th[1] + RCONST(-2.0) * th[0] * th[1] + RCONST(-2.0) * th[1] * th[2];
  phi_J[5] = RCONST(-1.0) * SUNRpowerI(th[1], 2);
  phi_J[6] = RCONST(2.0) * th[0] + RCONST(2.0) * th[1] + RCONST(-4.0) * th[0] * th[1] + RCONST(-2.0) * th[1] * th[2];
  phi_J[7] = RCONST(-2.0) * SUNRpowerI(th[0], 2) + RCONST(2.0) * th[0] + RCONST(-2.0) * th[0] * th[2];
  phi_J[8] = RCONST(-2.0) * th[0] * th[1];
}

void _theta_esdirk3_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = RCONST(0.0);
  A[1] = RCONST(0.0);
  A[2] = RCONST(0.0);
  A[3] = RCONST(-1.0) * th[0] + th[1];
  A[4] = th[0];
  A[5] = RCONST(0.0);
  A[6] = th[2];
  A[7] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  A[8] = th[0];
}

void _theta_esdirk3_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = th[2];
  b[1] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  b[2] = th[0];
}

void _theta_esdirk3_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = RCONST(0.0);
  c[1] = th[1];
  c[2] = RCONST(1.0);
}

/* theta_sdirk3 */

void _theta_sdirk3_lhs(sunrealtype* phi, const sunrealtype* th, const sunrealtype* rhs)
{
  phi[0] = RCONST(-1.0) * SUNRpowerI(th[0], 2) + RCONST(-1.0) * rhs[0] + RCONST(2.0) * th[0] + th[1] * th[2] + RCONST(-1.0) * th[0] * th[2];
  phi[1] = RCONST(-1.0) * SUNRpowerI(th[0], 3) + RCONST(-1.0) * rhs[1] + SUNRpowerI(th[1], 2) * th[2] + SUNRpowerI(th[0], 2) + RCONST(-1.0) * SUNRpowerI(th[0], 2) * th[2] + th[0];
  phi[2] = RCONST(-2.0) * SUNRpowerI(th[0], 3) + RCONST(-1.0) * rhs[2] + RCONST(3.0) * SUNRpowerI(th[0], 2) + RCONST(-3.0) * SUNRpowerI(th[0], 2) * th[2] + RCONST(3.0) * th[0] * th[1] * th[2];
}

void _theta_sdirk3_lhs_jac(sunrealtype* phi_J, const sunrealtype* th)
{
  phi_J[0] = RCONST(2.0) + RCONST(-1.0) * th[2] + RCONST(-2.0) * th[0];
  phi_J[1] = th[2];
  phi_J[2] = RCONST(-1.0) * th[0] + th[1];
  phi_J[3] = RCONST(1.0) + RCONST(2.0) * th[0] + RCONST(-3.0) * SUNRpowerI(th[0], 2) + RCONST(-2.0) * th[0] * th[2];
  phi_J[4] = RCONST(2.0) * th[1] * th[2];
  phi_J[5] = RCONST(-1.0) * SUNRpowerI(th[0], 2) + SUNRpowerI(th[1], 2);
  phi_J[6] = RCONST(-6.0) * SUNRpowerI(th[0], 2) + RCONST(6.0) * th[0] + RCONST(-6.0) * th[0] * th[2] + RCONST(3.0) * th[1] * th[2];
  phi_J[7] = RCONST(3.0) * th[0] * th[2];
  phi_J[8] = RCONST(-3.0) * SUNRpowerI(th[0], 2) + RCONST(3.0) * th[0] * th[1];
}

void _theta_sdirk3_btable_A(sunrealtype* A, const sunrealtype* th)
{
  A[0] = th[0];
  A[1] = RCONST(0.0);
  A[2] = RCONST(0.0);
  A[3] = RCONST(-1.0) * th[0] + th[1];
  A[4] = th[0];
  A[5] = RCONST(0.0);
  A[6] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  A[7] = th[2];
  A[8] = th[0];
}

void _theta_sdirk3_btable_b(sunrealtype* b, const sunrealtype* th)
{
  b[0] = RCONST(1.0) + RCONST(-1.0) * th[0] + RCONST(-1.0) * th[2];
  b[1] = th[2];
  b[2] = th[0];
}

void _theta_sdirk3_btable_c(sunrealtype* c, const sunrealtype* th)
{
  c[0] = th[0];
  c[1] = th[1];
  c[2] = RCONST(1.0);
}

