#include <stdio.h>
#include <tchar.h>
#include <vector>

#include "matrix.h"

#include <float.h>
//unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

#define X_DIM 2 // state dimension
#define U_DIM 2 // control input dimension
#define Z_DIM 2 // measurement dimension
#define M_DIM 2 // motion noise dimension
#define N_DIM 2 // measurement noise dimension

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM + S_DIM)

// Global constants
const double DT = 1;  // time step
const size_t LENGTH = 20;  // number of time steps
const double H1 = 0.0078125; // discretization for Jacobian and Hessian
const size_t NUM_ITER = 99; // number of iterations

// Global variables
Matrix<U_DIM,U_DIM> Rint;
Matrix<X_DIM,X_DIM> Qint;
Matrix<X_DIM,X_DIM> Qgoal;
Matrix<X_DIM> x0;
Matrix<X_DIM,X_DIM> SqrtSigma0;
Matrix<X_DIM> xGoal;

// Dynamics model
inline Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m) {  // m ~ N(0,I)
  Matrix<X_DIM> xNew;
  xNew[0] = x[0] + u[0] * DT + sqrt((0.1*0.1)*DT)*m[0];
  xNew[1] = x[1] + u[1] * DT + sqrt((0.1*0.1)*DT)*m[1];
  return xNew;
}

// Observation model
inline Matrix<Z_DIM> h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) {  // n ~ N(0,I)
  Matrix<Z_DIM> z;
  z[0] = x[0] + sqrt(0.5*(5-x[0])*(5-x[0]) + 0.1)*n[0];
  z[1] = x[1] + sqrt(0.5*(5-x[0])*(5-x[0]) + 0.1)*n[1];
  return z;
}

// Parameter initialization
inline void userInit() {
  Rint = identity<U_DIM>();
  Qgoal = 10*identity<X_DIM>();
  Qint = identity<X_DIM>();
  x0[0] = 2; x0[1] = 2;
  SqrtSigma0 = sqrt(5.0)*identity<X_DIM>();
  xGoal[0] = 0; xGoal[1] = 0;
}