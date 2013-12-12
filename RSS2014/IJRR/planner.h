#ifndef __PLANNER_H__
#define __PLANNER_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <queue>
#include <stack>
#include <algorithm>

#include "matrix.h"

#include <float.h>

#include "callisto.h"
#include "gammafunc.h"

#include "utils.h"

#define POINT2D
//#define CAR2D
//#define AIRCRAFT3D

//#define PROBLEM1
#define PROBLEM2

#define OBSTACLES

//#define MAXLIKELIHOOD

#if defined(CAR2D)

#define DIM 2 // workspace dimension

#if defined(PROBLEM1)
#define X_DIM 4 // state dimension
#define U_DIM 2 // control input dimension
#define Z_DIM 2 // measurement dimension
#define M_DIM 4 // motion noise dimension
#define N_DIM 2 // measurement noise dimension
#elif defined(PROBLEM2)
#define X_DIM 4 // state dimension
#define U_DIM 2 // control input dimension
#define Z_DIM 4 // measurement dimension
#define M_DIM 4 // motion noise dimension
#define N_DIM 4 // measurement noise dimension
#endif

const double DT = 0.25;  // time step
const double H = 0.0078125; // discretization for finite differencing
const size_t NUM_ITER = 100; // number of iterations
const double COST_INFTY = 1e10;
#endif

#if defined(POINT2D)
#define DIM 2 // workspace dimension

#if defined(PROBLEM1)
#define X_DIM 2 // state dimension
#define U_DIM 2 // control input dimension
#define Z_DIM 2 // measurement dimension
#define M_DIM 2 // motion noise dimension
#define N_DIM 2 // measurement noise dimension
#elif defined(PROBLEM2)
#define X_DIM 2 // state dimension
#define U_DIM 2 // control input dimension
#define M_DIM 2 // motion noise dimension

#define Z_DIM 2 // measurement dimension
#define N_DIM 2 // measurement noise dimension

#endif

const double DT = 1;  // time step
const double H = 0.0078125; // discretization for finite differencing
const size_t NUM_ITER = 50; // number of iterations
const double COST_INFTY = 1e10;
#endif


#if defined(AIRCRAFT3D)
#define DIM 3 // workspace dimension

#if defined(PROBLEM1)
#define X_DIM 6 // state dimension
#define U_DIM 3 // control input dimension
#define Z_DIM 3 // measurement dimension
#define M_DIM 3 // motion noise dimension
#define N_DIM 3 // measurement noise dimension
#elif defined(PROBLEM2)
#define X_DIM 6 // state dimension
#define U_DIM 3 // control input dimension
#define Z_DIM 3 // measurement dimension
#define M_DIM 3 // motion noise dimension
#define N_DIM 3 // measurement noise dimension
#endif

const double DT = 0.25;  // time step
const double H = 0.0078125; // discretization for finite differencing
const size_t NUM_ITER = 50; // number of iterations
const double COST_INFTY = 1e10;
#endif

#define S_DIM (((X_DIM+1)*X_DIM)/2)
#define B_DIM (X_DIM + S_DIM)

class Planner
{
public:
	virtual Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m, double step = DT) = 0;
	virtual Matrix<Z_DIM> h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) = 0;

	double cost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u);
	double finalCost(const Matrix<B_DIM>& b);

	virtual void initPlanner() = 0;
	
	void solveILQG();
	void solveDDP();

	void simulate(const int numParticles, bool display);
	
protected:
	void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma);
	void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b);
	Matrix<S_DIM,S_DIM> QvecS(const Matrix<X_DIM,X_DIM>& Q);

	void beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM>& g, Matrix<X_DIM,X_DIM>& SqrtW);
	
	// Obstacle cost
	double confLine(const Matrix<X_DIM>& p1, const Matrix<X_DIM>& p2, const Matrix<X_DIM, X_DIM>& Sigma);
	double confFunc(const double c);
	double confObs(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u);

	double value(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss);
	void initValue(const Matrix<B_DIM>& b, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss);

	// Jacobians
	Matrix<X_DIM,X_DIM> dfdx(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m);
	Matrix<X_DIM,U_DIM> dfdu(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m);
	Matrix<X_DIM,M_DIM> dfdm(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m);
	
	Matrix<Z_DIM,X_DIM> dhdx(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n);
	Matrix<Z_DIM,N_DIM> dhdn(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n);

	// Jacobians
	Matrix<1,U_DIM> dcdu(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss);
	Matrix<1,B_DIM> dcdb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss);
	Matrix<1,B_DIM> dvdb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<U_DIM,B_DIM>& L, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss);

	// Hessians
	Matrix<U_DIM,U_DIM> ddcdudu(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss);
	Matrix<B_DIM,B_DIM> ddcdbdb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss);
	Matrix<U_DIM,B_DIM> ddcdudb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss);
	Matrix<B_DIM,B_DIM> ddvdbdb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<U_DIM,B_DIM>& L, 
								const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss);

	// Obstacle cost Jacobians
	void obsApprox(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<U_DIM, U_DIM>& D, Matrix<B_DIM, B_DIM>& C, 
				   Matrix<U_DIM, B_DIM>& E, Matrix<1,B_DIM>& c, Matrix<1, U_DIM>& d, double& e);
	void obsApprox2(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<U_DIM, B_DIM>& L, Matrix<B_DIM, B_DIM>& C, Matrix<1,B_DIM>& c, double& e);

	virtual void RRT() = 0;
	virtual void initTrajectory() = 0;
	
	void quadratizeCost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& Q, Matrix<B_DIM>& q, 
						Matrix<U_DIM,B_DIM>& P, Matrix<U_DIM,U_DIM>& R, Matrix<U_DIM>& r, double& p);

	void linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM>& c, Matrix<B_DIM,B_DIM>& A, Matrix<B_DIM, U_DIM>& B, 
			 					 std::vector<Matrix<X_DIM> >& cn, std::vector<Matrix<X_DIM, B_DIM> >& An, std::vector<Matrix<X_DIM, U_DIM> >& Bn);

	void controlPolicyILQG(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, Matrix<U_DIM,B_DIM>& L, Matrix<U_DIM>& l);
	void controlPolicyDDP(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, Matrix<U_DIM,B_DIM>& L, Matrix<U_DIM>& l);

	void valueIteration(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, const Matrix<U_DIM,B_DIM>& L);

	void expectedCost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, const Matrix<U_DIM,B_DIM>& L);

	void backwardIterationILQG();
	void backwardIterationDDP();
	
	void forwardIterationILQG();
	void forwardIterationDDP();
	
	void displayMeanPath(const std::vector<Matrix<B_DIM> > & B);
	void displayTraces(std::vector<Matrix<X_DIM> > & xTrace);

	void saveControlPolicy(const std::vector<Matrix<U_DIM, B_DIM> >& L, const std::vector<Matrix<U_DIM> >& l);
	void loadControlPolicy(std::vector<Matrix<B_DIM> > & b, std::vector<Matrix<U_DIM> >& u, std::vector<Matrix<U_DIM, B_DIM> >& L, std::vector<Matrix<U_DIM> >& l);

protected:
	int pathlen;  // number of time steps

	// Global variables - Callisto (clearance queries and visualization)
	int cal_environment;
	int cal_obstacles;
	int cal_box;
	int cal_goal;
	int cal_point, cal_cylinder;
	int cal_path;
	int cal_rrt;

	// Parameters
	Matrix<U_DIM,U_DIM> Rint;
	Matrix<X_DIM,X_DIM> Qint;
	Matrix<X_DIM> x0, xGoal;
	Matrix<X_DIM,X_DIM> QGoal;
	Matrix<X_DIM,X_DIM> SqrtSigma0;

	// Variables
	std::vector<Matrix<B_DIM> > B;
	std::vector<Matrix<U_DIM> > U;
	
	std::vector<Matrix<U_DIM, B_DIM> > L;
	std::vector<Matrix<U_DIM> > l;

	Matrix<B_DIM, B_DIM> S;
	Matrix<B_DIM> s;
	double ss;
	
	double bestCost;
	double gSqr, Jp;
	double eps, lambda;
	
	int nsteps;
	double invnsteps;

	bool switchoffML;

	//std::ofstream resultsfptr;
};

#endif