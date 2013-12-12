#ifndef __CAR_PLANNER_H__
#define __CAR_PLANNER_H__

#include "Planner.h"

class CarPlanner: public Planner
{
public:
	Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m, double step = DT);
	Matrix<Z_DIM> h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n);

	CarPlanner();

	void initPlanner();
	
private:
	struct RRTNode {
		Matrix<X_DIM> x;
		Matrix<U_DIM> u;
		int bp;
	};

	struct PathNode {
		Matrix<X_DIM> x;
		Matrix<U_DIM> u;
	};

	void RRT();
	void initTrajectory();

private:
	// car-specific
	double car_l;
	double goalr;

	Matrix<DIM> b0, b1;
	int cal_obstacles2;

	Matrix<U_DIM> uMin, uMax;
	Matrix<X_DIM> xMin, xMax;
};

#endif