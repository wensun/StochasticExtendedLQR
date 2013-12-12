#ifndef __POINT_PLANNER_H__
#define __POINT_PLANNER_H__

#include "planner.h"

class PointPlanner: public Planner
{
public:
	Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m, double step = DT);
	Matrix<Z_DIM> h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n);

	void initPlanner();
	void initProblem();

	PointPlanner();
	
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
	double goalr;
	std::ifstream rrtfptr;
};

#endif