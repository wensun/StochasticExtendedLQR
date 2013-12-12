#ifndef __AIRCRAFT_PLANNER_H__
#define __AIRCRAFT_PLANNER_H__

#include "Planner.h"

class AircraftPlanner: public Planner
{
public:
	Matrix<X_DIM> f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m, double step = DT);
	Matrix<Z_DIM> h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n);

	AircraftPlanner();

	void initPlanner();
	
private:
	struct RRTNode {
		Matrix<X_DIM> x;
		Matrix<U_DIM> u;
		int bp;
		std::vector<int> children;
	};

	struct PathNode {
		Matrix<X_DIM> x;
		Matrix<U_DIM> u;

		void serialize(std::ofstream & fptr) {
			serializeMatrix(fptr, x);
			serializeMatrix(fptr, u);
		}

		void deserialize(std::ifstream & fptr) {
			deserializeMatrix(fptr, x);
			deserializeMatrix(fptr, u);
		}
	};

	void RRT();

	void initTrajectory();

private:
	// aircraft specific
	double goalr;
	
	Matrix<U_DIM> uMin, uMax;
	Matrix<X_DIM> xMin, xMax;
};

#endif