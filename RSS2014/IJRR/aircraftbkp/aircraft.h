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
		Matrix<4,4> T;
		Matrix<X_DIM> x;
		Matrix<U_DIM> u;
		
		int bp;
		double depth;
	
		std::vector<int> children;
	};

	struct PathNode {
		Matrix<4,4> T;
		Matrix<X_DIM> x;
		Matrix<U_DIM> u;

		void serialize(std::ofstream & fptr) {
			serializeMatrix(fptr, T);
			serializeMatrix(fptr, x);
			serializeMatrix(fptr, u);
		}

		void deserialize(std::ifstream & fptr) {
			deserializeMatrix(fptr, T);
			deserializeMatrix(fptr, x);
			deserializeMatrix(fptr, u);
		}
	};

	double dist(const Matrix<4,4>& T, const Matrix<3>& point);
	int nearestNeighbor(const Matrix<3>& point, const std::vector<RRTNode>& tree);

	Matrix<4,4> stepT(const Matrix<4,4>& T, const Matrix<3>& v, const Matrix<3>& w, double dt);
	void createInputs(Matrix<3>& v, Matrix<3>& w, double speed, double twist, double curvature);

	void RRT();

	void initTrajectory();

private:
	// aircraft specific
	double plane_l;
	double goalr;
	double factor;
	
	Matrix<U_DIM> uMin, uMax;
	Matrix<X_DIM> xGoal;
	double Xrange, Yrange, Zrange;
};

#endif