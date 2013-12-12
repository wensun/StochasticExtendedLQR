#ifndef _COSTFUNCTIONn_
#define _COSTFUNCTIONn_

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <list>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <ppl.h>
#include <omp.h>
#include <time.h>
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Obstacles.h"
#include "Dynamics.h" 

class Costfunction{

public:

	Matrix<XDIM> start;
	Matrix<XDIM> cgoal;
	std::vector<Obstacles::obs> obstacles;
	Matrix<3> bottomLeft;
	Matrix<3> topRight;
	Matrix<U_DIM> uNominal;

	double obstacleFactor; 
	double robotRadius;
	double scaleFactor;
	double rotCost;
	double dt;

	SymmetricMatrix<X_DIM> Q;
	SymmetricMatrix<U_DIM> R;

	Costfunction(const Matrix<XDIM>& _start, const Matrix<XDIM>& _cgoal, const double& _dt)
	{
		dt = _dt;
		start = _start;
		cgoal = _cgoal;
		Obstacles obs;
		obstacles = obs.obstacles;
		topRight = obs.Top;
		bottomLeft = obs.Bottom;
		Dynamics dyn(dt);
		uNominal = dyn.uNominal;

		Q = 500*identity<X_DIM>();
		R = 20*identity<U_DIM>();

		rotCost = 1.5; 
		obstacleFactor = 1.0;
		scaleFactor = 10;
		robotRadius = length + 0.1;  //length defined in the dynamics.h file, which contains all the robot parameters.
	}

	// Obstacle-cost term in local cost function
	inline double obstacleCost(const Matrix<X_DIM>& x) {
	
		double cost = 0;
		Matrix<DIM> pos = x.subMatrix<DIM>(0,0);

		for (size_t i = 0; i < obstacles.size(); ++i) {
			Matrix<DIM> d = pos - obstacles[i].pos;
			d[obstacles[i].dim] = 0;
			double dist = sqrt(scalar(~d*d)) - robotRadius - obstacles[i].radius;
			cost += obstacleFactor * exp(-scaleFactor*dist);
		}
		for (size_t i = 0; i < DIM; ++i) {
			double dist = (pos[i] - bottomLeft[i]) - robotRadius;
			cost += obstacleFactor * exp(-scaleFactor*dist);
		}
		for (size_t i = 0; i < DIM; ++i) {
			double dist = (topRight[i] - pos[i]) - robotRadius;
			cost += obstacleFactor * exp(-scaleFactor*dist);
		}
		return cost;	
	}


	inline void quadratizeObstacleCost(const Matrix<X_DIM>& x, SymmetricMatrix<X_DIM>& Q, Matrix<X_DIM>& q) {
		SymmetricMatrix<DIM> QObs = zeros<DIM>();
		Matrix<DIM> qObs = zero<DIM>();

		Matrix<DIM> pos = x.subMatrix<DIM>(0,0);

		for (size_t i = 0; i < obstacles.size(); ++i) {
			Matrix<DIM> d = pos - obstacles[i].pos;
			d[obstacles[i].dim] = 0;
			double distr = sqrt(tr(~d*d));
			d /= distr;
			double dist = distr - robotRadius - obstacles[i].radius;

			Matrix<DIM> n = zero<DIM>();
			n[obstacles[i].dim] = 1.0;
			Matrix<DIM> d_ortho = skewSymmetric(n)*d;
				
			double a0 = obstacleFactor * exp(-scaleFactor*dist);
			double a1 = -scaleFactor*a0;
			double a2 = -scaleFactor*a1;

			double b2 = a1 / distr;
				
			QObs += a2*SymProd(d,~d) + b2*SymProd(d_ortho,~d_ortho);
			qObs += a1*d;
		}
		for (size_t i = 0; i < DIM; ++i) {
			double dist = (pos[i] - bottomLeft[i]) - robotRadius;

			Matrix<DIM> d = zero<DIM>();
			d[i] = 1.0;

			double a0 = obstacleFactor * exp(-scaleFactor*dist);
			double a1 = -scaleFactor*a0;
			double a2 = -scaleFactor*a1;

			QObs += a2*SymProd(d,~d);
			qObs += a1*d;
		}
		for (size_t i = 0; i < DIM; ++i) {
			double dist = (topRight[i] - pos[i]) - robotRadius;

			Matrix<DIM> d = zero<DIM>();
			d[i] = -1.0;

			double a0 = obstacleFactor * exp(-scaleFactor*dist);
			double a1 = -scaleFactor*a0;
			double a2 = -scaleFactor*a1;

			QObs += a2*SymProd(d,~d);
			qObs += a1*d;
		}
		regularize(QObs);
		Q.insert(0, QObs + Q.subSymmetricMatrix<3>(0));
		q.insert(0,0, qObs - QObs*x.subMatrix<3>(0,0) + q.subMatrix<3>(0,0));
	}


	inline double ct(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const size_t& t) {
		double cost = 0;
		if (t == 0) {
			cost += 0.5 * tr(~(x - start)*Q*(x - start));
		} else {
			cost += obstacleCost(x);
		}
		cost += 0.5 * tr(~(u - uNominal)*R*(u - uNominal));
		return cost;
	}

	inline void quadratizeCost(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const int& t, const int& iter,
								Matrix<U_DIM,X_DIM>& Pt, SymmetricMatrix<X_DIM>& Qt, SymmetricMatrix<U_DIM>& Rt, 
								Matrix<X_DIM>& qt, Matrix<U_DIM>& rt, 
								double& qscalart) {
		if (t == 0) {
			Qt = Q;
			qt = -(Q*start);

		} else {

			Qt = zeros<X_DIM>(); 
			qt = zero<X_DIM>();

			if (iter < 1) {
				Qt.insert(6,rotCost*identity<3>());
			}

			quadratizeObstacleCost(x, Qt, qt);
		}

		Rt = R;
		rt = -(R*uNominal);
		Pt = zeros<U_DIM,X_DIM>();

		qscalart = ct(x, u, t) + 0.5 * tr(~x * Qt * x) + 0.5 * tr(~u * Rt * u) - tr(~x * (qt + Qt * x)) - tr(~u * (rt + Rt * u));

		//if(t == 0 && iter == 0)
		//{
		//	Qt = 100 * identity<XDIM>();
		//	qt = -(Qt * start);
		//	qscalart = 0.5 * tr(~(x - start)*Qt*(x - start)) +  0.5 * tr(~(u - uNominal)*R*(u - uNominal)) + 0.5 * tr(~x * Qt * x) + 0.5 * tr(~u * Rt * u) - tr(~x * (qt + Qt * x)) - tr(~u * (rt + Rt * u));
		//}
	}


	inline double cell(const Matrix<X_DIM>& x) {
		double cost = 0;
		cost += 0.5 * tr(~(x - cgoal)*Q*(x - cgoal));
		return cost;
	}	
	
	inline void quadratizeFinalCost(const Matrix<X_DIM>& x, SymmetricMatrix<X_DIM>& Qell, Matrix<X_DIM>& qell, double& qscalarell,
									const int& iter) 
	{
		/*Qell = hessian(x, cell); 
		qell = jacobian(x, cell) - Qell*x;*/

		Qell = Q;
		qell = -(Q*cgoal);

		qscalarell = cell(x) + 0.5 * tr(~x * Qell * x) - tr(~x * (qell + Qell * x));

		//if(iter == 0){
		//	Qell = 100 * identity<XDIM>();
		//	qell = -(Qell * cgoal);
		//	qscalarell = 0.5 * tr(~(x - cgoal)*Qell*(x - cgoal)) + 0.5 * tr(~x * Qell * x) - tr(~x * (qell + Qell * x));
		//}

	}
	
};


#endif