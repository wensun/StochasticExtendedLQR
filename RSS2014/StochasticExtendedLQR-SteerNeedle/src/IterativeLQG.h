#ifndef _ITERATIVELQG_
#define _ITERATIVELQG_

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
#include <time.h>
#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"
#include "CostFunction.h"

class IterativeLQG{

public:

	int T; //step number.
	Matrix<BDIM, 1> Start;
	Matrix<XDIM, 1> xGoal;
	
	int numTotalIter;
	double totalTime;
	double totalCost;

	double dt;

	Matrix<UDIM> uNorminal;

	std::vector<Matrix<UDIM, BDIM>> L;
	std::vector<Matrix<UDIM,1>> l;
	
	std::vector<std::pair<Matrix<BDIM, 1>, Matrix<UDIM, 1>>> nominalPlan;

	std::vector<Matrix<BDIM,1>> finalstate;
	std::vector<Matrix<UDIM,1>> finalu;
	std::vector<Matrix<UDIM,XDIM>> finalL;
	std::vector<Matrix<UDIM,1>> finall;

	int cal_p;
	int cal_env;
	int cal_obs;

	IterativeLQG(const Matrix<BDIM, 1>& s, const Matrix<XDIM,1>& cg,
				const int& numS,
				 const double& Deltat, const int& cal_point, const int& cal_environment, const int& cal_obstacles)
	{
		dt = Deltat;
		T = numS; //0, 1, ...., T (T+1 steps)
		Start = s;
		xGoal = cg;
		L.clear();
		l.clear();
		nominalPlan.clear();

		numTotalIter = 0;
		totalTime = 0.0;
		totalCost = 0.0;

		finalstate.clear();
		finalu.clear();
		finalL.clear();
		finall.clear();

		CostFunction cf(Start, xGoal, cal_p, cal_env, cal_obs);
		uNorminal = cf.uNominal;
	}

	inline void computeNominalPlanFromPolicy(std::vector<std::pair<Matrix<BDIM, 1>, Matrix<UDIM, 1>>>& plan, 
											const std::vector<Matrix<UDIM, BDIM>>& L,
											const std::vector<Matrix<UDIM,1>>& l)
	{
		Dynamics dyn(dt);
		Matrix<BDIM, 1> bnext;
		Matrix<BDIM, XDIM> M;

		plan.clear();
		plan.resize(T + 1);

		plan[0].first = Start;
		for(int i = 0; i < T; i++){
			//u = L * x + I.
			plan[i].second = L[i] * plan[i].first + l[i];
			dyn.BeliefDynamics(plan[i].first, plan[i].second, bnext, M);
			plan[i+1].first = bnext;
		}
		plan[T].second = zeros<UDIM, 1>();
		return;
	}

	//compute the expected cost by linearizing and quadratizing about the given plan. return the cost with costtogo.
	inline void BvalueIteration(const std::vector<std::pair<Matrix<BDIM>, Matrix<UDIM>>>& plan)
	{
		Dynamics dyn(dt);

		L.resize(T, zeros<UDIM,BDIM>());
		l.resize(T, zeros<UDIM,1>());

		CostFunction cf(Start, xGoal, cal_p, cal_env, cal_obs);
		
		SymmetricMatrix<BDIM> S = zeros<BDIM>();
		Matrix<BDIM, 1> svec = zeros<BDIM, 1>();
		double sscale = 0.0;

		Matrix<BDIM, BDIM> QT; Matrix<BDIM, 1> qvecT; double qT = 0;
		cf.quadratizefinalCost(plan[T].first, QT, qvecT, qT);
		S = SymProd(QT, identity<BDIM>()); svec = qvecT; sscale = qT;
		
		for(int t = T-1; t >= 0; t--){
			Matrix<BDIM, BDIM> At; Matrix<BDIM, UDIM> Bt; Matrix<BDIM, 1> avect;
			std::vector<Matrix<BDIM, BDIM>> Ft;
			std::vector<Matrix<BDIM, UDIM>> Gt;
			std::vector<Matrix<BDIM, 1>> evect;
			dyn.linearize_discrete_belief_dynamics(plan[t].first, plan[t].second, At, Bt, avect, Ft, Gt, evect);
			Matrix<BDIM, BDIM> Q; Matrix<UDIM, BDIM> Pt; Matrix<UDIM, UDIM> R;
			Matrix<BDIM, 1> qvect; Matrix<UDIM, 1> rvect;
			double qscalet = 0;
			cf.quadratizeCost_ov(plan[t].first, plan[t].second, Q, Pt, R, qvect, rvect, qscalet);
			SymmetricMatrix<BDIM> Qt = SymProd(Q, identity<BDIM>());
			SymmetricMatrix<UDIM> Rt = SymProd(R, identity<UDIM>());

			SymmetricMatrix<BDIM> Ct = Qt + SymProd(~At, S*At);
			SymmetricMatrix<UDIM> Dt = Rt + SymProd(~Bt, S*Bt);
			Matrix<UDIM, BDIM> Et = Pt + ~Bt * S * At;
			Matrix<BDIM, 1> cvect = qvect + ~At * S * avect + ~At * svec;
			Matrix<UDIM, 1> dvect = rvect + ~Bt * S * avect + ~Bt * svec;
			double escalet = qscalet + sscale + 0.5 * tr(~avect * S * avect) + tr(~avect * svec);

			for(int i = 0; i < XDIM; i++){
				Ct += SymProd(~Ft[i], S * Ft[i]);
				Dt += SymProd(~Gt[i], S * Gt[i]);
				Et += ~Gt[i] * S * Ft[i];
				cvect += ~Ft[i] * S * evect[i];
				dvect += ~Gt[i] * S * evect[i];
				escalet += 0.5 * tr(~evect[i] * S * evect[i]);
			}		

			regularize(Dt);
			L[t] = -(Dt % Et);
			l[t] = -(Dt % dvect);

			S = Ct + SymProd(~Et, L[t]);
			svec = cvect + ~Et * l[t];
			sscale = escalet + 0.5 * tr(~dvect * (l[t]));
		}
		return;
	}

	inline double ExpectedCost(const std::vector<std::pair<Matrix<BDIM>, Matrix<UDIM>>>& plan, 
							   const std::vector<Matrix<UDIM, BDIM>>& L)
	{
		double expcost = 0.0;

		Dynamics dyn(dt);
		CostFunction cf(Start, xGoal, cal_p, cal_env, cal_obs);

		SymmetricMatrix<BDIM> S = zeros<BDIM>();
		Matrix<BDIM, 1> svec = zeros<BDIM, 1>();
		double sscale = 0.0;

		Matrix<BDIM, BDIM> QT; Matrix<BDIM, 1> qvecT; double qT = 0;
		cf.quadratizefinalCost(plan[T].first, QT, qvecT, qT);
		S = SymProd(QT, identity<BDIM>()); svec = qvecT; sscale = qT;
		
		for(int t = T-1; t >= 0; t--){
			Matrix<BDIM, BDIM> At; Matrix<BDIM, UDIM> Bt; Matrix<BDIM, 1> avect;
			std::vector<Matrix<BDIM, BDIM>> Ft;
			std::vector<Matrix<BDIM, UDIM>> Gt;
			std::vector<Matrix<BDIM, 1>> evect;
			dyn.linearize_discrete_belief_dynamics(plan[t].first, plan[t].second, At, Bt, avect, Ft, Gt, evect);
			SymmetricMatrix<XDIM> Qt; Matrix<UDIM, XDIM> Pt; SymmetricMatrix<UDIM> Rt;
			Matrix<XDIM, 1> qvect; Matrix<UDIM, 1> rvect;
			double qscalet = 0;

			cf.quadratizeCost(plan[t].first, plan[t].second, t, 2, Pt, Qt, Rt, qvect, rvect, qscalet); 

			SymmetricMatrix<XDIM> Ct = Qt + SymProd(~At, S*At);
			SymmetricMatrix<UDIM> Dt = Rt + SymProd(~Bt, S*Bt);
			Matrix<UDIM, XDIM> Et = Pt + ~Bt * S * At;
			Matrix<XDIM, 1> cvect = qvect + ~At * S * avect + ~At * svec;
			Matrix<UDIM, 1> dvect = rvect + ~Bt * S * avect + ~Bt * svec;
			double escalet = qscalet + sscale + 0.5 * tr(~avect * S * avect) + tr(~avect * svec);

			for(int i = 0; i < XDIM; i++){
				Ct += SymProd(~Ft[i], S * Ft[i]);
				Dt += SymProd(~Gt[i], S * Gt[i]);
				Et += ~Gt[i] * S * Ft[i];
				cvect += ~Ft[i] * S * evect[i];
				dvect += ~Gt[i] * S * evect[i];
				escalet += 0.5 * tr(~evect[i] * S * evect[i]);
			}		

			Matrix<UDIM> kvect = L[t] * (-plan[t].first) + plan[t].second;
			S = Ct + SymProd(~L[t], Dt * L[t]) + SymSum(~L[t]*Et);
			svec = ~L[t] * Dt * kvect + ~Et * kvect + cvect + ~L[t] * dvect; 
			sscale = 0.5 * tr(~kvect * Dt * kvect) + tr(~kvect * dvect) + escalet;
		}
		
		return 0.5 * tr(~plan[0].first * S * plan[0].first) + tr(~svec * plan[0].first) + sscale;
	}
};


#endif