#ifndef _ITERATIVELQG_
#define _ITERATIVELQG_


#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"
#include "Costfunction.h"

class IterativeLQG{

public:
	int T; //step number.
	Matrix<XDIM, 1> Start;
	Matrix<2, 1> goal;
	Matrix<XDIM,1> cgoal;

	int numTotalIter;
	double totalTime;
	double totalCost;

	double dt;

	Matrix<UDIM> uNorminal;

	//std::vector<std::pair<Matrix<UDIM, XDIM>, Matrix<UDIM, 1>>> policy;  //policy: L and I.
	std::vector<Matrix<UDIM, XDIM>> L;
	std::vector<Matrix<UDIM,1>> l;
	
	std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>> nominalPlan;

	std::vector<Matrix<XDIM,1>> finalstate;
	std::vector<Matrix<UDIM,1>> finalu;
	std::vector<Matrix<UDIM,XDIM>> finalL;
	std::vector<Matrix<UDIM,1>> finall;

	IterativeLQG(const Matrix<XDIM, 1>& s, const Matrix<XDIM,1>& cg,
				const int& numS,
				 const double& Deltat)
	{
		dt = Deltat;
		T = numS; //0, 1, ...., T (T+1 steps)
		Start = s;
		cgoal = cg;
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

		Costfunction cf(Start, cgoal, dt);
		uNorminal = cf.uNominal;
	}

	
	inline void computeNominalPlanFromPolicy(std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>>& plan, 
											const std::vector<Matrix<UDIM, XDIM>>& L,
											const std::vector<Matrix<UDIM,1>>& l)
	{
		Dynamics dyn(dt);
		Matrix<XDIM, 1> g;
		Matrix<XDIM, XDIM> M;

		plan.clear();
		plan.resize(T + 1);

		plan[0].first = Start;
		for(int i = 0; i < T; i++){
			//u = L * x + I.
			plan[i].second = L[i] * plan[i].first + l[i];
			g.reset(); M.reset();
			dyn.discrete_dynamics(plan[i].first, plan[i].second, g, M);	
			plan[i+1].first = g;
		}
		plan[T].second = zeros<UDIM, 1>();

		return;
	}

	//compute cost given a plan.
	inline double computeCost(const std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>>& plan)
	{
		double cost = 0;
		Costfunction cf(Start, cgoal, dt);
		cost += cf.ct(plan[0].first, plan[0].second, 0);
		
		for(int i = 1; i < (int)plan.size() - 1; i++){
			cost += cf.ct(plan[i].first, plan[i].second, i);
			//std::cout<<cost<<std::endl;
		}
		
		cost += cf.cell(plan[(int)plan.size() - 1].first);
		return cost;
	}
	
	//compute the expected cost by linearizing and quadratizing about the given plan. return the cost with costtogo.
	inline void BvalueIteration(const std::vector<std::pair<Matrix<XDIM>, Matrix<UDIM>>>& plan)
	{

		Dynamics dyn(dt);

		L.resize(T, zeros<UDIM,XDIM>());
		l.resize(T, zeros<UDIM,1>());

		Costfunction cf(Start, cgoal, dt);

		SymmetricMatrix<XDIM> S = zeros<XDIM>();
		Matrix<XDIM, 1> svec = zeros<XDIM, 1>();
		double sscale = 0.0;

		SymmetricMatrix<XDIM> QT; Matrix<XDIM, 1> qvecT; double qT = 0;
		cf.quadratizeFinalCost(plan[T].first, QT, qvecT, qT, 1);
		S = QT; svec = qvecT; sscale = qT;
		
		for(int t = T-1; t >= 0; t--){
			Matrix<XDIM, XDIM> At; Matrix<XDIM, UDIM> Bt; Matrix<XDIM, 1> avect;
			std::vector<Matrix<XDIM, XDIM>> Ft;
			std::vector<Matrix<XDIM, UDIM>> Gt;
			std::vector<Matrix<XDIM, 1>> evect;
			dyn.linearize_discrete_dynamics(plan[t].first, plan[t].second, At, Bt, avect, Ft, Gt, evect);
			SymmetricMatrix<XDIM> Qt; Matrix<UDIM, XDIM> Pt; SymmetricMatrix<UDIM> Rt;
			Matrix<XDIM, 1> qvect; Matrix<UDIM, 1> rvect;
			double qscalet = 0;
			
			cf.quadratizeCost(plan[t].first, plan[t].second, t, 1, Pt, Qt, Rt, qvect, rvect, qscalet);

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

			//regularize(Dt);
			L[t] = -(Dt % Et);
			l[t] = -(Dt % dvect);

			S = Ct + SymProd(~Et, L[t]);
			svec = cvect + ~Et * l[t];
			sscale = escalet + 0.5 * tr(~dvect * (l[t]));
		}


		return;
	}

	inline double ExpectedCost(const std::vector<std::pair<Matrix<XDIM>, Matrix<UDIM>>>& plan, 
							   const std::vector<Matrix<UDIM, XDIM>>& L)
	{
		double expcost = 0.0;

		Dynamics dyn(dt);
		Costfunction cf(Start, cgoal, dt);

		SymmetricMatrix<XDIM> S = zeros<XDIM>();
		Matrix<XDIM, 1> svec = zeros<XDIM, 1>();
		double sscale = 0.0;

		SymmetricMatrix<XDIM> QT; Matrix<XDIM, 1> qvecT; double qT = 0;
		cf.quadratizeFinalCost(plan[T].first, QT, qvecT, qT, 2);
		S = QT; svec = qvecT; sscale = qT;
		
		for(int t = T-1; t >= 0; t--){
			Matrix<XDIM, XDIM> At; Matrix<XDIM, UDIM> Bt; Matrix<XDIM, 1> avect;
			std::vector<Matrix<XDIM, XDIM>> Ft;
			std::vector<Matrix<XDIM, UDIM>> Gt;
			std::vector<Matrix<XDIM, 1>> evect;
			dyn.linearize_discrete_dynamics(plan[t].first, plan[t].second, At, Bt, avect, Ft, Gt, evect);
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


	//line search
	inline void adjustPath(std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>>& plan, 
		const std::vector<std::pair<Matrix<XDIM,1>, Matrix<UDIM,1>>>& oldplan,
		const std::vector<Matrix<UDIM,XDIM>>& L,
		const std::vector<Matrix<UDIM, 1>>& l,
		const double& epsilon)
	{
		Dynamics dyn(dt);
		Matrix<XDIM, 1> g; Matrix<XDIM, XDIM> M;
		plan.resize(T+1);
		plan[0].first = Start;

		for(int i = 0; i < (int)plan.size() - 1; i++){
			Matrix<UDIM> u = (1 - epsilon) * oldplan[i].second + L[i]*(plan[i].first - (1 - epsilon) * oldplan[i].first) + epsilon * l[i];
			plan[i].second = u;
			g.reset(); M.reset();
			dyn.discrete_dynamics(plan[i].first, plan[i].second, g, M);
			plan[i+1].first = g;
			//std::cout<<g<<std::endl;
		}
		return;
	}

	inline double diffTwoPlans(const std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM,1>>>& p1, const std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM,1>>>& p2)
	{
		double diff = 0.0;
		double disp2 = 0.0;
		for(int i = 0; i <= T; i++){
			diff += sqrt(tr(~(p1[i].first - p2[i].first)*(p1[i].first - p2[i].first)));
			disp2 += sqrt(tr(~(p2[i].first)*(p2[i].first)));
		}

		return diff/disp2;
	}


	inline void iterate()
	{

		clock_t Tstart = clock();

		//CAL_SuspendVisualisation();

		std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>> plan;
		std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>> tmpplan;

		L.resize(T, zeros<UDIM,XDIM>());
		l.resize(T, uNorminal);
		
		computeNominalPlanFromPolicy(plan, L, l);
		
		double costnew = 0.0;
		double costold = ExpectedCost(plan, L);
		
		while(true){

			numTotalIter ++;
			BvalueIteration(plan); //update a new policy based on the nominalPlan and return the exepcted cost of this policy.
			computeNominalPlanFromPolicy(tmpplan, L, l);
			costnew = ExpectedCost(tmpplan,L);

			double epsilon = 1.0;
			while(costnew > costold)
			{
				epsilon *= 0.5;
				adjustPath(tmpplan, plan, L, l, epsilon);
				costnew = ExpectedCost(tmpplan, L);
			}
			std::cout<<costnew<<std::endl;

			if(abs(costnew - costold) / costold < 1e-4){
				totalCost = costnew;
				nominalPlan = plan;
				break;
			}
			else{
				plan = tmpplan;
				costold = costnew;
			}

		}


		totalTime = (double)(clock() - Tstart) / CLOCKS_PER_SEC;

		//CAL_ResumeVisualisation();

		return;  //the plan containts the final local optimal plan.
	}

	inline void nominalplan()
	{
		nominalPlan.resize(T+1);
		for(int i = 0; i <= T; i++){
			nominalPlan[i].first = finalstate[i];
			nominalPlan[i].second = finalu[i];
		}
		return;
	}
};

#endif