#ifndef _EXTENDEDLQR_
#define _EXTENDEDLQR_

#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"
#include "Costfunction.h"
#include "IterativeLQG.h"

class ExtendedLQR{

public:
	Matrix<XDIM, 1> start;
	Matrix<2, 1> goal;
	Matrix<XDIM,1> cgoal;
	int T;

	int numTotalIter;
	double totalTime;
	double totalCost;

	std::vector<std::pair<Matrix<UDIM, XDIM>, Matrix<UDIM,1>>> policy;
	std::vector<Matrix<XDIM,1>> stateplan;
	std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>> nominalPlan;

	std::vector<Matrix<XDIM,1>> minsequenceStates;

	int cal_l1;
	int cal_l2;
	int cal_obs;

	double lambda; //control the parameter of obstacle avoidance.

	ExtendedLQR(const Matrix<XDIM, 1>& s, const Matrix<2, 1>& g, const Matrix<XDIM,1>& cg, 
				const int& num,
				const int& cal_link1, const int& cal_link2, const int& cal_obstacles,
				const double& Lambda)
	{
		totalTime = 0.0;
		numTotalIter = 0;
		start = s;
		goal = g;
		cgoal = cg;
		T = num;
		totalCost = 0.0;

		cal_l1 = cal_link1;
		cal_l2 = cal_link2;
		cal_obs = cal_obstacles;

		lambda = Lambda;
		minsequenceStates.resize(T+1);
	}



	inline void computeNominalPlan()
	{
		nominalPlan.clear(); nominalPlan.resize(T+1);
		for(int i = 0; i <= T-1; i++){
			nominalPlan[i].first = stateplan[i];
			nominalPlan[i].second = policy[i].first * stateplan[i] + policy[i].second;
		}
		nominalPlan[T].first = stateplan[T];
		nominalPlan[T].second = zeros<UDIM, 1>();

		return;
	}

	//check if the total cost along the final trajecotry is consistent with each other
	inline bool checkTotalCost(const std::vector<Matrix<UDIM,XDIM> >& L, const std::vector<Matrix<UDIM> >& l,
							   const std::vector<SymmetricMatrix<XDIM>>& S, const std::vector<SymmetricMatrix<XDIM>>& SBar,
							   const Matrix<XDIM>& xHat,
							   const std::vector<Matrix<XDIM>>& s, const std::vector<Matrix<XDIM>>& sBar,
							   const std::vector<double>& sscale, const std::vector<double>& sbarscale)
	{
		Dynamics dyn;
		Matrix<XDIM> x(xHat);
		Matrix<UDIM> u;
		Matrix<XDIM,XDIM> M;
		double currentTotalCost = 0.5 * tr(~x * S[0] * x) + tr(~x * s[0]) + sscale[0];
		std::cout<<currentTotalCost<<std::endl;
		currentTotalCost += 0.5 * tr(~x * SBar[0] *x) + tr(~x * sBar[0]) + sbarscale[0];
		std::cout<<currentTotalCost<<std::endl;
		for(int i = 1; i <= T; i++){
			u = L[i-1] * x + l[i-1];
			Matrix<XDIM> tmpx;
			dyn.discrete_dynamics(x, u, tmpx, M);
			x = tmpx;
			double thisstepTotalCost = 0.5 * tr(~x * S[i] * x) + tr(~x * s[i]) + sscale[i]; //cost to go .
			thisstepTotalCost += 0.5 * tr(~x * SBar[i] * x) + tr(~x * sBar[i]) + sbarscale[i]; //cost so far.
			std::cout<<thisstepTotalCost<<std::endl;
		
		}

		return true;
	}

	inline bool checkTotalCost(const std::vector<SymmetricMatrix<XDIM>>& S, const std::vector<SymmetricMatrix<XDIM>>& SBar,
							   const std::vector<Matrix<XDIM>>& s, const std::vector<Matrix<XDIM>>& sBar,
							   const std::vector<double>& sscale, const std::vector<double>& sbarscale,
							   const std::vector<Matrix<XDIM>> msequenceStates){
	
		for(int i = 0; i <= T; i++){
			Matrix<XDIM> x = msequenceStates[i];
			double thisstepTotalCost = 0.5 * tr(~x * S[i] * x) + tr(~x * s[i]) + sscale[i]; //cost to go .
			thisstepTotalCost += 0.5 * tr(~x * SBar[i] * x) + tr(~x * sBar[i]) + sbarscale[i]; //cost so far.
			if(i%4 == 0)
				std::cout<<thisstepTotalCost<<std::endl;
		}
		return true;
	}


	inline double reldisTwoSequence(const std::vector<Matrix<XDIM>>& s1, const std::vector<Matrix<XDIM>>& s2)
	{
		double dis = 0.0;
		double diss2 = 0.0;
		for(int i = 0; i <= T; i++){
			dis += sqrt(tr(~(s1[i] - s2[i])*(s1[i] - s2[i])));
			diss2 += sqrt(tr(~(s2[i])*(s2[i])));
		}
		return dis / diss2;
	}


	inline double computeCumulativeDiffTotalCost(const std::vector<Matrix<XDIM>>& s1, 
												const std::vector<SymmetricMatrix<XDIM>>& S, const std::vector<SymmetricMatrix<XDIM>>& SBar,
												const std::vector<Matrix<XDIM>>&s, const std::vector<Matrix<XDIM>>& sBar,
												const std::vector<double>& sscale, const std::vector<double>& sbarscale)
	{
		double cdifftcost = 0.0;
		for(int i = 0; i <= T-1; i++){
			double tc1 = 0.0;
			double tc2 = 0.0;
			Matrix<XDIM> x(s1[i]); Matrix<XDIM> xn(s1[i+1]);
			tc1 = 0.5 * tr(~x * (S[i] + SBar[i]) * x) + tr(~x*(s[i] + sBar[i])) + sscale[i] + sbarscale[i];
			tc2 = 0.5 * tr(~xn * (S[i+1] + SBar[i+1]) * xn) + tr(~xn * (s[i+1] + sBar[i+1])) + sscale[i+1] + sbarscale[i+1];
			cdifftcost += abs(tc1 - tc2);
		}
		return cdifftcost;
	}


	inline double ExpectedCost(const std::vector<std::pair<Matrix<XDIM>, Matrix<UDIM>>>& plan, 
							   const std::vector<Matrix<UDIM, XDIM>>& L)
	{
		Dynamics dyn;
		Costfunction cf(start, goal, cgoal, lambda);

		SymmetricMatrix<XDIM> S = zeros<XDIM>();
		Matrix<XDIM, 1> svec = zeros<XDIM, 1>();
		double sscale = 0.0;

		SymmetricMatrix<XDIM> QT; Matrix<XDIM, 1> qvecT; double qT = 0;
		cf.QuadratizeCostcT(plan[T].first, plan[T].second, QT, qvecT, qT);
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
			
			if(t > 0)
				//cf.QuadratizeCostt_ov(plan[t].first, plan[t].second, cal_l1, cal_l2, cal_obs, Qt, Pt, Rt, qvect, rvect, qscalet);
				cf.QuadratizeCostt(Qt, Pt, Rt, qvect, rvect, qscalet, 1);
			else if(t == 0)
				cf.QuadratizeCost0(Qt, Pt, Rt, qvect, rvect, qscalet);

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

	inline void computePlan(std::vector<std::pair<Matrix<XDIM>, Matrix<UDIM>>>& plan, 
							const std::vector<Matrix<UDIM,XDIM>>& L,
							const std::vector<Matrix<UDIM,1>>& l,
							const Matrix<XDIM>& xHat)
	{
		plan.clear();
		plan.resize(T+1);
		Dynamics dyn;
		Matrix<XDIM,XDIM> M;
		Matrix<XDIM> x(xHat);  //using the current xHat instead of initial state.
		for(int i = 0; i <= T-1; i++){
			plan[i].first = x;
			plan[i].second = L[i] * plan[i].first + l[i];
			Matrix<XDIM> tmpx;
			dyn.discrete_dynamics(x, plan[i].second, tmpx, M);
			x = tmpx;
		}
		plan[T].first = x;
		plan[T].second = zeros<UDIM, 1>();
		
		return;
	}



	////////////////////Combined////////////////////////////////////////////////////////
	inline void extendedLQG(std::vector<Matrix<UDIM,XDIM> >& L, std::vector<Matrix<UDIM> >& l)
	{
		int maxIter = 1000;
		L.resize(T, zeros<UDIM,XDIM>());
		l.resize(T, zeros<UDIM,1>());

		//initialize a backward pass, so Set cost-to-go functions to zeros.
		std::vector<SymmetricMatrix<XDIM> > S(T + 1, zeros<XDIM>());
		std::vector<Matrix<XDIM> > s(T + 1, zero<XDIM>());
		std::vector<double> sscale(T+1, 0.0);

		std::vector<Matrix<XDIM>> tmpminsequenceStates(T+1, zeros<XDIM,1>());
		std::vector<std::pair<Matrix<XDIM>, Matrix<UDIM>>> plan; plan.clear();

		std::vector<SymmetricMatrix<XDIM> > SBar(T + 1);
		std::vector<Matrix<XDIM> > sBar(T + 1);
		std::vector<double> sbarscale(T+1);

		double oldCost = -log(0.0);
		Matrix<XDIM> xHat = start;
		//initialize the first step of the forward pass.
		SBar[0] = zeros<XDIM>();
		sBar[0] = zero<XDIM>();
		sbarscale[0] = 0.0;

		Dynamics dyn;
		Costfunction cf(start, goal, cgoal, lambda);

		Matrix<XDIM,XDIM> M;
		for (int iter = 0; iter < maxIter; ++iter) {
			// forward passs
			//xHat = start;
			for (int t = 0; t < T; ++t) {
				const Matrix<UDIM> uHat = L[t]*xHat + l[t];
				Matrix<XDIM> xHatPrime;
				dyn.discrete_dynamics(xHat, uHat, xHatPrime, M);  //xHatPrime: x_t+1

				Matrix<XDIM, XDIM> Abart; Matrix<XDIM, UDIM> Bbart; Matrix<XDIM, 1> abart;
				std::vector<Matrix<XDIM, XDIM>> Fbart; std::vector<Matrix<XDIM, UDIM>> Gbart; std::vector<Matrix<XDIM,1>> ebart;
				dyn.linearize_discrete_inverse_dynamics(xHatPrime, uHat, Abart, Bbart, abart, Fbart, Gbart, ebart);

				SymmetricMatrix<XDIM> Qt; Matrix<UDIM, XDIM> Pt; SymmetricMatrix<UDIM> Rt; Matrix<XDIM, 1> qvect; Matrix<UDIM,1> rvect; double qscalet;
				if(t == 0){
					cf.QuadratizeCost0(Qt, Pt, Rt, qvect, rvect, qscalet);
				}
				else
					cf.QuadratizeCostt(Qt, Pt, Rt, qvect, rvect, qscalet, iter); //this correponds to the non-obstacle case.
					//cf.QuadratizeCostt_ov(xHat, uHat, cal_l1, cal_l2, cal_obs, Qt, Pt, Rt, qvect, rvect, qscalet);

				SymmetricMatrix<XDIM> QtSbar = Qt + SBar[t];
				SymmetricMatrix<XDIM> Cbart = SymProd(~Abart, QtSbar * Abart);
				SymmetricMatrix<UDIM> Dbart = SymProd(~Bbart, QtSbar * Bbart) + SymSum(Pt*Bbart) + Rt;
				Matrix<UDIM, XDIM> Ebart = ~Bbart * QtSbar * Abart + Pt * Abart;
				Matrix<XDIM, 1> cbart = ~Abart * QtSbar *  abart + ~Abart * (qvect + sBar[t]);
				Matrix<UDIM, 1> dbart = ~Bbart * QtSbar * abart + Pt*abart + ~Bbart*(qvect + sBar[t]) + rvect;
				double escalebart = 0.5 * tr(~abart * QtSbar * abart) + tr(~abart*(qvect + sBar[t])) + qscalet + sbarscale[t];
				
				for(int i = 0; i < XDIM; i++){
					Cbart += SymProd(~Fbart[i], QtSbar	* Fbart[i]);
					Dbart += SymProd(~Gbart[i], QtSbar * Gbart[i]);
					Ebart += ~Gbart[i] * QtSbar * Fbart[i];
					cbart += ~Fbart[i] * QtSbar * ebart[i];
					dbart += ~Gbart[i] * QtSbar * ebart[i];
					escalebart += 0.5 * tr(~ebart[i] * QtSbar * ebart[i]);
				}
				
				L[t] = -(Dbart % Ebart);
				l[t] = -(Dbart % dbart);
				SBar[t+1] = Cbart + SymProd(~Ebart, L[t]);
				sBar[t+1] = cbart + ~Ebart * l[t];
				sbarscale[t+1] = escalebart + 0.5 * tr(~dbart * l[t]);
				xHat= -((S[t+1] + SBar[t+1]) % (s[t+1] + sBar[t+1]));    //check code helicopter.
				//std::cout<<xHat<<std::endl;
			}

			// backward pass
			cf.QuadratizeCostcT(xHat, zeros<UDIM,1>(), S[T], s[T], sscale[T]);
			xHat = -((S[T] + SBar[T])%(s[T] + sBar[T]));
			minsequenceStates[T] = xHat;
			//std::cout<<xHat<<std::endl;

			for (int t = T - 1; t != -1; --t) {
				const Matrix<UDIM> uHat = L[t]*xHat + l[t];
				Matrix<XDIM> xHatPrime;
				dyn.discrete_inverse_dynamics(xHat, uHat, xHatPrime, M);  //from x_t+1 to x_t

				Matrix<XDIM, XDIM> At; Matrix<XDIM, UDIM> Bt; Matrix<XDIM, 1> avect;
				std::vector<Matrix<XDIM, XDIM>> Ft;
				std::vector<Matrix<XDIM, UDIM>> Gt;
				std::vector<Matrix<XDIM, 1>> evect;
				dyn.linearize_discrete_dynamics(xHatPrime, uHat, At, Bt, avect, Ft, Gt, evect);

				SymmetricMatrix<XDIM> Qt; Matrix<UDIM, XDIM> Pt; SymmetricMatrix<UDIM> Rt;
				Matrix<XDIM, 1> qvect; Matrix<UDIM, 1> rvect;
				double qscalet = 0;

				if(t == 0){
					cf.QuadratizeCost0(Qt, Pt, Rt, qvect, rvect, qscalet);
				}
				else
					cf.QuadratizeCostt(Qt, Pt, Rt, qvect, rvect, qscalet, iter); //non obstacle case.
					//cf.QuadratizeCostt_ov(xHatPrime, uHat, cal_l1, cal_l2, cal_obs, Qt, Pt, Rt, qvect, rvect, qscalet); //obstacle case.
		
				SymmetricMatrix<XDIM> Ct = Qt + SymProd(~At, S[t+1] * At);
				SymmetricMatrix<UDIM> Dt = Rt + SymProd(~Bt, S[t+1] * Bt);
				Matrix<UDIM, XDIM> Et = Pt + ~Bt * S[t+1] * At;
				Matrix<XDIM, 1> cvect = qvect + ~At * S[t+1] * avect + ~At * s[t+1];
				Matrix<UDIM, 1> dvect = rvect + ~Bt * S[t+1] * avect + ~Bt * s[t+1];
				double escalet = qscalet + sscale[t+1] + 0.5 * tr(~avect * S[t+1] * avect) + tr(~avect * s[t+1]);
					
				for(int i = 0; i < XDIM; i++){
					Ct += SymProd(~Ft[i], S[t+1] * Ft[i]);
					Dt += SymProd(~Gt[i], S[t+1] * Gt[i]);
					Et += ~Gt[i] * S[t+1] * Ft[i];
					cvect += ~Ft[i] * S[t+1] * evect[i];
					dvect += ~Gt[i] * S[t+1] * evect[i];
					escalet += 0.5 * tr(~evect[i] * S[t+1] * evect[i]);
				}

				L[t] = -(Dt%Et);
				l[t] = -(Dt%dvect);

				S[t] = Ct + SymProd(~Et, L[t]);
				s[t] = cvect + ~Et * l[t];
				sscale[t] = escalet + 0.5 * tr(~dvect * l[t]);
				xHat = -((S[t] + SBar[t])%(s[t] + sBar[t]));
				minsequenceStates[t] = xHat;
				//std::cout<<xHat<<std::endl;
			}

			//Start terminiation criteria. Currently, only consider the approximated cost (not the exact expected cost)
			// compute cost
			computePlan(plan, L, l, xHat);
			double newCost = ExpectedCost(plan, L);
			if(abs(newCost - oldCost) / oldCost < 1e-4){
				numTotalIter = iter;
				totalCost = newCost;
				checkTotalCost(S, SBar, s, sBar, sscale, sbarscale, minsequenceStates);
				return;
			}
			else
				oldCost = newCost;

			//std::cout<<oldCost<<std::endl;
			//////////////////Third trial: convergence of the min-total cost sequence states.
			/*double dismin = 0.0;
			double diff = reldisTwoSequence(tmpminsequenceStates, minsequenceStates);
			//std::cout<<diff<<std::endl;
			if(abs(diff) < 1e-4){
				numTotalIter = iter;
				totalCost = 0.5 * tr(~xHat * S[0] * xHat) + tr(~xHat * s[0]) + sscale[0];
				checkTotalCost(S, SBar, s, sBar, sscale, sbarscale, minsequenceStates);
				return;
			}
			else{
				tmpminsequenceStates.clear(); tmpminsequenceStates = minsequenceStates;
			}*/
			/////////////////////////////////////////////////////////////////////////////////////

		}
		return;
	}

	inline void executeExtendedLQR()
	{
		CAL_SuspendVisualisation();

		std::vector<Matrix<UDIM,XDIM> > L; L.clear();
		std::vector<Matrix<UDIM> > l; l.clear();
		Costfunction cf(start, goal, cgoal, lambda);

		clock_t startT = clock();
		extendedLQG(L, l);
		totalTime = (double)(clock() - startT) / CLOCKS_PER_SEC;

		Dynamics dyn;
		Matrix<XDIM,1> x = start;
		Matrix<XDIM,XDIM> M;
		nominalPlan.clear(); nominalPlan.resize(T+1);

		//compute nominal plan, this time, start from the initial state.
		for(int i = 0; i <= T-1; i++){
			nominalPlan[i].first = x;
			nominalPlan[i].second = L[i] * nominalPlan[i].first + l[i];
			Matrix<XDIM> tmpx;
			dyn.discrete_dynamics(x, nominalPlan[i].second, tmpx, M);
			x = tmpx;
		}

		nominalPlan[T].first = x;
		nominalPlan[T].second = zeros<UDIM, 1>();

		CAL_ResumeVisualisation();

		return;
	}

};


#endif