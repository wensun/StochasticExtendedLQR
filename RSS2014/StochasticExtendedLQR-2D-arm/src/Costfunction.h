#ifndef _COST_
#define _COST_

#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"

class Costfunction{

public:
	Matrix<XDIM, 1> start; //configuration state.
	Matrix<2, 1> goal; //x,y 2d position;

	double K; //used in the exp function.
	double lambda; //used in the expfunction
	double alpha; //used in the exp
	double coeu; //used for u
	double coegoal; //used for goal.
	double coestart;
	double d;
	double rotCost;

	Matrix<XDIM> cgoal;

	Costfunction(const Matrix<XDIM,1>& s, const Matrix<2,1>& g, const Matrix<XDIM, 1>& cg,
				 const double& Lambda){
		start = s;
		goal = g;

		coeu = 0.01;
		coegoal = 10.0;
		coestart = 10.0;
		K = 1.0;
		//lambda = 0.015;
		lambda = Lambda;
		alpha = 0.1;

		rotCost = 0.1;

		d = 0.0078125;
		//d = 0.0009765625;

		cgoal = cg;

		//cgoal(0,0) = 3.14/2.0; cgoal(1,0) = 0.0; cgoal(2,0) = 0.0;
		//cgoal(3,0) = 0.0;
	}

	
	inline double setLambda(const double lam)
	{
		lambda = lam;
	}


	inline Matrix<2,1> ForwardK(const Matrix<XDIM, 1>& x)
	{
		Matrix<2> pos; pos.reset();
		Dynamics dyn;
		pos(0,0) = dyn.l(0,0)*cos(x(0,0)) + dyn.l(1,0) * cos(x(0,0) + x(1,0));
		pos(1,0) = dyn.l(0,0)*sin(x(0,0)) + dyn.l(1,0) * sin(x(0,0) + x(1,0));
		
		return pos;
	}



	//the cost for the final step, which is non quadratic.
	inline double costT(const Matrix<XDIM,1>& x)
	{
		double cost = 0;
		Matrix<2> pos; pos.reset();
		//forward kinematics;
		pos = ForwardK(x);
	
		cost += 0.5 * coegoal * tr(~(pos - goal)*(pos - goal)); //distance to goal
		cost += 0.5 * coeu * 10 * tr(~x.subMatrix<2,1>(2,0) * x.subMatrix<2,1>(2,0)); //need the final angular velocity to be zero.

		return cost;
	}


	//the exact final configuration position penality.
	inline double costcT(const Matrix<XDIM, 1>& x)
	{
		double cost = 0;
		cost += 0.5 * coegoal * tr(~(x - cgoal) * (x - cgoal));

		return cost;
	}

	//cost for the first term, needs to make sure that the start position 
	inline double cost0(const Matrix<XDIM,1>& x, const Matrix<UDIM,1>& u)
	{
		double cost = 0;
		cost += 0.5 * coestart * tr(~(start - x)*(start - x));
		cost += 0.5 * coeu * tr(~u*u);

		return cost;
	}

	//cost in the middle stage.
	inline double costt(const Matrix<XDIM,1>& x, const Matrix<UDIM,1>& u)
	{
		return 0.5 * coeu * tr(~u*u);
	}


	inline void QuadratizeCost0(SymmetricMatrix<XDIM>& Q0, Matrix<UDIM, XDIM>& P0, SymmetricMatrix<UDIM>& R0, 
								  Matrix<XDIM, 1>& qvec0, Matrix<UDIM, 1>& rvec0,
								  double&  q0)
	{
		Q0 = identity<XDIM>() * coestart;
		P0 = zeros<UDIM,XDIM>();
		R0 = identity<UDIM>() * coeu;

		qvec0 = -coestart * start;
		rvec0 = zeros<UDIM,1>();

		q0 = 0.5 * coestart * tr(~start * start);

		return;
	}

	inline void QuadratizeCostcT(const Matrix<XDIM,1>& xstar, const Matrix<UDIM, 1>& ustar, 
								SymmetricMatrix<XDIM>& QT,
								Matrix<XDIM, 1>& qvecT,
								double&  qT)
	{

		QT.reset();
		QT = coegoal * identity<XDIM>();
		qvecT = -coegoal * cgoal;
		qT = 0.5 * coegoal * tr(~cgoal * cgoal);

		return;
	}


	inline void QuadratizeCostt(SymmetricMatrix<XDIM>& Qt, Matrix<UDIM, XDIM>& Pt, SymmetricMatrix<UDIM>& Rt, 
								  Matrix<XDIM, 1>& qvect, Matrix<UDIM, 1>& rvect,
								  double&  qt,
								  const int& iter)
	{
		Qt = zeros<XDIM>();

		if(iter < 1)
			Qt.insert(2, rotCost * identity<2>());

		Pt = zeros<UDIM, XDIM>();
		Rt = coeu * identity<UDIM>();

		qvect = zeros<XDIM,1>();
		rvect = zeros<UDIM,1>();
		qt = 0;

		return;
	}


	inline void QuadratizeCostT_real(const Matrix<XDIM,1>& xstar, const Matrix<UDIM, 1>& ustar, 
								SymmetricMatrix<XDIM>& QT,
								Matrix<XDIM, 1>& qvecT,
								double&  qT)
	{
		QT = zeros<XDIM>();
		Matrix<XDIM, 1> jac = zeros<XDIM,1>();
		Dynamics dyn; double l1 = dyn.l(0,0); double l2 = dyn.l(1,0);
		double theta1 = xstar(0,0); double theta2 = xstar(1,0);
		jac(0,0) = coegoal * (l1*cos(theta1) + l2*cos(theta1 + theta2) - goal[0]) * (-l1*sin(theta1) - l2*sin(theta1+theta2))
					+ coegoal * (l1*sin(theta1) + l2*sin(theta1 + theta2) - goal[1]) * (l1*cos(theta1) + l2*cos(theta1+theta2));

		QT(0,0) = coegoal * (-l1*sin(theta1) - l2*sin(theta1 + theta2)) * (-l1*sin(theta1) - l2*sin(theta1+theta2))
				+ coegoal * (l1*cos(theta1) + l2*cos(theta1 + theta2) - goal[0]) * (-l1*cos(theta1) - l2*cos(theta1+theta2))
				+ coegoal * (l1*cos(theta1) + l2*cos(theta1 + theta2)) * (l1*cos(theta1) + l2*cos(theta1+theta2))
				+ coegoal * (l1*sin(theta1) + l2*sin(theta1 + theta2) - goal[1]) * (-l1*sin(theta1) - l2*sin(theta1+theta2));

		jac(1,0) = coegoal * (l1*cos(theta1) + l2*cos(theta1 + theta2) - goal[0]) * (-l2 * sin(theta1 + theta2)) 
				 + coegoal * (l1*sin(theta1) + l2*sin(theta1 + theta2) - goal[1]) * (l2 * cos(theta1 + theta2));

		QT(1,1) = coegoal * (-l2*sin(theta1 + theta2)) * (-l2 * sin(theta1 + theta2))
				 + coegoal * (l1*cos(theta1) + l2*cos(theta1 + theta2) - goal[0]) * (-l2 * cos(theta1 + theta2))
				 + coegoal * (l2*cos(theta1 + theta2)) * (l2 * cos(theta1 + theta2))
				 + coegoal * (l1*sin(theta1) + l2*sin(theta1 + theta2) - goal[1]) * (-l2 * sin(theta1 + theta2));

		QT(0,1) = coegoal * (-l2*sin(theta1 + theta2)) * (-l1*sin(theta1) - l2*sin(theta1+theta2))
				+ coegoal * (l1*cos(theta1) + l2*cos(theta1 + theta2) - goal[0]) * (-l2*cos(theta1+theta2))
				+ coegoal * (l2*cos(theta1 + theta2)) * (l1*cos(theta1) + l2*cos(theta1+theta2))
				+ coegoal * (l1*sin(theta1) + l2*sin(theta1 + theta2) - goal[1]) * (-l2*sin(theta1+theta2));


		QT(1,0) = QT(0,1);

		jac(2,0) = coeu * 10 * xstar(2,0);
		jac(3,0) = coeu * 10 * xstar(3,0);
		QT(2,2) = coeu * 10; QT(3,3) = coeu * 10;
		QT(2,3) = QT(3,2) = 0.0;

		regularize(QT);

		qvecT = jac - QT * xstar;

		qT = costT(xstar) - tr(~jac * xstar) + 0.5 * tr(~xstar * QT * xstar);		
	
		return;
	}


	inline void QuadratizeCostT(const Matrix<XDIM,1>& xstar, const Matrix<UDIM, 1>& ustar, 
								SymmetricMatrix<XDIM>& QT,
								Matrix<XDIM, 1>& qvecT,
								double&  qT)
	{
		//0.0078125
		//double d = 0.007765625;
		QT.reset();
		for(int i = 0; i < XDIM; i++){
			Matrix<XDIM, 1> augr = xstar;
			Matrix<XDIM, 1> augl = xstar;
			augr(i,0) += d;
			augl(i,0) -= d;
			QT(i,i) = (costT(augl) - 2.0*costT(xstar) + costT(augr)) / (d*d);
		}

		Matrix<XDIM> btr(xstar), btl(xstar), bbr(xstar), bbl(xstar);
		for (int i = 1; i < XDIM; ++i) {
			btr[i] += d; btl[i] -= d; bbr[i] += d; bbl[i] -= d;
			for (int j = 0; j < i; ++j) {
				btr[j] += d; btl[j] += d; bbr[j] -= d; bbl[j] -= d;
				QT(i,j) = QT(j,i) = (costT(bbl) + costT(btr) - costT(btl) - costT(bbr)) / (4.0*d*d); // O(n2 * n2 * n3) = O(n7)
				btr[j] = btl[j] = bbr[j] = bbl[j] = xstar[j];
			}
			btr[i] = btl[i] = bbr[i] = bbl[i] = xstar[i];
		}

		regularize(QT);

		Matrix<XDIM,1> jac; jac.reset();
		for(int i = 0; i < XDIM; i++){
			Matrix<XDIM,1> augrx(xstar), auglx(xstar);
			augrx(i, 0) += d; auglx(i,0) -= d;
			jac(i,0) = (costT(augrx) - costT(auglx)) / (2*d);
		}

		qvecT = jac - QT * xstar;

		qT = costT(xstar) - tr(~jac * xstar) + 0.5 * tr(~xstar * QT * xstar);		
	
		return;
	}


	/***********************************************************************************************/
	/****************************costt function with obstacles avoidance****************************/
	inline double costt_ov(const Matrix<XDIM>& state, const Matrix<UDIM>& u, const int& cal_link1, const int& cal_link2, const int& cal_obstacles){
		
		setTwoLinks(cal_link1, cal_link2, state);

		double distance1, distance2;
		double cost = 0.0;
		distance1 = 0.0; distance2 = 0.0;

		int col1 = -1;
		CAL_CheckGroupCollision(cal_link1, cal_obstacles, false, &col1);
		if(col1 == 0){ //no collision
			int num_pairs = 0;
			CAL_GetClosestPairs(cal_link1, cal_obstacles, &num_pairs);
			if(num_pairs > 0){
				SCALResult* results = new SCALResult[num_pairs];
				CAL_GetResults(results);
				distance1 = results[0].distance;
			}
		}
		else{
			int num_pairs = 0;
			CAL_GetPenetrationDepths(cal_link1, cal_obstacles, &num_pairs);
			//std::cout<<num_pairs<<std::endl;
			if(num_pairs > 0){
				SCALResult* results = new SCALResult[num_pairs];
				CAL_GetResults(results);
				distance1 = -results[0].distance; //negative distance for penetration.
			}
		}
		cost += exp(-K*distance1 + alpha);

		int col2 = -1;
		CAL_CheckGroupCollision(cal_link2, cal_obstacles, false, &col2);
		if(col2 == 0){ //no collision
			int num_pairs = 0;
			CAL_GetClosestPairs(cal_link2, cal_obstacles, &num_pairs);
			if(num_pairs > 0){
				SCALResult* results = new SCALResult[num_pairs];
				CAL_GetResults(results);
				distance2 = results[0].distance;
			}
		}
		else{
			int num_pairs = 0;
			CAL_GetPenetrationDepths(cal_link2, cal_obstacles, &num_pairs);
			//std::cout<<num_pairs<<std::endl;
			if(num_pairs > 0){
				SCALResult* results = new SCALResult[num_pairs];
				CAL_GetResults(results);
				distance2 = -results[0].distance; //negative distance for penetration.
			}
		}
		cost += exp(-K*distance2 + alpha);
		cost *= lambda;

		cost += costt(state, u);

		return cost;
	}

	inline void QuadratizeCostt_ov(const Matrix<XDIM>& xstar, const Matrix<UDIM>& u, const int& cal_link1, const int& cal_link2, const int& cal_obstacles,
								  SymmetricMatrix<XDIM>& Qt, Matrix<UDIM, XDIM>& Pt, SymmetricMatrix<UDIM>& Rt, 
								  Matrix<XDIM, 1>& qvect, Matrix<UDIM, 1>& rvect,
								  double&  qt)
	{
		setTwoLinks(cal_link1, cal_link2, xstar);
		Pt = zeros<UDIM, XDIM>();
		Rt = coeu * identity<UDIM>();
		rvect = zeros<UDIM,1>();

		//compute contact points and norms.
		double distance1, distance2;
		Matrix<2> p1l, p1o, p2l, p2o, n1, n2;

		int col1 = -1;
		CAL_CheckGroupCollision(cal_link1, cal_obstacles, false, &col1);
		if(col1 == 0){ //no collision
			int num_pairs = 0;
			CAL_GetClosestPairs(cal_link1, cal_obstacles, &num_pairs);
			if(num_pairs > 0){
				SCALResult* results = new SCALResult[num_pairs];
				CAL_GetResults(results);
				distance1 = results[0].distance;
				p1l[0] = results[0].vector0[0]; p1l[1] = results[0].vector0[1];
				p1o[0] = results[0].vector1[0]; p1o[1] = results[0].vector1[1];
				n1 = (p1l - p1o) / sqrt(tr(~(p1l-p1o)*(p1l-p1o)) + 0.00000001);
				//std::cout<<tr(~n1 * (p1l - p1o))<<std::endl<<distance1<<std::endl;
			}
		}
		else{
			int num_pairs = 0;
			CAL_GetPenetrationDepths(cal_link1, cal_obstacles, &num_pairs);
			if(num_pairs > 0){
				SCALResult* results = new SCALResult[num_pairs];
				CAL_GetResults(results);
				distance1 = -results[0].distance; //negative distance for penetration.
				p1l[0] = results[0].vector0[0]; p1l[1] = results[0].vector0[1];
				p1o[0] = results[0].vector1[0]; p1o[1] = results[0].vector1[1];
				n1 = (-p1l + p1o) / sqrt(tr(~(p1l-p1o)*(p1l-p1o)) + 0.0000001);
				//std::cout<<tr(~n1 * (p1l - p1o))<<std::endl<<distance1<<std::endl;
			}
		}

		int col2 = -1;
		CAL_CheckGroupCollision(cal_link2, cal_obstacles, false, &col2);
		if(col2 == 0){ //no collision
			int num_pairs = 0;
			CAL_GetClosestPairs(cal_link2, cal_obstacles, &num_pairs);
			if(num_pairs > 0){
				SCALResult* results = new SCALResult[num_pairs];
				CAL_GetResults(results);
				distance2 = results[0].distance;
				p2l[0] = results[0].vector0[0]; p2l[1] = results[0].vector0[1];
				p2o[0] = results[0].vector1[0]; p2o[1] = results[0].vector1[1];
				n2 = (p2l - p2o) / sqrt(tr(~(p2l-p2o)*(p2l-p2o)) + 0.00000001);
				//std::cout<<tr(~n2 * (p2l - p2o))<<std::endl<<distance2<<std::endl;
			}
		}
		else{
			int num_pairs = 0;
			CAL_GetPenetrationDepths(cal_link2, cal_obstacles, &num_pairs);
			//std::cout<<num_pairs<<std::endl;
			if(num_pairs > 0){
				SCALResult* results = new SCALResult[num_pairs];
				CAL_GetResults(results);
				distance2 = -results[0].distance; //negative distance for penetration.
				p2l[0] = results[0].vector0[0]; p2l[1] = results[0].vector0[1];
				p2o[0] = results[0].vector1[0]; p2o[1] = results[0].vector1[1];
				n2 = (-p2l + p2o) / sqrt(tr(~(p2l-p2o)*(p2l-p2o)) + 0.000001);
				//std::cout<<tr(~n2 * (p2l - p2o))<<std::endl<<distance2<<std::endl;
			}
		}

		//compute FK jacobian
		Matrix<2, XDIM> J1 = zeros<2,XDIM>();
		Dynamics dyn;
		J1 = dyn.Jacobian_FK(xstar, 1, sqrt(tr(~p1l*p1l)));
		Matrix<2, XDIM> J2 = zeros<2,XDIM>();
		Matrix<2> middle = dyn.FK(xstar, 1, dyn.l(0,0));
		J2 = dyn.Jacobian_FK(xstar, 2, sqrt(tr(~(p2l - middle)*(p2l - middle))));

		//compute Qt
		Qt = lambda * K*K * (exp(-K*distance1 + alpha) * SymProd(~J1 * n1, ~n1 * J1) + exp(-K*distance2 + alpha) * SymProd(~J2 * n2, ~n2 * J2));
		//regularize(Qt);
		
		Matrix<XDIM> tmpja = -lambda * K * ~(exp(-K*distance1 + alpha) * ~n1 * J1 + exp(-K*distance2 + alpha) * ~n2 * J2);

		qvect = tmpja - Qt * xstar;

		qt = costt_ov(xstar, u, cal_link1, cal_link2, cal_obstacles) - tr(~tmpja * xstar) + 0.5 * tr(~xstar * Qt * xstar);

		return;
	}


	inline void QuadratizeCostt_ov_finite_d(const Matrix<XDIM>& xstar, const Matrix<UDIM>& u, const int& cal_link1, const int& cal_link2, const int& cal_obstacles,
								  Matrix<XDIM, XDIM>& Qt, Matrix<UDIM, XDIM>& Pt, Matrix<UDIM, UDIM>& Rt, 
								  Matrix<XDIM, 1>& qvect, Matrix<UDIM, 1>& rvect,
								  double&  qt)
	{
		//double d = 0.0009765625;
		Pt = zeros<UDIM, XDIM>();
		Rt = coeu * identity<UDIM>();
		rvect = zeros<UDIM,1>();

		for(int i = 0; i < XDIM; i++){
			for(int j = i; j < XDIM; j++){
				Matrix<XDIM,1> augx1 = xstar;
				Matrix<XDIM,1> augx2 = xstar;
				Matrix<XDIM,1> augx3 = xstar;
				Matrix<XDIM,1> augx4 = xstar;
				augx1(i,0) += d;
				augx1(j,0) += d;
				augx2(j,0) += d;
				augx3(i,0) += d;
				Qt(i,j) = (costt_ov(augx1,u, cal_link1, cal_link2, cal_obstacles) - costt_ov(augx2, u, cal_link1, cal_link2, cal_obstacles) - costt_ov(augx3,u,cal_link1, cal_link2, cal_obstacles) + costt_ov(augx4,u,cal_link1, cal_link2, cal_obstacles)) / d / d;
				Qt(j,i) = Qt(i,j);
			}
		}

		Matrix<XDIM,1> jac; jac.reset();
		for(int i = 0; i < XDIM; i++){
			Matrix<XDIM,1> augx = xstar;
			augx(i, 0) += d;
			jac(i,0) = (costt_ov(augx,u, cal_link1, cal_link2, cal_obstacles) - costt_ov(xstar, u, cal_link1, cal_link2, cal_obstacles)) / d;
		}
		qvect = jac - Qt * xstar;

		qt = costt_ov(xstar, u, cal_link1, cal_link2, cal_obstacles) - tr(~jac * xstar) + 0.5 * tr(~xstar * Qt * xstar);

		return;
	}


};


#endif