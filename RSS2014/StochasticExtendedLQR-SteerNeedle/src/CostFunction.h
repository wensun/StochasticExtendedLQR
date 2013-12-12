#ifndef _COSTFUNCTION_
#define _COSTFUNCTION_

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

class CostFunction{

public :

	Matrix<BDIM> bstart;
	Matrix<XDIM> xGoal;
	double kmax;
	Matrix<UDIM, UDIM> Rint;
	Matrix<XDIM, XDIM> Qint;
	Matrix<XDIM, XDIM> QGoal;
	Matrix<UDIM> uNominal;

	double obsscale;
	double k;

	double H;

	int cal_p;
	int cal_env;
	int cal_obs;

	CostFunction(const Matrix<BDIM>& _bstart, const Matrix<XDIM>& _goal,
				const int& _cal_point, const int& _cal_environment, const int& _cal_obstacle)
	{
		bstart = _bstart;
		xGoal = _goal;

		kmax = 5.0;  //0.2

		Rint = identity<UDIM>() * 1.0;
		Qint = identity<XDIM>() * 1.0;
		QGoal = identity<XDIM>();
		QGoal.insert<3,3>(0,0, identity<3>() * 50); //penalize the 3 d position. 

		H = 0.0078;

		uNominal = zeros<UDIM,1>();
		uNominal[2] = 0.7 * kmax;
		
		k = 1.0;
		obsscale = 5.0;

		cal_p = _cal_point;
		cal_env = _cal_environment;
		cal_obs = _cal_obstacle;

	}

	//compute number of standard deviation.
	inline double obstacleConfid(const Matrix<BDIM>& b, const int& cal_environment, const int& cal_obstacles, const int& cal_point)
	{
		Matrix<XDIM> x; Matrix<XDIM, XDIM> sqrtSigma;
		unVec(b, x, sqrtSigma);
		Matrix<3> p = x.subMatrix<3,1>(0,0);
		Matrix<3,3> Sigma = (sqrtSigma*sqrtSigma).subMatrix<3,3>(0,0); //x, y, z
		regularize(Sigma);
		SymmetricMatrix<3> symSigma = SymProd(Sigma, identity<3>());
		Matrix<3,3> V;
		SymmetricMatrix<3> D;
		jacobi(symSigma, V, D);
		Matrix<3> scale;
		scale[0] = sqrt(D(0,0)); scale[1] = sqrt(D(1,1)); scale[2] = sqrt(D(2,2));

		Matrix<3,3> invScale = zeros<3,3>();
		invScale(0,0) = 1.0/scale[0];
		invScale(1,1) = 1.0/scale[1];
		invScale(2,2) = 1.0/scale[2];

		Matrix<4,1> q1 = quatFromRot(~V);
		CAL_SetGroupQuaternion(cal_obstacles, (float)q1[0], (float)q1[1], (float)q1[2], (float)q1[3]);
		CAL_SetGroupScaling(cal_environment, (float)invScale(0,0), (float)invScale(1,1), (float)invScale(2,2));

		Matrix<3> tp;
		tp = invScale * ~V * p; //Trans p.
		CAL_SetGroupPosition(cal_point, tp[0], tp[1], tp[2]);

		double distance = 0.0;
		int col = 0;
		CAL_CheckGroupCollision(cal_point, cal_obstacles, false, &col);
		if(col != 0) //in collision
		{
			int num_pairs;
			CAL_GetPenetrationDepths(cal_point, cal_environment, &num_pairs);
			SCALResult* results = new SCALResult[num_pairs];
			CAL_GetResults(results);
			distance = -results[0].distance;
			delete[] results;
		}
		else // non-collision
		{
			int num_pairs;
			CAL_GetClosestPairs(cal_point, cal_environment, &num_pairs);
			SCALResult* results = new SCALResult[num_pairs];
			CAL_GetResults(results);
			distance = results[0].distance;
			delete[] results;
		}

		CAL_SetGroupQuaternion(cal_obstacles, 0.0f, 0.0f, 0.0f, 1.0f);
		CAL_SetGroupScaling(cal_environment, 1.0f, 1.0f, 1.0f);

		return distance;
	}

	inline double obstacleCost(const Matrix<BDIM>& b, const int& cal_environment, const int& cal_obstacle, const int& cal_point)
	{
		double signdis = obstacleConfid(b, cal_environment, cal_obstacle, cal_point);
		return k * exp(-obsscale * signdis);
	}

	inline void quadratizeObsCost(const Matrix<BDIM>& b, const int& cal_environment, const int& cal_obstacle, const int& cal_point,
								Matrix<BDIM, BDIM>& Hessian, Matrix<BDIM>& jac)
	{
		double d0 = obstacleConfid(b, cal_environment, cal_obstacle, cal_point);
		double a = obsscale * obsscale * k * exp(-obsscale * d0);
		double c = -obsscale * k * exp(-obsscale * d0);

		Matrix<BDIM> obsJacob;
		Matrix<BDIM> br(b); Matrix<BDIM> bl(b);
		for(int i = 0; i < BDIM; i++){
			br = b; bl = b;
			br[i] += H; bl[i] -= H;
			double dr = obstacleConfid(br, cal_environment, cal_obstacle, cal_point);
			double dl = obstacleConfid(bl, cal_environment, cal_obstacle, cal_point);
			obsJacob[i] = (dr - dl) / (2.0 * H);
		}

		Hessian = a * obsJacob * ~obsJacob;
		jac = c * obsJacob;

	}

	//quadratize the cost part that is not related to the obstacle avoidance.
	inline void quadratizeCost(const Matrix<BDIM>& b, const Matrix<UDIM>& u,
								Matrix<BDIM, BDIM>& bH, Matrix<UDIM, BDIM>& P, Matrix<UDIM,UDIM>& uH,
								Matrix<BDIM, 1>& bjac, Matrix<UDIM,1>& ujac)
	{
		// diag(Q). q
		double p = cost(b,u);
		Matrix<BDIM> br(b), bl(b);
		for (int i = 0; i < BDIM; ++i) {
			br[i] += H; bl[i] -= H; 
			bjac[i] = (cost(br,u) - cost(bl,u)) / (2.0*H);      // O(n2 * n3) = O(n5)
			bH(i,i) = ( cost(bl,u) - 2.0*p + cost(br, u) ) / (H*H);
			br[i] = bl[i] = b[i];
		}

		Matrix<BDIM> btr(b), btl(b), bbr(b), bbl(b);
		for (int i = 1; i < BDIM; ++i) {
			btr[i] += H; btl[i] -= H; bbr[i] += H; bbl[i] -= H;
			for (int j = 0; j < i; ++j) {
				btr[j] += H; btl[j] += H; bbr[j] -= H; bbl[j] -= H;
				bH(i,j) = bH(j,i) = (cost(bbl,u) + cost(btr,u) - cost(btl,u) - cost(bbr,u)) / (4.0*H*H); // O(n2 * n2 * n3) = O(n7)
				btr[j] = btl[j] = bbr[j] = bbl[j] = b[j];
			}
			btr[i] = btl[i] = bbr[i] = bbl[i] = b[i];
		}

		uH = Rint;
		ujac = Rint * (u - uNominal);

		P = zeros<UDIM,BDIM>();
		return;
	}

	//no obstacle part.
	inline double cost(const Matrix<BDIM>& b, const Matrix<UDIM>& u)
	{
		Matrix<XDIM> x; Matrix<XDIM, XDIM> sqrtSigma;
		unVec(b, x, sqrtSigma);
		return 0.5*tr(~(u - uNominal)*Rint*(u - uNominal)) + 0.5*tr(~sqrtSigma*Qint*sqrtSigma);
	}

	//combined cost function non-obstacle part + obstacle_avoidance part.
	inline double cost_ov(const Matrix<BDIM>& b, const Matrix<UDIM>& u)
	{
		double c = 0.0;
		c += cost(b, u);
		c += obstacleCost(b, cal_env, cal_obs, cal_p);
		return c;
	}

	inline void quadratizeCost_ov(const Matrix<BDIM>& b, const Matrix<UDIM>& u,
									Matrix<BDIM, BDIM>& Q, Matrix<UDIM, BDIM>& P, Matrix<UDIM,UDIM>& R,
									Matrix<BDIM>& q, Matrix<UDIM>& r, double& qscalar)
	{
		Matrix<BDIM, BDIM> obsHessian; Matrix<BDIM> obsJacob;
		quadratizeObsCost(b, cal_env, cal_obs, cal_p, obsHessian, obsJacob);

		Matrix<BDIM, BDIM> H; Matrix<BDIM> jac; Matrix<UDIM,UDIM> uH; Matrix<UDIM> ujac;
		quadratizeCost(b, u, H, P, uH, jac, ujac);

		Q = H + obsHessian;
		R = uH;
		regularize(Q);
		q = (obsJacob + jac) - Q * b;
		r = (ujac) - R * u;
		
		double c0 = cost_ov(b, u);
		qscalar = c0 + 0.5 * tr(~b * Q * b) + 0.5 * tr(~u * R * u) - tr(~b * (q + Q * b)) - tr(~u * (r + R * u));
		P = zeros<UDIM, BDIM>();

		return;
	}


	//deal with the final cost step.
	inline double finalcost(const Matrix<BDIM>& b)
	{
		Matrix<XDIM> x; Matrix<XDIM, XDIM> SqrtSigma;
		unVec(b, x, SqrtSigma);
		return 0.5*tr(~(x - xGoal)*QGoal*(x - xGoal)) + 0.5*tr(~SqrtSigma*QGoal*SqrtSigma);
	}

	inline void quadratizefinalCost(const Matrix<BDIM>& b, Matrix<BDIM, BDIM>& Q, Matrix<BDIM>& q, double& qscalar)
	{
		double p = finalcost(b);

		//q here: jacobian.
		Matrix<BDIM> br(b), bl(b);
		for (int i = 0; i < BDIM; ++i) {
			br[i] += H; bl[i] -= H; 
			q[i] = (finalcost(br) - finalcost(bl)) / (2.0*H);
			Q(i,i) = ( finalcost(bl) - 2.0*p + finalcost(br) ) / (H*H);
			br[i] = bl[i] = b[i];
		}

		Matrix<BDIM> btr(b), btl(b), bbr(b), bbl(b);
		for(int i = 1; i < BDIM; ++i) {
			btr[i] += H; btl[i] -= H; bbr[i] += H; bbl[i] -= H;
			for (int j = 0; j < i; ++j) {
				btr[j] += H; btl[j] += H; bbr[j] -= H; bbl[j] -= H;
				Q(i,j) = Q(j,i) = (finalcost(bbl) + finalcost(btr) - finalcost(btl) - finalcost(bbr)) / (4.0*H*H);
				btr[j] = btl[j] = bbr[j] = bbl[j] = b[j];
			}
			btr[i] = btl[i] = bbr[i] = bbl[i] = b[i];
		}

		//jacobian - Hessian * bstar.
		q = q - Q * b;
		qscalar = p - tr(~(q + Q*b)*b) + 0.5 * tr(~b * Q * b);

		return;
	}


};



#endif