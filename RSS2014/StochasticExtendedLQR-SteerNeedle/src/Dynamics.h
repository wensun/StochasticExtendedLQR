#ifndef _DYNAMICS_
#define _DYNAMICS_

#define XDIM 6  //p and r
#define UDIM 3  // v w, k
#define BDIM 27
#define ZDIM 3
const double percent = 0.2; //noise. 20% of the norm of the control input.


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


class Dynamics{

public:

	double d;
	double dt;

	Matrix<UDIM, UDIM> MotionNoise;
	Matrix<3> sensorP;
	

	Dynamics(const double& tau){

		dt = tau;
		//d = 0.005;
		//d = 0.0078125;
		d = 0.0009765625; //for finite difference.
		sensorP[0] = 0.0; sensorP[1] = 4.0;  sensorP[2] = 0.0;  //at the top middle.
		MotionNoise = identity<UDIM>() * 0.05;
	}

	//vector version.
	inline Matrix<6> seDynamics(const Matrix<6>& x, const Matrix<3>& u)
	{
		Matrix<4,4> U;
		Matrix<3,1> v, w;
		w[0] = u[0]*u[2];
		w[1] = 0.0;
		w[2] = u[1];
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = u[0];
		U = zeros<4,4>();
		U.insert(0,0, cpMatrix(w));
		U.insert(0,3, v);

		Matrix<4,4> SEX = zeros<4,4>();
		SEX.insert<3,3>(0,0, exp(cpMatrix(x.subMatrix<3,1>(3,0))));  //exp([r]);
		SEX.insert<3,1>(0,3, x.subMatrix<3,1>(0,0));
		SEX(3,3) = 1.0;

		Matrix<4,4> SEXN = SEX * exp(dt * U);
		Matrix<6> xnext(zeros<6,1>());

		xnext.insert<3,1>(0,0, SEXN.subMatrix<3,1>(0,3));
		xnext.insert<3,1>(3,0, down(logm(SEXN.subMatrix<3,3>(0,0))));  //[r] = log(R).

		return xnext;
	}

	inline void Linearize_A_B_V(const Matrix<6>& x, const Matrix<3>& u,
								Matrix<6,6>& At, Matrix<6,3>& Bt, Matrix<6,3>& Vt)
	{
		At.reset(); Bt.reset(); Vt.reset();
		Matrix<6> oldValue = seDynamics(x, u);
		Matrix<6, 1> augr, augl; augr.reset(); augl.reset();
		for(int c = 0; c < 6; c++){
			augr.reset(); augr = x;
			augl.reset(); augl = x;
			augr(c,0) += d;
			augl(c,0) -= d;
			Matrix<6, 1> tmpvaluer = seDynamics(augr, u);
			Matrix<6, 1> tmpvaluel = seDynamics(augl, u);
			At.insert<6, 1>(0, c, (tmpvaluer - tmpvaluel)/(2.0*d));
		}

		Matrix<3,1> ur, ul; ur.reset(); ul.reset();
		for(int c = 0; c < 3; c++){
			ur = u; ul = u;
			ur(c,0) += d;
			ul(c,0) -= d;
			Matrix<6,1> tmpr = seDynamics(x, ur);
			Matrix<6,1> tmpl = seDynamics(x, ul);
			Bt.insert<6,1>(0,c, (tmpr - tmpl)/(2.0*d));
		}
		Vt = Bt;
	}

	inline Matrix<ZDIM> obs(const Matrix<XDIM,1>& x)
	{
	
		Matrix<3> p = x.subMatrix<3,1>(0,0);
		
		double strenght = sqrt(tr(~(p-sensorP)*(p-sensorP))) * 0.05; //proportial to the distance to sensors.

		return p + sampleGaussian(zeros<3,1>(), identity<3>() * strenght);
	}

	inline void Linear_H_N(const Matrix<XDIM,1>& x,
							Matrix<ZDIM, XDIM>& H, Matrix<ZDIM, ZDIM>& N)
	{
		Matrix<3> p = x.subMatrix<3,1>(0,0);
		
		double strenght = sqrt(tr(~(p-sensorP)*(p-sensorP))) * 0.05; //proportial to the distance to sensors.

		H.insert<ZDIM,ZDIM>(0,0, identity<ZDIM>());
		N = identity<3>() * strenght;
		
	}

	//stochastic belie dynamics.
	inline void BeliefDynamics(const Matrix<BDIM>& belief, const Matrix<UDIM>& u, 
								Matrix<BDIM>& nbelief, Matrix<BDIM, XDIM>& M)
	{
		Matrix<XDIM, XDIM> At; Matrix<XDIM, UDIM> Vt; Matrix<XDIM, UDIM> Bt;
		Matrix<XDIM> x; Matrix<XDIM, XDIM> sqrtS;
		unVec(belief, x, sqrtS);
		Linearize_A_B_V(x, u, At, Bt, Vt);

		Matrix<XDIM, XDIM> Tau = At * (sqrtS * sqrtS) * ~At + Vt * MotionNoise * ~Vt;
		Matrix<ZDIM,XDIM> H; Matrix<ZDIM, ZDIM> N;
		Linear_H_N(x, H, N);
		Matrix<XDIM, ZDIM> K = Tau * ~H * !(H * Tau * ~H + N);
		Matrix<XDIM, XDIM> Sigma = Tau - K * H * Tau;
		Matrix<XDIM, XDIM> EVec;
		SymmetricMatrix<XDIM> EVal;
		SymmetricMatrix<XDIM> SymSigma = SymProd(Sigma, identity<XDIM>());
		jacobi(SymSigma, EVec, EVal);
		for(int i = 0; i < XDIM; i++){
			if(EVal(i,i) > 0)
				EVal(i,i) = sqrt(EVal(i,i));
			else
				EVal(i,i) = 0.0;
		}
		Matrix<XDIM, XDIM> SqrtSigma = EVec * EVal * ~EVec;
		Matrix<XDIM> xnext = seDynamics(x, u);
		vec(xnext, SqrtSigma, nbelief);

		M = zeros<BDIM, XDIM>();
		Matrix<XDIM, XDIM> tmpS = K * H * Tau;
		SymmetricMatrix<XDIM> SymtmpS = SymProd(tmpS, identity<XDIM>());
		jacobi(SymtmpS, EVec, EVal);
		for(int i = 0; i < XDIM; i++){
			if(EVal(i,i) > 0)
				EVal(i,i) = sqrt(EVal(i,i));
			else
				EVal(i,i) = 0.0;
		}
		Matrix<XDIM,XDIM> sqrttmpS = EVec * EVal * ~EVec;
		M.insert<XDIM, XDIM>(0,0, sqrttmpS);
	}


	inline Matrix<4,4> SEDynamics(const Matrix<4,4>& T, const Matrix<3>& u)
	{
		Matrix<4,4> U;
		Matrix<3,1> v, w;
		w[0] = u[0]*u[2];
		w[1] = 0.0;
		w[2] = u[1];
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = u[0];
		U = zeros<4,4>();
		U.insert(0,0, cpMatrix(w));
		U.insert(0,3, v);

		return T*exp(dt*U);
	}



	//linearization for discrete forward dynamics
	//given the pair of state and control (xstar, ustar), return At, Bt, ct and F, G,e
	inline void linearize_discrete_belief_dynamics(const Matrix<BDIM,1>& bstar, const Matrix<UDIM, 1>& ustar, 
											Matrix<BDIM, BDIM>& At, Matrix<BDIM,UDIM>& Bt, Matrix<BDIM,1>& ct, 
											std::vector<Matrix<BDIM, BDIM>>& F, std::vector<Matrix<BDIM, UDIM>>& G, std::vector<Matrix<BDIM,1>>& e)
	{
		At = zeros<BDIM,BDIM>(); Bt = zeros<BDIM, UDIM>(); ct = zeros<BDIM,1>(); F.clear(); G.clear(); e.clear();
		F.resize(XDIM); G.resize(XDIM); e.resize(XDIM);

		Matrix<BDIM,1> originalbn; originalbn.reset();
		Matrix<BDIM, XDIM> originalM; originalM.reset();
		BeliefDynamics(bstar, ustar, originalbn, originalM);

		//linearize g to get At, Bt, and ct.
		Matrix<BDIM, XDIM> M = zeros<BDIM, XDIM>();
		for(int c = 0; c < BDIM; c++){
			Matrix<BDIM,1> augr = bstar; Matrix<BDIM, 1> augl = bstar;
			augr(c,0) += d; augl(c,0) -= d;
			Matrix<BDIM,1> augrxn; augrxn.reset();
			Matrix<BDIM,1> auglxn; auglxn.reset();
			Matrix<BDIM,XDIM> augrM, auglM; augrM.reset(); auglM.reset(); 
			BeliefDynamics(augr, ustar, augrxn, augrM);
			BeliefDynamics(augl, ustar, auglxn, auglM);
		
			At.insert<BDIM, 1>(0, c, (augrxn-auglxn)/(2.0*d));
			for(int i = 0; i < XDIM; i++){
				F[i].insert<BDIM,1>(0, c, (augrM.subMatrix<BDIM, 1>(0,i) - auglM.subMatrix<BDIM,1>(0,i))/(2.0*d));
			}
		}
		for(int c = 0; c < UDIM; c++){
			Matrix<UDIM, 1> augru(ustar);
			Matrix<UDIM, 1> auglu(ustar);
			augru(c,0)+= d; auglu(c,0) -=d;
			Matrix<BDIM, 1> augrxn, auglxn; augrxn.reset(); auglxn.reset();
			Matrix<BDIM, XDIM> augrM, auglM; augrM.reset(); auglM.reset();
			BeliefDynamics(bstar, augru, augrxn, augrM);
			BeliefDynamics(bstar, auglu, auglxn, auglM);
			
			Bt.insert<BDIM,1>(0, c, (augrxn - auglxn)/(2.0*d));
			for(int i = 0; i < XDIM; i++){
				G[i].insert<BDIM, 1>(0, c, (augrM.subMatrix<BDIM,1>(0,i) - auglM.subMatrix<BDIM,1>(0,i))/(2.0*d));
			}
		}

		ct = originalbn - At*bstar - Bt*ustar;
		for(int i = 0; i < XDIM; i++){
			e[i] = originalM.subMatrix<BDIM,1>(0,i) - F[i]*bstar - G[i]*ustar;
		}

		return;
	}

};



#endif