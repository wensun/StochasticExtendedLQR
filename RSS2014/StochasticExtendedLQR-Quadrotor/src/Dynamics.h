#ifndef _DYNAMICSS_
#define _DYNAMICSS_

#define XDIM 12  //dimension of state .
#define UDIM 4  //Diemons of control input
#define DIM 3

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
#include "Obstacles.h"

/////Robot physical parameters.

const double percent = 0.1;
const double gravity = 9.80665;
const double mass = 0.5;             // mass, kg  (source: paper)
const SymmetricMatrix<3> inertia = 0.05 * identity<3>(); // moment of inertia matrix 
const double momentConst = 1.5e-9 / 6.11e-8;  // ratio between force and moment of rotor
const double dragConst = 0.15; //0.15;
const double length = 0.3429/2;          // distance between center and rotor, m
const SymmetricMatrix<3> invInertia  = !inertia;
Matrix<3> eX, eY, eZ;


class Dynamics{

public:
	double dt;
	double d; //d is used for finite difference computing.
	Matrix<U_DIM> uNominal;

	Dynamics(const double& Deltat){

		dt = Deltat;
		d = 0.0009765625;

		eX[0] = 1; eX[1] = 0; eX[2] = 0;
		eY[0] = 0; eY[1] = 1; eY[2] = 0;
		eZ[0] = 0; eZ[1] = 0; eZ[2] = 1;

		uNominal[0] = uNominal[1] = uNominal[2] = uNominal[3] = gravity*mass/4;

	}

	//deterministic continuous dynamics. function: f.
	inline Matrix<XDIM> continuous_dynamics_determin(const Matrix<XDIM, 1>& x, const Matrix<UDIM,1>& u)
	{
		Matrix<XDIM> xDot;

		Matrix<3> p = x.subMatrix<3,1>(0,0);
		Matrix<3> v = x.subMatrix<3,1>(3,0);
		Matrix<3> r = x.subMatrix<3,1>(6,0);
		Matrix<3> w = x.subMatrix<3,1>(9,0);

		// \dot{p} = v
		xDot.insert(0, 0, v);

		// \dot{v} = [0,0,-g]^T + R*exp([r])*[0,0,(f_1 + f_2 + f_3 + f_4) / m]^T; 
		xDot.insert(3, 0, -gravity*eZ + exp(skewSymmetric(r))*eZ*((u[0]+u[1]+u[2]+u[3])/mass) - v*(dragConst/mass)); 

		// \dot{r} = w + 0.5*skewSymmetric(r)*w + (1.0/tr(~r*r))*(1.0 - 0.5*sqrt(tr(~r*r))/tan(0.5*sqrt(tr(~r*r))))*skewSymmetric(r)*(skewSymmetric(r)*w)
		double l = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
		if (0.5*l > 0.0) {
			xDot.insert(6, 0, w + 0.5*skewSymmetric(r)*w + (1.0 - 0.5*l/tan(0.5*l))*skewSymmetric(r / l)*(skewSymmetric(r / l)*w));
		} else {
			xDot.insert(6, 0, w);
		}

		// \dot{w} = J^{-1}*([l*(f_2 - f_4), l*(f_3 - f_1), (f_1 - f_2 + f_3 - f_4)*k_M]^T - [w]*J*w)
		xDot.insert(9, 0, invInertia*( length*(u[1] - u[3])*eX + length*(u[2] - u[0])*eY + (u[0] - u[1] + u[2] - u[3])*momentConst*eZ - skewSymmetric(w)*inertia*w));		

		return xDot;
	
	}

	//the Jacobian of f with respect to state x.
	inline Matrix<XDIM,XDIM> Jacobian_cont(const Matrix<XDIM,1>& x, const Matrix<UDIM,1>& control)
	{
		
		Matrix<XDIM,XDIM> Jacob(zeros<XDIM, XDIM>());
		Matrix<XDIM,1> oldvalue = continuous_dynamics_determin(x, control);

		Matrix<XDIM, 1> augr, augl; augr.reset(); augl.reset();
		for(int c = 0; c < XDIM; c++){
			augr.reset(); augr = x;
			augl.reset(); augl = x;
			augr(c,0) += d;
			augl(c,0) -= d;
			Matrix<XDIM, 1> tmpvaluer = continuous_dynamics_determin(augr, control);
			Matrix<XDIM, 1> tmpvaluel = continuous_dynamics_determin(augl, control);
			
			Jacob.insert<XDIM, 1>(0, c, (tmpvaluer - tmpvaluel)/(2.0*d));
		}
		return Jacob;
	}

	//return dot Sigma.
	inline Matrix<XDIM, XDIM> continuous_noise(const Matrix<XDIM, 1>& x, const Matrix<UDIM, 1>& u, const Matrix<XDIM, XDIM>& sigma)
	{
		Matrix<XDIM, XDIM> dotSigma =  zeros<XDIM, XDIM>();
		Matrix<XDIM, XDIM> Jacob = Jacobian_cont(x, u);

		double normu = sqrt(tr(~u * u));
		Matrix<XDIM,XDIM> N = identity<XDIM>() * percent * normu;   //N: 0.2||u|| * I

		//Matrix<XDIM, 1> N(zeros<XDIM,1>());
		//N.insert<3,1>(3,0, 4.0 / mass * percent * normu * exp(skewSymmetric(r)) * eZ);
		//N.insert<3,1>(9,0, invInertia*(2*length*percent*normu * eX + 2*length * percent*normu*eY + 4*momentConst*eZ*percent*normu));


		dotSigma = Jacob * sigma + sigma * (~Jacob) + N * (~N);
		return dotSigma;
	}

	inline Matrix<XDIM,XDIM> continuous_noise_inverse(const Matrix<XDIM, 1>& x, const Matrix<UDIM, 1>& u, const Matrix<XDIM, XDIM>& sigma)
	{
		Matrix<XDIM, XDIM> dotSigma =  zeros<XDIM, XDIM>();
		Matrix<XDIM, XDIM> Jacob = Jacobian_cont(x, u);

		double normu = sqrt(tr(~u * u));
		Matrix<XDIM,XDIM> N = identity<XDIM>() * percent * normu;
	
		dotSigma = -Jacob * sigma - sigma * ~Jacob + N * ~N;
		return dotSigma;
	}



	/***************************************************************************************/
	/********************the following are discrete dynamics implementation*****************/
	//discrete_dynamics (forward, computes x_n+1 and M_n+1 given x_n and u_n)
	inline void discrete_dynamics(const Matrix<XDIM, 1>& x, const Matrix<UDIM, 1>& u, Matrix<XDIM, 1>& g, Matrix<XDIM,XDIM>& M)
	{
		g = zeros<XDIM, 1>();
		M = zeros<XDIM, XDIM>();

		Matrix<XDIM,XDIM> sigma = zeros<XDIM,XDIM>();

		//RK4 for x.
		Matrix<XDIM, 1> k1 = continuous_dynamics_determin(x, u);
		Matrix<XDIM, 1> k2 = continuous_dynamics_determin(x + 0.5 * dt * k1, u);
		Matrix<XDIM, 1> k3 = continuous_dynamics_determin(x + 0.5 * dt * k2, u);
		Matrix<XDIM, 1> k4 = continuous_dynamics_determin(x + dt * k3, u);

		g = x + 1.0/6.0 * dt * (k1 + 2*k2 + 2*k3 + k4);

		//RK4 for sigma
		Matrix<XDIM, XDIM> K1 = continuous_noise(x, u, sigma);
		Matrix<XDIM, XDIM> K2 = continuous_noise(x + 0.5 * dt * k1, u, sigma + 0.5 * dt * K1);
		Matrix<XDIM, XDIM> K3 = continuous_noise(x + 0.5 * dt * k2, u, sigma + 0.5 * dt * K2);
		Matrix<XDIM, XDIM> K4 = continuous_noise(x + dt * k3, u, sigma + dt * K3);

		Matrix<XDIM, XDIM> Sigma = sigma + 1.0/6.0 * dt * (K1 + 2*K2 + 2*K3 + K4);
		Matrix<XDIM,XDIM> EVec = zeros<XDIM,XDIM>();
		SymmetricMatrix<XDIM> EVal = zeros<XDIM>();
		SymmetricMatrix<XDIM> Sigmap = SymProd(Sigma, identity<XDIM>());
		jacobi(Sigmap, EVec, EVal);
		
		for(int i = 0; i < XDIM; i++){
			if(EVal(i,i) < 0)
				EVal(i,i) = 0.0;

			EVal(i,i) = sqrt(EVal(i,i));
		}
		M = EVec * EVal * ~EVec;  //M = U * EVal^0.5.

		return;
	}
	
	//discrete inverse dynamics
	//discrete inverse dynamics(backward, computes x_n and bar_M_n+1 given x_n+1 and u_n)
	inline void discrete_inverse_dynamics(const Matrix<XDIM, 1>& xn, const Matrix<UDIM, 1>&u, Matrix<XDIM, 1>& gbar, Matrix<XDIM, XDIM>& Mbar)
	{
		gbar = zeros<XDIM, 1>(); 
		Mbar = zeros<XDIM, XDIM>();

		Matrix<XDIM, XDIM> sigman = zeros<XDIM, XDIM>();
		
		//inverse RK4 for x
		Matrix<XDIM, 1> k1 = continuous_dynamics_determin(xn, u);
		Matrix<XDIM, 1> k2 = continuous_dynamics_determin(xn - 0.5 * dt * k1, u);
		Matrix<XDIM, 1> k3 = continuous_dynamics_determin(xn - 0.5 * dt * k2, u);
		Matrix<XDIM, 1> k4 = continuous_dynamics_determin(xn - dt * k3, u);

		gbar = xn - 1.0/6.0 * dt * (k1 + 2*k2 + 2*k3 + k4);

		//inverse RK4 for sigma;
		Matrix<XDIM, XDIM> K1 = continuous_noise_inverse(xn, u, sigman);
		Matrix<XDIM, XDIM> K2 = continuous_noise_inverse(xn - 0.5 * dt * k1, u, sigman + 0.5 * dt * K1);
		Matrix<XDIM, XDIM> K3 = continuous_noise_inverse(xn - 0.5 * dt * k2, u, sigman + 0.5 * dt * K2);
		Matrix<XDIM, XDIM> K4 = continuous_noise_inverse(xn - dt * k3, u, sigman + dt * K3);

		Matrix<XDIM, XDIM> Sigma = sigman + 1.0/6.0 * dt * (K1 + 2*K2 + 2*K3 + K4);
		Matrix<XDIM,XDIM> EVec = zeros<XDIM,XDIM>();
		SymmetricMatrix<XDIM> EVal = zeros<XDIM>();
		SymmetricMatrix<XDIM> Sigmap = SymProd(Sigma, identity<XDIM>());
		jacobi(Sigmap, EVec, EVal);

		for(int i = 0; i < XDIM; i++){
			EVal(i,i) = sqrt(EVal(i,i));
		}
		Mbar = EVec * EVal * ~EVec;  //M = U * EVal^0.5.

		return;
	}

	//linearization for discrete forward dynamics
	//given the pair of state and control (xstar, ustar), return At, Bt, ct and F, G,e
	inline void linearize_discrete_dynamics(const Matrix<XDIM,1>& xstar, const Matrix<UDIM, 1>& ustar, 
											Matrix<XDIM, XDIM>& At, Matrix<XDIM,UDIM>& Bt, Matrix<XDIM,1>& ct, 
											std::vector<Matrix<XDIM, XDIM>>& F, std::vector<Matrix<XDIM, UDIM>>& G, std::vector<Matrix<XDIM,1>>& e)
	{
		//double d = 0.0009765625;
		At = zeros<XDIM,XDIM>(); Bt = zeros<XDIM, UDIM>(); ct = zeros<XDIM,1>(); F.clear(); G.clear(); e.clear();
		F.resize(XDIM); G.resize(XDIM); e.resize(XDIM);

		Matrix<XDIM,1> originalxn; originalxn.reset();
		Matrix<XDIM, XDIM> originalM; originalM.reset();
		discrete_dynamics(xstar, ustar, originalxn, originalM);

		//linearize g to get At, Bt, and ct.
		Matrix<XDIM,1> augrxn; augrxn.reset();
		Matrix<XDIM,1> auglxn; auglxn.reset();
		Matrix<XDIM,XDIM> augrM, auglM; augrM.reset(); auglM.reset(); 
		for(int c = 0; c < XDIM; c++){
			Matrix<XDIM,1> augr = xstar; Matrix<XDIM, 1> augl = xstar;
			augr(c,0) += d; augl(c,0) -= d;
			discrete_dynamics(augr, ustar, augrxn, augrM);
			discrete_dynamics(augl, ustar, auglxn, auglM);
			At.insert<XDIM, 1>(0, c, (augrxn-auglxn)/(2.0*d));
			for(int i =0; i < XDIM; i++){
				F[i].insert<XDIM,1>(0, c, (augrM.subMatrix<XDIM, 1>(0,i) - auglM.subMatrix<XDIM,1>(0,i))/(2.0*d));
			}
		}
		for(int c = 0; c < UDIM; c++){
			Matrix<UDIM, 1> augru(ustar);
			Matrix<UDIM, 1> auglu(ustar);
			augru(c,0)+= d; auglu(c,0) -=d;
			discrete_dynamics(xstar, augru, augrxn, augrM);
			discrete_dynamics(xstar, auglu, auglxn, auglM);
			Bt.insert<XDIM,1>(0, c, (augrxn - auglxn)/(2.0*d));
			for(int i = 0; i < XDIM; i++){
				G[i].insert<XDIM, 1>(0, c, (augrM.subMatrix<XDIM,1>(0,i) - auglM.subMatrix<XDIM,1>(0,i))/(2.0*d));
			}
		}
		ct = originalxn - At*xstar - Bt*ustar;
		for(int i = 0; i < XDIM; i++){
			e[i] = originalM.subMatrix<XDIM,1>(0,i) - F[i]*xstar - G[i]*ustar;
		}

		return;
	}

	//linearization for discrete backward dynamics
	inline void linearize_discrete_inverse_dynamics(const Matrix<XDIM,1>& xnstar, const Matrix<UDIM, 1>& ustar, 
											Matrix<XDIM, XDIM>& Abart, Matrix<XDIM,UDIM>& Bbart, Matrix<XDIM,1>& cbart, 
											std::vector<Matrix<XDIM, XDIM>>& Fbar, std::vector<Matrix<XDIM, UDIM>>& Gbar, std::vector<Matrix<XDIM,1>>& ebar)
	{
		Abart = zeros<XDIM,XDIM>(); Bbart = zeros<XDIM, UDIM>(); cbart = zeros<XDIM,1>(); Fbar.clear(); Gbar.clear(); ebar.clear();
		Fbar.resize(XDIM, zeros<XDIM,XDIM>()); Gbar.resize(XDIM, zeros<XDIM, UDIM>()); ebar.resize(XDIM, zeros<XDIM,1>());

		Matrix<XDIM, 1> originalx; originalx.reset();
		Matrix<XDIM, XDIM> originalMbar; originalMbar.reset();
		discrete_inverse_dynamics(xnstar, ustar, originalx, originalMbar);

		//linearize gbar to get At, Bt, and ct.
		Matrix<XDIM, 1> grbar, glbar; grbar.reset(); glbar.reset();
		Matrix<XDIM, XDIM> augrMbar, auglMbar; augrMbar.reset(); auglMbar.reset();
		for(int c = 0; c < XDIM; c++){
			Matrix<XDIM, 1> augrxn(xnstar);
			Matrix<XDIM, 1> auglxn(xnstar);
			augrxn(c,0) += d; auglxn(c,0) -= d;
			discrete_inverse_dynamics(augrxn, ustar, grbar, augrMbar);
			discrete_inverse_dynamics(auglxn, ustar, glbar, auglMbar);
			Abart.insert<XDIM,1>(0,c, (grbar - glbar)/(2*d));
			for(int i = 0; i < XDIM; i++){
				Fbar[i].insert<XDIM, 1>(0, c, (augrMbar.subMatrix<XDIM,1>(0,i) - auglMbar.subMatrix<XDIM,1>(0,i))/(2.0*d));
			}
		}
		for(int c = 0; c < UDIM; c++){
			Matrix<UDIM, 1> augru(ustar);
			Matrix<UDIM, 1> auglu(ustar);
			augru(c, 0) += d; auglu(c,0) -= d;
			discrete_inverse_dynamics(xnstar, augru, grbar, augrMbar);
			discrete_inverse_dynamics(xnstar, auglu, glbar, auglMbar);
			Bbart.insert<XDIM,1>(0,c, (grbar - glbar)/(2.0*d));
			for(int i = 0; i < XDIM; i++){
				Gbar[i].insert<XDIM,1>(0, c, (augrMbar.subMatrix<XDIM,1>(0,i) - auglMbar.subMatrix<XDIM,1>(0,i))/(2.0*d));
			}
		}
		cbart = originalx - Abart * xnstar - Bbart * ustar;

		for(int i = 0; i < XDIM; i++){
			ebar[i] = originalMbar.subMatrix<XDIM,1>(0,i) - Fbar[i]*xnstar - Gbar[i]*ustar;
		}

		return;
	}


	/**********************************************************************************************/
	/**********************************************************************************************/
	/**********************************************************************************************/
	//implement of Euler method, in order to check if consistent with the RK4
	//checked with the RK4, works well.
	inline Matrix<XDIM,1> Euler_integral_forward_noise(const Matrix<XDIM,1>& x, const Matrix<UDIM,1>&u)
	{
		Matrix<XDIM,1> xn = x;

		int set = 100;
		double tau = dt / set;
		for(int i = 0; i < set; i++){
			Matrix<XDIM, XDIM> N  = zeros<XDIM, XDIM>();
			double normu = sqrt(tr(~u * u));
			N = identity<XDIM>() * percent * normu;
			xn += tau * continuous_dynamics_determin(xn, u); 
			xn += sqrt(tau) * N * sampleGaussian(zeros<XDIM,1>(), identity<XDIM>());
		}
		return xn;
	}


	inline Matrix<XDIM,1> Euler_integral_inverse_noise(const Matrix<XDIM,1> xn, const Matrix<UDIM,1>& u)
	{
		Matrix<XDIM,1> x = xn;
		int seg = 100;
		double tau = dt / seg;

		for(int i = 0; i < seg;i++){
			Matrix<XDIM, XDIM> N  = zeros<XDIM, XDIM>();
			double normu = sqrt(tr(~u * u));
			N = identity<XDIM>() * percent * normu;
			x -= tau * continuous_dynamics_determin(x, u);
			x -= sqrt(tau) * N * sampleGaussian(zeros<XDIM,1>(), identity<XDIM>());
		}
		return x;
	}

	inline Matrix<XDIM,1> Euler_integral_forward_determin(const Matrix<XDIM,1>& x, const Matrix<UDIM,1>&u)
	{
		Matrix<XDIM,1> xn = x;

		int set = 10000;
		double tau = dt / set;
		for(int i = 0; i < set; i++)
			xn += tau * continuous_dynamics_determin(xn, u); 

		return xn;
	}

	inline Matrix<XDIM,1> Euler_integral_inverse_determin(const Matrix<XDIM,1> xn, const Matrix<UDIM,1>& u)
	{
		Matrix<XDIM,1> x = xn;
		int seg = 10000;
		double tau = dt / seg;
		for(int i = 0; i < seg;i++)
			x -= tau * continuous_dynamics_determin(x, u);

		return x;
	}




};



#endif