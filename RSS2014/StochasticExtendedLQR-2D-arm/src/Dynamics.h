#ifndef _DYNAMICS_
#define _DYNAMICS_

#define XDIM 4  //dimension of state theta1 theta2, angular velocity 1 and 2.
#define UDIM 2  //torques u1 and u2

const double percent = 0.2; //noise. 20% of the norm of the control input.
const double dt = 0.02;

#include "matrix.h"
#include "utils.h"
#include <vector>

class Dynamics{

public:
	
	Matrix<2,2> B;
	Matrix<2,1> m; //mass of two links
	Matrix<2,1> l; //length of two links.
	Matrix<2,1> s; //distance from joint center to the center of the mass for two links.
	Matrix<2,1> Inertial; //moment of interial for two links.
	double a1, a2, a3;
	double d;


	Dynamics(){

		//define constant terms.
		m(0,0) = 1.4;  m(1,0) = 1.0;
		l(0,0) = 0.3;  l(1,0) = 0.33;
		s(0,0) = 0.11; s(1,0) = 0.16;
		Inertial(0,0) = 0.025; Inertial(1,0) = 0.045;
		B(0,0) = 0.05; B(0,1) = 0.025; B(1,0) = 0.025; B(1,1) = 0.05;
		a1 = Inertial(0,0) + Inertial(1,0) + m(1,0)*l(0,0)*l(0,0);
		a2 = m(1,0)*l(0,0)*s(1,0);
		a3 = Inertial(1,0);
		
		//d = 0.005;
		//d = 0.0078125;
		d = 0.0009765625; //for finite difference.
	}

	//deterministic continuous dynamics. function: f.
	inline Matrix<XDIM> continuous_dynamics_determin(const Matrix<XDIM, 1>& x, const Matrix<UDIM,1>& control)
	{
		Matrix<2,2> M; Matrix<2,1> C;
		M(0,0) = a1 + 2 * a2 * cos(x(1,0)); M(0,1) = a3 + a2 * cos(x(1,0));
		M(1,0) = a3 + a2 * cos(x(1,0)); M(1,1) = a3;
		C(0,0) = -x(3,0) * (2 * x(2,0) + x(3,0)); C(1,0) = x(2,0) * x(2,0); C = C * a2 * sin(x(1,0));

		Matrix<4> dotx = zeros<4,1>();
		dotx(0,0) = x(2,0);  //w1;
		dotx(1,0) = x(3,0);  //w2;

		Matrix<2> dotw = zeros<2,1>();
		dotw = !M * (control - C - B * x.subMatrix<2,1>(2,0));

		dotx.insert<2,1>(2,0, dotw);

		return dotx;
	}

	//the Jacobian of f with respect to state x.
	inline Matrix<XDIM,XDIM> Jacobian_cont(const Matrix<XDIM,1>& x, const Matrix<UDIM,1>& control)
	{

		Matrix<XDIM,XDIM> Jacob = zeros<XDIM,XDIM>();
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

	//return dot Sigma: Sigma * Jacobian + Jacobian * Sigma + N N. 
	inline Matrix<XDIM, XDIM> continuous_noise(const Matrix<XDIM, 1>& x, const Matrix<UDIM, 1>& u, const Matrix<XDIM, XDIM>& sigma)
	{
		Matrix<2,2> M; M.reset();
		M(0,0) = a1 + 2 * a2 * cos(x(1,0)); M(0,1) = a3 + a2 * cos(x(1,0));
		M(1,0) = a3 + a2 * cos(x(1,0)); M(1,1) = a3;
		Matrix<XDIM, XDIM> dotSigma =  zeros<XDIM, XDIM>();
		Matrix<XDIM, XDIM> Jacob = Jacobian_cont(x, u);

		Matrix<XDIM, XDIM> N  = zeros<XDIM, XDIM>();
		double normu = sqrt(tr(~u * u));
		//Matrix<2,2> noise_torque = identity<2>() * percent * normu;
		//N.insert<2,2>(2,2, !M * noise_torque);

		///different N trial (the above did not make sense theoretically)
		N = identity<XDIM>() * percent * normu;
		
		dotSigma = Jacob * sigma + sigma * (~Jacob) + N * (~N);
		return dotSigma;
	}

	inline Matrix<XDIM,XDIM> continuous_noise_inverse(const Matrix<XDIM, 1>& x, const Matrix<UDIM, 1>& u, const Matrix<XDIM, XDIM>& sigma)
	{
		Matrix<2,2> M; M.reset();
		M(0,0) = a1 + 2 * a2 * cos(x(1,0)); M(0,1) = a3 + a2 * cos(x(1,0));
		M(1,0) = a3 + a2 * cos(x(1,0)); M(1,1) = a3;
		Matrix<XDIM, XDIM> dotSigma =  zeros<XDIM, XDIM>();
		Matrix<XDIM, XDIM> Jacob = Jacobian_cont(x, u);

		Matrix<XDIM, XDIM> N  = zeros<XDIM, XDIM>();
		double normu = sqrt(tr(~u * u));
		//Matrix<2,2> noise_torque = identity<2>() * percent * normu;
		//N.insert<2,2>(2,2, !M * noise_torque);
		N = identity<XDIM>() * percent * normu;
	
		dotSigma = -Jacob * sigma - sigma * ~Jacob + N * ~N;
		return dotSigma;
	}



	/***************************************************************************************/
	/********************the following are discrete dynamics implementation*****************/
	//discrete_dynamics (forward, computes x_n+1 and M_n+1 given x_n and u_n)
	inline void discrete_dynamics(const Matrix<XDIM, 1>& x, const Matrix<UDIM, 1>& u, Matrix<XDIM, 1>& g, Matrix<XDIM,XDIM>& M)
	{
		g.reset();
		M.reset();

		Matrix<XDIM,XDIM> A = identity<XDIM>();
		A.insert<2,2>(0,2, identity<2>() * dt);
		Matrix<XDIM, UDIM> B = zeros<XDIM,UDIM>();
		B.insert<2,2>(2,0, identity<2>()*dt);

		g = A * x + B * u;

		M = identity<XDIM>() * sqrt(0.1);

		/*g = zeros<XDIM, 1>();
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
		
		//std::cout<<EVal<<std::endl;

		//std::cout<<EVec<<std::endl;

		for(int i = 0; i < XDIM; i++){
			EVal(i,i) = sqrt(EVal(i,i));
		}
		M = EVec * EVal * ~EVec;  //M = U * EVal^0.5.*/

		return;
	}
	
	//discrete inverse dynamics
	//discrete inverse dynamics(backward, computes x_n and bar_M_n+1 given x_n+1 and u_n)
	inline void discrete_inverse_dynamics(const Matrix<XDIM, 1>& xn, const Matrix<UDIM, 1>&u, Matrix<XDIM, 1>& gbar, Matrix<XDIM, XDIM>& Mbar)
	{
		gbar.reset();
		Mbar.reset();

		Matrix<XDIM,XDIM> A = identity<XDIM>();
		A.insert<2,2>(0,2, identity<2>() * dt);
		Matrix<XDIM, UDIM> B = zeros<XDIM,UDIM>();
		B.insert<2,2>(2,0, identity<2>()*dt);

		Matrix<XDIM, XDIM> Ain = !A;
		gbar = Ain * xn - Ain * B * u;

		SymmetricMatrix<XDIM> MM = SymProd(Ain, identity<XDIM>() * 0.1 * ~Ain);
		Matrix<XDIM,XDIM> EVec;
		SymmetricMatrix<XDIM> EVal;
		jacobi(MM, EVec, EVal);

		for(int i = 0; i < XDIM; i++){
			EVal(i,i) = sqrt(EVal(i,i));
		}
		Mbar = EVec * EVal * ~EVec;  //M = U * EVal^0.5.*/


		/*gbar = zeros<XDIM, 1>(); 
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
		Mbar = EVec * EVal * ~EVec;  //M = U * EVal^0.5.*/

		return;
	}

	//linearization for discrete forward dynamics
	//given the pair of state and control (xstar, ustar), return At, Bt, ct and F, G,e
	inline void linearize_discrete_dynamics(const Matrix<XDIM,1>& xstar, const Matrix<UDIM, 1>& ustar, 
											Matrix<XDIM, XDIM>& At, Matrix<XDIM,UDIM>& Bt, Matrix<XDIM,1>& ct, 
											std::vector<Matrix<XDIM, XDIM>>& F, std::vector<Matrix<XDIM, UDIM>>& G, std::vector<Matrix<XDIM,1>>& e)
	{
		
		At = zeros<XDIM,XDIM>(); Bt = zeros<XDIM, UDIM>(); ct = zeros<XDIM,1>(); F.clear(); G.clear(); e.clear();
		F.resize(XDIM); G.resize(XDIM); e.resize(XDIM);

		Matrix<XDIM,1> originalxn; originalxn.reset();
		Matrix<XDIM, XDIM> originalM; originalM.reset();
		discrete_dynamics(xstar, ustar, originalxn, originalM);

		//linearize g to get At, Bt, and ct.
		Matrix<XDIM, XDIM> M = zeros<XDIM, XDIM>();
		for(int c = 0; c < XDIM; c++){
			Matrix<XDIM,1> augr = xstar; Matrix<XDIM, 1> augl = xstar;
			augr(c,0) += d; augl(c,0) -= d;
			Matrix<XDIM,1> augrxn; augrxn.reset();
			Matrix<XDIM,1> auglxn; auglxn.reset();
			discrete_dynamics(augr, ustar, augrxn, M);
			discrete_dynamics(augl, ustar, auglxn, M);
			At.insert<XDIM, 1>(0, c, (augrxn-auglxn)/(2.0*d));
		}
		for(int c = 0; c < UDIM; c++){
			Matrix<UDIM, 1> augru(ustar);
			Matrix<UDIM, 1> auglu(ustar);
			augru(c,0)+= d; auglu(c,0) -=d;
			Matrix<XDIM, 1> augrxn, auglxn; augrxn.reset(); auglxn.reset();
			discrete_dynamics(xstar, augru, augrxn, M);
			discrete_dynamics(xstar, auglu, auglxn, M);
			Bt.insert<XDIM,1>(0, c, (augrxn - auglxn)/(2.0*d));
		}
		ct = originalxn - At*xstar - Bt*ustar;

		//linearize M to get F, G and e.
		Matrix<XDIM,1> x; x.reset();
		for(int c = 0; c < XDIM; c++){
			Matrix<XDIM,1> augrx(xstar), auglx(xstar);
			augrx(c,0) += d;
			auglx(c,0) -= d;
			Matrix<XDIM,XDIM> augrM, auglM; augrM.reset(); auglM.reset(); 
			discrete_dynamics(augrx, ustar, x, augrM);
			discrete_dynamics(auglx, ustar, x, auglM);
			for(int i = 0; i < XDIM; i++){
				F[i].insert<XDIM,1>(0, c, (augrM.subMatrix<XDIM, 1>(0,i) - auglM.subMatrix<XDIM,1>(0,i))/(2.0*d));
			}
		}
		for(int c = 0; c < UDIM; c++){
			Matrix<UDIM,1> augru(ustar), auglu(ustar);
			augru(c,0) += d; auglu(c,0) -= d;
			Matrix<XDIM, XDIM> augrM, auglM; augrM.reset(); auglM.reset();
			discrete_dynamics(xstar, augru, x, augrM);
			discrete_dynamics(xstar, auglu, x, auglM);
			for(int i = 0; i < XDIM; i++){
				G[i].insert<XDIM, 1>(0, c, (augrM.subMatrix<XDIM,1>(0,i) - auglM.subMatrix<XDIM,1>(0,i))/(2.0*d));
			}
		}
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
		Matrix<XDIM, XDIM> M = zeros<XDIM, XDIM>();
		for(int c = 0; c < XDIM; c++){
			Matrix<XDIM, 1> augrxn(xnstar);
			Matrix<XDIM, 1> auglxn(xnstar);
			augrxn(c,0) += d; auglxn(c,0) -= d;
			Matrix<XDIM, 1> grbar, glbar; grbar.reset(); glbar.reset();
			discrete_inverse_dynamics(augrxn, ustar, grbar, M);
			discrete_inverse_dynamics(auglxn, ustar, glbar, M);
			Abart.insert<XDIM,1>(0,c, (grbar - glbar)/(2*d));
		}
		for(int c = 0; c < UDIM; c++){
			Matrix<UDIM, 1> augru(ustar);
			Matrix<UDIM, 1> auglu(ustar);
			augru(c, 0) += d; auglu(c,0) -= d;
			Matrix<XDIM, 1> grbar, glbar; grbar.reset(); glbar.reset();
			discrete_inverse_dynamics(xnstar, augru, grbar, M);
			discrete_inverse_dynamics(xnstar, auglu, glbar, M);
			Bbart.insert<XDIM,1>(0,c, (grbar - glbar)/(2.0*d));
		}
		cbart = originalx - Abart * xnstar - Bbart * ustar;

		//linearize M to get F, G and e
		Matrix<XDIM, 1> x; x.reset();
		for(int c = 0; c < XDIM; c++){
			Matrix<XDIM, 1> augrxn(xnstar);
			Matrix<XDIM, 1> auglxn(xnstar);
			augrxn(c,0) += d; auglxn(c,0) -= d;
			Matrix<XDIM, XDIM> augrMbar, auglMbar; augrMbar.reset(); auglMbar.reset();
			discrete_inverse_dynamics(augrxn, ustar, x, augrMbar);
			discrete_inverse_dynamics(auglxn, ustar, x, auglMbar);
			for(int i = 0; i < XDIM; i++){
				Fbar[i].insert<XDIM, 1>(0, c, (augrMbar.subMatrix<XDIM,1>(0,i) - auglMbar.subMatrix<XDIM,1>(0,i))/(2.0*d));
			}
		}
		for(int c = 0; c < UDIM; c++){
			Matrix<UDIM,1> augru(ustar);
			Matrix<UDIM, 1> auglu(ustar);
			augru(c,0) += d; auglu(c,0) -= d;
			Matrix<XDIM, XDIM> augrMbar, auglMbar; augrMbar.reset(); auglMbar.reset();
			discrete_inverse_dynamics(xnstar, augru, x, augrMbar);
			discrete_inverse_dynamics(xnstar, auglu, x, auglMbar);
			for(int i = 0; i < XDIM; i++){
				Gbar[i].insert<XDIM, 1>(0, c, (augrMbar.subMatrix<XDIM,1>(0,i) - auglMbar.subMatrix<XDIM,1>(0,i))/(2.0*d));
			}
		}
		for(int i = 0; i < XDIM; i++){
			ebar[i] = originalMbar.subMatrix<XDIM,1>(0,i) - Fbar[i]*xnstar - Gbar[i]*ustar;
		}

		return;
	}



	/*************************Forward Kinematics***********************************************/
	inline Matrix<2> FK(const Matrix<4,1>& state, const int& link, const double& length)
	{
		Matrix<2> pos = zeros<2,1>();
		double theta1 = state(0,0);
		double theta2 = state(1,0);

		if(link == 1) //first link
		{
			pos[0] = length * cos(theta1);
			pos[1] = length * sin(theta1);
		}
		else if(link == 2)
		{
			pos[0] = l(0,0) * cos(theta1) + length * cos(theta1 + theta2);
			pos[1] = l(0,0) * sin(theta1) + length * sin(theta1 + theta2);
		}
		
		return pos;
	}


	inline Matrix<2, 4> Jacobian_FK(const Matrix<4,1>& state, const int& link, const double& length)
	{
		Matrix<2,4> J = zeros<2,4>();
		if(link == 1){
			J(0,0) = -length * sin(state[0]);
			J(1,0) = length * cos(state[0]);
		}
		else{
			J(0,0) = -l(0,0) * sin(state[0]) - length * sin(state[0] + state[1]);  J(0,1) = -length * sin(state[0]+state[1]);
			J(1,0) = l(0,0) * cos(state[0]) + length * cos(state[0] + state[1]);   J(1,1) = length * cos(state[0] + state[1]);
		}

		return J;
	}



	/**********************************************************************************************/
	/**********************************************************************************************/
	/**********************************************************************************************/
	//implement of Euler method, in order to check if consistent with the RK4.
	inline Matrix<XDIM,1> Euler_integral_forward_noise(const Matrix<XDIM,1>& x, const Matrix<UDIM,1>&u)
	{
		Matrix<XDIM,1> xn = x;

		int set = 100;
		double tau = dt / set;
		for(int i = 0; i < set; i++){
			Matrix<2,2> M = zeros<2,2>();
			M(0,0) = a1 + 2 * a2 * cos(xn(1,0)); M(0,1) = a3 + a2 * cos(xn(1,0));
			M(1,0) = a3 + a2 * cos(xn(1,0)); M(1,1) = a3;
			Matrix<XDIM, XDIM> N  = zeros<XDIM, XDIM>();
			double normu = sqrt(tr(~u * u));
			//Matrix<2,2> noise_torque = identity<2>() * percent * normu;
			//N.insert<2,2>(2,2, !M * noise_torque);
			N = identity<XDIM>() * percent * normu;
		
			xn += tau * continuous_dynamics_determin(xn, u); 
			xn += sqrt(tau) * N * sampleGaussian(zeros<4,1>(), identity<4>());
		}
		return xn;
	}


	inline Matrix<XDIM,1> Euler_integral_inverse_noise(const Matrix<XDIM,1> xn, const Matrix<UDIM,1>& u)
	{
		Matrix<XDIM,1> x = xn;
		int seg = 100;
		double tau = dt / seg;

		for(int i = 0; i < seg;i++){
			Matrix<2,2> M = zeros<2,2>();
			M(0,0) = a1 + 2 * a2 * cos(xn(1,0)); M(0,1) = a3 + a2 * cos(xn(1,0));
			M(1,0) = a3 + a2 * cos(xn(1,0)); M(1,1) = a3;
			Matrix<XDIM, XDIM> N  = zeros<XDIM, XDIM>();
			double normu = sqrt(tr(~u * u));
			//Matrix<2,2> noise_torque = identity<2>() * percent * normu;
			//N.insert<2,2>(2,2, !M * noise_torque);
			N = identity<XDIM>() * percent * normu;
	
			x -= tau * continuous_dynamics_determin(x, u);
			x -= sqrt(tau) * N * sampleGaussian(zeros<4,1>(), identity<4>());
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