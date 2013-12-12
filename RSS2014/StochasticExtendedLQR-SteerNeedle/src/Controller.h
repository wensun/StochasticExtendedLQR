#ifndef _CONTROLLER_
#define _CONTROLLER_

#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>

#include "RRTNode.h"
#include "Dynamics.h"
#include "callisto.h"
#include "matrix.h"
#include "utils.h"

class Controller
{
public:

	int path_size;
	double dt;
	double car_l;
	double goal_radius;
	std::vector<RRTNode> path;
	std::vector<Matrix<X_DIM, X_DIM>> A;
	std::vector<Matrix<X_DIM, U_DIM>> B;
	std::vector<Matrix<X_DIM, U_DIM>> V;
	std::vector<Matrix<X_DIM, Z_DIM>> K;
	std::vector<Matrix<U_DIM, X_DIM>> L;
	Matrix<Z_DIM, X_DIM> H;
	Matrix<Z_DIM, Z_DIM> W;
	Matrix<X_DIM,X_DIM> C; //penalty for state deviation
	Matrix<U_DIM,U_DIM> D; //penalty for control deviation
	Matrix<U_DIM, U_DIM> M; //process noise;
	Matrix<Z_DIM, Z_DIM> N; //sensing noise;


	Controller();

	Controller(const std::vector<RRTNode>& Path, const double& d_t, const double& CAR_L, const double& g_r)
	{
		path.clear();
		path = Path;
		dt = d_t;
		car_l = CAR_L;
		goal_radius = g_r;
		path_size = path.size();
		A.resize(path_size);
		B.resize(path_size);
		V.resize(path_size);
		K.resize(path_size);
		L.resize(path_size);

		H(0,0) = 1; H(0,1) = 0; H(0,2) = 0; H(0,3) = 0;
		H(1,0) = 0; H(1,1) = 1; H(1,2) = 0; H(1,3) = 0;
		H(2,0) = 0; H(2,1) = 0; H(2,2) = 0; H(2,3) = 1;

		W = identity<Z_DIM>();
		M = identity<U_DIM>() * 0.08*0.08;
		N = identity<Z_DIM>() * 0.001*0.001;
		C = identity<X_DIM>();
		D = identity<U_DIM>();
	}

	void compute_A_B_V()
	{
		for(int i = 1; i < path_size; i++){
			RRTNode tmp = path[i-1];
			Matrix<U_DIM> u = tmp.u;
			Matrix<X_DIM> X = tmp.x;
			A[i](0,0) = 1; A[i](0,1) = 0; A[i](0,2) = -dt*X[3]*sin(X[2]); A[i](0,3) = dt*cos(X[2]);
			A[i](1,0) = 0; A[i](1,1) = 1; A[i](1,2) = dt*X[3]*cos(X[2]); A[i](1,3) = dt*sin(X[2]);
			A[i](2,0) = 0; A[i](2,1) = 0; A[i](2,2) = 1; A[i](2,3) = dt*tan(u[1]) / car_l;
			A[i](3,0) = 0; A[i](3,1) = 0; A[i](3,2) = 0; A[i](3,3) = 1;

			B[i](0,0) = 0; B[i](0,1) = 0;
			B[i](1,0) = 0; B[i](1,1) = 0;
			B[i](2,0) = 0; B[i](2,1) = dt*X[3]*(1 + tan(u[1])*tan(u[1])) / car_l;
			B[i](3,0) = dt; B[i](3,1) = 0;

			V[i] = B[i];
		}
	}


	//kalman Filter
	void compute_K(const Matrix<X_DIM, X_DIM>& P0) //given initial covariance of the robot's initial distribution.
	{
		Matrix<X_DIM, X_DIM> P;
		P = P0;
		K[0] = zeros<X_DIM, Z_DIM>();
		for(int i = 1; i < path_size; i++){
			P = A[i]*P*~A[i] + V[i]*M*~V[i];
			K[i] = P*~H*!(H*P*~H + W*N*~W);
			P = (identity<X_DIM>() - K[i]*H)*P;
		}
	}
	void compute_L()
	{
		Matrix<X_DIM, X_DIM> S;
		S = C;
		L[path_size - 1] = zeros<U_DIM, X_DIM>();
		for(int k = path_size - 2; k >= 0; k--){
			L[k] = - !(~B[k+1]*S*B[k+1] + D) * ~B[k+1]*S*A[k+1];
			S = C + ~A[k+1] * S * A[k+1] + ~A[k+1]*S*B[k+1]*L[k];
		}
	}

	bool execute(int& cal_rrt)
	{
		Matrix<X_DIM> x,x_old;
		Matrix<U_DIM> u;
		x = path[0].x;
		
		for(int k = 1; k < path.size(); k++){
			u = path[k-1].u;
			x_old = x;
			Dynamics dy(dt/5, car_l);
			for(int m = 1; m <= 5; m++){
				x = dy.f(x_old, u, zeros<U_DIM,1>());

				int np[1] = {2};
				float p[6] = {x_old[0], x_old[1], 0.0, x[0], x[1], 0.0};
				CAL_CreatePolyline(cal_rrt, 1, np, p);
				x_old = x;
			}
		}
		return true;

	}

	bool simulateLQG(const Matrix<X_DIM>& start, const Matrix<X_DIM, X_DIM>& P0, const int& cal_rrt, const int& cal_environment, const Matrix<DIM>& goal, double& dis)
	{
	
		compute_A_B_V();
		compute_L();
		compute_K(P0);

		Matrix<U_DIM> u, u_d; //u_d is deviation u - u*.
		Matrix<X_DIM, X_DIM> P;
		Matrix<U_DIM> m; // gaussian noise for control input.
		Matrix<Z_DIM> n;
		Matrix<X_DIM> x_true, x_est, x_true_old, x_d_est;
		Matrix<Z_DIM> z_d; //observation z - h(x*,0)
		P = P0;
		//at the starting position, x_est is the mean, P0 is the covairance,x_true is the simulated true position
		x_est = start;
		x_d_est = x_est - path[0].x; //x_0 - x_0*
		
		x_true = start + sampleGaussian(zeros<X_DIM,1>(), P0);
		//x_true = start;

		for(int k = 1; k < path_size; k++){
			u_d = L[k-1] * x_d_est;
			u = path[k-1].u + u_d;

			//check if the control u violate the constraints.
		/*	if(u[0] < -0.4)
				u[0] = -0.4;
			if(u[0] > 0.8)
				u[0] = 0.8;
			if(u[1] < -M_PI*0.35)
				u[1] = -M_PI*0.35;
			if(u[1] > M_PI*0.35)
				u[1] = M_PI*0.35;*/

			//divide the segment to five segment to draw
			Dynamics dy(dt/5, car_l);
			m = sampleGaussian(zeros<U_DIM, 1>(), M);
				
			for(int iter = 1; iter <=5; iter++){
				x_true_old = x_true;
				x_true = dy.f(x_true_old, u, m);
				//x_true = dy.f(x_true_old, u, zeros<U_DIM,1>());
				dis += sqrt((x_true[0]-x_true_old[0])*(x_true[0]-x_true_old[0])+(x_true[1]-x_true_old[1])*(x_true[1]-x_true_old[1]));
				//check collision
				int col;
				CAL_CheckLineCollision(cal_environment, x_true_old[0], x_true_old[1], 0.0, x_true[0], x_true[1],0.0, false, &col);
				if(col != 0){ //if collision happen
					std::cout<<"LQG Controller fail: Collision happen"<<std::endl;
					return false;
				}

				//if no collision, then draw this small segment.
				int np[1] = {2};
				float p[6] = {x_true_old[0], x_true_old[1], 0.0, x_true[0], x_true[1], 0.0};
				CAL_CreatePolyline(cal_rrt, 1, np, p);


			}
			
			//kalman Filter
			x_d_est = A[k]*x_d_est + B[k]*u_d;
			//Dynamics dy1(dt, car_l);
			//x_d_est = dy1.f_d(x_d_est + path[k-1].x, u_d + path[k-1].u, zeros<U_DIM,1>()) - path[k].x;
			n = sampleGaussian(zeros<Z_DIM,1>(),N);
			
			z_d = dy.h(x_true, n) - dy.h(path[k].x, zeros<Z_DIM,1>());
			//z_d = dy.h(x_true, zeros<Z_DIM,1>()) - dy.h(path[k].x, zeros<Z_DIM,1>());
			x_d_est = x_d_est + K[k]*(z_d - H * x_d_est);

		}
		if((x_true[0]-goal[0])*(x_true[0]-goal[0]) + (x_true[1]-goal[1])*(x_true[1]-goal[1]) > goal_radius*goal_radius){
			std::cout<<"LQG Controller fail: Out of Goal Region" << std::endl;
			return false;
		}
		return true;
	}

	//for replanning, compute one step ahead
	bool simulateOneStep(Matrix<X_DIM>& start_true, Matrix<X_DIM>& start, Matrix<X_DIM, X_DIM>& P0, const int& cal_rrt, const int& cal_environment, const Matrix<DIM>& goal, double& dis){
	
		Matrix<U_DIM> u;
		Matrix<X_DIM> x_est = start;
		Matrix<X_DIM> x_true = start_true;
		Matrix<U_DIM> m;
		Matrix<X_DIM> x_true_old = start_true;
		Matrix<X_DIM> x_d_est = x_est - path[0].x;
		Matrix<X_DIM, X_DIM> P = P0;
		compute_A_B_V();
		compute_K(P0);
		u = path[0].u;

		//simulate one step using the above control input u.
		Dynamics dy(dt/5, car_l);
		for(int iter = 1; iter <=5; iter++){
			m = sampleGaussian(zeros<U_DIM, 1>(), M);
			x_true_old = x_true;
			x_true = dy.f(x_true_old, u, m);
			dis += sqrt((x_true[0]-x_true_old[0])*(x_true[0]-x_true_old[0])+(x_true[1]-x_true_old[1])*(x_true[1]-x_true_old[1]));
			//check collision
			int col;
			CAL_CheckLineCollision(cal_environment, x_true_old[0], x_true_old[1], 0.0, x_true[0], x_true[1],0.0, false, &col);
			if(col != 0){ //if collision happen
				std::cout<<"LQG Controller fail: Collision happen"<<std::endl;
				return false;
			}

			//if no collision, then draw this small segment.
			int np[1] = {2};
			float p[6] = {x_true_old[0], x_true_old[1], 0.0, x_true[0], x_true[1], 0.0};
			CAL_CreatePolyline(cal_rrt, 1, np, p);

		//    clock_t start = clock();
		//	while((double)((clock() - start) / CLOCKS_PER_SEC) < 0.1){;}
		}
		//return the true position of the robot after applying u.
		start_true = x_true;
		
		//kalman Filter
		x_d_est = A[1]*x_d_est + B[1]*(u - path[0].u);
		P = A[1]*P0*~A[1] + V[1]*M*~V[1];
		Matrix<Z_DIM> n;
		n = sampleGaussian(zeros<Z_DIM,1>(),N);
		Matrix<Z_DIM> z_d = dy.h(x_true, n) - dy.h(path[1].x, zeros<Z_DIM,1>());
		x_d_est = x_d_est + K[1]*(z_d - H * x_d_est);
		
		//return the estimate position of the robot and its covariance.
		start = path[1].x + x_d_est;
		P0 = (identity<X_DIM>() - K[1]*H)*P;
	
		return true;
	}





};


#endif