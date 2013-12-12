#define _USE_MATH_DEFINES
#include "planner.h"
#include <windows.h>

// Switch between belief vector and matrices
void Planner::unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma) {
	x = b.subMatrix<X_DIM,1>(0,0);
	int idx = X_DIM;
	for (int j = 0; j < X_DIM; ++j) {
		for (int i = j; i < X_DIM; ++i) {
			SqrtSigma(i,j) = SqrtSigma(j,i) = b[idx];
			++idx;
		}
	}
}

void Planner::vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b) {
	b.insert(0,0,x);
	int idx = X_DIM;
	for (int j = 0; j < X_DIM; ++j) {
		for (int i = j; i < X_DIM; ++i) {
			b[idx] = 0.5 * (SqrtSigma(i,j) + SqrtSigma(j,i));
			++idx;
		}
	}
}

// returns QvecS such that tr(SqrtS*Q*SqrtS) == ~vec(SqrtS) * QvecS * vec(SqrtS)
Matrix<S_DIM,S_DIM> Planner::QvecS (const Matrix<X_DIM,X_DIM>& Q) {
	Matrix<S_DIM,S_DIM> Qvec;

	int row_id = 0;
	for (int i = 0; i < X_DIM; ++i) {
		for (int j = i; j < X_DIM; ++j) {
			int col_id = 0;
			for (int k = 0; k < X_DIM; ++k) {
				for (int l = k; l < X_DIM; ++l) {
					if (i == j && i == k && i == l) { // i == j == k == l
						Qvec(row_id, col_id) = Q(i,i);
					} else if (i == k && j == l) {
						Qvec(row_id, col_id) = Q(i,i) + Q(j,j);
					} else if (i == k) {
						Qvec(row_id, col_id) = Q(j,l);
					} else if (j == l) {
						Qvec(row_id, col_id) = Q(i,k);
					} else {
						Qvec(row_id, col_id) = 0;
					}
					++col_id;
				}
			}
			++row_id;
		}
	}
	return Qvec;
}

// Belief dynamics
void Planner::beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM>& g, Matrix<X_DIM,X_DIM>& SqrtW) 
{
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;
	unVec(b, x, SqrtSigma);

	Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*~SqrtSigma;

	Matrix<X_DIM,X_DIM> A = dfdx(x, u, zeros<M_DIM,1>());
	Matrix<X_DIM,M_DIM> M = dfdm(x, u, zeros<M_DIM,1>());

	//x = f(x, u, zeros<M_DIM,1>());
	for(int i = 1; i <= nsteps; ++i) {
		x = f(x, u, zeros<M_DIM,1>(), DT*invnsteps);
	}

	Sigma = A*Sigma*~A + M*~M;

	Matrix<Z_DIM,X_DIM> H = dhdx(x, zeros<N_DIM,1>());
	Matrix<Z_DIM,N_DIM> N = dhdn(x, zeros<N_DIM,1>());
	Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*~N);

	Matrix<X_DIM, X_DIM> W = K*(H*Sigma);
	Sigma -= W;

	Matrix<X_DIM, X_DIM> V, D;

	jacobi2(Sigma, V, D);
	for (int i = 0; i < X_DIM; ++i) {
		if (D(i,i) > 0) {
			D(i,i) = sqrt(D(i,i));
		} else {
			D(i,i) = 0;
		}
	}
	SqrtSigma = V * D * ~V;
	vec(x, SqrtSigma, g);

	jacobi2(W, V, D);
	for (int i = 0; i < X_DIM; ++i) {
		if (D(i,i) > 0) {
			D(i,i) = sqrt(D(i,i));
		} else {
			D(i,i) = 0;
		}
	}
	SqrtW = V * D * ~V;
}

double Planner::confLine(const Matrix<X_DIM>& p1, const Matrix<X_DIM>& p2, const Matrix<X_DIM, X_DIM>& Sigma)
{
#if (DIM == 2)
	Matrix<DIM,DIM> V, E;
	jacobi2(Sigma.subMatrix<DIM,DIM>(0,0), V, E);

	Matrix<DIM> scale;
	scale[0] = sqrt(E(0,0)); scale[1] = sqrt(E(1,1));

	Matrix<DIM,DIM> invScale = zeros<DIM,DIM>();
	invScale(0,0) = 1.0/scale[0];
	invScale(1,1) = 1.0/scale[1];

	CAL_SetGroupOrientation(cal_obstacles, 0.0f, 0.0f, (float)mod2pi(atan2(V(0,1), V(0,0))) );
	CAL_SetGroupScaling(cal_environment, (float)invScale(0,0), (float)invScale(1,1), 0.0f);

	Matrix<DIM> tp1, tp2;
	tp1 = invScale * ~V * p1.subMatrix<DIM,1>(0,0);
	tp2 = invScale * ~V * p2.subMatrix<DIM,1>(0,0);
	double h = sqrt(tr(~(tp2 - tp1)*(tp2 - tp1)));
	Matrix<DIM> tp =  (tp1+tp2)*0.5;
	CAL_SetGroupScaling(cal_cylinder, 1.0f, (float)h, 1.0f);
	CAL_SetGroupPosition(cal_cylinder, (float)tp[0], (float)tp[1], 0.0f);
	CAL_SetGroupOrientation(cal_cylinder, 0.0f, 0.0f, (float)mod2pi(atan2(tp2[1] - tp1[1], tp2[0] - tp1[0]) - M_PI*0.5));

	int num_pairs;
	CAL_GetClosestPairs (cal_cylinder, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);
	double distance = results[0].distance;
	delete[] results;

	CAL_SetGroupOrientation(cal_obstacles, 0.0f, 0.0f, 0.0f);
	CAL_SetGroupScaling(cal_environment, 1.0f, 1.0f, 0.0f);

	CAL_SetGroupOrientation(cal_cylinder, 0.0f, 0.0f, 0.0f);
	CAL_SetGroupPosition(cal_cylinder, 0.0f, 0.0f, 0.0f);
	CAL_SetGroupScaling(cal_cylinder, 1.0f, 1.0f, 1.0f);

	// Display scaled ellipses (debug only)
	//int steps = 20;
	//double invsteps = 1.0/(double)steps;
	//for(int i = 0; i <= steps; ++i) 
	//{
	//	Matrix<DIM> pos = p1 + (p2 - p1)*invsteps*i;
	//	int id;
	//	CAL_CreateSphere(cal_environment, 1.0f, 0.0f, 0.0f, 0.0f, &id);
	//	CAL_SetObjectColor(id, 0.0f, 0.0f, 0.8f, 0.5f);
	//	CAL_SetObjectPosition(id, (float)pos[0], (float)pos[1], 0.0f);
	//	CAL_SetObjectOrientation(id, 0.0f, 0.0f, (float)mod2pi(atan2(V(1,0), V(0,0))));
	//	CAL_SetObjectScaling(id, (float)(sqrt(E(0,0))*distance), (float)(sqrt(E(1,1))*distance), 0.0f);
	//}

	return distance;
#elif (DIM == 3)
	Matrix<DIM,DIM> V, E;
	jacobi2(Sigma.subMatrix<DIM,DIM>(0,0), V, E);

	Matrix<DIM> scale;
	scale[0] = sqrt(E(0,0)); scale[1] = sqrt(E(1,1)); scale[2] = sqrt(E(2,2));

	Matrix<DIM,DIM> invScale = zeros<DIM,DIM>();
	invScale(0,0) = 1.0/scale[0];
	invScale(1,1) = 1.0/scale[1];
	invScale(2,2) = 1.0/scale[2];

	Matrix<4,1> q1 = quatFromRot(~V), q2;

	CAL_SetGroupQuaternion(cal_obstacles, (float)q1[0], (float)q1[1], (float)q1[2], (float)q1[3]);
	CAL_SetGroupScaling(cal_environment, (float)invScale(0,0), (float)invScale(1,1), (float)invScale(2,2));
	
	Matrix<DIM> tp1, tp2;
	tp1 = invScale * ~V * p1.subMatrix<DIM,1>(0,0);
	tp2 = invScale * ~V * p2.subMatrix<DIM,1>(0,0);

	Matrix<DIM> m, t, c;
	float len1, len2;
	
	m = (tp1 + tp2)*0.5;
	len1 = sqrt(tr((tp2 - tp1)*~(tp2 - tp1)));
	t = (tp2 - tp1)/len1;
	c[0] = -t[2]; c[1] = 0.0; c[2] = t[0];
	len2 = sqrt(t[0]*t[0] + t[2]*t[2]);
	c = c/len2;

	q2 = quatFromAA(c, acos(t[1]));
	CAL_SetGroupScaling(cal_cylinder, 1.0f, (float)len1, 1.0f);
	CAL_SetGroupPosition(cal_cylinder, (float)m[0], (float)m[1], (float)m[2]);
	CAL_SetGroupQuaternion(cal_cylinder, (float) q2[0], (float) q2[1], (float) q2[2], (float) q2[3]);
		
	int num_pairs;
	CAL_GetClosestPairs (cal_cylinder, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);
	double distance = results[0].distance;
	delete[] results;

	CAL_SetGroupQuaternion(cal_cylinder, 0.0f, 0.0f, 0.0f, 1.0f);
	CAL_SetGroupPosition(cal_cylinder, 0.0f, 0.0f, 0.0f);
	CAL_SetGroupScaling(cal_cylinder, 1.0f, 1.0f, 1.0f);

	CAL_SetGroupQuaternion(cal_obstacles, 0.0f, 0.0f, 0.0f, 1.0f);
	CAL_SetGroupScaling(cal_environment, 1.0f, 1.0f, 1.0f);
	
	return distance;
#endif
}

double Planner::confFunc(const double c)
{
	return -5.0*log(incompletegamma(DIM*0.5, 0.5*c*c));
}

double Planner::confObs(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u)
{
	// Penalize obstacle cost (line segment)
	Matrix<X_DIM> p1, p2;
	Matrix<X_DIM, X_DIM> SqrtSigma;

	Matrix<M_DIM> m = zeros<M_DIM,1>();

	unVec(b, p1, SqrtSigma);
	Matrix<X_DIM, X_DIM> Sigma = SqrtSigma * ~SqrtSigma;

	p2 = f(p1, u, m);

	int cols = 0;

#if (DIM == 2)
	CAL_CheckLineCollision(cal_environment, (float) p1[0], (float) p1[1], 0, (float) p2[0], (float) p2[1], 0, false, &cols);
#elif (DIM == 3)
	CAL_CheckLineCollision(cal_environment, (float) p1[0], (float) p1[1], (float) p1[2], (float) p2[0], (float) p2[1], (float) p2[2], false, &cols);
#endif
	if (cols != 0) {
		return 0.0;
	}

	double cl = confLine(p1, p2, Sigma);
	return cl;
}

// cost function
double Planner::cost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;

	unVec(b, x, SqrtSigma);
	return 0.5*tr(~u*Rint*u) + 0.5*tr(~SqrtSigma*Qint*SqrtSigma);
}

// final cost function
double Planner::finalCost(const Matrix<B_DIM>& b) {
	Matrix<X_DIM> x;
	Matrix<X_DIM,X_DIM> SqrtSigma;

	unVec(b, x, SqrtSigma);

	return 0.5*tr(~(x - xGoal)*QGoal*(x - xGoal)) + 0.5*tr(~SqrtSigma*QGoal*SqrtSigma);

	// log-likelihood
	//Matrix<X_DIM, X_DIM> Sigma = SqrtSigma * ~SqrtSigma;
	//Matrix<2> term = Sigma.subMatrix<2,2>(0,0)%(xGoal.subMatrix<2,1>(0,0) - x.subMatrix<2,1>(0,0));
	//return -1.0*(-1.0*log(det(Sigma.subMatrix<2,2>(0,0))) - 0.5*tr(~(xGoal.subMatrix<2,1>(0,0) - x.subMatrix<2,1>(0,0))*term));
}

// value function
double Planner::value(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss) 
{
	Matrix<B_DIM> g;

	Matrix<X_DIM,X_DIM> W;
	beliefDynamics(b, u, g, W);

	if (!switchoffML) {
		return cost(b, u) + 0.5*tr(~(g - bt1)*S*(g - bt1)) + tr(~(g - bt1)*s) + ss;
	} else {
		return cost(b, u) + 0.5*tr(~(g - bt1)*S*(g - bt1)) + tr(~(g - bt1)*s) + ss + 0.5*tr(S.subMatrix<X_DIM,X_DIM>(0,0) * W);
	}
}

// Compute second-order expansion of finalCost around b
void Planner::initValue(const Matrix<B_DIM>& b, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss) {
	ss = finalCost(b);

	Matrix<B_DIM> br(b), bl(b);
	for (int i = 0; i < B_DIM; ++i) {
		br[i] += H; bl[i] -= H; 
		s[i] = (finalCost(br) - finalCost(bl)) / (2.0*H);
		S(i,i) = ( finalCost(bl) - 2.0*ss + finalCost(br) ) / (H*H);
		br[i] = bl[i] = b[i];
	}

	Matrix<B_DIM> btr(b), btl(b), bbr(b), bbl(b);
	for (int i = 1; i < B_DIM; ++i) {
		btr[i] += H; btl[i] -= H; bbr[i] += H; bbl[i] -= H;
		for (int j = 0; j < i; ++j) {
			btr[j] += H; btl[j] += H; bbr[j] -= H; bbl[j] -= H;
			S(i,j) = S(j,i) = (finalCost(bbl) + finalCost(btr) - finalCost(btl) - finalCost(bbr)) / (4.0*H*H);
			btr[j] = btl[j] = bbr[j] = bbl[j] = b[j];
		}
		btr[i] = btl[i] = bbr[i] = bbl[i] = b[i];
	}

	//Matrix<B_DIM> bGoal = zeros<B_DIM,1>();
	//bGoal.insert(0,0, xGoal);
	//S = zeros<B_DIM,B_DIM>();
	//S.insert(0,0, QGoal);
	//S.insert(X_DIM,X_DIM, QvecS(QGoal));
	//s = S*(b - bGoal);
	//ss = 0.5*tr(~(b - bGoal)*S*(b - bGoal));
}

// Compute second-order expansion of cost around b and u
void Planner::quadratizeCost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& Q, Matrix<B_DIM>& q, 
							   Matrix<U_DIM,B_DIM>& P, Matrix<U_DIM,U_DIM>& R, Matrix<U_DIM>& r, double& p) {
	// p
	p = cost(b,u);

	// q, diag(Q)
	Matrix<B_DIM> br(b), bl(b);
	for (int i = 0; i < B_DIM; ++i) {
		br[i] += H; bl[i] -= H; 
		q[i] = (cost(br,u) - cost(bl,u)) / (2.0*H);      // O(n2 * n3) = O(n5)
		Q(i,i) = ( cost(bl,u) - 2.0*p + cost(br, u) ) / (H*H);
		br[i] = bl[i] = b[i];
	}

	// r, diag(R), P
	Matrix<U_DIM> ur(u), ul(u);
	for (int i = 0; i < U_DIM; ++i) {
		ur[i] += H; ul[i] -= H; 
		r[i] = (cost(b,ur) - cost(b,ul)) / (2.0*H);
		R(i,i) = ( cost(b,ul) - 2.0*p + cost(b, ur) ) / (H*H);
		for (int j = 0; j < B_DIM; ++j) {
			br[j] += H; bl[j] -= H;
			P(i,j) = (cost(bl, ul) + cost(br, ur) - cost(br, ul) - cost(bl, ur)) / (4.0*H*H);
			br[j] = bl[j] = b[j];
		}
		ur[i] = ul[i] = u[i];
	}

	// Q
	Matrix<B_DIM> btr(b), btl(b), bbr(b), bbl(b);
	for (int i = 1; i < B_DIM; ++i) {
		btr[i] += H; btl[i] -= H; bbr[i] += H; bbl[i] -= H;
		for (int j = 0; j < i; ++j) {
			btr[j] += H; btl[j] += H; bbr[j] -= H; bbl[j] -= H;
			Q(i,j) = Q(j,i) = (cost(bbl,u) + cost(btr,u) - cost(btl,u) - cost(bbr,u)) / (4.0*H*H); // O(n2 * n2 * n3) = O(n7)
			btr[j] = btl[j] = bbr[j] = bbl[j] = b[j];
		}
		btr[i] = btl[i] = bbr[i] = bbl[i] = b[i];
	}

	// R
	Matrix<U_DIM> utr(u), utl(u), ubr(u), ubl(u);
	for (int i = 1; i < U_DIM; ++i) {
		utr[i] += H; utl[i] -= H; ubr[i] += H; ubl[i] -= H;
		for (int j = 0; j < i; ++j) {
			utr[j] += H; utl[j] += H; ubr[j] -= H; ubl[j] -= H;
			R(i,j) = R(j,i) = (cost(b,ubl) + cost(b,utr) - cost(b,utl) - cost(b,ubr)) / (4.0*H*H);
			utr[j] = utl[j] = ubr[j] = ubl[j] = u[j];
		}
		utr[i] = utl[i] = ubr[i] = ubl[i] = u[i];
	}

	// Test positive-definiteness
	Matrix<B_DIM+U_DIM,B_DIM+U_DIM> CombinedHessian;
	CombinedHessian.insert(0,0,Q);
	CombinedHessian.insert(0,B_DIM,~P);
	CombinedHessian.insert(B_DIM,0,P);
	CombinedHessian.insert(B_DIM,B_DIM,R);

	Matrix<B_DIM+U_DIM, B_DIM+U_DIM> V, E;
	jacobi2(CombinedHessian, V, E);
	for (int i = 0; i < B_DIM+U_DIM; ++i) {
		if (E(i,i) < -1e-10) {
			std::cerr << "CombinedHessian is not positive-definite!" << std::endl;
			std::exit(-1);
		}
	}

	/*
	double delta = 0.001;
	Matrix<B_DIM+U_DIM, B_DIM+U_DIM> chvec, chval;
	jacobi2(CombinedHessian, chvec, chval);
	double minchval = chval(0,0);
	for (int i = 1; i < U_DIM; ++i) {
		if (chval(i,i) < minchval) {
			minchval = chval(i,i);
		}
	}
	double CHeps = 0.0;
	if (minchval < delta)
		CHeps = delta - minchval;
	CombinedHessian = CombinedHessian + CHeps*identity<B_DIM+U_DIM>();

	forcePD(CombinedHessian);

	P = CombinedHessian.subMatrix<U_DIM,B_DIM>(B_DIM,0);
	R = CombinedHessian.subMatrix<U_DIM,U_DIM>(B_DIM,B_DIM);
	Q = CombinedHessian.subMatrix<B_DIM,B_DIM>(0,0);
	*/
}

/*
// Compute closed form quadratic finalCost function around around b and u
void Planner::initValue(const Matrix<B_DIM>& b, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss) 
{
Matrix<B_DIM> bGoal = zeros<B_DIM,1>();
bGoal.insert(0,0, xGoal);

S = zeros<B_DIM,B_DIM>();
S.insert(0,0, QGoal);
S.insert(X_DIM,X_DIM, QvecS(QGoal));
s = S*(b - bGoal);
ss = 0.5*tr(~(b - bGoal)*S*(b - bGoal));
}

// Compute closed form quadratic cost function around around b and u
void Planner::quadratizeCost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& Q, Matrix<B_DIM>& q, 
Matrix<U_DIM,B_DIM>& P, Matrix<U_DIM,U_DIM>& R, Matrix<U_DIM>& r, double& p) 
{
Q = zeros<B_DIM,B_DIM>();
Q.insert(X_DIM,X_DIM, QvecS(Qint));  // penalize uncertainty
q = Q*b;
R = Rint;
r = R*u;
P = zeros<U_DIM,B_DIM>();
p = 0.5*tr(~u*R*u) + 0.5*tr(~b*Q*b);
}
*/

// Jacobian df/dx(x,u,m)
Matrix<X_DIM,X_DIM> Planner::dfdx(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m) {
	Matrix<X_DIM,X_DIM> A;
	Matrix<X_DIM> xr(x), xl(x);
	for (int i = 0; i < X_DIM; ++i) {
		xr[i] += H; xl[i] -= H;
		A.insert(0,i, (f(xr, u, m) - f(xl, u, m)) / (2.0*H));
		xr[i] = xl[i] = x[i];
	}
	return A;
}

// Jacobian df/du(x,u,m)
Matrix<X_DIM,U_DIM> Planner::dfdu(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m) {
	Matrix<X_DIM,U_DIM> B;
	Matrix<U_DIM> ur(u), ul(u);
	for (int i = 0; i < U_DIM; ++i) {
		ur[i] += H; ul[i] -= H;
		B.insert(0,i, (f(x, ur, m) - f(x, ul, m)) / (2.0*H));
		ur[i] = ul[i] = u[i];
	}
	return B;
}

// Jacobian df/dm(x,u,m)
Matrix<X_DIM,M_DIM> Planner::dfdm(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m) {
	Matrix<X_DIM,M_DIM> M;
	Matrix<M_DIM> mr(m), ml(m);
	for (int i = 0; i < M_DIM; ++i) {
		mr[i] += H; ml[i] -= H;
		M.insert(0,i, (f(x, u, mr) - f(x, u, ml)) / (2.0*H));
		mr[i] = ml[i] = m[i];
	}
	return M;
}

// Jacobian dh/dx(x,n)
Matrix<Z_DIM,X_DIM> Planner::dhdx(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) {
	Matrix<Z_DIM,X_DIM> Hx;
	Matrix<X_DIM> xr(x), xl(x);
	for (int i = 0; i < X_DIM; ++i) {
		xr[i] += H; xl[i] -= H;
		Hx.insert(0,i, (h(xr, n) - h(xl, n)) / (2.0*H));
		xr[i] = xl[i] = x[i];
	}
	return Hx;
}

// Jacobian dh/dn(x,n)
Matrix<Z_DIM,N_DIM> Planner::dhdn(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) {
	Matrix<Z_DIM,N_DIM> N;
	Matrix<N_DIM> nr(n), nl(n);
	for (int i = 0; i < N_DIM; ++i) {
		nr[i] += H; nl[i] -= H;
		N.insert(0,i, (h(x, nr) - h(x, nl)) / (2.0*H));
		nr[i] = nl[i] = n[i];
	}
	return N;
}

// Jacobian dc/db(b,u)
Matrix<1,B_DIM> Planner::dcdb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss) {
	Matrix<1,B_DIM> c;
	for (int i = 0; i < B_DIM; ++i) {
		Matrix<B_DIM> br(b);
		br[i] += H;
		Matrix<B_DIM> bl(b);
		bl[i] -= H; 

		c(0,i) = (value(br, u, bt1, S, s, ss) - value(bl, u, bt1, S, s, ss)) / (2*H);
	}
	return c;
}

// Jacobian dc/du(b,u)
Matrix<1,U_DIM> Planner::dcdu(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss) {
	Matrix<1,U_DIM> d;
	for (int i = 0; i < U_DIM; ++i) {
		Matrix<U_DIM> ur(u);
		ur[i] += H;
		Matrix<U_DIM> ul(u);
		ul[i] -= H; 

		d(0,i) = (value(b, ur, bt1, S, s, ss) - value(b, ul, bt1, S, s, ss)) / (2*H);
	}
	return d;
}

// Jacobian dv/db(b)
Matrix<1,B_DIM> Planner::dvdb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<U_DIM,B_DIM>& L, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss) {
	Matrix<1,B_DIM> d;
	for (int i = 0; i < B_DIM; ++i) {
		Matrix<B_DIM> br(b);
		br[i] += H;
		Matrix<B_DIM> bl(b);
		bl[i] -= H; 

		d(0,i) = (value(br, u + L*(br - b), bt1, S, s, ss) - value(bl, u + L*(bl - b), bt1, S, s, ss)) / (2*H);
	}
	return d;
}

// Hessian d^2c/dudu(b,u)
Matrix<U_DIM,U_DIM> Planner::ddcdudu(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss) {
	Matrix<U_DIM,U_DIM> D;

	double midcost = value(b, u, bt1, S, s, ss);

	for (int i = 0; i < U_DIM; ++i) {
		Matrix<U_DIM> ur(u);
		ur[i] += H;
		Matrix<U_DIM> ul(u);
		ul[i] -= H;
		D(i,i) = ( value(b, ul, bt1, S, s, ss) - 2*midcost + value(b, ur, bt1, S, s, ss) ) / (H*H);
	}

	for (int i = 1; i < U_DIM; ++i) {
		for (int j = 0; j < i; ++j) {
			Matrix<U_DIM> utr(u);
			utr[i] += H; utr[j] += H;
			Matrix<U_DIM> utl(u);
			utl[i] -= H; utl[j] += H;
			Matrix<U_DIM> ubr(u);
			ubr[i] += H; ubr[j] -= H;
			Matrix<U_DIM> ubl(u);
			ubl[i] -= H; ubl[j] -= H;

			D(i,j) = D(j,i) = (value(b, ubl, bt1, S, s, ss) + value(b, utr, bt1, S, s, ss) - value(b, utl, bt1, S, s, ss) - value(b, ubr, bt1, S, s, ss)) / (4*H*H);
		}
	}
	return D;
}

// Hessian d^2c/dbdb(b,u)
Matrix<B_DIM,B_DIM> Planner::ddcdbdb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss) {
	Matrix<B_DIM,B_DIM> C;

	double midcost = value(b, u, bt1, S, s, ss);

	for (int i = 0; i < B_DIM; ++i) {
		Matrix<B_DIM> br(b);
		br[i] += H;
		Matrix<B_DIM> bl(b);
		bl[i] -= H;
		C(i,i) = ( value(bl, u, bt1, S, s, ss) - 2*midcost + value(br, u, bt1, S, s, ss) ) / (H*H);
	}

	for (int i = 1; i < B_DIM; ++i) {
		for (int j = 0; j < i; ++j) {
			Matrix<B_DIM> btr(b);
			btr[i] += H; btr[j] += H;
			Matrix<B_DIM> btl(b);
			btl[i] -= H; btl[j] += H;
			Matrix<B_DIM> bbr(b);
			bbr[i] += H; bbr[j] -= H;
			Matrix<B_DIM> bbl(b);
			bbl[i] -= H; bbl[j] -= H;

			C(i,j) = C(j,i) = (value(bbl, u, bt1, S, s, ss) + value(btr, u, bt1, S, s, ss) - value(btl, u, bt1, S, s, ss) - value(bbr, u, bt1, S, s, ss)) / (4*H*H);
		}
	}
	return C;
}

// Hessian d^2c/dudb(b,u)
Matrix<U_DIM,B_DIM> Planner::ddcdudb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss) {
	Matrix<U_DIM,B_DIM> E;

	for (int i = 0; i < U_DIM; ++i) {
		Matrix<U_DIM> ur(u);
		ur[i] += H;
		Matrix<U_DIM> ul(u);
		ul[i] -= H; 
		for (int j = 0; j < B_DIM; ++j) {
			Matrix<B_DIM> br(b);
			br[j] += H;
			Matrix<B_DIM> bl(b);
			bl[j] -= H;

			E(i,j) = (value(bl, ul, bt1, S, s, ss) + value(br, ur, bt1, S, s, ss) - value(br, ul, bt1, S, s, ss) - value(bl, ur, bt1, S, s, ss)) / (4*H*H);
		}
	}
	return E;
}

Matrix<B_DIM,B_DIM> Planner::ddvdbdb(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, const Matrix<U_DIM,B_DIM>& L, 
									  const Matrix<B_DIM,B_DIM>& S, const Matrix<B_DIM>& s, double ss) 
{
	Matrix<B_DIM,B_DIM> C;

	double midval = value(b, u, bt1, S, s, ss);

	for (int i = 0; i < B_DIM; ++i) {
		Matrix<B_DIM> br(b);
		br[i] += H;
		Matrix<B_DIM> bl(b);
		bl[i] -= H;
		C(i,i) = ( value(bl, u + L*(bl - b), bt1, S, s, ss) - 2*midval + value(br, u + L*(br - b), bt1, S, s, ss) ) / (H*H);
	}

	for (int i = 1; i < B_DIM; ++i) {
		for (int j = 0; j < i; ++j) {
			Matrix<B_DIM> btr(b);
			btr[i] += H; btr[j] += H;
			Matrix<B_DIM> btl(b);
			btl[i] -= H; btl[j] += H;
			Matrix<B_DIM> bbr(b);
			bbr[i] += H; bbr[j] -= H;
			Matrix<B_DIM> bbl(b);
			bbl[i] -= H; bbl[j] -= H;

			C(i,j) = C(j,i) = (value(bbl, u + L*(bbl - b), bt1, S, s, ss) + value(btr, u + L*(btr - b), bt1, S, s, ss) - 
							   value(btl, u + L*(btl - b), bt1, S, s, ss) - value(bbr, u + L*(bbr - b), bt1, S, s, ss)) / (4*H*H);
		}
	}
	return C;
}


void Planner::linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM>& c, Matrix<B_DIM,B_DIM>& A, Matrix<B_DIM, U_DIM>& B, 
										std::vector<Matrix<X_DIM> >& cn, std::vector<Matrix<X_DIM, B_DIM> >& An, std::vector<Matrix<X_DIM, U_DIM> >& Bn) 
{
	cn.resize(X_DIM);
	An.resize(X_DIM);
	Bn.resize(X_DIM);

	Matrix<X_DIM, X_DIM> sqrtW;
	// set c
	beliefDynamics(b, u, c, sqrtW);  // O(n3)
	// set cn's
	for (int i = 0; i < X_DIM; ++i) {
		cn[i] = sqrtW.subMatrix<X_DIM, 1>(0,i);
	}

	// Compute Jacobian's
	Matrix<B_DIM> br(b), bl(b);
	for (int i = 0; i < B_DIM; ++i) {
		br[i] += H; bl[i] -= H;

		Matrix<B_DIM> cr, cl;
		Matrix<X_DIM,X_DIM> sqrtWr, sqrtWl;

		beliefDynamics(br, u, cr, sqrtWr);  // O(n2 * n3) = O(n5)
		beliefDynamics(bl, u, cl, sqrtWl);

		A.insert(0,i, (cr - cl) / (2.0*H));

		for(int j = 0; j < X_DIM; ++j) {
			An[j].insert(0,i, (sqrtWr.subMatrix<X_DIM,1>(0,j) - sqrtWl.subMatrix<X_DIM,1>(0,j)) / (2.0*H));
		}

		br[i] = bl[i] = b[i];
	}

	Matrix<U_DIM> ur(u), ul(u);
	for (int i = 0; i < U_DIM; ++i) {
		ur[i] += H; ul[i] -= H;

		Matrix<B_DIM> cr, cl;
		Matrix<X_DIM,X_DIM> sqrtWr, sqrtWl;

		beliefDynamics(b, ur, cr, sqrtWr);
		beliefDynamics(b, ul, cl, sqrtWl);

		B.insert(0,i, (cr - cl) / (2.0*H));

		for(int j = 0; j < X_DIM; ++j) {
			Bn[j].insert(0,i, (sqrtWr.subMatrix<X_DIM,1>(0,j) - sqrtWl.subMatrix<X_DIM,1>(0,j)) / (2.0*H));
		}

		ur[i] = ul[i] = u[i];
	}
}

void Planner::obsApprox(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<U_DIM, U_DIM>& D, Matrix<B_DIM, B_DIM>& C, 
						  Matrix<U_DIM, B_DIM>& E, Matrix<1,B_DIM>& c, Matrix<1, U_DIM>& d, double& e)
{
	Matrix<1,B_DIM> A;
	for (int i = 0; i < B_DIM; ++i) {
		Matrix<B_DIM> br(b);
		br[i] += H;
		Matrix<B_DIM> bl(b);
		bl[i] -= H; 

		double delta = confObs(br, u) - confObs(bl, u);
		A(0,i) = delta / (2*H);
	}
	Matrix<1,U_DIM> B;
	for (int i = 0; i < U_DIM; ++i) {
		Matrix<U_DIM> ur(u);
		ur[i] += H;
		Matrix<U_DIM> ul(u);
		ul[i] -= H; 

		double delta = confObs(b, ur) - confObs(b, ul);
		B(0,i) = delta / (2*H);
	}

	double cObs = confObs(b, u);
	double kappa = confFunc(cObs);
	double beta = (confFunc(cObs + H) - confFunc(cObs - H)) / (2*H);
	double H2 = H;
	double alpha = (confFunc(cObs + H2) - 2.0*kappa + confFunc(cObs - H2)) / (H2*H2);

	// Original
	e = kappa;
	c = beta*A;
	d = beta*B;
	C = alpha*~A*A;
	D = alpha*~B*B;
	E = alpha*~B*A;
}

void Planner::obsApprox2(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<U_DIM, B_DIM>& L, Matrix<B_DIM, B_DIM>& C, Matrix<1,B_DIM>& c, double& e)
{
	double h = H;

	double cObs = confObs(b, u);
	if (cObs == 0.0) {
		e = COST_INFTY;
		return;
	}

	Matrix<1,B_DIM> A;
	for (int i = 0; i < B_DIM; ++i) {
		Matrix<B_DIM> br(b);
		br[i] += h;
		Matrix<B_DIM> bl(b);
		bl[i] -= h; 

		double delta = confObs(br, u + L*(br - b)) - confObs(bl, u + L*(bl - b));
		A(0,i) = delta / (2*h);
	}

	// dydx, ddydxdx, constant
	double kappa = confFunc(cObs);
	double beta = (confFunc(cObs + H) - confFunc(cObs - H)) / (2*H);
	double H2 = H;
	double alpha = (confFunc(cObs + H2) - 2.0*kappa + confFunc(cObs - H2)) / (H2*H2);

	// Original
	e = kappa;
	c = beta*A;
	C = alpha*~A*A;
}

void Planner::controlPolicyDDP(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, Matrix<B_DIM,B_DIM>& S, 
							 Matrix<B_DIM>& s, double& ss, Matrix<U_DIM,B_DIM>& L, Matrix<U_DIM>& l) 
{
#if defined(OBSTACLES)
	Matrix<U_DIM, U_DIM> DObs;
	Matrix<B_DIM, B_DIM> CObs;
	Matrix<1, B_DIM> cObs;
	Matrix<1, U_DIM> dObs;
	Matrix<U_DIM, B_DIM> EObs;
	double eObs;

	obsApprox(b, u, DObs, CObs, EObs, cObs, dObs, eObs);

	Matrix<U_DIM,U_DIM> D = ddcdudu(b, u, bt1, S, s, ss) + DObs;
	Matrix<U_DIM,B_DIM> E = ddcdudb(b, u, bt1, S, s, ss) + EObs;
	Matrix<U_DIM> d = ~dcdu(b, u, bt1, S, s, ss) + ~dObs;
	Matrix<B_DIM,B_DIM> C = ddcdbdb(b, u, bt1, S, s, ss) + CObs;
	Matrix<B_DIM> c = ~dcdb(b, u, bt1, S, s, ss) + ~cObs;
	double e = value(b, u, bt1, S, s, ss) + eObs;
#else
	Matrix<U_DIM,U_DIM> D = ddcdudu(b, u, bt1, S, s, ss);
	Matrix<U_DIM,B_DIM> E = ddcdudb(b, u, bt1, S, s, ss);
	Matrix<U_DIM> d = ~dcdu(b, u, bt1, S, s, ss);
	Matrix<B_DIM,B_DIM> C = ddcdbdb(b, u, bt1, S, s, ss);
	Matrix<B_DIM> c = ~dcdb(b, u, bt1, S, s, ss);
	double e = value(b, u, bt1, S, s, ss);
#endif

	Matrix<B_DIM+U_DIM, B_DIM+U_DIM> CombinedHessian;
	CombinedHessian.insert(0,0, C);
	CombinedHessian.insert(0,B_DIM, ~E);
	CombinedHessian.insert(B_DIM,0, E);
	CombinedHessian.insert(B_DIM,B_DIM, D - Rint);
	
	double delta = 0.001;
	Matrix<B_DIM+U_DIM, B_DIM+U_DIM> chvec, chval;
	jacobi2(CombinedHessian, chvec, chval);
	double minchval = chval(0,0);
	for (int i = 1; i < U_DIM; ++i) {
		if (chval(i,i) < minchval) {
			minchval = chval(i,i);
		}
	}
	double CHeps = 0.0;
	if (minchval < delta)
		CHeps = delta - minchval;
	CombinedHessian = CombinedHessian + CHeps*identity<B_DIM+U_DIM>();

	forcePD(CombinedHessian);

	E = CombinedHessian.subMatrix<U_DIM,B_DIM>(B_DIM,0);
	D = CombinedHessian.subMatrix<U_DIM,U_DIM>(B_DIM,B_DIM) + Rint;
	C = CombinedHessian.subMatrix<B_DIM,B_DIM>(0,0);

	Matrix<U_DIM, U_DIM> invD = !D;

	L = -invD*E;
	l = -invD*d;

	// Compute S, s, and ss using valueIteration
	S = C + ~E*L;
	s = c + ~E*l;
	ss = e + 0.5*tr(~d*l);

	gSqr += tr(~d*d);
	Jp += 0.5*tr(~d*l);
}

void Planner::controlPolicyILQG(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, Matrix<U_DIM,B_DIM>& L, Matrix<U_DIM>& l) 
{
	Matrix<B_DIM, B_DIM> A;
	Matrix<B_DIM, U_DIM> B;
	Matrix<B_DIM> g;
	std::vector<Matrix<X_DIM> > cn(X_DIM);
	std::vector<Matrix<X_DIM, B_DIM> > An(X_DIM);
	std::vector<Matrix<X_DIM, U_DIM> > Bn(X_DIM);

	double p;
	Matrix<U_DIM,B_DIM> P;
	Matrix<B_DIM,B_DIM> Q;
	Matrix<U_DIM,U_DIM> R;
	Matrix<B_DIM> q;
	Matrix<U_DIM> r;

	linearizeBeliefDynamics(b, u, g, A, B, cn, An, Bn); // O(n5)
	quadratizeCost(b, u, Q, q, P, R, r, p);

#if defined(OBSTACLES)
	Matrix<U_DIM, U_DIM> DObs;
	Matrix<B_DIM, B_DIM> CObs;
	Matrix<1, B_DIM> cObs;
	Matrix<1, U_DIM> dObs;
	Matrix<U_DIM, B_DIM> EObs;
	double eObs;

	obsApprox(b, u, DObs, CObs, EObs, cObs, dObs, eObs);

	Matrix<U_DIM,U_DIM> D = R + ~B*S*B + DObs;
	Matrix<U_DIM,B_DIM> E = P + ~B*S*A + EObs;
	Matrix<U_DIM> d = r + ~B*s + ~dObs;
	Matrix<B_DIM,B_DIM> C = Q + ~A*S*A + CObs;  // O(n6)
	Matrix<B_DIM> c = q + ~A*s + ~cObs;
	double e = p + ss + eObs;
#else
	Matrix<U_DIM,U_DIM> D = R + ~B*S*B;
	Matrix<U_DIM,B_DIM> E = P + ~B*S*A;
	Matrix<U_DIM> d = r + ~B*s;
	Matrix<B_DIM,B_DIM> C = Q + ~A*S*A;  // O(n6)
	Matrix<B_DIM> c = q + ~A*s;
	double e = p + ss;
#endif

	Matrix<X_DIM,X_DIM> Sx = S.subMatrix<X_DIM,X_DIM>(0,0);
	for (int i = 0; i < X_DIM; ++i) {
		D += ~Bn[i] * Sx * Bn[i];
		E += ~Bn[i] * Sx * An[i];
		d += ~Bn[i] * (Sx * cn[i]);
		C += ~An[i] * Sx * An[i];     // O(n*n5) = O(n6)
		c += ~An[i] * (Sx * cn[i]);
		e += 0.5*tr(~cn[i] * Sx * cn[i]);
	}

	Matrix<U_DIM,U_DIM> invD = !D;

	// control policy du = L dx + l
	L = -invD*E;
	l = -invD*d;

	//if (tr(~l*l) > 0.01*0.01) {
	//	l = l*0.01/sqrt(tr(~l*l));
	//}

	// update value function
	S = C + ~E*L;
	s = c + ~E*l;
	ss = e + 0.5*tr(~d*l);

	//S = C + ~L*D*L + ~L*E + ~E*L;
	//s = c + ~L*D*l + ~L*d + ~E*l;
	//ss = e + 0.5*tr(~l*D*l) + tr(~l*d);

	gSqr += tr(~d*d);
	Jp += 0.5*tr(~d*l);
}

// DDP value iteration
void Planner::valueIteration(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, const Matrix<B_DIM>& bt1, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, const Matrix<U_DIM,B_DIM>& L) 
{
#if defined(OBSTACLES)
	Matrix<B_DIM, B_DIM> CObs;
	Matrix<1,B_DIM> cObs;
	double eObs;

	obsApprox2(b, u, L, CObs, cObs, eObs);

	if (eObs > COST_INFTY) {
		ss = eObs;
		return;
	}

	Matrix<B_DIM,B_DIM> SNew = ddvdbdb(b, u, bt1, L, S, s, ss) + CObs;

	// Check if SNew is nonnegative-definite?
	forcePD(SNew);
	
	Matrix<B_DIM> sNew = ~dvdb(b, u, bt1, L, S, s, ss) + ~cObs;
	double ssNew = value(b, u, bt1, S, s, ss) + eObs;
#else
	Matrix<B_DIM,B_DIM> SNew = ddvdbdb(b, u, bt1, L, S, s, ss);

	// Check if SNew is nonnegative-definite?
	forcePD(SNew);
	
	Matrix<B_DIM> sNew = ~dvdb(b, u, bt1, L, S, s, ss);
	double ssNew = value(b, u, bt1, S, s, ss);
#endif

	S = SNew; s = sNew; ss = ssNew;
}

void Planner::expectedCost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, const Matrix<U_DIM,B_DIM>& L) 
{
	// u_dev = L * x_dev

	Matrix<B_DIM, B_DIM> A;
	Matrix<B_DIM, U_DIM> B;
	Matrix<B_DIM> g;
	std::vector<Matrix<X_DIM> > cn(X_DIM);
	std::vector<Matrix<X_DIM, B_DIM> > An(X_DIM);
	std::vector<Matrix<X_DIM, U_DIM> > Bn(X_DIM);

	double p;
	Matrix<U_DIM,B_DIM> P;
	Matrix<B_DIM,B_DIM> Q;
	Matrix<U_DIM,U_DIM> R;
	Matrix<B_DIM> q;
	Matrix<U_DIM> r;

	linearizeBeliefDynamics(b, u, g, A, B, cn, An, Bn);
	quadratizeCost(b, u, Q, q, P, R, r, p);

#if defined(OBSTACLES)
	Matrix<B_DIM, B_DIM> CObs;
	Matrix<1,B_DIM> cObs;
	double eObs;

	obsApprox2(b, u, L, CObs, cObs, eObs);

	if (eObs > COST_INFTY) {
		ss = eObs;
		return;
	}

	Matrix<B_DIM,B_DIM> C = Q + ~(A + B*L)*S*(A + B*L) + ~L*R*L + ~L*P + ~P*L + CObs;
	Matrix<B_DIM> c = q + ~(A + B*L)*s + ~L*r + ~cObs;
	double e = p + ss + eObs;
#else
	Matrix<B_DIM,B_DIM> C = Q + ~(A + B*L)*S*(A + B*L) + ~L*R*L + ~L*P + ~P*L;
	Matrix<B_DIM> c = q + ~(A + B*L)*s + ~L*r;
	double e = p + ss;
#endif

	Matrix<X_DIM,X_DIM> Sx = S.subMatrix<X_DIM,X_DIM>(0,0);
	for (int i = 0; i < X_DIM; ++i) {
		C += ~(An[i] + Bn[i]*L) * Sx * (An[i] + Bn[i]*L);
		c += ~(An[i] + Bn[i]*L) * (Sx * cn[i]);
		e += 0.5*tr(~cn[i] * Sx * cn[i]);
	}

	// update value function
	S = C;
	s = c;
	ss = e;
}

void Planner::displayMeanPath(const std::vector<Matrix<B_DIM> > & B)
{
	// Clear out group for visualization
	CAL_DestroyGroup(cal_path);

	// Display final path
	CAL_CreateGroup(&cal_path, 0, false, "Path");
	//CAL_SetGroupColor(cal_path, 1.0f, 0.8f, 0.0f);
	CAL_SetGroupColor(cal_path, 0.1f, 0.1f, 0.1f);

#if (DIM == 2)
	int id;
	Matrix<X_DIM> xt = x0;
	CAL_CreateSphere(cal_path, 0.04f, (float)x0[0], (float)x0[1], 0.1f, &id);
	//CAL_SetObjectColor(id, 1.0f, 0.7f, 0.0f);
	CAL_SetObjectColor(id, 0.1f, 0.1f, 0.1f);

	//nsteps = 10;
	//invnsteps = 0.1;
	// Integrated path
	for (int i = 0; i < pathlen; ++i) 
	{
		int np[1] = {nsteps+1};
		float* side = new float[3*(nsteps+1)];
		side[0] = (float)xt[0]; side[1] = (float)xt[1]; side[2] = 0.1f;

		for(int j = 1; j <= nsteps; ++j) {
			xt = f(xt, U[i], zeros<M_DIM,1>(), DT*invnsteps);
			side[3*j] = (float)xt[0]; side[3*j+1] = (float)xt[1]; side[3*j+2] = 0.1f;
		}

		CAL_CreatePolyline(cal_path, 1, np, side, &id);
		delete[] side;
	}
	//nsteps = 1;
	//invnsteps = 1.0;
	
	// Beliefs
	for (int i = 0; i <= pathlen; ++i) 
	{
		Matrix<X_DIM> x;
		Matrix<X_DIM, X_DIM> SqrtSigma;
		unVec(B[i], x, SqrtSigma);
		Matrix<X_DIM, X_DIM> Sigma = SqrtSigma * ~SqrtSigma;

		CAL_CreateSphere(cal_path, 0.04f, (float)x[0], (float)x[1], 0.1f, &id);
		//CAL_SetObjectColor(id, 1.0f, 0.7f, 0.0f);
		CAL_SetObjectColor(id, 0.1f, 0.1f, 0.1f);

		//if (i == 0) {
			//CAL_CreateSphere(cal_path, 1.0f, 0.0f, 0.0f, 0.0f, &id);
			//CAL_SetObjectColor(id, 1.0f, 0.0f, 0.0f);
			CAL_CreateUserDrawn(cal_path, drawUnitCircle, NULL, 0.0f, 0.0f, 0.0f, &id);
			Matrix<DIM,DIM> V, E;
			jacobi(Sigma.subMatrix<DIM,DIM>(0,0), V, E);
			CAL_SetObjectPosition(id, (float)x[0], (float)x[1], (float)0.1);
			CAL_SetObjectOrientation(id, 0.0f, 0.0f, (float)mod2pi(atan2(V(1,0), V(0,0))) );
			CAL_SetObjectScaling(id, (float)(sqrt(E(0,0))), (float)(sqrt(E(1,1))), 1.0f);
		//}
	}
#elif (DIM == 3)
	int id;
	Matrix<X_DIM> xt = x0;
	CAL_CreateSphere(cal_path, 0.04f, (float)x0[0], (float)x0[1], (float)x0[2], &id);
	//CAL_SetObjectColor(id, 1.0f, 0.7f, 0.0f);
	CAL_SetObjectColor(id, 0.1f, 0.1f, 0.1f);

	nsteps = 10;
	invnsteps = 0.1;
	// Integrated path
	for (int i = 0; i < pathlen; ++i) 
	{
		int np[1] = {nsteps+1};
		float* side = new float[3*(nsteps+1)];
		side[0] = (float)xt[0]; side[1] = (float)xt[1]; side[2] = (float)xt[2];

		for(int j = 1; j <= nsteps; ++j) {
			xt = f(xt, U[i], zeros<M_DIM,1>(), DT*invnsteps);
			side[3*j] = (float)xt[0]; side[3*j+1] = (float)xt[1]; side[3*j+2] = (float)xt[2];
		}

		CAL_CreatePolyline(cal_path, 1, np, side, &id);
		delete[] side;
	}
	nsteps = 1;
	invnsteps = 1.0;

	//Matrix<3> p1, p2, m, t, c;
	//Matrix<4,1> q;
	//float len1, len2;
	//int obj;

	//p1[0] = (float)path[0].T(0,3); p1[1] = (float)path[0].T(1,3); p1[2] = (float)path[0].T(2,3);
	//for (int i = 1; i < l; ++i) {
	//	p2[0] = (float)path[i].T(0,3); p2[1] = (float)path[i].T(1,3); p2[2] = (float)path[i].T(2,3);
	//	m = (p1 + p2)*0.5;
	//	len1 = sqrt(tr((p2 - p1)*~(p2 - p1)));
	//	t = (p2 - p1)/len1;
	//	c[0] = -t[2]; c[1] = 0.0; c[2] = t[0];
	//	len2 = sqrt(t[0]*t[0] + t[2]*t[2]);
	//	c = c/len2;
	//	CAL_CreateCylinder(groupid, 0.04f, len1, 0.0f, 0.0f, 0.0f, &obj);
	//	q = quatFromAA(c, acos(t[1]));
	//	CAL_SetObjectQuaternion(obj, (float) q[0], (float) q[1], (float) q[2], (float) q[3]);
	//	CAL_SetObjectPosition(obj, (float)m[0], (float)m[1], (float)m[2]);
	//	p1 = p2;
	//}

	// Beliefs
	for (int i = 0; i <= pathlen; ++i) 
	{
		Matrix<X_DIM> x;
		Matrix<X_DIM, X_DIM> SqrtSigma;
		unVec(B[i], x, SqrtSigma);
		Matrix<X_DIM, X_DIM> Sigma = SqrtSigma * ~SqrtSigma;

		//CAL_CreateSphere(cal_path, 0.04f, (float)x[0], (float)x[1], (float)x[2], &id);
		//CAL_SetObjectColor(id, 1.0f, 0.7f, 0.0f);
		//CAL_SetObjectColor(id, 0.1f, 0.1f, 0.1f);

		if (i == pathlen || i%2 == 0) {
			CAL_CreateSphere(cal_path, 1.0f, 0.0f, 0.0f, 0.0f, &id);
			CAL_SetObjectColor(id, 0.0f, 0.0f, 0.7f);
			
			//CAL_CreateUserDrawn(cal_path, drawUnitSphere, NULL, 0.0f, 0.0f, 0.0f, &id);
			
			Matrix<DIM,DIM> V, E;
			jacobi(Sigma.subMatrix<DIM,DIM>(0,0), V, E);
			CAL_SetObjectPosition(id, (float)x[0], (float)x[1], (float)x[2]);

			Matrix<3,3> R = identity<3>();
			R.insert(0,0, V);
			Matrix<4,1> q = quatFromRot(R);
			CAL_SetObjectQuaternion(id, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));

			double scale = 1.0;
			(i == 0)? scale = 0.5:scale = 2.0;
			CAL_SetObjectScaling(id, (float) (scale*sqrt(E(0,0))), (float) (scale*sqrt(E(1,1))), (float) (scale*sqrt(E(2,2))));
		}
	}
		
	//Matrix<3,3> V, E;
	//jacobi(S, V, E);
	//Matrix<3,3> R = identity<3>();
	//R.insert(0,0, V);
	//Matrix<4,1> q = quatFromRot(R);

	//int obj;
	//if (contour) {
	//	CAL_CreateUserDrawn(groupid, drawUnitSphere, NULL, 0.0f, 0.0f, 0.0f, &obj);
	//} else {
	//	CAL_CreateSphere(groupid, 1.0f, 0.0f, 0.0f, 0.0f, &obj);
	//}
	//CAL_SetObjectPosition(obj, (float) x[0], (float) x[1], (float) x[2]);
	//CAL_SetObjectQuaternion(obj, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	//
	//double scale = 3.0;
	//CAL_SetObjectScaling(obj, (float) (scale*sqrt(E(0,0))), (float) (scale*sqrt(E(1,1))), (float) (scale*sqrt(E(2,2))));

#endif
}

void Planner::displayTraces(std::vector<Matrix<X_DIM> > & xTrace)
{
	// Display final path
	CAL_CreateGroup(&cal_path, 0, false, "Path");
	//CAL_SetGroupColor(cal_path, 1.0f, 0.8f, 0.0f);
	CAL_SetGroupColor(cal_path, 0.1f, 0.1f, 0.1f);

	CAL_CreateUserDrawn(cal_path, drawUnitCircle, NULL, 110.0f, 110.0f, 110.0f);
	
	int id;
	float red = (float)random(), green = (float)random(), blue = (float)random();
	
	CAL_CreateSphere(cal_path, 0.05f, (float)xTrace[0][0], (float)xTrace[0][1], 0.0f, &id);
	CAL_SetObjectColor(id, 1.0f, 0.0f, 0.0f);

	int ns = (int)xTrace.size()-1;

	for (int i = 0; i < ns; ++i) 
	{
		Matrix<X_DIM> & pi = xTrace[i];
		
		if (i < ns-1) {
			int np[1] = {2};
			Matrix<X_DIM> & pi1 = xTrace[i+1];
			float side[6] = {(float)pi[0], (float)pi[1], 0.001f, (float)pi1[0], (float)pi1[1], 0.001f};
			//CAL_CreateSphere(cal_path, 0.05f, (float)pi[0], (float)pi[1], 0.01f, &id);
			//CAL_SetObjectColor(id, 1.0f, 0.0f, 0.0f);
			CAL_CreatePolyline(cal_path, 1, np, side, &id);
			CAL_SetObjectColor(id, 1.0f, 0.0f, 0.0f);
		}
	}

	CAL_CreateSphere(cal_path, 0.05f, (float)xTrace[ns-1][0], (float)xTrace[ns-1][1], 0.0f, &id);
	CAL_SetObjectColor(id, 1.0f, 0.0f, 0.0f);
}

void Planner::simulate(const int numParticles, bool display) 
{
	std::vector<Matrix<B_DIM> > Bc(pathlen+1);
	std::vector<Matrix<U_DIM> > Uc(pathlen);

	//nsteps = 10;
	//invnsteps = 0.1;

	std::vector<Matrix<X_DIM> > xTrace(nsteps*pathlen+1);
	
	Matrix<X_DIM> xTrue, xEst, xPrev;
	Matrix<X_DIM, X_DIM> SqrtSigmaEst, SigmaEst, V, D;

	double trajectoryCost, avgCost = 0.0;
	std::vector<double> costs;

	int numCollisions = 0;

	for(int run = 0; run < numParticles; ++run) 
	{
		unVec(B[0], xEst, SqrtSigmaEst);
		
		//x0[0] = -2.0 + random()*4.0;
		//x0[1] = -2.25 + random()*4.5;
		
		SigmaEst = SqrtSigmaEst*~SqrtSigmaEst;
		
		//x0[0] = 2.5; x0[1] = 0.0;
		//xEst = sampleGaussian(x0, SigmaEst);
		//xTrue = x0;
		
		xTrue = sampleGaussian(xEst, 0.2*SigmaEst);
		//std::cout << xEst[0] << " " << xEst[1] << " " << xTrue[0] << " " << xTrue[1] << std::endl;

		xTrace[0] = xTrue;

		trajectoryCost = 0.0;

		bool collision = false;

		for (int t = 0; t < pathlen && !collision; ++t) 
		{
			vec(xEst, SqrtSigmaEst, Bc[t]);
			Uc[t] = U[t] + L[t]*(Bc[t] - B[t]);

			trajectoryCost += cost(Bc[t], Uc[t]);

			xPrev = xTrue;

			// Propagate true state
			for(int i = 1; i <= nsteps; ++i) {
				xTrace[t*nsteps+i] = f(xTrue, Uc[t], 0.75*sampleGaussian(zeros<M_DIM,1>(), identity<M_DIM>()), DT*invnsteps);
				xTrue = xTrace[t*nsteps+i];
			}

			int col = 0;
			CAL_CheckLineCollision(cal_environment, (float)xPrev[0], (float)xPrev[1], 0.0f, (float)xTrue[0], (float)xTrue[1], 0.0f, false, &col);
			if (col != 0) {
				collision = true;
				++numCollisions;
			}
			
			// Kalman filter step 1
			Matrix<X_DIM,X_DIM> A = dfdx(xEst, Uc[t], zeros<M_DIM,1>());
			Matrix<X_DIM,M_DIM> M = dfdm(xEst, Uc[t], zeros<M_DIM,1>());

			SigmaEst = SqrtSigmaEst*~SqrtSigmaEst;

			xEst = f(xEst, Uc[t], zeros<M_DIM,1>());
			SigmaEst = A*SigmaEst*~A + M*~M;

			// Measurement
			Matrix<Z_DIM,1> z = h(xTrue, 0.75*sampleGaussian(zeros<N_DIM,1>(), identity<N_DIM>()));

			// Kalman filter step 2
			Matrix<Z_DIM,X_DIM> H = dhdx(xEst, zeros<N_DIM,1>());
			Matrix<Z_DIM,N_DIM> N = dhdn(xEst, zeros<N_DIM,1>());
			Matrix<X_DIM,Z_DIM> K = SigmaEst*~H/(H*SigmaEst*~H + N*~N);

			xEst += K*(z - h(xEst, zeros<N_DIM,1>()));
			SigmaEst = (identity<X_DIM>() - K*H)*SigmaEst;

			jacobi2(SigmaEst, V, D);
			for (int i = 0; i < X_DIM; ++i) {
				if (D(i,i) > 0) {
					D(i,i) = sqrt(D(i,i));
				} else {
					D(i,i) = 0;
				}
			}
			SqrtSigmaEst = V * D * ~V;
		} // for each time step

		vec(xEst, SqrtSigmaEst, Bc[pathlen]);
		trajectoryCost += finalCost(Bc[pathlen]);

		//avgCost += trajectoryCost;

		if (trajectoryCost < 10000) {
			costs.push_back(trajectoryCost);
		}

		if (display && !collision) {
			displayTraces(xTrace);
		}
	} // for each particle

	std::sort(costs.begin(), costs.end());

	int initialIndex = (int)(0.05*(double)costs.size());
	int finalIndex = (int)(0.95*(double)costs.size());

	std::ofstream fptr("data\\costs.txt", std::ios::out);
	for(int i = initialIndex; i < finalIndex; ++i) {
		avgCost += costs[i];
	    fptr << costs[i] << std::endl;
	}
	fptr.close();

	//std::cout << costs[0] << std::endl;
	//std::cout << costs[finalIndex] << std::endl;
	
	std::cout << "Expected cost (by sampling): " << avgCost/((double)finalIndex-initialIndex) << std::endl;
	
	//resultsfptr << numCollisions << std::endl;
}

void Planner::saveControlPolicy(const std::vector<Matrix<U_DIM, B_DIM> >& L, const std::vector<Matrix<U_DIM> >& l) 
{
	std::ofstream fptr("data\\policy.txt", std::ios::out);
	fptr << pathlen << std::endl;
	for (int t = 0; t < pathlen; ++t) 
	{
		serializeMatrix(fptr, L[t]);
		serializeMatrix(fptr, l[t]);
		if (t < pathlen-1) {
			fptr << std::endl;
		}
	}
	if (fptr.is_open()) { fptr.close(); }
}

void Planner::loadControlPolicy(std::vector<Matrix<B_DIM> > & b, std::vector<Matrix<U_DIM> >& u, std::vector<Matrix<U_DIM, B_DIM> >& L, std::vector<Matrix<U_DIM> >& l)
{
	std::ifstream fptr("data\\policy.txt", std::ios::in);
	fptr >> pathlen;

	L.assign(pathlen, zeros<U_DIM,B_DIM>());
	l.assign(pathlen, zeros<U_DIM, 1>());

	for (int t = 0; t < pathlen; ++t) 
	{
		deserializeMatrix(fptr, L[t]);
		deserializeMatrix(fptr, l[t]);
	}
	if (fptr.is_open()) { fptr.close(); }

	// display mean path
	u.resize(pathlen);
	b.resize(pathlen+1);
}

void Planner::backwardIterationILQG()
{
	// backward iteration: compute new control policy around nominal beliefs and control inputs
	initValue(B[pathlen], S, s, ss);
	gSqr = 0.0;
	Jp = 0.0;
	for (int t = pathlen - 1; t >= 0; --t) {
		controlPolicyILQG(B[t], U[t], S, s, ss, L[t], l[t]);
	}
}

void Planner::backwardIterationDDP()
{
	// backward iteration: compute new control policy around nominal beliefs and control inputs
	initValue(B[pathlen], S, s, ss);
	gSqr = 0.0;
	Jp = 0.0;
	for (int t = pathlen - 1; t >= 0; --t) {
		controlPolicyDDP(B[t], U[t], B[t+1], S, s, ss, L[t], l[t]);
	}
}

void Planner::forwardIterationILQG()
{
	std::vector<Matrix<B_DIM> > Bnxt(pathlen+1);
	std::vector<Matrix<U_DIM> > Unxt(pathlen);

	double curCost = DBL_MAX;
	
	Matrix<X_DIM,X_DIM> SqrtW;
	bool found = false;

	while(eps > 1e-6 && !found) 
	{
		Bnxt[0] = B[0];
		for (int t = 0; t < pathlen; ++t) 
		{
			Unxt[t] = (eps*l[t] + L[t]*(Bnxt[t] - B[t])) + U[t];
			beliefDynamics(Bnxt[t], Unxt[t], Bnxt[t+1], SqrtW); 
		}

		// Compute expected cost
		initValue(Bnxt[pathlen], S, s, ss);
		for (int t = pathlen - 1; t >= 0; --t) {
			expectedCost(Bnxt[t], Unxt[t], S, s, ss, L[t]);
		}
		curCost = ss;

		//std::cout << std::setprecision(12) << "Eps: " << eps << " curCost: " << curCost << " bestCost: " << bestCost << "\n\n";

		/*
		//if ( (curCost - bestCost) < eps*Jp ) {
		if (curCost < bestCost*1.1) {
			bestCost = curCost;
			eps *= 2.0;
			if (eps > 1.0) { eps = 1.0; }
			U = Unxt;
			B = Bnxt;
			found = true;
		} else {
			eps *= 0.5;
		}
		*/

		if ((bestCost < 500 && ((curCost - bestCost) < eps*Jp)) || (bestCost >= 500 && (curCost < bestCost*1.01))) {
			bestCost = curCost;
			eps *= 2.0;
			if (eps > 1.0) { eps = 1.0; }
			U = Unxt;
			B = Bnxt;
			found = true;
		} else {
			eps *= 0.5;
		}
	}
}

void Planner::forwardIterationDDP()
{
	std::vector<Matrix<B_DIM> > Bnxt(pathlen+1);
	std::vector<Matrix<U_DIM> > Unxt(pathlen);

	double curCost = DBL_MAX;
	double eps = 1.0;

	Matrix<X_DIM,X_DIM> SqrtW;
	bool found = false;

	while(eps > 1e-07 && !found) 
	{
		Bnxt[0] = B[0];
		for (int t = 0; t < pathlen; ++t) 
		{
			Unxt[t] = (eps*l[t] + L[t]*(Bnxt[t] - B[t])) + U[t];
			beliefDynamics(Bnxt[t], Unxt[t], Bnxt[t+1], SqrtW); 
		}

		initValue(Bnxt[pathlen], S, s, ss);
		for (int t = pathlen - 1; t >= 0; --t) {
			valueIteration(Bnxt[t], Unxt[t], Bnxt[t+1], S, s, ss, L[t]);
			if (ss > COST_INFTY) {
				break;
			}
		}
		curCost = ss;

		//std::cout << std::setprecision(12) << "Eps: " << eps << " curCost: " << curCost << " bestCost: " << bestCost << std::endl;

		/*
		//if ( (curCost - bestCost) < eps*Jp ) {
		if (curCost < bestCost*1.1) {
			//if (curCost < bestCost) {
			bestCost = curCost;
			eps *= 2.0;
			if (eps > 1.0) { eps = 1.0; }
			U = Unxt;
			B = Bnxt;
			found = true;
		} else {
			eps *= 0.5;
		}
		*/

		if ((bestCost < 20 && ((curCost - bestCost) < eps*Jp)) || (bestCost >= 20 && (curCost < bestCost*1.01))) {
			bestCost = curCost;
			eps *= 2.0;
			if (eps > 1.0) { eps = 1.0; }
			U = Unxt;
			B = Bnxt;
			found = true;
		} else {
			eps *= 0.5;
		}
	}
}

void Planner::solveILQG()
{
#if defined(MAXLIKELIHOOD)
	switchoffML = false;
#else 
	switchoffML = true;
#endif

	//bestCost = COST_INFTY; 
	initValue(B[pathlen], S, s, ss);
	for (int t = pathlen - 1; t >= 0; --t) {
		expectedCost(B[t], U[t], S, s, ss, L[t]);
	}
	bestCost = ss; 

	std::cout << "Initial expected cost: " << bestCost << std::endl;
	//std::cout << std::setprecision(10) << bestCost << std::endl;

	// For stepping through individual ilqr iterations
	int it = 0;
	CAL_SuspendVisualisation();

	clock_t tstart = clock();
	bool terminate = false;
	double prevBestCost;

	eps = 1.0;
	lambda = 1.0;

	while(it < NUM_ITER && !terminate)
	{
		//CAL_SuspendVisualisation();
		backwardIterationILQG();

		prevBestCost = bestCost;

		forwardIterationILQG();

		terminate = (abs(prevBestCost - bestCost) < 1e-07);

		++it; 

		//CAL_ResumeVisualisation();
		//displayMeanPath(B);

		//Sleep(100);
		//int num;
		//std::cin >> num;

		std::cout << "Iteration: " << it << ", Cost: " << std::setprecision(10) << bestCost << " Gradient sqr: " << gSqr << std::endl;
		//std::cout << std::setprecision(10) << bestCost << std::endl;
	}
	std::cout << "Computation time: " << (clock() - tstart) / (double) (CLOCKS_PER_SEC) << std::endl;

	CAL_ResumeVisualisation();
	displayMeanPath(B);

	//std::cout << std::setprecision(6);
	//for(int t = 0; t < pathlen; ++t) {
	//	std::cout << L[t] << std::endl;
	//}

	//saveControlPolicy(L, l);
} 

void Planner::solveDDP()
{
#if defined(MAXLIKELIHOOD)
	switchoffML = false;
#else 
	switchoffML = true;
#endif

	initValue(B[pathlen], S, s, ss);
	for (int t = pathlen - 1; t >= 0; --t) {
		valueIteration(B[t], U[t], B[t+1], S, s, ss, L[t]);
		if (ss > COST_INFTY) {
			break;
		}
	}
	bestCost = ss;
	
	std::cout << "Initial expected cost: " << bestCost << std::endl;
	//std::cout << std::setprecision(10) << bestCost << std::endl;

	//resultsfptr << bestCost << " ";

	// For stepping through individual ilqr iterations
	int it = 0;
	CAL_SuspendVisualisation();

	clock_t tstart = clock();
	bool terminate = false;
	double prevBestCost;

	eps = 1.0;
	lambda = 1.0;

	while(it < NUM_ITER && !terminate)
	{
		//CAL_SuspendVisualisation();
		gSqr = 0.0;
		Jp = 0.0;
 
		backwardIterationDDP();

		prevBestCost = bestCost;

		forwardIterationDDP();

		terminate = (abs(prevBestCost - bestCost) < 1e-07);
		//terminate = (gSqr < 0.1);

		++it; 

		//CAL_ResumeVisualisation();
		//displayMeanPath(B);

		//Sleep(100);

		std::cout << "Iteration: " << it << ", Cost: " << std::setprecision(10) << bestCost << " Gradient sqr: " << gSqr << std::endl;
		//std::cout << std::setprecision(10) << bestCost << std::endl;
	}
	//std::cout << "Computation time: " << (clock() - tstart) / (double) (CLOCKS_PER_SEC) << std::endl;

	CAL_ResumeVisualisation();
	displayMeanPath(B);

	switchoffML = true;

	initValue(B[pathlen], S, s, ss);
	for (int t = pathlen - 1; t >= 0; --t) {
		valueIteration(B[t], U[t], B[t+1], S, s, ss, L[t]);
		if (ss > COST_INFTY) {
			break;
		}
	}
	std::cout << "Cost: " << bestCost << " Final expected cost: " << ss << std::endl;

	//resultsfptr << bestCost << " " << ss << " ";

	//std::cout << std::setprecision(6);
	//for(int t = 0; t < pathlen; ++t) {
	//	std::cout << L[t] << std::endl;
	//}

	//std::string filename = "data\\Path1-ML.txt";

	//std::ofstream fptr(filename.c_str(),std::ios::out);
	//fptr << U.size()  << std::endl;
	//for(int i = 0; i < (int)U.size(); ++i) {
	//	fptr << ~U[i];
	//}
	//fptr.close();

	//std::cout << "Finished writing to file" << std::endl;
}