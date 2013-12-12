#include "userinput.h"

// Switch between belief vector and matrices
inline void unVec(const Matrix<B_DIM>& b, Matrix<X_DIM>& x, Matrix<X_DIM,X_DIM>& SqrtSigma) {
  x = b.subMatrix<X_DIM,1>(0,0);
  size_t idx = X_DIM;
  for (size_t j = 0; j < X_DIM; ++j) {
    for (size_t i = j; i < X_DIM; ++i) {
      SqrtSigma(i,j) = SqrtSigma(j,i) = b[idx];
      ++idx;
    }
  }
}

inline void vec(const Matrix<X_DIM>& x, const Matrix<X_DIM,X_DIM>& SqrtSigma, Matrix<B_DIM>& b) {
  b.insert(0,0,x);
  size_t idx = X_DIM;
  for (size_t j = 0; j < X_DIM; ++j) {
    for (size_t i = j; i < X_DIM; ++i) {
      b[idx] = 0.5 * (SqrtSigma(i,j) + SqrtSigma(j,i));
      ++idx;
    }
  }
}

// final cost function
double finalCost(const Matrix<B_DIM>& b) {
  Matrix<X_DIM> x;
  Matrix<X_DIM,X_DIM> SqrtSigma;

  unVec(b, x, SqrtSigma);

  return 0.5*tr(~(x - xGoal)*Qgoal*(x - xGoal)) + 0.5*tr(~SqrtSigma*Qgoal*SqrtSigma);
}

// cost function
double cost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u) {
  Matrix<X_DIM> x;
  Matrix<X_DIM,X_DIM> SqrtSigma;

  unVec(b, x, SqrtSigma);

  return  0.5*tr(~u*Rint*u) + 0.5*tr(~SqrtSigma*Qint*SqrtSigma);
}

// Compute second-order expansion of finalCost around b
inline void initValue(const Matrix<B_DIM>& b, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss) {
	ss = finalCost(b);

  Matrix<B_DIM> br(b), bl(b);
  for (size_t i = 0; i < B_DIM; ++i) {
		br[i] += H1; bl[i] -= H1; 
    s[i] = (finalCost(br) - finalCost(bl)) / (2.0*H1);
    S(i,i) = ( finalCost(bl) - 2.0*ss + finalCost(br) ) / (H1*H1);
    br[i] = bl[i] = b[i];
	}

  Matrix<B_DIM> btr(b), btl(b), bbr(b), bbl(b);
	for (size_t i = 1; i < B_DIM; ++i) {
    btr[i] += H1; btl[i] -= H1; bbr[i] += H1; bbl[i] -= H1;
		for (size_t j = 0; j < i; ++j) {
			btr[j] += H1; btl[j] += H1; bbr[j] -= H1; bbl[j] -= H1;
      S(i,j) = S(j,i) = (finalCost(bbl) + finalCost(btr) - finalCost(btl) - finalCost(bbr)) / (4.0*H1*H1);
      btr[j] = btl[j] = bbr[j] = bbl[j] = b[j];
		}
    btr[i] = btl[i] = bbr[i] = bbl[i] = b[i];
	}
}

// Compute second-order expansion of cost around b and u
inline void quadratizeCost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& Q, Matrix<B_DIM>& q, Matrix<U_DIM,B_DIM>& P, Matrix<U_DIM,U_DIM>& R, Matrix<U_DIM>& r, double& p) {
	// p
  p = cost(b,u);

  // q, diag(Q)
  Matrix<B_DIM> br(b), bl(b);
  for (size_t i = 0; i < B_DIM; ++i) {
		br[i] += H1; bl[i] -= H1; 
    q[i] = (cost(br,u) - cost(bl,u)) / (2.0*H1);
    Q(i,i) = ( cost(bl,u) - 2.0*p + cost(br, u) ) / (H1*H1);
    br[i] = bl[i] = b[i];
	}

  // r, diag(R), P
  Matrix<U_DIM> ur(u), ul(u);
  for (size_t i = 0; i < U_DIM; ++i) {
		ur[i] += H1; ul[i] -= H1; 
    r[i] = (cost(b,ur) - cost(b,ul)) / (2.0*H1);
    R(i,i) = ( cost(b,ul) - 2.0*p + cost(b, ur) ) / (H1*H1);
    for (size_t j = 0; j < B_DIM; ++j) {
			br[j] += H1; bl[j] -= H1;
      P(i,j) = (cost(bl, ul) + cost(br, ur) - cost(br, ul) - cost(bl, ur)) / (4.0*H1*H1);
      br[j] = bl[j] = b[j];
		}
    ur[i] = ul[i] = u[i];
	}

  // Q
  Matrix<B_DIM> btr(b), btl(b), bbr(b), bbl(b);
	for (size_t i = 1; i < B_DIM; ++i) {
    btr[i] += H1; btl[i] -= H1; bbr[i] += H1; bbl[i] -= H1;
		for (size_t j = 0; j < i; ++j) {
			btr[j] += H1; btl[j] += H1; bbr[j] -= H1; bbl[j] -= H1;
      Q(i,j) = Q(j,i) = (cost(bbl,u) + cost(btr,u) - cost(btl,u) - cost(bbr,u)) / (4.0*H1*H1);
      btr[j] = btl[j] = bbr[j] = bbl[j] = b[j];
		}
    btr[i] = btl[i] = bbr[i] = bbl[i] = b[i];
	}

  // R
  Matrix<U_DIM> utr(u), utl(u), ubr(u), ubl(u);
  for (size_t i = 1; i < U_DIM; ++i) {
    utr[i] += H1; utl[i] -= H1; ubr[i] += H1; ubl[i] -= H1;
		for (size_t j = 0; j < i; ++j) {
			utr[j] += H1; utl[j] += H1; ubr[j] -= H1; ubl[j] -= H1;
      R(i,j) = R(j,i) = (cost(b,ubl) + cost(b,utr) - cost(b,utl) - cost(b,ubr)) / (4.0*H1*H1);
      utr[j] = utl[j] = ubr[j] = ubl[j] = u[j];
		}
    utr[i] = utl[i] = ubr[i] = ubl[i] = u[i];
	}
}

// Jacobian df/dx(x,u,m)
inline Matrix<X_DIM,X_DIM> dfdx(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m) {
  Matrix<X_DIM,X_DIM> A;
  Matrix<X_DIM> xr(x), xl(x);
  for (size_t i = 0; i < X_DIM; ++i) {
    xr[i] += H1; xl[i] -= H1;
    A.insert(0,i, (f(xr, u, m) - f(xl, u, m)) / (2.0*H1));
    xr[i] = xl[i] = x[i];
  }
  return A;
}

// Jacobian df/du(x,u,m)
inline Matrix<X_DIM,U_DIM> dfdu(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m) {
  Matrix<X_DIM,U_DIM> B;
  Matrix<U_DIM> ur(u), ul(u);
  for (size_t i = 0; i < U_DIM; ++i) {
    ur[i] += H1; ul[i] -= H1;
    B.insert(0,i, (f(x, ur, m) - f(x, ul, m)) / (2.0*H1));
    ur[i] = ul[i] = u[i];
  }
  return B;
}

// Jacobian df/dm(x,u,m)
inline Matrix<X_DIM,M_DIM> dfdm(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m) {
  Matrix<X_DIM,M_DIM> M;
  Matrix<M_DIM> mr(m), ml(m);
  for (size_t i = 0; i < M_DIM; ++i) {
    mr[i] += H1; ml[i] -= H1;
    M.insert(0,i, (f(x, u, mr) - f(x, u, ml)) / (2.0*H1));
    mr[i] = ml[i] = m[i];
  }
  return M;
}

// Jacobian dh/dx(x,n)
inline Matrix<Z_DIM,X_DIM> dhdx(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) {
  Matrix<Z_DIM,X_DIM> H;
  Matrix<X_DIM> xr(x), xl(x);
  for (size_t i = 0; i < X_DIM; ++i) {
    xr[i] += H1; xl[i] -= H1;
    H.insert(0,i, (h(xr, n) - h(xl, n)) / (2.0*H1));
    xr[i] = xl[i] = x[i];
  }
  return H;
}

// Jacobian dh/dn(x,n)
inline Matrix<Z_DIM,N_DIM> dhdn(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) {
  Matrix<Z_DIM,N_DIM> N;
  Matrix<N_DIM> nr(n), nl(n);
  for (size_t i = 0; i < N_DIM; ++i) {
    nr[i] += H1; nl[i] -= H1;
    N.insert(0,i, (h(x, nr) - h(x, nl)) / (2.0*H1));
    nr[i] = nl[i] = n[i];
  }
  return N;
}


// Belief dynamics
inline void beliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM>& g, Matrix<X_DIM,X_DIM>& SqrtW) {
  Matrix<X_DIM> x;
  Matrix<X_DIM,X_DIM> SqrtSigma;
  unVec(b, x, SqrtSigma);

  Matrix<X_DIM,X_DIM> Sigma = SqrtSigma*~SqrtSigma;
  
  Matrix<X_DIM,X_DIM> A = dfdx(x, u, zeros<M_DIM,1>());
  Matrix<X_DIM,M_DIM> M = dfdm(x, u, zeros<M_DIM,1>());
  
  x = f(x, u, zeros<M_DIM,1>());
  Sigma = A*Sigma*~A + M*~M;
  
  Matrix<Z_DIM,X_DIM> H = dhdx(x, zeros<N_DIM,1>());
  Matrix<Z_DIM,N_DIM> N = dhdn(x, zeros<N_DIM,1>());
  Matrix<X_DIM,Z_DIM> K = Sigma*~H/(H*Sigma*~H + N*~N);

  Matrix<X_DIM, X_DIM> W = K*H*Sigma;
  Sigma -= W;


  Matrix<X_DIM, X_DIM> V, D;

  jacobi(Sigma, V, D);
  for (size_t i = 0; i < X_DIM; ++i) {
    if (D(i,i) > 0) {
      D(i,i) = sqrt(D(i,i));
    } else {
      D(i,i) = 0;
    }
  }
  SqrtSigma = V * D * ~V;
  vec(x, SqrtSigma, g);
  
  jacobi(W, V, D);
  for (size_t i = 0; i < X_DIM; ++i) {
    if (D(i,i) > 0) {
      D(i,i) = sqrt(D(i,i));
    } else {
      D(i,i) = 0;
    }
  }
  SqrtW = V * D * ~V;
}

inline void linearizeBeliefDynamics(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM>& c, Matrix<B_DIM,B_DIM>& A, Matrix<B_DIM, U_DIM>& B, std::vector<Matrix<X_DIM> >& cn, std::vector<Matrix<X_DIM, B_DIM> >& An, std::vector<Matrix<X_DIM, U_DIM> >& Bn) {
  cn.resize(X_DIM);
  An.resize(X_DIM);
  Bn.resize(X_DIM);

  Matrix<X_DIM, X_DIM> sqrtW;
  // set c
  beliefDynamics(b, u, c, sqrtW);
  // set cn's
  for (size_t i = 0; i < X_DIM; ++i) {
    cn[i] = sqrtW.subMatrix<X_DIM, 1>(0,i);
  }

  // Compute Jacobian's
  Matrix<B_DIM> br(b), bl(b);
  for (size_t i = 0; i < B_DIM; ++i) {
    br[i] += H1; bl[i] -= H1;

    Matrix<B_DIM> cr, cl;
    Matrix<X_DIM,X_DIM> sqrtWr, sqrtWl;

    beliefDynamics(br, u, cr, sqrtWr);
    beliefDynamics(bl, u, cl, sqrtWl);

    A.insert(0,i, (cr - cl) / (2.0*H1));
    
    for(size_t j = 0; j < X_DIM; ++j) {
      An[j].insert(0,i, (sqrtWr.subMatrix<X_DIM,1>(0,j) - sqrtWl.subMatrix<X_DIM,1>(0,j)) / (2.0*H1));
    }

    br[i] = bl[i] = b[i];
  }

  Matrix<U_DIM> ur(u), ul(u);
  for (size_t i = 0; i < U_DIM; ++i) {
    ur[i] += H1; ul[i] -= H1;

    Matrix<B_DIM> cr, cl;
    Matrix<X_DIM,X_DIM> sqrtWr, sqrtWl;

    beliefDynamics(b, ur, cr, sqrtWr);
    beliefDynamics(b, ul, cl, sqrtWl);

    B.insert(0,i, (cr - cl) / (2.0*H1));
    
    for(size_t j = 0; j < X_DIM; ++j) {
      Bn[j].insert(0,i, (sqrtWr.subMatrix<X_DIM,1>(0,j) - sqrtWl.subMatrix<X_DIM,1>(0,j)) / (2.0*H1));
    }

    ur[i] = ul[i] = u[i];
  }
}

inline void controlPolicy(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, Matrix<U_DIM,B_DIM>& L, Matrix<U_DIM>& l) {
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

  Matrix<U_DIM,U_DIM> D = R + ~B*S*B;
  Matrix<U_DIM,B_DIM> E = P + ~B*S*A;
  Matrix<U_DIM> d = r + ~B*s;
  Matrix<B_DIM,B_DIM> C = Q + ~A*S*A;
  Matrix<B_DIM> c = q + ~A*s;
  double e = p + ss;

  for (size_t i = 0; i < X_DIM; ++i) {
    D += ~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * Bn[i];
    E += ~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * An[i];
    d += ~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * cn[i];
    C += ~An[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * An[i];
    c += ~An[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * cn[i];
    e += 0.5*tr(~cn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * cn[i]);
  }
  
  Matrix<U_DIM,U_DIM> invD = !D;

  // control policy du = L dx + l
  L = -invD*E;
  l = -invD*d;

  // update value function
  S = C - ~E*invD*E;
  s = c - ~E*invD*d;
  ss = e - 0.5*tr(~d*invD*d);
}

inline void expectedCost(const Matrix<B_DIM>& b, const Matrix<U_DIM>& u, Matrix<B_DIM,B_DIM>& S, Matrix<B_DIM>& s, double& ss, const Matrix<U_DIM,B_DIM>& L) {
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

  Matrix<B_DIM,B_DIM> C = Q + ~A*S*A + ~L*(R + ~B*S*B)*L + ~L*(P + ~B*S*A) + ~(P + ~B*S*A)*L;
  Matrix<B_DIM> c = q + ~A*s + ~L*(r + ~B*s); //+ ~(P + ~B*S*A)*l + ~L*(R + ~B*S*B)*l;
  double e = p + ss; // + tr(~l*(r + ~B*s)) + 0.5*tr(~l*(R + ~B*S*B)*l);

  for (size_t i = 0; i < X_DIM; ++i) {
    C += ~An[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * An[i] + 
         ~L*(~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * Bn[i])*L + 
         ~L*(~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * An[i]) + 
         ~(~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * An[i])*L;
    c += ~An[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * cn[i] 
       + ~L*(~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * cn[i]); 
    // + ~(~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * An[i])*l 
    // + ~L*(~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * Bn[i])*l;
    e += 0.5*tr(~cn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * cn[i]);
    // + tr(~l*(~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * cn[i])) 
    // + 0.5*tr(~l*(~Bn[i] * S.subMatrix<X_DIM,X_DIM>(0,0) * Bn[i])*l);
  }
  
  // update value function
  S = C;
  s = c;
  ss = e;
}

int _tmain(int argc, _TCHAR* argv[])
{
  userInit();

  Matrix<B_DIM,B_DIM> S;
  Matrix<B_DIM> s;
  double ss;
  
  Matrix<X_DIM,X_DIM> SqrtW;

  std::vector<Matrix<B_DIM> > b(LENGTH+1), bCopy(LENGTH+1);
  std::vector<Matrix<U_DIM> > u(LENGTH), uCopy(LENGTH);
  
  // initialize trajectory
  vec(x0, SqrtSigma0, b[0]);
  for (size_t t = 0; t < LENGTH; ++t) {
    u[t] = zeros<U_DIM, 1>();
    beliefDynamics(b[t], u[t], b[t+1], SqrtW);
  }

  std::vector<Matrix<U_DIM, B_DIM> > L(LENGTH, zeros<U_DIM,B_DIM>());
  std::vector<Matrix<U_DIM> > l(LENGTH, zeros<U_DIM, 1>());

  double bestCost = -log(0.0);
  double eps = 1.0;

  int i = 0;
  while (eps > 1.0e-05) {
    // backward iteration: compute new control policy around nominal belief and control input
    initValue(b[LENGTH], S, s, ss);
    for (int t = LENGTH - 1; t >= 0; --t) {
      controlPolicy(b[t], u[t], S, s, ss, L[t], l[t]);
    }

    double estCost = ss;

    bCopy = b;
    uCopy = u;

    // forward iteration: compute nominal belief and control input with current control policy
    while (eps > 1.0e-05) {
      Matrix<B_DIM> b_curr = b[0];
      for (size_t t = 0; t < LENGTH; ++t) {
        Matrix<U_DIM> u_dev = L[t]*(b_curr - b[t]) + eps*l[t];

        b[t] = b_curr;
        u[t] = u[t] + u_dev;

        beliefDynamics(b[t], u[t], b_curr, SqrtW);
      }
      b[LENGTH] = b_curr;

      // Compute expected cost
      initValue(b[LENGTH], S, s, ss);
      for (int t = LENGTH - 1; t >= 0; --t) {
        expectedCost(b[t], u[t], S, s, ss, L[t]);
      }

      std::cout << "Iteration: " << i << ", LQG cost: " << estCost << ", Expected cost: " << ss << ", Eps: " << eps << std::endl;

      // Accept new solution?
      if (ss < bestCost) {
        eps = 1.0;
        bestCost = ss;
        ++i;
        break;
      } else {
        eps *= 0.5;
        b = bCopy;
        u = uCopy;
      }
    }
  }

  for (size_t i = 0; i <= LENGTH; ++i) {
    std::cout << ~b[i].subMatrix<2,1>(0,0);
  }

  int k;
  std::cin >> k;
	return 0;
}

