#ifndef __UTILS_H__
#define __UTILS_H__

#include <math.h>

inline double random() {
	unsigned int i;
	rand_s(&i);
	return ((double) i) / ((double)UINT_MAX + 1.0);
}

inline double normal() {
	double u_1 = 0;
	while (u_1 == 0) {
		u_1 = random();
	}
	double u_2 = 0;
	while (u_2 == 0) {
		u_2 = random();
	}
	return sqrt(-2*log(u_1)) * sin(2*M_PI*u_2);
}

inline unsigned int randomInt(unsigned int upper_bound) {
	unsigned int i;
	rand_s(&i);
	return i % upper_bound;
}

template <size_t size>
inline Matrix<size> sampleGaussian(const Matrix<size>& mean, const SymmetricMatrix<size>& var) {
	Matrix<size> sample;
	for (int j = 0; j < size; ++j) {
		sample[j] = normal();
	}
	Matrix<size, size> SVec; 
	SymmetricMatrix<size> SVal;
	jacobi(var, SVec, SVal);
	for (int i = 0; i < size; ++i) {
		SVal(i,i) = sqrt(SVal(i,i));
	}
	return SVec * SVal * sample + mean;
}

template <size_t size>
inline Matrix<size> sampleUniform(const Matrix<size>& mean, const Matrix<size, size>& var) {
	Matrix<size> sample;
	for (int j = 0; j < size; ++j) {
		sample[j] = 2*(random()-0.5) * sqrt(3.0 * var(j,j));
	}
	return sample + mean;
}

inline Matrix<3,3> cpMatrix(const Matrix<3,1>& a) 
{
	Matrix<3,3> A;
	A(0,0) = 0;       A(0,1) = -a(2,0); A(0,2) = a(1,0);
	A(1,0) = a(2,0);  A(1,1) = 0;       A(1,2) = -a(0,0);
	A(2,0) = -a(1,0); A(2,1) = a(0,0);  A(2,2) = 0;

	return A;
}

inline Matrix<3,3> rotFromErr(const Matrix<3,1>& q) {
	double rr = q(0,0)*q(0,0)+q(1,0)*q(1,0)+q(2,0)*q(2,0);
	if (rr == 0) {
		return identity<3>();
	} else {
		double r = sqrt(rr);
		//printf_s("%24.24g ", cos(r) - 1);
		return cpMatrix(q * (sin(r) / r)) + identity<3>() * cos(r) + (q*~q) * ((1 - cos(r)) / rr);
		//return cpMatrix(q * (sin(r) / r)) + identity<3>() * cos(r) + (q*~q)/rr - (q*~q)*(cos(r)/rr);
	}
}

inline Matrix<3,1> errFromRot(const Matrix<3,3>& R) {
	Matrix<3,1> q;
	q(0,0) = R(2,1) - R(1,2);
	q(1,0) = R(0,2) - R(2,0);
	q(2,0) = R(1,0) - R(0,1);

	double r = sqrt(q(0,0)*q(0,0)+q(1,0)*q(1,0)+q(2,0)*q(2,0));
	double t = R(0,0) + R(1,1) + R(2,2) - 1;

	if (r == 0) {
		return zeros<3,1>();
	} else {
		return q * (atan2(r, t) / r);
	}
}

inline Matrix<4,4> transFromErr(const Matrix<6,1>& x) {
	Matrix<4,4> X;
	X = identity<4>();
	X.insert(0,0, rotFromErr(x.subMatrix<3,1>(3,0)));
	X.insert(0,3, x.subMatrix<3,1>(0,0));
	return X;
}

inline Matrix<4,1> randQuat() {
	double s = random();
	double s1 = sqrt(1 - s);
	double s2 = sqrt(s);
	double t1 = 2 * M_PI * random();
	double t2 = 2 * M_PI * random();

	Matrix<4,1> q;

	q(0,0) = sin(t1) * s1; // x
	q(1,0) = cos(t1) * s1; // y
	q(2,0) = sin(t2) * s2; // z
	q(3,0) = cos(t2) * s2; // w

	return q;
}

inline Matrix<3,3> rotFromQuat(const Matrix<4,1>& q) {
	double x = q(0,0);
	double y = q(1,0);
	double z = q(2,0);
	double w = q(3,0);
	Matrix<3,3> R;
	R(0,0) =  1 - 2*y*y - 2*z*z; R(0,1) = 2*x*y - 2*z*w;     R(0,2) = 2*x*z + 2*y*w;
	R(1,0) =  2*x*y + 2*z*w;     R(1,1) = 1 - 2*x*x - 2*z*z; R(1,2) = 2*y*z - 2*x*w;
	R(2,0) =  2*x*z - 2*y*w;     R(2,1) = 2*y*z + 2*x*w;     R(2,2) = 1 - 2*x*x - 2*y*y;
	return R;
}

inline Matrix<4,1> quatFromRot(const Matrix<3,3>& R) {
	double x = R(2,1) - R(1,2);
	double y = R(0,2) - R(2,0);
	double z = R(1,0) - R(0,1);
	double r = sqrt(x*x+y*y+z*z);
	double t = R(0,0) + R(1,1) + R(2,2);
	double angle = atan2(r,t-1);
	if (angle != 0) {
		x /= r;
		y /= r;
		z /= r;
	} else {
		x = 0;
		y = 0;
		z = 0;
	}
	Matrix<4,1> q;
	q(0,0) = sin(angle/2)*x;
	q(1,0) = sin(angle/2)*y;
	q(2,0) = sin(angle/2)*z;
	q(3,0) = cos(angle/2);

	return q;
}

#include "glut.h"
const double DEG2RAD = M_PI/180;

inline double mod2pi(double x) {
	// Returns a value 0 <= x < 2*PI
	double result = fmod(x, 2*M_PI);
	if (result < 0) {
		result += (2*M_PI);
	}
	return result;
}


inline void linewidth(void *userdef, float time){
	glLineWidth(4.0f);
}


inline void drawUnitCircle(void *userdef, float time) {
	glLineWidth(4.0f);
	glDisable(GL_LIGHTING);
	glColor3f(0.0f, 0.0f, 0.0f);
	glBegin(GL_LINE_LOOP);
	//glColor3f(0.5f, 0.0f, 1.0f);

	for (int i = 0; i < 72; i++)
	{
		float theta = (float)(i*10.0*DEG2RAD);
		glVertex3f(cosf(theta), sinf(theta), 0.0f);
	}
	glEnd();
	glEnable(GL_LIGHTING);
	glLineWidth(1.0f);
}


/*inline void drawEllipse2d(const Matrix<2>& x, const Matrix<2, 2>& S, int groupid, bool contour = false)
{
	Matrix<2, 2> V, E;
	jacobi(S, V, E);
	Matrix<3,3> R = identity<3>();
	R.insert(0,0, V);
	Matrix<3,3> T = identity<3>();
	if (!contour) {
		T(1,1) = T(2,2) = 0.0;
		T(1,2) = -1.0; T(2,1) = 1.0;
	}
	Matrix<4,1> q = quatFromRot(R*T);

	int obj;
	if (contour) {
		CAL_CreateUserDrawn(groupid, drawUnitCircle, NULL, 0.0f, 0.0f, 0.0f, &obj);
		CAL_SetObjectPosition(obj, (float) x[0], (float) x[1], 0.05f);
	} else {
		CAL_CreateCylinder(groupid, 1.0f, 0.02f, 0.0f, 0.0f, 0.0f, &obj);
		CAL_SetObjectPosition(obj, (float) x[0], (float) x[1], 0.01f);
	}
	CAL_SetObjectQuaternion(obj, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	if (contour) {
		CAL_SetObjectScaling(obj, (float) (3*sqrt(E(0,0))), (float) (3*sqrt(E(1,1))), 1.0f);
	} else {
		CAL_SetObjectScaling(obj, (float) (3*sqrt(E(0,0))), 1.0f, (float) (3*sqrt(E(1,1))));
	}
}*/

/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to pseudocode.

See subroutines comments for additional copyrights.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

/*************************************************************************
Natural logarithm of gamma function

Input parameters:
    X       -   argument

Result:
    logarithm of the absolute value of the Gamma(X).

Output parameters:
    SgnGam  -   sign(Gamma(X))

Domain:
    0 < X < 2.55e305
    -2.55e305 < X < 0, X is not an integer.

ACCURACY:
arithmetic      domain        # trials     peak         rms
   IEEE    0, 3                 28000     5.4e-16     1.1e-16
   IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
The error criterion was relative when the function magnitude
was greater than one but absolute when it was less than one.

The following test used the relative error criterion, though
at certain points the relative error could be much higher than
indicated.
   IEEE    -200, -4             10000     4.8e-16     1.3e-16
*************************************************************************/
inline double lngamma(double x, double& sgngam)
{
    double result;
    double a;
    double b;
    double c;
    double p;
    double q;
    double u;
    double w;
    double z;
    int i;
    double logpi;
    double ls2pi;
    double tmp;

    sgngam = 1;
    logpi = 1.14472988584940017414;
    ls2pi = 0.91893853320467274178;
    if( x < -34.0 )
    {
        q = -x;
        w = lngamma(q, tmp);
        p = floor(q);
        i = (int) floor(p+0.5);
        if( i%2==0 )
        {
            sgngam = -1;
        }
        else
        {
            sgngam = 1;
        }
        z = q-p;
        if( z > 0.5 )
        {
            p = p+1;
            z = p-q;
        }
        z = q*sin(M_PI*z);
        result = logpi-log(z)-w;
        return result;
    }
    if( x < 13 )
    {
        z = 1;
        p = 0;
        u = x;
        while( u >= 3 )
        {
            p = p-1;
            u = x+p;
            z = z*u;
        }
        while( u < 2 )
        {
            z = z/u;
            p = p+1;
            u = x+p;
        }
        if( z < 0 )
        {
            sgngam = -1;
            z = -z;
        }
        else
        {
            sgngam = 1;
        }
        if( u == 2 )
        {
            result = log(z);
            return result;
        }
        p = p-2;
        x = x+p;
        b = -1378.25152569120859100;
        b = -38801.6315134637840924+x*b;
        b = -331612.992738871184744+x*b;
        b = -1162370.97492762307383+x*b;
        b = -1721737.00820839662146+x*b;
        b = -853555.664245765465627+x*b;
        c = 1;
        c = -351.815701436523470549+x*c;
        c = -17064.2106651881159223+x*c;
        c = -220528.590553854454839+x*c;
        c = -1139334.44367982507207+x*c;
        c = -2532523.07177582951285+x*c;
        c = -2018891.41433532773231+x*c;
        p = x*b/c;
        result = log(z)+p;
        return result;
    }
    q = (x-0.5)*log(x)-x+ls2pi;
    if( x > 100000000 )
    {
        result = q;
        return result;
    }
    p = 1/(x*x);
    if( x >= 1000.0 )
    {
        q = q+((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
    }
    else
    {
        a = 8.11614167470508450300*0.0001;
        a = -5.95061904284301438324*0.0001+p*a;
        a = 7.93650340457716943945*0.0001+p*a;
        a = -2.77777777730099687205*0.001+p*a;
        a = 8.33333333333331927722*0.01+p*a;
        q = q+a/x;
    }
    result = q;
    return result;
}

/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15
*************************************************************************/
// forward declarations
inline double incompletegamma(double a, double x);
inline double incompletegammac(double a, double x);

inline double incompletegammac(double a, double x)
{
    double result;
    double igammaepsilon;
    double igammabignumber;
    double igammabignumberinv;
    double ans;
    double ax;
    double c;
    double yc;
    double r;
    double t;
    double y;
    double z;
    double pk;
    double pkm1;
    double pkm2;
    double qk;
    double qkm1;
    double qkm2;
    double tmp;

    igammaepsilon = 0.000000000000001;
    igammabignumber = 4503599627370496.0;
    igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
    if( x <= 0 || a <= 0 )
    {
        result = 1;
        return result;
    }
    if( x < 1 || x < a )
    {
        result = 1-incompletegamma(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ax < -709.78271289338399 )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    y = 1-a;
    z = x+y+1;
    c = 0;
    pkm2 = 1;
    qkm2 = x;
    pkm1 = x+1;
    qkm1 = z*x;
    ans = pkm1/qkm1;
    do
    {
        c = c+1;
        y = y+1;
        z = z+2;
        yc = y*c;
        pk = pkm1*z-pkm2*yc;
        qk = qkm1*z-qkm2*yc;
        if( qk != 0 )
        {
            r = pk/qk;
            t = fabs((ans-r)/r);
            ans = r;
        }
        else
        {
            t = 1;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( fabs(pk) > igammabignumber )
        {
            pkm2 = pkm2*igammabignumberinv;
            pkm1 = pkm1*igammabignumberinv;
            qkm2 = qkm2*igammabignumberinv;
            qkm1 = qkm1*igammabignumberinv;
        }
    }
    while( t > igammaepsilon );
    result = ans*ax;
    return result;
}

/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14
*************************************************************************/
inline double incompletegamma(double a, double x)
{
	//std::cout << "a: " << a << " x: " << x << std::endl;

    double result;
    double igammaepsilon;
    double ans;
    double ax;
    double c;
    double r;
    double tmp;

    igammaepsilon = 0.000000000000001;
    if( x <= 0 || a <= 0 )
    {
        result = 0;
        return result;
    }
    if( x > 1 && x > a )
    {
        result = 1-incompletegammac(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ax < -709.78271289338399 )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    r = a;
    c = 1;
    ans = 1;
    do
    {
        r = r+1;
        c = c*x/r;
        ans = ans+c;
    }
    while( c/ans > igammaepsilon );
    result = ans*ax/a;
    return result;
}



inline double pdf(double x) {
	return exp(-0.5*x*x) * (1.0/sqrt(M_PI*2.0));
}

inline double cdf(double x) {
	if (x < 0) {
		return 0.5 - 0.5*incompletegamma(0.5, 0.5*x*x);
	} else {
		return 0.5 + 0.5*incompletegamma(0.5, 0.5*x*x);
	}

	//return 0.5*(1.0 + erf(x*M_SQRT1_2));
}

inline void truncate(double x, double mean, double var, double& newMean, double& newVar) {
	double stddev = sqrt(var);
	double y = (x - mean) / stddev;
	double z = pdf(y) / cdf(y);
	//std::cout << "y: " << y << " z: " << z << " (y+z): " << (y+z) << " pdf(y): " << pdf(y) << " cdf(y): " << cdf(y) << std::endl;
	newMean = mean - z*stddev;
	//newVar = var*(1.0 - y*z - z*z);
	newVar = var*(1.0 - z*(y + z));
	//std::cout << "newVar: " << newVar << std::endl;
}

inline double Constraint(double theta)
{
    if(theta >= 0 && theta <= 2 * 3.1416)
		;

	else if(theta < 0){
		while(theta < 0){
			theta += 2 * 3.1416;
		}
	}
	else if(theta > 2*3.1416){
		while(theta > 2 * 3.1416)
			theta -= 2 * 3.1416;
	}

	return theta;
}


inline void setTwoLinks(const int& cal_link1, const int& cal_link2, const Matrix<4>& state)
{

	Matrix<2> poslink1 = zeros<2,1>();
	Matrix<2> poslink2 = zeros<2,1>();

	poslink1[0] = cos(state[0]) * 0.3/2.0; poslink1[1] = sin(state[0]) * 0.3 / 2.0;
	
	poslink2[0] = cos(state[0]) * 0.3 + cos(state[0] + state[1]) * 0.33/2.0;
	poslink2[1] = sin(state[0]) * 0.3 + sin(state[0] + state[1]) * 0.33/2.0;

	CAL_SetGroupOrientation(cal_link1, 0.0, 0.0, Constraint(state[0]));
	CAL_SetGroupPosition(cal_link1, poslink1[0], poslink1[1], 0.0);
	CAL_SetGroupPosition(cal_link2, poslink2[0], poslink2[1], 0.0);
	CAL_SetGroupOrientation(cal_link2, 0.0, 0.0, Constraint(state[1] + state[0]));


	return;
}

template <size_t xDim>
void regularize(SymmetricMatrix<xDim>& Q) {
	SymmetricMatrix<xDim> D;
	Matrix<xDim,xDim> V;
	jacobi(Q, V, D);
	for (size_t i = 0; i < xDim; ++i) {
		if (D(i,i) < 0) {
			D(i,i) = 0;
		}
	}
	Q = SymProd(V,D*~V);
}

template <size_t xDim>
void regularize(SymmetricMatrix<xDim>& Q, const double& lambda)
{
	SymmetricMatrix<xDim> D;
	Matrix<xDim,xDim> V;
	jacobi(Q, V, D);
	for (size_t i = 0; i < xDim; ++i) {
		if (D(i,i) < 1e-10) {
			D(i,i) = 1e-10;
		}
		D(i,i) += lambda; //add a postive number. if this posivtive approachs to infinity, then Q^-1 will approaches to zero.
	}
	Q = SymProd(V,D*~V);
}


inline Matrix<3,3> skewSymmetric(const Matrix<3>& vector) {
	Matrix<3,3> result = zeros<3,3>();
	result(0,1) = -vector[2]; result(0,2) = vector[1];
	result(1,0) = vector[2];    result(1,2) = -vector[0];
	result(2,0) = -vector[1]; result(2,1) = vector[0];

	return result;
}


template <size_t xDim>
void computeMeanCov(const std::vector<Matrix<xDim>>& samples, Matrix<xDim>& mean, Matrix<xDim, xDim>& Cov)
{
	mean.reset(); Cov.reset();
	int Num = (int)samples.size();

	for(int i = 0; i < Num; i++){
		mean += samples[i];
	}
	mean /= (Num * 1.0);

	for(int i = 0; i < Num; i++){
		Cov += (samples[i] - mean) * ~(samples[i] - mean);
	}
	Cov /= (Num * 1.0);
}

inline void drawHelicopter(const int& cal_quadrotor)
{
	// Visualization parameters
	double beamWidth     = 0.015; // m
	double beamHeight    = 0.0065; // m
	double beamRadius    = 0.02; // m
	double motorRadius   = 0.015; // m
	double motorHeight   = 0.02; // m
	double rotorRadius   = 0.10; // m
	double rotorHeight   = 0.005; // m
	double centerSide    = 0.0889; // m
	double centerHeight  = 0.0365; // m
	double centerTopSide = 0.03; // m
	double flagLength    = 0.0508; // m
	double tileSize      = 1;  // m
	double length      = 0.3429/2; 
	// visualization
	//CAL_Initialisation(true, true, true);
	int obj;

	// Quadrotor
	CAL_SetGroupColor(cal_quadrotor, 0.05, 0.05, 0.05);
	CAL_CreateBox(cal_quadrotor, 2*length, beamWidth, beamHeight, 0, 0, 0);
	CAL_CreateBox(cal_quadrotor, beamWidth, 2*length, beamHeight, 0, 0, 0);
	CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, length, 0, beamHeight / 2 + motorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, -length, 0, beamHeight / 2 + motorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, 0, length, beamHeight / 2 + motorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, motorRadius, motorHeight, 0, -length, beamHeight / 2 + motorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, length, 0, 0, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, -length, 0, 0, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, 0, length, 0, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, beamRadius, beamHeight, 0, -length, 0, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, length, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
	CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, -length, 0, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
	CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, 0, length, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
	CAL_CreateCylinder(cal_quadrotor, rotorRadius, rotorHeight, 0, -length, beamHeight / 2 + motorHeight + rotorHeight / 2, &obj);
	CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);
	CAL_SetObjectColor(obj, 0, 0, 0, 0.1);
	CAL_CreateBox(cal_quadrotor, centerSide, centerSide, beamHeight, 0, 0, 0, &obj);
	CAL_SetObjectOrientation(obj, 0, 0, (float) (M_PI*0.25));
	CAL_CreateBox(cal_quadrotor, flagLength, beamWidth + 0.001, beamHeight + 0.001, length / 1.65, 0, 0, &obj);
	CAL_SetObjectColor(obj, 1, 0.15, 0);

	float flagTriangle[18] = {length / 1.65 - flagLength / 2, 0, -beamHeight / 2,
		length / 1.65, 0, -beamHeight / 2 - flagLength / 2,
		length / 1.65 + flagLength / 2, 0, -beamHeight / 2,
		length / 1.65 + flagLength / 2, 0, -beamHeight / 2,
		length / 1.65, 0, -beamHeight / 2 - flagLength / 2,
		length / 1.65 - flagLength / 2, 0, -beamHeight / 2};
	CAL_CreateTriangles(cal_quadrotor, 2, flagTriangle, &obj);
	CAL_SetObjectColor(obj, 1, 0.15, 0);

	float polygon1[18] = {-sqrt(2.0)*centerSide/2, 0, 0,
		-sqrt(2.0)*centerSide/2+centerHeight, 0, centerHeight,
		sqrt(2.0)*centerSide/2-centerHeight, 0, centerHeight,
		sqrt(2.0)*centerSide/2, 0, 0,
		sqrt(2.0)*centerSide/2-centerHeight, 0, -centerHeight,
		-sqrt(2.0)*centerSide/2+centerHeight, 0, -centerHeight};
	CAL_CreatePolygon(cal_quadrotor, 6, polygon1, &obj);
	CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
	float polygon2[18] = {-sqrt(2.0)*centerSide/2, 0, 0,
		-sqrt(2.0)*centerSide/2+centerHeight, 0, -centerHeight,
		sqrt(2.0)*centerSide/2-centerHeight, 0, -centerHeight,
		sqrt(2.0)*centerSide/2, 0, 0,
		sqrt(2.0)*centerSide/2-centerHeight, 0, centerHeight,
		-sqrt(2.0)*centerSide/2+centerHeight, 0, centerHeight};
	CAL_CreatePolygon(cal_quadrotor, 6, polygon2, &obj);
	CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
	float polygon3[18] = {0, -sqrt(2.0)*centerSide/2, 0,
		0, -sqrt(2.0)*centerSide/2+centerHeight, centerHeight,
		0, sqrt(2.0)*centerSide/2-centerHeight, centerHeight,
		0, sqrt(2.0)*centerSide/2, 0,
		0, sqrt(2.0)*centerSide/2-centerHeight, -centerHeight,
		0, -sqrt(2.0)*centerSide/2+centerHeight, -centerHeight};
	CAL_CreatePolygon(cal_quadrotor, 6, polygon3, &obj);
	CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);
	float polygon4[18] = {0, -sqrt(2.0)*centerSide/2, 0,
		0, -sqrt(2.0)*centerSide/2+centerHeight, -centerHeight,
		0, sqrt(2.0)*centerSide/2-centerHeight, -centerHeight,
		0, sqrt(2.0)*centerSide/2, 0,
		0, sqrt(2.0)*centerSide/2-centerHeight, centerHeight,
		0, -sqrt(2.0)*centerSide/2+centerHeight, centerHeight};
	CAL_CreatePolygon(cal_quadrotor, 6, polygon4, &obj);
	CAL_SetObjectColor(obj, 0.15, 0.15, 0.15);

}

inline void visualize(const Matrix<3>& xTrue, const Matrix<3,3>& RTrue, double t, const int& cal_quadrotor) {
	
	float p[3] = {(float) xTrue[0], (float) xTrue[1], (float) xTrue[2]};

	Matrix<3,3> rx(zeros<3,3>());
	rx(0,0) = 1.0; rx(1,1) = 0.0; rx(2,2) = 0.0; rx(1,2) = 1.0; rx(2,1) = -1.0;
	Matrix<4> q = quatFromRot(RTrue * rx);

	float o[4] = {(float) q[0], (float) q[1], (float) q[2], (float) q[3]};
	//CAL_SetGroupOrientation(cal_quadrotor, M_PI*1.5, 0.0, 0.0);
	CAL_AddGroupKeyState(cal_quadrotor, (float) t, p, o);

	return;
}

#endif