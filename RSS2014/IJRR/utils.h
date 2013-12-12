#ifndef __UTILS_H__
#define __UTILS_H__

inline double sqr(double x) {
	return x*x;
}

// Noise model
inline double sigmoid(double x, double mean)
{
	double y = (x - mean);
	double s = (y/sqrt(1+y*y))+1.0;
	//double s = (y/(1.0 + abs(y)))+1.0;

	if (x < mean) 
		return s*0.1;
	else
		return s*10.0;
}

// Noise model
inline double sigmoid2(double x, double mean)
{
	double y = (mean - x);
	double s = (y/sqrt(1+y*y))+1.0;
	//double s = (y/(1.0 + abs(y)))+1.0;

	if (x >= mean) 
		return 0.1; //0.1/(x - mean);
	else
		return s*10.0;
}

inline double random() {
	return ((double) rand()) / RAND_MAX;
}

inline double random(double low, double high) { 
	return low + ((double)rand()/(double)RAND_MAX)*(high-low); 
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

template <int size>
inline Matrix<size> sampleGaussian(const Matrix<size>& mean, const Matrix<size, size>& var) {
	Matrix<size> sample;
	for (int j = 0; j < size; ++j) {
		sample[j] = normal();
	}
	Matrix<size, size> SVec, SVal;
	jacobi(var, SVec, SVal);
	for (int i = 0; i < size; ++i) {
		SVal(i,i) = sqrt(SVal(i,i));
	}
	return SVec * SVal * sample + mean;
}

#include "glut.h"
const double DEG2RAD = M_PI/180;

inline double mod2pi(double x) 
{
	// Returns a value 0 <= x < 2*PI
	double result = fmod(x, 2*M_PI);
	if (result < 0) {
		result += (2*M_PI);
	}
	return result;
}

inline void drawUnitCircle(void *userdef, float time)
{
	glLineWidth(4.0f);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINE_LOOP);
	glColor3f(0.0f, 0.0f, 0.0f);
	//glColor3f(0.5f, 0.0f, 1.0f);

	for (int i = 0; i < 36; i++)
	{
		float theta = (float)(i*10.0*DEG2RAD);
		glVertex3f(cosf(theta), sinf(theta), 0.0f);
	}
	glEnd();
	glEnable(GL_LIGHTING);
	glLineWidth(4.0f);
}

inline void drawUnitSphere(void *userdef, float time)
{
	glDisable(GL_LIGHTING);
	glLineWidth(4.0f);
	glColor3f(0.0f, 0.0f, 0.0f);
	
	glutWireSphere(1.0, 15, 10);

	glLineWidth(4.0f);
	glEnable(GL_LIGHTING);
}

inline Matrix<4,1> quatFromRot(const Matrix<3,3>& R) 
{
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

inline Matrix<3,3> rotFromQuat(const Matrix<4,1>& q)
{
	Matrix<3,3> R;

	R(0,0) = q[3]*q[3] + q[0]*q[0] - q[1]*q[1] - q[2]*q[2];
	R(0,1) = 2.0*q[0]*q[1] - 2.0*q[3]*q[2];
	R(0,2) = 2.0*q[0]*q[2] + 2.0*q[3]*q[1];

	R(1,0) = 2.0*q[0]*q[1] + 2.0*q[3]*q[2];
	R(1,1) = q[3]*q[3] - q[0]*q[0] + q[1]*q[1] - q[2]*q[2];
	R(1,2) = 2.0*q[1]*q[2] - 2.0*q[3]*q[0];

	R(2,0) = 2.0*q[0]*q[2] - 2.0*q[3]*q[1];
	R(2,1) = 2.0*q[1]*q[2] + 2.0*q[3]*q[0];
	R(2,2) = q[3]*q[3] - q[0]*q[0] - q[1]*q[1] + q[2]*q[2];
	
	return R;
}

inline Matrix<4,1> quatFromAA(const Matrix<3>& axis, double angle)
{
	Matrix<4,1> q;
	double sa = sin(-angle*0.5);
	double ca = cos(angle*0.5);
	q(0,0) = sa*axis[0];
	q(1,0) = sa*axis[1];
	q(2,0) = sa*axis[2];
	q(3,0) = ca;

	return q;
}

inline Matrix<3,3> cpMatrix(const Matrix<3>& a) {
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

#endif