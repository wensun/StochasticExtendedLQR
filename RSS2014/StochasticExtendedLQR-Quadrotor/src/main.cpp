///////////////////
//This code is used to implement the idea of real-time replanning. 
//namely, replanning ahead. 
//assume the time for replanning is also dt
///////////////////////////////////

#define _CRT_RAND_S

#define DIM 3
#define X_DIM 12  //X = (x,y,theta, v)  theta: orientation, v speed.
#define U_DIM 4  //a, phi, a: acceleration, phi: steering wheel angel
#define Z_DIM 3  //it can observe position x and position y.  
#define INFTY 9e9

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
#include <ppl.h>
#include <omp.h>
#include <time.h>
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Obstacles.h"
#include "Dynamics.h" 
#include "Costfunction.h"
#include "IterativeLQG.h"
#include "ExtendedLQR.h"

using namespace Concurrency;
//using namespace std;

static Matrix<3> goal;
static Matrix<X_DIM> cgoal; //goal in configuraiton space.
static Matrix<X_DIM> start;
static Matrix<X_DIM, X_DIM> P0; //covariance

static Matrix<3> Top;
static Matrix<3> Bottom;


double lambda;

std::ofstream foutextendedlqr;
std::ofstream foutiterativelqg;
/***********************Environment settings***************************/
int cal_environment;
int cal_obstacles;
int cal_goal;
int cal_line;
int cal_rrt;
int cal_paths;
int cal_ellipse;
int cal_ellipse_trunc;
int cal_point;
int cal_cienvironment, cal_ciobstacles, cal_cipoint;
int cal_plan;
int cal_execute;
int cal_robot;
int cal_link1;
int cal_link2;
int cal_box;
int cal_quadrotor;
static std::vector<Matrix<2,1>> obstaclesSet; //ax + by = b
/******************End of Environment settings************************/


double Random()
{
	double a;
	a = rand()%1000;
	return a / 1000;
}

void initEnvironment() 
{
	CAL_SetViewParams(0, 0, 0, 14, 0, 0, 0, 0.0); 
	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 0, 0, 1, 1);
//	CAL_SetGroupVisibility(cal_obstacles, 0, false);

	Obstacles obss;
	for(int i = 0; i < (int)obss.obstacles.size(); i++){
		Obstacles::obs tmpobs = obss.obstacles[i];
		int obj;
		CAL_CreateCylinder(cal_obstacles, tmpobs.radius, 6.0, 0.0, 0.0, 0.0, &obj);
		CAL_SetObjectPosition(obj, tmpobs.pos[0], tmpobs.pos[1], tmpobs.pos[2]);
		if(tmpobs.dim == 2)
			CAL_SetObjectOrientation(obj, M_PI/2.0, 0.0, 0.0);
		else if(tmpobs.dim == 1)
			CAL_SetObjectOrientation(obj, 0.0, 0.0, 0.0);
		else
			CAL_SetObjectOrientation(obj, 0.0, 0.0, M_PI/2.0);
	}


	CAL_CreateGroup(&cal_line, 0, false, "line");
	CAL_SetGroupColor(cal_line, 0, 0, 0, 1);
	int np = 2;

	Obstacles obs;
	Top = obs.Top;
	Bottom = obs.Bottom;

	float line1[6] = {Top[0], Top[1], Top[2], Top[0], Top[1], -Top[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line1);
	float line2[6] = {Top[0], Top[1], Top[2], Top[0], -Top[1], Top[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line2);
	float line3[6] = {Top[0], Top[1], Top[2], -Top[0], Top[1], Top[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line3);
	float line4[6] = {Bottom[0], Bottom[1], Bottom[2], Bottom[0], Bottom[1], -Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line4);
	float line5[6] = {Bottom[0], Bottom[1], Bottom[2], Bottom[0], -Bottom[1], Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line5);
	float line6[6] = {Bottom[0], Bottom[1], Bottom[2], -Bottom[0], Bottom[1], Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line6);
	float line7[6] = {Bottom[0], Bottom[1], -Bottom[2], Bottom[0], -Bottom[1], -Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line7);
	float line8[6] = {Bottom[0], Bottom[1], -Bottom[2], -Bottom[0], Bottom[1], -Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line8);
	float line9[6] = {Bottom[0], -Bottom[1], Bottom[2], Bottom[0], -Bottom[1], -Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line9);
	float line10[6] = {Bottom[0], -Bottom[1], Bottom[2], -Bottom[0], -Bottom[1], Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line10);
	float line11[6] = {-Bottom[0], Bottom[1], Bottom[2], -Bottom[0], -Bottom[1], Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line11);
	float line12[6] = {-Bottom[0], Bottom[1], Bottom[2], -Bottom[0], Bottom[1], -Bottom[2]};
	CAL_CreatePolyline(cal_line, 1, &np,  line12);


	CAL_CreateGroup(&cal_paths, 0, true, "Paths");
	CAL_SetGroupColor(cal_paths, 1, 0, 1);

	CAL_CreateGroup(&cal_link1, 0, true, "link1");
	CAL_SetGroupColor(cal_link1, 0, 1, 0);
	//CAL_SetGroupCastShadows(cal_link1, false);

	CAL_CreateGroup(&cal_link2, 0, true, "link2");
	CAL_SetGroupColor(cal_link2, 1, 1, 0);
	//CAL_SetGroupCastShadows(cal_link2, false);

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 0, 0.0, 0.6);

	CAL_CreateGroup(&cal_ellipse_trunc, 0, false, "Ellipse_trunc");
	CAL_SetGroupColor(cal_ellipse_trunc, 1, 0, 0);

	CAL_CreateGroup(&cal_execute, 0, false, "execute");
	CAL_SetGroupColor(cal_execute, 1, 0, 0);

	CAL_CreateGroup(&cal_quadrotor, 0, false, "QuadRotor");



	//CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
	//CAL_SetGroupCastShadows(cal_goal, false);
	//CAL_SetGroupColor(cal_goal, 0, 1, 1, 1);
	//CAL_CreateCylinder(cal_goal, (float) goal_radius, 0.01f, 0, 0, 0);
	//CAL_SetGroupPosition(cal_goal, (float) goal[0], (float) goal[1], -0.025f);
	//CAL_SetGroupOrientation(cal_goal, (float) (M_PI*0.5), 0, 0);
	//CAL_CreateSphere(cal_point, 0.2, 7.0/4.0*3.1416, 0.0, 0.0);
	
}


void showplan(const std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>>& nominalPlan)
{
	for(int i = 0; i < (int)nominalPlan.size(); i++){
		CAL_CreateSphere(cal_paths, 0.05, nominalPlan[i].first[0], nominalPlan[i].first[1], nominalPlan[i].first[2]);
	}
	return;
}


////////////////////End of replanning////////////////////////////////
int main()
{

	srand(1372474623);
	CAL_Initialisation(true, true, true);

	initEnvironment();
	//std::cout.precision(4);

	/*int obj;
	CAL_CreateSphere(cal_paths, 0.05, 0.0,0.0,0.0, &obj);

	int T = 150;
	double dt = 0.032;

	start = zeros<XDIM,1>(); start[0] = -3; start[1] = -3; start[2] = -3;
	cgoal = -start;

	drawHelicopter(cal_quadrotor);
	Matrix<3,3> RotTrue = identity<3>();
	visualize(start.subMatrix<3,1>(0,0), RotTrue, 0, cal_quadrotor);

	std::ofstream fout("comparisonNew1.txt", std::ofstream::app);
	while(dt <= 0.05){
	
		for(int i = 1; i <= 20; i++){

			/*int r = rand()%4;
			start[1] = -3.0;
			if(r == 0){
				start[0] = -3.0 + (double)rand()/RAND_MAX * 6.0;
				start[2] = 3.0;
			}
			else if(r == 1){
				start[0] = -3.0 + (double)rand()/RAND_MAX * 6.0;
				start[2] = -3.0;
			}
			else if(r == 2){
				start[0] = -3.0;
				start[2] = -3.0 + (double)rand()/RAND_MAX * 6.0;
			}
			else{
				start[0] = 3.0;
				start[2] = -3.0 + (double)rand()/RAND_MAX * 6.0;
			}
		
			std::cout<<start.subMatrix<3,1>(0,0)<<std::endl;
			//start[0] = 3; start[1] = -3.0; start[2] = 2.72533;

			int random = rand() % 12;
			cgoal[random / 4] = ((double) rand() / RAND_MAX) * 4.5 - 2.25; //M_PI; //0; //
			cgoal[((random / 4) + 1) % 3] = ((random % 4) / 2 == 0 ? 2.25 : -2.25) + ((double) rand() / RAND_MAX) * 0.02 - 0.01;
			cgoal[((random / 4) + 2) % 3] = ((random % 4) % 2 == 0 ? 2.25 : -2.25) + ((double) rand() / RAND_MAX) * 0.02 - 0.01;
			start = -cgoal;

			std::cout<<start.subMatrix<3,1>(0,0)<<std::endl;
			CAL_SetObjectPosition(obj, start[0], start[1], start[2]);

			IterativeLQG ilqg(start, cgoal, T, dt);
			ilqg.iterate();
			std::cout<<dt<<"	"<<"Iterative LQG:	"<<ilqg.totalCost<<"	"<<ilqg.totalTime<<"	"<<ilqg.numTotalIter<<std::endl;

			ExtendedLQR elqr(start, cgoal, T, dt, cal_paths);
			elqr.executeExtendedLQR();
			std::cout<<dt<<"	"<<"ExtendedLQG:	"<<elqr.totalCost<<"	"<<elqr.totalTime<<"	"<<elqr.numTotalIter<<std::endl;

			fout<<dt<<"	"<<elqr.totalCost<<"	"<<elqr.totalTime<<"	"<<elqr.numTotalIter<<"	"<<ilqg.totalCost<<"	"<<ilqg.totalTime<<"	"<<ilqg.numTotalIter<<std::endl;
		}
		dt += 0.002;
	}
	fout.close();*/
	
	/*Matrix<XDIM> x;
	for(int i = 1; i < (int)elqr.nominalPlan.size(); i++){
		
		x = elqr.nominalPlan[i].first;
		RotTrue = exp(skewSymmetric(x.subMatrix<3,1>(6,0)));
		visualize(x.subMatrix<3,1>(0,0), RotTrue, i*dt, cal_quadrotor);
	}*/

	
	//showplan(elqr.nominalPlan);
	//std::cout<<"end"<<std::endl;
	

	int num;
	std::cin>>num;

	
	// end Callisto
	CAL_End();

	return 0;
}