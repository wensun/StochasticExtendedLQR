///////////////////
//This code is used to implement the idea of real-time replanning. 
//namely, replanning ahead. 
//assume the time for replanning is also dt
///////////////////////////////////

#define _CRT_RAND_S

//define constratins:
const static double goal_radius = 0.02;
const static double plan_goal_radius = 0.05;

const static double car_l = 1.0;

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
#include "CostFunction.h"
#include "IterativeLQG.h"
//#include "ExtendedLQR.h"

using namespace Concurrency;
//using namespace std;


double lambda;
static Matrix<3> Top;
static Matrix<3> Bottom;
static Matrix<XDIM> xGoal;
static Matrix<BDIM> bstart;
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

	Obstacles obss;
	for(int i = 0; i < (int)obss.obstacles.size(); i++){
		Obstacles::obs tmpobs = obss.obstacles[i];
		int obj;
		if(tmpobs.dim == 2)
			CAL_CreateCylinder(cal_obstacles, tmpobs.radius, 8.0, 0.0, 0.0, 0.0, &obj);
		else
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

	CAL_CreateGroup(&cal_link2, 0, true, "link2");
	CAL_SetGroupColor(cal_link2, 1, 1, 0);

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 0, 0.0, 0.6);

	CAL_CreateGroup(&cal_ellipse_trunc, 0, false, "Ellipse_trunc");
	CAL_SetGroupColor(cal_ellipse_trunc, 1, 0, 0);

	CAL_CreateGroup(&cal_execute, 0, false, "execute");
	CAL_SetGroupColor(cal_execute, 1, 0, 0);

	CAL_CreateGroup(&cal_point, 0.0,true, "Point");
	CAL_CreateSphere(cal_point, 0.001, 0.0, 0.0, 0.0);

//	CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
//	CAL_SetGroupColor(cal_goal, 0, 1, 1, 1);
//	CAL_CreateSphere(cal_goal, 0.2, xGoal[0], xGoal[1], xGoal[2]);
}

////////////////////End of replanning////////////////////////////////
int main()
{

	srand(1000);
	CAL_Initialisation (true);
	initEnvironment();

	double dt = 0.1;
	xGoal = zeros<XDIM,1>();
	xGoal[0] = 0.0; xGoal[1] = 0.0; xGoal[2] = 3.0;

	



	int num;
	std::cin>>num;

	
	// end Callisto
	CAL_End();

	return 0;
}