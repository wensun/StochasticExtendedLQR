///////////////////
//This code is used to implement the idea of real-time replanning. 
//namely, replanning ahead. 
//assume the time for replanning is also dt
///////////////////////////////////

#define _CRT_RAND_S

#define DIM 2
#define X_DIM 4  //X = (x,y,theta, v)  theta: orientation, v speed.
#define U_DIM 2  //a, phi, a: acceleration, phi: steering wheel angel
#define Z_DIM 3  //it can observe position x and position y.  
#define INFTY 9e9


//define constratins:


const static double goal_radius = 0.02;
const static double plan_goal_radius = 0.05;

const static double car_l = 1.0;
const static double MaxPlanTime = 1.0;


#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>
#include <conio.h>
#include <windows.h>
#include <algorithm>
#include <ppl.h>
#include <omp.h>
#include "callisto.h"
//#include "callisto41.h"
#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"
#include "Costfunction.h"
#include "IterativeLQG.h"
#include "ExtendedLQR.h"

using namespace Concurrency;
//using namespace std;

static Matrix<2> goal;
static Matrix<X_DIM> cgoal; //goal in configuraiton space.
static Matrix<X_DIM> start;
static Matrix<X_DIM, X_DIM> P0; //covariance

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
	goal[0] = 0.0;
	goal[1] = 0.5;

	start[0] = 0.01; start[1] = 0.01; start[2] = 0.01; start[3] = 0.01;
	//cgoal[0] = 1.8;  cgoal[1] = -0.5; cgoal[3] = 0.0; cgoal[4] = 0.0;
	// callisto params
	//CAL_ShowView (0, "Callisto");
	//CAL_SetViewOptions (0, CAL_HIDEGROUNDPLANE);
	CAL_SetViewParams(0, 0, 0, 2, 0, 0, 0, 0.0); 
	//CAL_SetViewNavigationSpeed (0, 1000.0);
	//CAL_SetViewNearClippingDistance (0, 2.0);
	//CAL_SetViewFogDistance(10);
	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 0, 0, 1, 1);
	CAL_SetGroupVisibility(cal_obstacles, 0, false);

	CAL_CreateGroup(&cal_box, 0, false, "box");
	//CAL_SetGroupCastShadows(cal_box, false);
	CAL_CreateBox(cal_obstacles, 0.2, 0.3, 5.0, 0.3, 0.40, 0.0);
	CAL_CreateBox(cal_box, 0.2, 0.3, 0.02, 0.3, 0.40, 0.0);

	CAL_CreateGroup(&cal_line, 0, false, "line");
	CAL_SetGroupColor(cal_line, 0, 0, 0, 1);
	int np = 2;
	float line1[6] = {0.8, -0.8, 0.0, 0.8, 0.8, 0.0};
	CAL_CreatePolyline(cal_line, 1, &np,  line1);
	float line2[6] = {-0.8, -0.8, 0.0, -0.8, 0.8, 0.0};
	CAL_CreatePolyline(cal_line, 1, &np,  line2);
	float line3[6] = {0.8, 0.8, 0.0, -0.8, 0.8, 0.0};
	CAL_CreatePolyline(cal_line, 1, &np,  line3);
	float line4[6] = {0.8, -0.8, 0.0, -0.8, -0.8, 0.0};
	CAL_CreatePolyline(cal_line, 1, &np,  line4);

	CAL_CreateGroup(&cal_paths, 0, true, "Paths");
	CAL_SetGroupColor(cal_paths, 0, 0, 1);

	CAL_CreateGroup(&cal_link1, 0, true, "link1");
	CAL_SetGroupColor(cal_link1, 0, 1, 0);
	//CAL_SetGroupCastShadows(cal_link1, false);

	CAL_CreateGroup(&cal_link2, 0, true, "link2");
	CAL_SetGroupColor(cal_link2, 1, 1, 0);
	//CAL_SetGroupCastShadows(cal_link2, false);
	
	Dynamics dyn;
	Matrix<2> poslink1 = dyn.FK(start, 1, dyn.l(0,0)/2.0);
	Matrix<2> poslink2 = dyn.FK(start, 2, dyn.l(1,0)/2.0);

	CAL_CreateBox(cal_link1, dyn.l(0,0), 0.05, 0.02, 0.0, 0.0, 0.0);
	CAL_CreateBox(cal_link2, dyn.l(1,0), 0.02, 0.02, 0.0,0.0, 0.0);
	setTwoLinks(cal_link1, cal_link2, start);


	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 0, 0.0, 0.6);

	CAL_CreateGroup(&cal_ellipse_trunc, 0, false, "Ellipse_trunc");
	CAL_SetGroupColor(cal_ellipse_trunc, 1, 0, 0);

	CAL_CreateGroup(&cal_execute, 0, false, "execute");
	CAL_SetGroupColor(cal_execute, 1, 0, 0);
	
	CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
	//CAL_SetGroupCastShadows(cal_goal, false);
	CAL_SetGroupColor(cal_goal, 0, 1, 1, 1);
	CAL_CreateCylinder(cal_goal, (float) goal_radius, 0.01f, 0, 0, 0);
	CAL_SetGroupPosition(cal_goal, (float) goal[0], (float) goal[1], -0.025f);
	CAL_SetGroupOrientation(cal_goal, (float) (M_PI*0.5), 0, 0);

	//CAL_CreateSphere(cal_point, 0.2, 7.0/4.0*3.1416, 0.0, 0.0);
	
}


void showaplan(const std::vector<std::pair<Matrix<XDIM, 1>, Matrix<UDIM, 1>>>& plan)
{
	
	double l = (double)plan.size();
	for(int i = 1; i < (int)plan.size(); i+=1){

		Matrix<XDIM> state = plan[i].first;
		setTwoLinks(cal_link1, cal_link2, state);
		std::cout<<state<<std::endl;
		
		clock_t start_t = clock();
		while(true){
			if ((double)(clock() - start_t)/ CLOCKS_PER_SEC < 0.01)
				;
			else
				break;
		}

/*		int obj1, obj2;
		Dynamics dyn;
		CAL_CreateBox(cal_link1, dyn.l(0,0), 0.05, 0.02, 0.0, 0.0, 0.0, &obj1);
		CAL_CreateBox(cal_link2, dyn.l(1,0), 0.02, 0.02, 0.0,0.0, 0.0, &obj2);
		CAL_SetObjectColor(obj1, 0, (l-i)/l*1.0, 0);
		CAL_SetObjectColor(obj2, (l-i)/l*1.0, (l-i)/l*1.0, 0);
		Matrix<XDIM> state = plan[i].first;
		Matrix<2> poslink1 = dyn.FK(state, 1, dyn.l(0,0)/2.0);
		Matrix<2> poslink2 = dyn.FK(state, 2, dyn.l(1,0)/2.0);
		CAL_SetObjectPosition(obj1, poslink1[0], poslink1[1], 0.0);
		CAL_SetObjectPosition(obj2, poslink2[0], poslink2[1], 0.0);
		CAL_SetObjectOrientation(obj1, 0.0, 0.0, state[0]);
		CAL_SetObjectOrientation(obj2, 0.0, 0.0, state[1]+state[0]);*/
	}

	setTwoLinks(cal_link1, cal_link2, plan[(int)plan.size() - 1].first);
	return;
}

void runIterativeLQG(const int& T)
{
	IterativeLQG ilqg(start, goal, cgoal, T, cal_link1, cal_link2, cal_obstacles, lambda);
	std::vector<std::pair<Matrix<UDIM, XDIM>, Matrix<UDIM, 1>>> policy;
	policy.clear(); policy.resize(T+1);

	for(int i = 0; i < (int)policy.size(); i++){
		policy[i].first = zeros<UDIM, XDIM>();
		policy[i].second = zeros<UDIM,1>();
	}

	//ilqg.initializePolicy(policy);
	ilqg.iterate();
	double cost = ilqg.computeCost(ilqg.nominalPlan);

	std::cout<<ilqg.totalTime<<"	"<<cost<<std::endl;
	//foutiterativelqg<<ilqg.totalTime<<"	"<<cost<<std::endl;
	showaplan(ilqg.nominalPlan);
}

void runStoExtendedLQR(const int& T)
{
	ExtendedLQR elqr(start, goal, cgoal, T, cal_link1, cal_link2, cal_obstacles, lambda);

	elqr.executeExtendedLQR();

	std::cout<<elqr.totalTime<<"	"<<elqr.totalCost<<std::endl;
	//foutextendedlqr<<elqr.totalTime<<"	"<<elqr.totalCost<<std::endl;
	showaplan(elqr.nominalPlan);

}

////////////////////End of replanning////////////////////////////////
int main()
{

	//foutextendedlqr.open("ExtendedLQRdata.txt");
	//foutiterativelqg.open("IterativeLQGdata.txt");
	srand(1000);
	CAL_Initialisation (true);
	initEnvironment();
	int T = 200;
	std::cout.precision(5);

	//std::ofstream fout("Comparison.txt", std::ofstream::app);
	for(int i = 0; i < 1; i++){
		
		cgoal[0] = -3.14 + 2 * 3.14 * random();
		//cgoal[0] = 2.7196;
		cgoal[1] = -3.14 + 2 * 3.14 * random();
		//cgoal[1] = -2.7347;
		cgoal[2] = 0.0;
		cgoal[3] = 0.0;

		//lambda = 0.01 + i * 0.00015; 
		lambda = 0.00;
		ExtendedLQR elqr(start, goal, cgoal, T, cal_link1, cal_link2, cal_obstacles, lambda);
		elqr.executeExtendedLQR();
		//showaplan(elqr.nominalPlan);
		std::cout<<"ExtendedLQR: "<<elqr.totalCost<<"	"<<elqr.totalTime<<"	"<<elqr.numTotalIter<<std::endl;

		//IterativeLQG ilqg(start, goal, cgoal, T, cal_link1, cal_link2, cal_obstacles, lambda);
		//ilqg.iterate();
		//showaplan(ilqg.nominalPlan);
		//std::cout<<"iterativeLQG: "<<ilqg.totalCost<<"	"<<ilqg.totalTime<<"	"<<ilqg.numTotalIter<<std::endl;
		//fout<<elqr.totalCost<<"	"<<elqr.totalTime<<"	"<<elqr.numTotalIter<<"	"<<ilqg.totalCost<<"	"<<ilqg.totalTime<<"	"<<ilqg.numTotalIter<<std::endl;
	}
	//fout.close();

	/*Matrix<XDIM> state = start;
	state(0,0) += 0.1; state(1,0) += 3.1; state(2,0) -= 0.24; state(3,0) += 0.04;

	Matrix<XDIM> xnext; 
	Matrix<XDIM, XDIM> M;
	Matrix<UDIM, 1> u; u(0,0) = 0.3; u(1,0) = -0.5;

	Dynamics dyn;
	dyn.discrete_dynamics(state, u, xnext, M);

	std::cout<<xnext<<std::endl<<M<<std::endl;

	std::vector<Matrix<XDIM, XDIM>> F;
	std::vector<Matrix<XDIM, UDIM>> G;
	std::vector<Matrix<XDIM,1>> e;
	Matrix<XDIM, XDIM> At; Matrix<XDIM, UDIM> Bt; Matrix<XDIM,1> ct;
	dyn.linearize_discrete_dynamics(state, u, At, Bt, ct, F, G, e);

	Matrix<XDIM> augstate = state;
	augstate[0] += 0.01; augstate[1] += 0.01; augstate[2] += 0.01; augstate[3] += 0.01;
	Matrix<XDIM> augnext; Matrix<XDIM, XDIM> augM;
	dyn.discrete_dynamics(augstate, u, augnext, augM);
	std::cout<<augnext<<std::endl;
	std::cout<<augM<<std::endl;
	

	//std::vector<Matrix<XDIM>> samples; samples.clear();
	//for(int i = 0; i < 4000; i++){
	//	Matrix<XDIM> s = dyn.Euler_integral_forward_noise(augstate, u);
	//	samples.push_back(s);
	//}
	//Matrix<XDIM> mean; Matrix<XDIM, XDIM> Cov;
	//computeMeanCov(samples, mean, Cov);
	//std::cout<<mean<<std::endl;
	//std::cout<<Cov<<std::endl;


	std::cout<<At*augstate + Bt * u + ct<<std::endl;
	std::cout<<e[0] + F[0] * augstate + G[0] * u<<std::endl;*/


	//Matrix<XDIM> xold;
	//dyn.discrete_inverse_dynamics(xnext, u, xold, M);

	//std::cout<<state<<std::endl<<xold<<std::endl;


	int num;
	std::cin>>num;

	
	// end Callisto
	CAL_End();

	return 0;
}