#include "point.h"

PointPlanner::PointPlanner() {
#ifndef POINT2D
	std::cerr << "Improper setup. Define POINT2D in planner.h" << std::endl;
	std::exit(-1);
#endif
}

// Dynamics model
Matrix<X_DIM> PointPlanner::f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m, double tstep) // m ~ N(0,I)
{  
	Matrix<X_DIM> xNew;

#if defined(PROBLEM1)
	//xNew[0] = x[0] + x[2]*DT + 0.5 * (u[0] + 0.1*u[0]*m[0]) * DT*DT;
	//xNew[1] = x[1] + x[3]*DT + 0.5 * (u[1] + 0.1*u[1]*m[1]) * DT*DT;
	//xNew[2] = x[2] + (u[0] + 0.1*u[0]*m[0]) * DT;
	//xNew[3] = x[3] + (u[1] + 0.1*u[1]*m[1]) * DT;

	xNew[0] = x[0] + u[0] * DT + sqrt(0.01*u[0]*u[0]*DT + 0.001)*m[0];
	xNew[1] = x[1] + u[1] * DT + sqrt(0.01*u[1]*u[1]*DT + 0.001)*m[1];

#elif defined(PROBLEM2)
	//xNew[0] = x[0] + u[0] * DT + 0.1*sqrt(DT)*m[0];
	//xNew[1] = x[1] + u[1] * DT + 0.1*sqrt(DT)*m[1];

	xNew[0] = x[0] + u[0] * DT + 0.1*abs(u[0])*sqrt(DT)*m[0];
	xNew[1] = x[1] + u[1] * DT + 0.1*abs(u[1])*sqrt(DT)*m[1];

	//xNew[0] = x[0] + u[0] * DT + sqrt(0.01*u[0]*u[0]*DT + 0.01)*m[0];
	//xNew[1] = x[1] + u[1] * DT + sqrt(0.01*u[1]*u[1]*DT + 0.01)*m[1];
#endif

	return xNew;
}

// Observation model
Matrix<Z_DIM> PointPlanner::h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) // n ~ N(0,I)
{
	Matrix<Z_DIM> z;

#if defined(PROBLEM1)
	z[0] = x[0] + sqrt(0.5*0.5*(5-x[0])*(5-x[0]) + 0.001)*n[0];
	z[1] = x[1] + sqrt(0.5*0.5*(5-x[0])*(5-x[0]) + 0.001)*n[1];
#elif defined(PROBLEM2)
	//double scale = sigmoid(x[0], -2);
	//z[0] = x[0] + sqrt(scale)*n[0];
	//z[1] = x[1] + sqrt(scale)*n[1];

	//z[0] = x[0] + sqrt(sqr(0.5*(-0.5-x[0])) + 0.1)*n[0];
	//z[1] = x[1] + sqrt(sqr(0.5*(-0.5-x[0])) + 0.1)*n[1];

	z[0] = x[0] + sqrt(sqr(0.5*(-x[0])) + 0.2)*n[0];
	z[1] = x[1] + sqrt(sqr(0.5*(-x[0])) + 0.01)*n[1];

	//z[0] = 1.0/(sqr(x[0]) + sqr(x[1]) + 1.0) + 0.075*n[0];
#endif

	return z;
}

void PointPlanner::RRT() 
{
	std::vector<RRTNode> rrtTree;

	RRTNode startNode;
	startNode.x = x0;

	rrtTree.push_back(startNode);

	double goaldist = tr(~(startNode.x - xGoal) * (startNode.x - xGoal));

	while (goaldist > goalr * goalr) 
	{
		Matrix<X_DIM> sample;
		sample[0] = random() * 10.0 - 5.0;
		sample[1] = random() * 6.0 - 3.0;
		
		int col;

		CAL_CheckPointCollision(cal_environment, (float) sample[0], (float) sample[1], 0.0f, false, &col);
		if (col != 0) {
			continue;
		}

		int node = -1;
		double mindist = 9e99;
		for (int j = 0; j < (int) rrtTree.size(); ++j) {
			double ddist = tr(~(rrtTree[j].x - sample) * (rrtTree[j].x - sample));
			if (ddist < mindist) {
				mindist = ddist;
				node = j;
			}
		}
		if (node == -1) {
			continue;
		}

		Matrix<U_DIM> input;
		input[0] = random() * 1.4 - 0.7;
		input[1] = random() * 1.4 - 0.7;
		
		Matrix<X_DIM> xNew = f(rrtTree[node].x, input, zeros<M_DIM,1>());

		bool valid = true;
		if (xNew[0] < -5 || xNew[0] > 5 || xNew[1] < -3 || xNew[1] > 3) {
			valid = false;
		}
		if (!valid) {
			continue;
		}

		CAL_CheckLineCollision(cal_environment, (float) rrtTree[node].x[0], (float) rrtTree[node].x[1], 0.0f, (float) xNew[0], (float) xNew[1], 0.0f, false, &col);
		if (col != 0) {
			continue;
		}

		RRTNode newnode;
		newnode.x = xNew;
		newnode.u = input;
		newnode.bp = node;

		rrtTree.push_back(newnode);

		int np[1] = {2};
		float p[6] = {(float) rrtTree[node].x[0], (float) rrtTree[node].x[1], 0.0f, (float) xNew[0], (float) xNew[1], 0.0f};
		CAL_CreatePolyline(cal_rrt, 1, np, p);

		goaldist = tr(~(newnode.x - xGoal) * (newnode.x - xGoal));
	}

	std::deque<PathNode> path;

	int i = (int) rrtTree.size() - 1;
	PathNode node;
	node.x = rrtTree[i].x;
	node.u = zeros<U_DIM, 1>();
	path.push_front(node);

	while (i != 0) 
	{
		node.u = rrtTree[i].u;
		i = rrtTree[i].bp;
		node.x = rrtTree[i].x;

		path.push_front(node);
	}

	for (int i = 0; i < (int)path.size()-1; ++i) {
		U.push_back(path[i].u);
	}

	// Set pathlen based on solution size
	pathlen = path.size()-1;

	// Clear rrt tree visualization
	CAL_DestroyGroup(cal_rrt);
}

void PointPlanner::initTrajectory()
{
	//rrtfptr >> pathlen;
	//std::cout << "Path len: " << pathlen << std::endl;
	//U.resize(pathlen);
	//for(int i = 0; i < (int)U.size(); ++i) {
	//	Matrix<U_DIM,1>& Ui = U[i];
	//	rrtfptr >> Ui[0] >> Ui[1];
	//}

//#if defined(PROBLEM1)
//	std::string filename = "data\\rrtPoint1.txt";
//#elif defined(PROBLEM2)
//	std::string filename = "data\\rrtPointPath1.txt";
//#endif

	//RRT();

	//std::ofstream fptr(filename.c_str(),std::ios::out);
	//fptr << U.size()  << std::endl;
	//for(int i = 0; i < (int)U.size(); ++i) {
	//	fptr << ~U[i];
	//}
	//fptr.close();

	std::string filename = "data\\Path1.txt";

	std::ifstream fptr(filename.c_str(), std::ios::in);
	fptr >> pathlen;
	std::cout << "Path len: " << pathlen << std::endl;
	U.resize(pathlen);
	for(int i = 0; i < (int)U.size(); ++i) {
		Matrix<U_DIM,1>& Ui = U[i];
		fptr >> Ui[0] >> Ui[1];
	}
	fptr.close();

	// 150-200 RRT paths
	/*
	std::string filename = "data\\RRT-Paths-ML-NoML.txt";
	std::ofstream fptr(filename.c_str(),std::ios::out);
	int pathcount = 0;
	srand(time(0));
	while(pathcount < 200) 
	{
		RRT();
		
		if ((int) U.size() > 25 && (int) U.size() < 35) {
			fptr << U.size();
			for(int j = 0; j < (int)U.size(); ++j) {
				Matrix<U_DIM,1>& Uj = U[j];
				fptr << " " << Uj[0] << " " << Uj[1];
			}
			fptr << std::endl;	
			std::cout << ++pathcount << " ";
		}
		U.clear();
	}
	fptr.close();

	std::cout << "Dumped RRT paths" << std::endl;
	
	int num;
	std::cin >> num;
	*/

	//pathlen = 30;
	//U.resize(pathlen);
	//Matrix<U_DIM> uStep;

	//uStep[0] = -5.85/12.0; uStep[1] = 1.2/12.0;
	//for(int i = 0; i < 12; ++i) {
	//	U[i] = uStep;
	//}

	//uStep[0] = 0.0; uStep[1] = 1.6/6.0;
	//for(int i = 12; i < 18; ++i) {
	//	U[i] = uStep;
	//}

	//uStep[0] = 5.85/12.0; uStep[1] = 1.2/12.0;
	//for(int i = 18; i < 30; ++i) {
	//	U[i] = uStep;
	//}

	/*
	pathlen = 30;
	U.resize(pathlen);
	Matrix<U_DIM> uStep;

	// path1
	static bool toggle = true;
	if (toggle) {
		uStep[0] = 0.0; uStep[1] = 2.5/(double)20.0;
		for(int i = 0; i < 20; ++i) {
			U[i] = uStep;
		}
		uStep[0] = 0.0; uStep[1] = 1.5/(double)10.0;
		for(int i = 20; i < 30; ++i) {
			U[i] = uStep;
		}

		toggle = false;
	}
	else {
		// path2
		uStep[0] = -5.75/12.0; uStep[1] = 1.2/12.0;
		for(int i = 0; i < 12; ++i) {
			U[i] = uStep;
		}

		uStep[0] = 0.0; uStep[1] = 1.6/6.0;
		for(int i = 12; i < 18; ++i) {
			U[i] = uStep;
		}

		uStep[0] = 5.75/12.0; uStep[1] = 1.2/12.0;
		for(int i = 18; i < 30; ++i) {
			U[i] = uStep;
		}
	}
	*/

	B.resize(pathlen+1);
	Matrix<X_DIM,X_DIM> SqrtW;

	vec(x0, SqrtSigma0, B[0]);
	for (int t = 0; t < pathlen; ++t) {
		beliefDynamics(B[t], U[t], B[t+1], SqrtW);
	}

	/*
#if defined(PROBLEM1)
	//pathlen = 20;
	//U.assign(pathlen, zeros<U_DIM,1>());

	pathlen = 20;
	U.resize(pathlen);
	Matrix<U_DIM> uStep;
	uStep[0] = -0.1; uStep[1] = -0.1;
	for(int i = 0; i < 20; ++i) {
		U[i] = uStep;
	}
#elif defined(PROBLEM2)
	//RRT();

	pathlen = 30;
	U.resize(pathlen);
	Matrix<U_DIM> uStep;
	uStep[0] = 0.0; uStep[1] = -0.18;
	for(int i = 0; i < 10; ++i) {
		U[i] = uStep;
	}
	uStep[0] = 0.34; uStep[1] = 0.0;
	for(int i = 10; i < 20; ++i) {
		U[i] = uStep;
	}
	uStep[0] = 0.0; uStep[1] = -0.26;
	for(int i = 20; i < 30; ++i) {
		U[i] = uStep;
	}
#endif

	B.resize(pathlen+1);
	Matrix<X_DIM,X_DIM> SqrtW;

	vec(x0, SqrtSigma0, B[0]);
	for (int t = 0; t < pathlen; ++t) {
		beliefDynamics(B[t], U[t], B[t+1], SqrtW);
	}
	*/
}

void PointPlanner::initProblem()
{
	// Initialize with collision-free trajectory
	initTrajectory();
	displayMeanPath(B);

	L.assign(pathlen, zeros<U_DIM,B_DIM>());
	l.assign(pathlen, zeros<U_DIM, 1>());

	S = zeros<B_DIM, B_DIM>();
	s = zeros<B_DIM, 1>();
	ss = 0; 

	bestCost = COST_INFTY; 

	int num;
	std::cin >> num;
}

void PointPlanner::initPlanner()
{
#if defined(PROBLEM1)
	//Rint = identity<U_DIM>();
	//Qint = identity<X_DIM>();
	//QGoal = 300.0*identity<X_DIM>();
	//SqrtSigma0 = sqrt(0.75)*identity<X_DIM>();

	Rint = identity<U_DIM>();
	Qint = identity<X_DIM>();
	QGoal = 250.0*identity<X_DIM>();
	SqrtSigma0 = identity<X_DIM>();

	x0[0] = 2; x0[1] = 2;
	xGoal[0] = 0; xGoal[1] = 0;
	goalr = 0.25;

	nsteps = 1;
	invnsteps = 1.0/(double)nsteps;

	// Init environment
	CAL_Initialisation(true, true, true);

	// Set view params
	CAL_SetViewParams(0, 1.5, 2, 10.0, 1.5, 2, 0, 0, 1, 0);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");

	// Create bounding box
	CAL_CreateGroup(&cal_box, 0, false, "Box");
	//CAL_SetGroupColor(cal_box, 0.0f, 0.0f, 0.0f);
	//int np[1] = {2};
	//float side1[6] = {-2.0f, -5.0f, 0.0f, 8.0f, -5.0f, 0.0f};
	//float side2[6] = {8.0f, -5.0f, 0.0f, 8.0f, 5.0f, 0.0f};
	//float side3[6] = {8.0f, 5.0f, 0.0f, -2.0f, 5.0f, 0.0f};
	//float side4[6] = {-2.0f, 5.0f, 0.0f, -2.0f, -5.0f, 0.0f};
	//CAL_CreatePolyline(cal_box, 1, np, side1);
	//CAL_CreatePolyline(cal_box, 1, np, side2);
	//CAL_CreatePolyline(cal_box, 1, np, side3);
	//CAL_CreatePolyline(cal_box, 1, np, side4);
	
	// Visualize goal region
	CAL_CreateGroup(&cal_goal, 0, false, "Goal");
	CAL_CreateSphere(cal_goal, (float)goalr, (float)xGoal[0], (float)xGoal[1], 0.01f);
	CAL_SetGroupColor(cal_goal, 0.0f, 1.0f, 0.0f, 0.5f);
	CAL_SetGroupScaling(cal_goal, 1.0f, 1.0f, 0.01f);

	// Visualize sensing region
	int cal_sensor;
	CAL_CreateGroup(&cal_sensor, 0, false, "Sensor");
	
	int boxid;
	CAL_LoadTexture(0, "data\\platt-env.png", 0.8f);
	CAL_CreateBox(cal_sensor, 17.73f, 10.63f, 0.0f, 1.5f, 2.0f, -0.05f, &boxid);
	CAL_SetObjectTexture (boxid, 0, 1.0f, 1.0f);

	// RRT planning routines
	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 0, 1, 0);

	// Initialize with collision-free trajectory
	initTrajectory();
	displayMeanPath(B);

	L.assign(pathlen, zeros<U_DIM,B_DIM>());
	l.assign(pathlen, zeros<U_DIM, 1>());

	S = zeros<B_DIM, B_DIM>();
	s = zeros<B_DIM, 1>();
	ss = 0; 

	bestCost = COST_INFTY; 
#elif defined(PROBLEM2)
	/*
	nsteps = 1;
	invnsteps = 1.0/(double)nsteps;

	// obstacle
	Rint = identity<U_DIM>();
	//Qint = identity<X_DIM>();
	Qint = 5.0*identity<X_DIM>();
	QGoal = 250.0*identity<X_DIM>();
	x0[0] = 1.2; x0[1] = 1.8;
	xGoal[0] = 4.6; xGoal[1] = -2.6;
	goalr = 0.25;
	SqrtSigma0 = identity<X_DIM>(); //sqrt(0.75)*identity<X_DIM>();

	// Init environment
	CAL_Initialisation(true, true, true);

	// Set view params
	CAL_SetViewParams(0, 0, 0, 12, 0, 0, 0, 0, 1, 0);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");

	// Create bounding box
	CAL_CreateGroup(&cal_box, 0, false, "Box");
	CAL_SetGroupColor(cal_box, 0.0f, 0.0f, 0.0f);
	int np[1] = {2};
	float side1[6] = {-5.0f, -3.0f, 0.0f, 5.0f, -3.0f, 0.0f};
	float side2[6] = {5.0f, -3.0f, 0.0f, 5.0f, 3.0f, 0.0f};
	float side3[6] = {5.0f, 3.0f, 0.0f, -5.0f, 3.0f, 0.0f};
	float side4[6] = {-5.0f, 3.0f, 0.0f, -5.0f, -3.0f, 0.0f};
	CAL_CreatePolyline(cal_box, 1, np, side1);
	CAL_CreatePolyline(cal_box, 1, np, side2);
	CAL_CreatePolyline(cal_box, 1, np, side3);
	CAL_CreatePolyline(cal_box, 1, np, side4);

	// Create obstacles
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 0.93f, 0.11f, 0.14f);
	// Blocks
	CAL_CreateBox(cal_obstacles, 0.75f, 2.6f, 0.0f, 2.5f, 1.7f, 0.0f);
	CAL_CreateBox(cal_obstacles, 0.75f, 2.6f, 0.0f, 2.5f, -1.7f, 0.0f);
	// Bounding box proxy
	int boxid;
	CAL_CreateBox(cal_obstacles, 10.0f, 0.01f, 0.1f, 0.0f, 3.0f, 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 10.0f, 0.01f, 0.1f, 0.0f, -3.0f, 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 0.01f, 6.0f, 0.1f, 5.0f, 0.0f, 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 0.01f, 6.0f, 0.1f, -5.0f, 0.0f, 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.0f, 0.0f);

	// Visualize goal region
	CAL_CreateGroup(&cal_goal, 0, false, "Goal");
	CAL_CreateSphere(cal_goal, (float)goalr, (float)xGoal[0], (float)xGoal[1], 0.01f);
	CAL_SetGroupColor(cal_goal, 0.0f, 1.0f, 0.0f, 0.5f);
	CAL_SetGroupScaling(cal_goal, 1.0f, 1.0f, 0.01f);

	// Visualize sensing region
	int cal_sensor;
	CAL_CreateGroup(&cal_sensor, 0, false, "Sensor");
	//CAL_SetGroupColor(cal_sensor, 1.0f, 0.625f, 0.4765f, 0.5f);
	//CAL_CreateBox(cal_sensor, 0.1f, 6.0f, 0.0f, -2.0f, 0.0f, 0.0f);
	
	// Sensing environment texture?
	CAL_LoadTexture(0, "data\\sensing.png", 0.5f);
	CAL_CreateBox(cal_sensor, 10.0f, 6.0f, 0.0f, 0.0f, 0.0f, -0.05f, &boxid);
	CAL_SetObjectTexture (boxid, 0, 1.0f, 1.0f);

	// Create Point/Cylinder in Callisto for distance calculation
	CAL_CreateGroup(&cal_cylinder, 0, true, "Cylinder", false, CAL_FALSE);
	CAL_CreateCylinder(cal_cylinder, 0.001f, 1.0f, 0.0f, 0.0f, 0.0f);
	CAL_SetGroupColor(cal_cylinder, 0.0f, 0.0f, 1.0f);

	CAL_CreateGroup(&cal_point, 0, true, "Point", false, CAL_FALSE);
	CAL_CreateSphere(cal_point, 0.0f, 0.0f, 0.0f, 0.0f);

	// RRT planning routines
	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 0, 1, 0);

	// Initialize with collision-free trajectory
	initTrajectory();
	displayMeanPath(B);

	L.assign(pathlen, zeros<U_DIM,B_DIM>());
	l.assign(pathlen, zeros<U_DIM, 1>());

	S = zeros<B_DIM, B_DIM>();
	s = zeros<B_DIM, 1>();
	ss = 0; 

	bestCost = COST_INFTY; 
	*/

	// ML assumption example
	nsteps = 1;
	invnsteps = 1.0/(double)nsteps;

	// obstacle
	Rint = 3.0*identity<U_DIM>();
	//Qint = identity<X_DIM>();
	Qint = 1.0*identity<X_DIM>();
	QGoal = 5.0*identity<X_DIM>();
	x0[0] = 3.0; x0[1] = -2.0;
	xGoal[0] = 3.0; xGoal[1] = 2.0;
	goalr = 0.25;
	SqrtSigma0 = identity<X_DIM>(); //sqrt(0.75)*identity<X_DIM>();

	// Init environment
	CAL_Initialisation(true, true, true);

	// Set view params
	CAL_SetViewParams(0, 0, 0, 12, 0, 0, 0, 0, 1, 0);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");

	// Create bounding box
	CAL_CreateGroup(&cal_box, 0, false, "Box");
	CAL_SetGroupColor(cal_box, 0.0f, 0.0f, 0.0f);
	int np[1] = {2};
	float side1[6] = {-5.0f, -3.0f, 0.0f, 5.0f, -3.0f, 0.0f};
	float side2[6] = {5.0f, -3.0f, 0.0f, 5.0f, 3.0f, 0.0f};
	float side3[6] = {5.0f, 3.0f, 0.0f, -5.0f, 3.0f, 0.0f};
	float side4[6] = {-5.0f, 3.0f, 0.0f, -5.0f, -3.0f, 0.0f};
	CAL_CreatePolyline(cal_box, 1, np, side1);
	CAL_CreatePolyline(cal_box, 1, np, side2);
	CAL_CreatePolyline(cal_box, 1, np, side3);
	CAL_CreatePolyline(cal_box, 1, np, side4);

	// Create obstacles
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 0.93f, 0.11f, 0.14f);
	// Blocks
	//CAL_CreateBox(cal_obstacles, 2.0f, 0.5f, 0.0f, -4.0f, 0.0f, 0.0f);
	//CAL_CreateBox(cal_obstacles, 3.8f, 0.5f, 0.0f, 0.6f, 0.0f, 0.0f);

	//CAL_CreateBox(cal_obstacles, 2.2f, 0.5f, 0.0f, -3.9f, 0.0f, 0.0f);
	//CAL_CreateBox(cal_obstacles, 4.0f, 0.5f, 0.0f, 0.6f, 0.0f, 0.0f);
	//CAL_CreateBox(cal_obstacles, 1.5f, 0.5f, 0.0f, 4.25f, 0.0f, 0.0f);

	//CAL_CreateBox(cal_obstacles, 1.7f, 0.4f, 0.0f, -4.15f, 0.0f, 0.0f);
	//CAL_CreateBox(cal_obstacles, 5.0f, 0.4f, 0.0f, 0.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 7.5f, 0.4f, 0.0f, -1.25f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 1.7f, 0.4f, 0.0f, 4.15f, 0.0f, 0.0f);

	//CAL_CreateBox(cal_obstacles, 3.0f, 0.5f, 0.0f, -3.0f, 0.0f, 0.0f);

	// Bounding box proxy
	int boxid;
	CAL_CreateBox(cal_obstacles, 10.0f, 0.01f, 0.1f, 0.0f, 3.0f, 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 10.0f, 0.01f, 0.1f, 0.0f, -3.0f, 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 0.01f, 6.0f, 0.1f, 5.0f, 0.0f, 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 0.01f, 6.0f, 0.1f, -5.0f, 0.0f, 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.0f, 0.0f);

	// Visualize goal region
	CAL_CreateGroup(&cal_goal, 0, false, "Goal");
	CAL_CreateSphere(cal_goal, (float)goalr, (float)xGoal[0], (float)xGoal[1], 0.01f);
	CAL_SetGroupColor(cal_goal, 0.0f, 1.0f, 0.0f, 0.5f);
	CAL_SetGroupScaling(cal_goal, 1.0f, 1.0f, 0.01f);

	// Visualize sensing region
	int cal_sensor;
	CAL_CreateGroup(&cal_sensor, 0, false, "Sensor");
	//CAL_SetGroupColor(cal_sensor, 1.0f, 0.625f, 0.4765f, 0.5f);
	//CAL_CreateBox(cal_sensor, 0.1f, 6.0f, 0.0f, -2.0f, 0.0f, 0.0f);
	
	// Sensing environment texture?
	// Sensing environment texture?
	CAL_LoadTexture(0, "data\\sensingML.png", 0.6f);
	CAL_CreateBox(cal_sensor, 10.0f, 6.0f, 0.0f, 0.0f, 0.0f, -0.01f, &boxid);
	CAL_SetObjectTexture (boxid, 0, 1.0f, 1.0f);

	// Create Point/Cylinder in Callisto for distance calculation
	CAL_CreateGroup(&cal_cylinder, 0, true, "Cylinder", false, CAL_FALSE);
	CAL_CreateCylinder(cal_cylinder, 0.001f, 1.0f, 0.0f, 0.0f, 0.0f);
	CAL_SetGroupColor(cal_cylinder, 0.0f, 0.0f, 1.0f);

	CAL_CreateGroup(&cal_point, 0, true, "Point", false, CAL_FALSE);
	CAL_CreateSphere(cal_point, 0.0f, 0.0f, 0.0f, 0.0f);

	// RRT planning routines
	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 0, 1, 0);

	rrtfptr.open("data\\RRT-Paths-ML-NoML.txt", std::ios::in);
	if (!rrtfptr.is_open()) {
		std::cout << "RRT path file ptr not opened" << std::endl;
		std::exit(-1);
	}

	//resultsfptr.open("data\\Paths-NoML-collision.txt", std::ios::out);
	//if (!resultsfptr.is_open()) {
	//	std::cout << "Results file ptr not opened" << std::endl;
	//	std::exit(-1);
	//}

#endif
}