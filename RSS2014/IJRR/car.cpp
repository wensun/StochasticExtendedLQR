#include "car.h"

CarPlanner::CarPlanner() {
#ifndef CAR2D
	std::cerr << "Improper setup. Define CAR2D in planner.h" << std::endl;
	std::exit(-1);
#endif
}

// Dynamics model
Matrix<X_DIM> CarPlanner::f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m, double tstep) // m ~ N(0,I)
{  
	Matrix<X_DIM> xNew;

	/* Euler integration
	xNew = x;
	xNew[0] += tstep*x[3]*cos(x[2]);
	xNew[1] += tstep*x[3]*sin(x[2]);
	xNew[2] += tstep*x[3]*tan(u[1])/car_l;
	xNew[3] += tstep*(u[0]);
	*/

	// RK4 integration
	Matrix<X_DIM> x1, x2, x3, x4, xtmp;

	xtmp = x;

	x1[0] = xtmp[3]*cos(xtmp[2]); 
	x1[1] = xtmp[3]*sin(xtmp[2]); 
	x1[2] = xtmp[3]*tan(u[1])/car_l; 
	x1[3] = u[0];

	xtmp = x + 0.5*tstep*x1;

	x2[0] = xtmp[3]*cos(xtmp[2]); 
	x2[1] = xtmp[3]*sin(xtmp[2]); 
	x2[2] = xtmp[3]*tan(u[1])/car_l; 
	x2[3] = u[0];

	xtmp = x + 0.5*tstep*x2;

	x3[0] = xtmp[3]*cos(xtmp[2]); 
	x3[1] = xtmp[3]*sin(xtmp[2]); 
	x3[2] = xtmp[3]*tan(u[1])/car_l; 
	x3[3] = u[0];

	xtmp = x + tstep*x3;

	x4[0] = xtmp[3]*cos(xtmp[2]); 
	x4[1] = xtmp[3]*sin(xtmp[2]); 
	x4[2] = xtmp[3]*tan(u[1])/car_l; 
	x4[3] = u[0];

	xNew = x + (tstep/6.0)*(x1 + 2.0*(x2 + x3) + x4);

#if defined(PROBLEM1)
	xNew[0] += 0.01*m[0];
	xNew[1] += 0.01*m[1];
	xNew[2] += 0.02*m[2];
	xNew[3] += 0.02*m[3];
#elif defined(PROBLEM2)
	xNew[0] += 0.1*abs(u[1])*sqrt(DT)*m[0];
	xNew[1] += 0.1*abs(u[1])*sqrt(DT)*m[1];

	xNew[2] += 0.1*abs(u[1])*sqrt(DT)*m[2];
	xNew[3] += 0.1*abs(u[0])*sqrt(DT)*m[3];

	//xNew[0] += 0.05*m[0];
	//xNew[1] += 0.05*m[1];
	//xNew[2] += 0.05*m[2];
	//xNew[3] += 0.05*m[3];
#endif

	return xNew;
}

// Observation model
Matrix<Z_DIM> CarPlanner::h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) // n ~ N(0,I)
{
	Matrix<Z_DIM> z;

#if defined(PROBLEM1)
	double noise = sqrt(0.5*0.5*(5-x[0])*(5-x[0]) + 0.001);
	z[0] = x[0] + noise*n[0];
	z[1] = x[1] + noise*n[1];
	//z[2] = x[3] + 0.01*n[2];
#elif defined(PROBLEM2)
	z[0] = 0.5/(sqr(x[0] - b0[0]) + sqr(x[1] - b0[1]) + 1.0) + 0.075*n[0];
	z[1] = 0.5/(sqr(x[0] - b1[0]) + sqr(x[1] - b1[1]) + 1.0) + 0.075*n[1];
	
	//z[2] = x[3] + 0.01*n[2];

	z[2] = x[2] + 0.025*n[2];
	z[3] = x[3] + 0.025*n[3];
#endif

	return z;
}

void CarPlanner::RRT() 
{
	srand(time(0));

	std::vector<RRTNode> rrtTree;

	RRTNode startNode;
	startNode.x = x0;

	rrtTree.push_back(startNode);

	double goaldist = tr(~(startNode.x.subMatrix<DIM,1>(0,0) - xGoal.subMatrix<DIM,1>(0,0)) * (startNode.x.subMatrix<DIM,1>(0,0) - xGoal.subMatrix<DIM,1>(0,0)));

	while (goaldist > goalr * goalr * 0.25) 
	{
		Matrix<X_DIM> sample;
		for (int i = 0; i < X_DIM; ++i) {
			sample[i] = random() * (xMax[i] - xMin[i]) + xMin[i];
		}
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
		for (int i = 0; i < U_DIM; ++i) {
			input[i] = random() * (uMax[i] - uMin[i]) + uMin[i];
		}

		Matrix<X_DIM> xNew = f(rrtTree[node].x, input, zeros<M_DIM,1>());

		bool valid = true;
		for (int i = 0; i < X_DIM; ++i) {
			if (xNew[i] < xMin[i] || xNew[i] > xMax[i]) {
				valid = false;
				break;
			}
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

		goaldist = tr(~(newnode.x.subMatrix<DIM,1>(0,0) - xGoal.subMatrix<DIM,1>(0,0)) * (newnode.x.subMatrix<DIM,1>(0,0) - xGoal.subMatrix<DIM,1>(0,0)));
	}

	std::deque<PathNode> path;

	int i = (int) rrtTree.size() - 1;
	PathNode node;
	node.x = rrtTree[i].x;
	node.u = zeros<U_DIM, 1>();
	path.push_front(node);

	while (i != 0) {
		node.u = rrtTree[i].u;
		i = rrtTree[i].bp;
		node.x = rrtTree[i].x;

		path.push_front(node);
	}

	//std::cout << path.size() << std::endl;
	for (int i = 0; i < (int)path.size()-1; ++i) {
		U.push_back(path[i].u);
		//std::cout << ~path[i].x;
		//std::cout << ~path[i].u;
	}

	// Set pathlen based on solution size
	pathlen = path.size()-1;

	// Clear rrt tree visualization
	CAL_DestroyGroup(cal_rrt);
}

void CarPlanner::initTrajectory()
{

#if defined(PROBLEM1)
	std::string filename = "data\\rrtCar1.txt";
#elif defined(PROBLEM2)
	std::string filename = "data\\rrtCar2Bottom.txt";
#endif

	//RRT();

	//std::ofstream fptr(filename.c_str(),std::ios::out);
	//fptr << U.size()  << std::endl;
	//for(int i = 0; i < (int)U.size(); ++i) {
	//	fptr << ~U[i];
	//}
	//fptr.close();

	std::ifstream fptr(filename.c_str(), std::ios::in);
	fptr >> pathlen;
	std::cout << "Path len: " << pathlen << std::endl;
	U.resize(pathlen);
	for(int i = 0; i < (int)U.size(); ++i) {
		Matrix<U_DIM,1>& Ui = U[i];
		fptr >> Ui[0] >> Ui[1];
	}
	fptr.close();

	B.resize(pathlen+1);
	Matrix<X_DIM,X_DIM> SqrtW;

	vec(x0, SqrtSigma0, B[0]);
	for (int t = 0; t < pathlen; ++t) {
		beliefDynamics(B[t], U[t], B[t+1], SqrtW);
	}
}

void CarPlanner::initPlanner()
{
#if defined(PROBLEM1)
	// Init variables
	// rrtCar1Good
	//Rint = identity<U_DIM>();
	//Qint = 3.0*identity<X_DIM>();
	//QGoal = 10.0*identity<X_DIM>();
	//QGoal(2,2) = QGoal(3,3) = 0.0;

	Rint = 0.5*identity<U_DIM>();
	Qint = 7.0*identity<X_DIM>();
	QGoal = 10.0*identity<X_DIM>();
	QGoal(2,2) = QGoal(3,3) = 0.0;

	x0[0] = 2.0; x0[1] = 2.0; x0[2] = 1.75*M_PI; x0[3] = 0.0;
	
	SqrtSigma0 = identity<X_DIM>();
	SqrtSigma0(0,0) = SqrtSigma0(1,1) = sqrt(0.025);
	SqrtSigma0(2,2) = SqrtSigma0(3,3) = sqrt(0.01);

	xGoal[0] = 0; xGoal[1] = 0; xGoal[2] = 0; xGoal[3] = 0;
	
	// car specific
	goalr = 0.25;
	car_l = 0.5;

	xMin[0] = -3; xMin[1] = -0; xMin[2] = 0.0; xMin[3] = 0;
	xMax[0] = 5;  xMax[1] = 2;  xMax[2] = 2.0*M_PI;  xMax[3] = 2;
	
	uMin[0] = -1.0; uMin[1] = -M_PI*0.25;
	uMax[0] = 1.0;  uMax[1] = M_PI*0.25;

	nsteps = 1;
	invnsteps = 1.0/(double)nsteps;

	CAL_Initialisation(true, true, true);

	// Set view params
	CAL_SetViewParams(0, 1.0, 1.0, 7.75, 1.0, 1.0, 0, 0, 1, 0);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");

	// Visualize goal region
	CAL_CreateGroup(&cal_goal, 0, false, "Goal");
	CAL_CreateSphere(cal_goal, (float)goalr, (float)xGoal[0], (float)xGoal[1], 0.0f);
	CAL_SetGroupColor(cal_goal, 0.0f, 0.8f, 0.0f, 0.7f);
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

	xGoal = B[pathlen].subMatrix<X_DIM,1>(0,0);

	L.assign(pathlen, zeros<U_DIM,B_DIM>());
	l.assign(pathlen, zeros<U_DIM, 1>());

	S = zeros<B_DIM, B_DIM>();
	s = zeros<B_DIM, 1>();
	ss = 0; 
#elif defined(PROBLEM2)
	// Init variables
	// bottom?
	//Rint = 0.5*identity<U_DIM>();
	//Qint = 3.0*identity<X_DIM>();
	//QGoal = 100.0*identity<X_DIM>();
	//QGoal(2,2) = QGoal(3,3) = 0.0;

	//bottom
	Rint = 0.5*identity<U_DIM>();
	Qint = 3.0*identity<X_DIM>();
	QGoal = 100.0*identity<X_DIM>();
	QGoal(2,2) = QGoal(3,3) = 0.0;

	//Rint = identity<U_DIM>();
	//Qint = 3.0*identity<X_DIM>();
	//QGoal = 200.0*identity<X_DIM>();
	//QGoal(2,2) = QGoal(3,3) = 0.0;

	x0[0] = -4.0; x0[1] = 0.0; x0[2] = 0.0; x0[3] = 0.2;

	SqrtSigma0 = identity<X_DIM>();
	//SqrtSigma0(0,0) = SqrtSigma0(1,1) = sqrt(0.025);
	//SqrtSigma0(2,2) = SqrtSigma0(3,3) = sqrt(0.01);

	SqrtSigma0(0,0) = SqrtSigma0(1,1) = sqrt(0.5);
	SqrtSigma0(2,2) = SqrtSigma0(3,3) = sqrt(0.25);

	xGoal[0] = 4.0; xGoal[1] = 0.0; xGoal[2] = 0.0; xGoal[3] = 0.0;

	// car specific
	goalr = 0.25;
	car_l = 0.5;

	xMin[0] = -5; xMin[1] = -3; xMin[2] = -M_PI; xMin[3] = 0;
	xMax[0] = 5;  xMax[1] = 3;  xMax[2] = M_PI;  xMax[3] = 2;

	// beacons
	b0[0] = -2.25; b0[1] = 2.25;
	b1[0] = -2.25; b1[1] = -2.25;
	
	uMin[0] = -1.0; uMin[1] = -M_PI*0.25;
	uMax[0] = 1.0;  uMax[1] = M_PI*0.25;

	nsteps = 1;
	invnsteps = 1.0/(double)nsteps;

	CAL_Initialisation(true, true, true);

	// Set view params
	CAL_SetViewParams(0, 0, 0, 7.75, 0, 0, 0, 0, 1, 0);

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

	// Blocks
	CAL_CreateBox(cal_obstacles, 2.75f, 1.5f, 0.1f, 1.0f, 1.25f, 0.0f);
	CAL_CreateBox(cal_obstacles, 0.8f, 1.8f, 0.1f, 1.0f, -1.25f, 0.0f);
	
	// For RRT only
	CAL_CreateGroup(&cal_obstacles2, cal_environment, true, "Obstacles2");
	CAL_CreateBox(cal_obstacles2, 3.25f, 1.6f, 0.1f, 1.0f, 1.25f, 0.0f);
	CAL_CreateBox(cal_obstacles2, 1.0f, 1.9f, 0.1f, 1.0f, -1.25f, 0.0f);

	// Visualize goal region
	CAL_CreateGroup(&cal_goal, 0, false, "Goal");
	CAL_CreateSphere(cal_goal, (float)goalr, (float)xGoal[0], (float)xGoal[1], 0.0f);
	CAL_SetGroupColor(cal_goal, 0.0f, 1.0f, 0.0f, 0.7f);
	CAL_SetGroupScaling(cal_goal, 1.0f, 1.0f, 0.01f);
	
	// Visualize sensing region
	int cal_sensor;
	CAL_CreateGroup(&cal_sensor, 0, false, "Sensor");

	CAL_CreateBox(cal_sensor, 0.2f, 0.2f, 0.0f, (float)b0[0], (float)b0[1], 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.44f, 1.0f);
	
	CAL_CreateBox(cal_sensor, 0.2f, 0.2f, 0.0f, (float)b1[0], (float)b1[1], 0.0f, &boxid);
	CAL_SetObjectColor(boxid, 0.0f, 0.44f, 1.0f);
	
	// Sensing environment texture?
	CAL_LoadTexture(0, "data\\sensing2.png", 0.6f);
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

	xGoal = B[pathlen].subMatrix<X_DIM,1>(0,0);

	L.assign(pathlen, zeros<U_DIM,B_DIM>());
	l.assign(pathlen, zeros<U_DIM, 1>());

	S = zeros<B_DIM, B_DIM>();
	s = zeros<B_DIM, 1>();
	ss = 0; 
#endif
}