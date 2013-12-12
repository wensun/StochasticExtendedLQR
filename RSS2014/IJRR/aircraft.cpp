#include "aircraft.h"
const double DIST_INFTY = 1e10;

AircraftPlanner::AircraftPlanner() {
#ifndef AIRCRAFT3D
	std::cerr << "Improper setup. Define AIRCRAFT3D in planner.h" << std::endl;
	std::exit(-1);
#endif
}

// Dynamics model
Matrix<X_DIM> AircraftPlanner::f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m, double tstep) // m ~ N(0,I)
{  
	Matrix<X_DIM> xNew;
	
	double scale = sigmoid2(x[1], 8);

	//xNew[0] = x[0] + x[3]*tstep + 0.5*(u[0] + 0.25*u[0]*m[0])*tstep*tstep;
	//xNew[1] = x[1] + x[4]*tstep + 0.5*(u[1] + 0.25*u[1]*m[1])*tstep*tstep;
	//xNew[2] = x[2] + x[5]*tstep + 0.5*(u[2] + 0.25*u[2]*m[2])*tstep*tstep;

	xNew[0] = x[0] + x[3]*tstep + 0.5*(u[0] + scale*m[0])*tstep*tstep;
	xNew[1] = x[1] + x[4]*tstep + 0.5*(u[1] + scale*m[1])*tstep*tstep;
	xNew[2] = x[2] + x[5]*tstep + 0.5*(u[2] + scale*m[2])*tstep*tstep;
		
	xNew[3] = x[3] + (u[0] + scale*m[0])*tstep;
	xNew[4] = x[4] + (u[1] + scale*m[1])*tstep;
	xNew[5] = x[5] + (u[2] + scale*m[2])*tstep;

	return xNew;
}

// Observation model
Matrix<Z_DIM> AircraftPlanner::h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) // n ~ N(0,I)
{
	Matrix<Z_DIM> z;
		
	z[0] = x[0] + 0.2*n[0];
	z[1] = x[1] + 0.1*n[1];
	z[2] = x[2] + 0.2*n[2];

	return z;
}

void AircraftPlanner::RRT() 
{
	srand(time(0));

	std::vector<RRTNode> rrtTree;

	RRTNode startNode;
	startNode.x = x0;

	rrtTree.push_back(startNode);

	double goaldist = tr(~(startNode.x - xGoal) * (startNode.x - xGoal));

	while (goaldist > goalr * goalr * 0.5) 
	{
		Matrix<X_DIM> sample;
		if (random() < 0.2) { 
			sample = xGoal;
		} else {
			for (int i = 0; i < X_DIM; ++i) {
				sample[i] = random() * (xMax[i] - xMin[i]) + xMin[i];
			}
		}
		
		int col;

		CAL_CheckPointCollision(cal_environment, (float) sample[0], (float) sample[1], (float) sample[2], false, &col);
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
		if (xNew[0] < xMin[0] || xNew[0] > xMax[0] || xNew[1] < xMin[1] || xNew[1] > xMax[1] || xNew[2] < xMin[2] || xNew[2] > xMax[2]) {
			valid = false;
		}
		if (!valid) {
			continue;
		}

		CAL_CheckLineCollision(cal_environment, (float) rrtTree[node].x[0], (float) rrtTree[node].x[1], (float) rrtTree[node].x[2], 
												(float) xNew[0], (float) xNew[1], (float) xNew[2], false, &col);
		if (col != 0) {
			continue;
		}

		RRTNode newnode;
		newnode.x = xNew;
		newnode.u = input;
		newnode.bp = node;

		rrtTree.push_back(newnode);

		int np[1] = {2};
		float p[6] = {(float) rrtTree[node].x[0], (float) rrtTree[node].x[1], (float) rrtTree[node].x[2], (float) xNew[0], (float) xNew[1], (float) xNew[2]};
		CAL_CreatePolyline(cal_rrt, 1, np, p);

		goaldist = tr(~(newnode.x.subMatrix<DIM,1>(0,0) - xGoal.subMatrix<DIM,1>(0,0)) * (newnode.x.subMatrix<DIM,1>(0,0) - xGoal.subMatrix<DIM,1>(0,0)));
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
	
	//std::cout << "Finished RRT" << std::endl;
}

void AircraftPlanner::initTrajectory()
{

#if defined(PROBLEM1)
	std::string filename = "data\\rrtAircraft1.txt";
#elif defined(PROBLEM2)
	std::string filename = "data\\rrtAircraft2.txt";
#endif

	//RRT();

	//std::ofstream fptr(filename.c_str(),std::ios::out);
	//fptr << U.size() << std::endl;
	//for(int i = 0; i < (int)U.size(); ++i) {
	//	fptr << ~U[i];
	//}
	//fptr.close();

	std::ifstream fptr(filename.c_str(), std::ios::in);
	fptr >> pathlen;
	std::cout << "Path len: " << pathlen << std::endl;
	U.resize(pathlen);
	for(int i = 0; i < pathlen; ++i) {
		Matrix<U_DIM,1>& Ui = U[i];
		for(int j = 0; j < U_DIM; ++j) { fptr >> Ui[j]; }
	}
	fptr.close();

	B.resize(pathlen+1);
	Matrix<X_DIM,X_DIM> SqrtW;

	vec(x0, SqrtSigma0, B[0]);
	for (int t = 0; t < pathlen; ++t) {
		beliefDynamics(B[t], U[t], B[t+1], SqrtW);
	}
}

void AircraftPlanner::initPlanner()
{
#if defined(PROBLEM1)
	xMin[0] = 0.0; xMax[0] = 10.0;
	xMin[1] = 0.0; xMax[1] = 10.0;
	xMin[2] = 0.0; xMax[2] = 10.0;

	xMin[3] = 0.25; xMax[3] = 2.0;
	xMin[4] = 0.25; xMax[4] = 2.0;
	xMin[5] = 0.25; xMax[5] = 2.0;
	
	uMin[0] = -1.0; uMax[0] = 1.0;
	uMin[1] = -1.0; uMax[1] = 1.0;
	uMin[2] = -1.0; uMax[2] = 1.0;

	x0[0] = 5.0; x0[1] = 5.0; x0[2] = 1.0; x0[3] = 0.25; x0[4] = 0.25; x0[5] = 0.5;
	
	goalr = 0.4;
	
	nsteps = 1;
	invnsteps = 1.0/(double)nsteps;

	xGoal[0] = 5.0; xGoal[1] = 5.0; xGoal[2] = 9.0;
	xGoal[3] = xGoal[4] = xGoal[5] = 0.0;
	
	//Rint = identity<U_DIM>();
	//Qint = 10.0*identity<X_DIM>();
	//QGoal = 100.0*identity<X_DIM>();

	Rint = identity<U_DIM>();
	Qint = 10.0*identity<X_DIM>();
	QGoal = 200.0*identity<X_DIM>();
	
	QGoal(3,3) = QGoal(4,4) = QGoal(5,5) = 0.0;
	
	SqrtSigma0 = identity<X_DIM>() * sqrt(0.5);
	SqrtSigma0(3,3) = SqrtSigma0(4,4) = SqrtSigma0(5,5) = sqrt(0.25);
	
	CAL_Initialisation(true, true, true);

	// Set view params
	//CAL_SetViewParams(0, 23.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 0.0f, 1.0f, 0.0f);
	CAL_SetViewParams(0, 18.9f, 10.9f, 21.9f, 5.0f, 5.0f, 5.0f, 0.0f, 1.0f, 0.0f);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	
	CAL_CreateGroup (&cal_box, 0, false, "Box");
	CAL_SetGroupColor(cal_box, 0, 0, 0);
	float w = 0.05f;
	CAL_CreateBox(cal_box, (float)xMax[0], w, w, 0.5*(float)xMax[0], 0.0f, 0.0f);
	CAL_CreateBox(cal_box, (float)xMax[0], w, w, 0.5*(float)xMax[0], (float)xMax[1], 0.0f);
	CAL_CreateBox(cal_box, (float)xMax[0], w, w, 0.5*(float)xMax[0], (float)xMax[1], (float)xMax[2]);
	CAL_CreateBox(cal_box, (float)xMax[0], w, w, 0.5*(float)xMax[0], 0.0f, (float)xMax[2]);

	CAL_CreateBox(cal_box, w, (float)xMax[1], w, 0.0f, 0.5*(float)xMax[1], 0.0f);
	CAL_CreateBox(cal_box, w, (float)xMax[1], w, (float)xMax[0], 0.5*(float)xMax[1], 0.0f);
	CAL_CreateBox(cal_box, w, (float)xMax[1], w, (float)xMax[0], 0.5*(float)xMax[1], (float)xMax[2]);
	CAL_CreateBox(cal_box, w, (float)xMax[1], w, 0.0f, 0.5*(float)xMax[1], (float)xMax[2]);

	CAL_CreateBox(cal_box, w, w, (float)xMax[2], 0.0f, 0.0f, 0.5*(float)xMax[2]);
	CAL_CreateBox(cal_box, w, w, (float)xMax[2], (float)xMax[0], 0.0f, 0.5*(float)xMax[2]);
	CAL_CreateBox(cal_box, w, w, (float)xMax[2], (float)xMax[0], (float)xMax[1], 0.5*(float)xMax[2]);
	CAL_CreateBox(cal_box, w, w, (float)xMax[2], 0.0f, (float)xMax[1], 0.5*(float)xMax[2]);
	
	// Obstacles in the scene
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 0.93f, 0.11f, 0.14f);

	//CAL_CreateBox(cal_obstacles, (float)Xrange, (float)(Yrange * 0.25), 0.5f, (float)(0.5*Xrange), (float)(0.125*Yrange), (float)(0.75*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)Xrange, (float)(Yrange * 0.25), 0.5f, (float)(0.5*Xrange), (float)(0.875*Yrange), (float)(0.75*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.25), (float)Yrange, 0.5f, (float)(0.125*Xrange), (float)(0.5*Yrange), (float)(0.75*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.25), (float)Yrange, 0.5f, (float)(0.875*Xrange), (float)(0.5*Yrange), (float)(0.75*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.4), (float)(Yrange * 0.4), 0.5f, (float)(0.55*Xrange), (float)(0.55*Yrange), (float)(0.75*Zrange));

	//CAL_CreateBox(cal_obstacles, 1.0f, 2.5f, 1.0f, 2.0f, 4.0f, 5.0f);
	//CAL_CreateBox(cal_obstacles, 2.5f, 1.0f, 1.0f, 5.0f, 5.0f, 8.0f);

	//CAL_CreateCylinder(cal_obstacles, 0.7f, 10.0f, 5.0f, 5.0f, 7.0f);
	CAL_CreateCylinder(cal_obstacles, 0.7f, 10.0f, 5.0f, 5.0f, 3.5f);

	//CAL_CreateBox(cal_obstacles, 2.0f, 1.0f, 1.0f, 5.0f, 5.0f, 4.0f);
	int obsid;
	CAL_CreateCylinder(cal_obstacles, 0.7f, 10.0f, 5.0f, 5.0f, 7.0f, &obsid);
	CAL_SetObjectOrientation(obsid, 0.0f, 0.0f, (float)(0.5*M_PI));

	int cal_obstacles2;
	CAL_CreateGroup(&cal_obstacles2, cal_obstacles, true, "", false, CAL_FALSE);
	CAL_CreateBox(cal_obstacles2, 10.0f, 1.0f, 10.0f, 5.0f, 10.0f, 5.0f);
	CAL_CreateBox(cal_obstacles2, 10.0f, 0.01f, 10.0f, 5.0f, 0.0f, 5.0f);
	CAL_CreateBox(cal_obstacles2, 0.01f, 10.0f, 10.0f, 10.0f, 5.0f, 5.0f);
	CAL_CreateBox(cal_obstacles2, 0.01f, 10.0f, 10.0f, 0.0f, 5.0f, 5.0f);
	CAL_CreateBox(cal_obstacles2, 10.0f, 10.0f, 0.01f, 5.0f, 5.0f, 10.0f);
	CAL_CreateBox(cal_obstacles2, 10.0f, 10.0f, 0.01f, 5.0f, 5.0f, 0.0f);
	
	// Visualize goal region
	CAL_CreateGroup(&cal_goal, 0, false, "Goal");
	//CAL_CreateSphere(cal_goal, (float)goalr, (float)xGoal[0], (float)xGoal[1], (float)xGoal[2]);
	int id;
	CAL_CreateUserDrawn(cal_goal, drawUnitSphere, NULL, (float)xGoal[0], (float)xGoal[1], (float)xGoal[2], &id);
	CAL_SetObjectScaling(id, (float)goalr, (float)goalr, (float)goalr);
		
	CAL_SetGroupColor(cal_goal, 0.0f, 1.0f, 0.0f, 0.5f);
	
	// Visualize sensing region
	int cal_sensor;
	CAL_CreateGroup(&cal_sensor, 0, false, "Sensor");
	CAL_CreateBox(cal_sensor, 10.0f, 2.0f, 10.0f, 5.0f, 9.0f, 5.0f, &id);
	//CAL_CreateBox(cal_sensor, 10.0f, 2.0f, 10.0f, 5.0f, 1.0f, 5.0f, &id);
	//CAL_SetObjectColor(id, 1.0f, 0.76f, 0.6f, 0.5f);
	CAL_SetObjectColor(id, 1.0f, 0.8f, 0.2f, 0.3f);

	// RRT planning routines
	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 0, 1, 0);

	// Create Point/Cylinder in Callisto for distance calculation
	CAL_CreateGroup(&cal_cylinder, 0, true, "Cylinder", false, CAL_FALSE);
	CAL_CreateCylinder(cal_cylinder, 0.001f, 1.0f, 0.0f, 0.0f, 0.0f);
	CAL_SetGroupColor(cal_cylinder, 0.0f, 0.0f, 1.0f);

	CAL_CreateGroup(&cal_point, 0, true, "Point");
	CAL_CreateSphere(cal_point, 0, 0, 0, 0);

	//CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	//CAL_SetGroupColor(cal_ellipse, 0.2, 0.2, 0.2, 0.6f);

	// Initialize with collision-free trajectory
	initTrajectory();

	displayMeanPath(B);

	L.assign(pathlen, zeros<U_DIM,B_DIM>());
	l.assign(pathlen, zeros<U_DIM, 1>());

	S = zeros<B_DIM, B_DIM>();
	s = zeros<B_DIM, 1>();
	ss = 0; 
#elif defined(PROBLEM2)
	std::cerr << "PROBLEM2 not defined for Aircraft3D." << std::endl;
	std::exit(-1);
#endif
}