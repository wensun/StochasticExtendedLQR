#include "aircraft.h"
const double DIST_INFTY = 1e10;

AircraftPlanner::AircraftPlanner() {
#ifndef AIRCRAFT3D
	std::cerr << "Improper setup. Define AIRCRAFT3D in planner.h" << std::endl;
	std::exit(-1);
#endif
}

Matrix<4,4> AircraftPlanner::stepT(const Matrix<4,4>& T, const Matrix<3>& v, const Matrix<3>& w, double dt) {
	Matrix<4,4> U = zeros<4,4>();
	U.insert(0,0, cpMatrix(w));
	U.insert(0,3, v);
	return T*exp(dt*U);
}

void AircraftPlanner::createInputs(Matrix<3>& v, Matrix<3>& w, double speed, double twist, double curvature) 
{
	v(0,0) = 0;
	v(1,0) = 0;
	v(2,0) = speed; 

	w(0,0) = speed * curvature;
	w(1,0) = 0;
	w(2,0) = twist;  
}

// Dynamics model
Matrix<X_DIM> AircraftPlanner::f(const Matrix<X_DIM>& x, const Matrix<U_DIM>& u, const Matrix<M_DIM>& m, double tstep) // m ~ N(0,I)
{  
	Matrix<X_DIM> xNew;

	Matrix<3,1> v;
	Matrix<3,1> w;
	//std::cout << "Controls: " << u[0] << " " << u[1] << " " << u[2] << std::endl;

	createInputs(v, w, u[0], u[1], u[2]);

	v += m.subMatrix<3,1>(0,0);
	w += m.subMatrix<3,1>(3,0);

	Matrix<4,4> T;
	T.insert(0,3, x.subMatrix<3,1>(0,0));
	T.insert(0,0, rotFromQuat(x.subMatrix<4,1>(3,0)));
	T(3,0) = T(3,1) = T(3,2) = 0.0; T(3,3) = 1.0;

	T = stepT(T, v, w, tstep);

	xNew.insert(0,0, T.subMatrix<3,1>(0,3));
	xNew.insert(3,0, quatFromRot(T.subMatrix<3,3>(0,0)));

	return xNew;
}

// Observation model
Matrix<Z_DIM> AircraftPlanner::h(const Matrix<X_DIM>& x, const Matrix<N_DIM>& n) // n ~ N(0,I)
{
	Matrix<Z_DIM> z;

	z[0] = x[0] + 0.25*n[0];
	z[1] = x[1] + 0.15*n[1];
	z[2] = x[2] + 0.2*n[2];

	return z;
}

double AircraftPlanner::dist(const Matrix<4,4>& T, const Matrix<3>& point) 
{
	Matrix<3> proj_point = ~T.subMatrix<3,3>(0,0) * (point - T.subMatrix<3,1>(0,3));
	double y = proj_point(2,0);
	if (y > 0) { 
		return DIST_INFTY;
	}
	double x = sqrt(proj_point(0,0)*proj_point(0,0) + proj_point(1,0)*proj_point(1,0));
	if (x == 0) {
		return -y;
	}
	double r = (x*x + y*y) / (2*x);
	if (r < (1.0/uMax[2])) {
		return DIST_INFTY;
	}
	//return sqrt(x*x + y*y);
	return atan2(-y, r-x) * r;
}

int AircraftPlanner::nearestNeighbor(const Matrix<3>& point, const std::vector<RRTNode>& rrtTree) 
{
	int closest = -1;
	double mindist = DIST_INFTY;
	int col = 1;
	CAL_CheckLineCollision(cal_obstacles, (float)point(0,0), (float)point(1,0), (float)point(2,0), 
						  (float)rrtTree[0].T(0,3), (float)rrtTree[0].T(1,3), (float)rrtTree[0].T(2,3), false, &col);
	if (col == 0) {
		closest = 0;
		mindist = factor*sqrt(tr(~(rrtTree[0].T.subMatrix<3,1>(0,3) - point) * (rrtTree[0].T.subMatrix<3,1>(0,3) - point)));
	} 

	std::stack<int> st;
	for (int i = 0; i < (int) rrtTree[0].children.size(); ++i) {
		st.push(rrtTree[0].children[i]);
	}

	while (!st.empty()) {
		int i = st.top();
		st.pop();
		double d = dist(rrtTree[i].T, point);
		if (rrtTree[i].T(2,3) < 0) { // do not expand nodes outside entry plane
			d = DIST_INFTY;
		}
		if (d < DIST_INFTY) {
			if (factor*d + rrtTree[i].depth * uMax[0] * DT < mindist) {
				col = 1;
				CAL_CheckLineCollision(cal_obstacles, (float)point(0,0), (float)point(1,0), (float)point(2,0), 
					(float)rrtTree[i].T(0,3), (float)rrtTree[i].T(1,3), (float)rrtTree[i].T(2,3), 
					false, &col);
				if (col == 0) {
					closest = i;
					mindist = factor*d + rrtTree[i].depth * uMax[0] * DT;
				}
			}
			for (int c = 0; c < (int) rrtTree[i].children.size(); ++c) {
				st.push(rrtTree[i].children[c]);
			}
		}
	}

	return closest;
}

void AircraftPlanner::RRT() 
{
	srand(time(0));

	std::vector<RRTNode> rrtTree;

	RRTNode n;
	n.T = identity<4>();
	n.T.insert(0,3, xGoal.subMatrix<3,1>(0,0));
	n.x.insert(0,0, xGoal.subMatrix<3,1>(0,0));
	n.x[3] = n.x[4] = n.x[5] = 0.0; n.x[6] = 1.0;
	n.u[0] = 0.0; n.u[1] = 0.0; n.u[2] = uMax[2];
	n.bp = -1;
	n.depth = -1;
	rrtTree.push_back(n);

	bool found = false;

	int node;
	Matrix<3> v; 
	Matrix<3> w;

	while (!found)  
	{
		Matrix<3> point; 
		do {
			if (random() < 0.2) {
				point(0,0) = random(4.0, 6.0); point(1,0) = random(4.0, 6.0); point(2,0) = -1.0;
			} else {
				point(0,0) = Xrange * random(); point(1,0) = Yrange * random(); point(2,0) = Zrange * random();
			}
			node = nearestNeighbor(point, rrtTree);
		} while (node == -1);
	
		if (node == 0) {
			RRTNode goalnode;
			goalnode.bp = 0;
			goalnode.depth = 0;
			goalnode.T = rrtTree[0].T;
			goalnode.u[1] = 0;

			Matrix<3> z = (rrtTree[0].T.subMatrix<3,1>(0,3) - point); 
			z = z / sqrt(tr(~z*z)); 
			Matrix<3> random_point; 
			do {
				random_point(0,0) = random(); random_point(1,0) = random(); random_point(2,0) = random();
			} while (tr(~random_point * random_point) > 1);
			Matrix<3> y = random_point - z * tr(~random_point * z);
			y = y / sqrt(tr(~y*y));
			Matrix<3> x = cpMatrix(y) * z;

			goalnode.T.insert(0,0,x);
			goalnode.T.insert(0,1,y);
			goalnode.T.insert(0,2,z);

			goalnode.x.insert(0,0, goalnode.T.subMatrix<3,1>(0,3));
			goalnode.x.insert(3,0, quatFromRot(goalnode.T.subMatrix<3,3>(0,0)));

			node = (int)rrtTree.size();
			rrtTree[0].children.push_back(node);
			rrtTree.push_back(goalnode);
		}
		
		createInputs(v, w, uMax[0], (2.0 * uMax[1] * random()) - uMax[1], uMax[2]);
		
		RRTNode newnode;
		newnode.bp = node;
		newnode.depth = rrtTree[node].depth + uMax[0]*DT;
		newnode.T = stepT(rrtTree[node].T, v, w, -DT);

		newnode.x.insert(0,0, newnode.T.subMatrix<3,1>(0,3));
		newnode.x.insert(3,0, quatFromRot(newnode.T.subMatrix<3,3>(0,0)));
		
		newnode.u[0] = uMax[0];
		newnode.u[1] = w(2,0);
		newnode.u[2] = w(0,0)/uMax[0];

		float line[6] = {(float) newnode.T(0,3), (float) newnode.T(1,3), (float) newnode.T(2,3), 
						 (float)rrtTree[newnode.bp].T(0,3), (float)rrtTree[newnode.bp].T(1,3), (float)rrtTree[newnode.bp].T(2,3)};
		int col = 1;
		CAL_CheckLineCollision(cal_obstacles, line[0], line[1], line[2], line[3], line[4], line[5], false, &col);
		if (col != 0) { // collision
			if (rrtTree[newnode.bp].bp == 0) {
				rrtTree.pop_back();
				rrtTree[0].children.pop_back();
			}
		} else {
			int np[1] = {2};
			CAL_CreatePolyline(cal_rrt, 1, np, line);

			rrtTree[node].children.push_back((int)rrtTree.size());
			rrtTree.push_back(newnode);

			if (newnode.T(2,3) < 0.0 && (newnode.T(1,3) < 5.5 && newnode.T(1,3) > 4.5) && (newnode.T(0,3) < 5.5 && newnode.T(0,3) > 4.5)) {
				found = true;
			}
		}
	}

	std::vector<PathNode> path;

	int i = (int) rrtTree.size() - 1;
	x0 = rrtTree[i].x;

	while (i != 0) {
		PathNode stage;

		stage.T = rrtTree[i].T;
		stage.x = rrtTree[i].x;
		stage.u = rrtTree[i].u;

		path.push_back(stage);

		i = rrtTree[i].bp;
	}

	// Set pathlen based on solution size
	pathlen = (int)(path.size()-1);

	//std::ofstream fptr("data\\aircraft.txt", std::ios::out);
	//for (int i = 0; i < pathlen; ++i) {
	//	fptr << path[i].T;
	//}
	//fptr.close();

	//std::cout << path.size() << std::endl;
	for (int i = 0; i < pathlen; ++i) {
		U.push_back(path[i].u);
	}

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
	//fptr << ~x0;
	//for(int i = 0; i < (int)U.size(); ++i) {
	//	fptr << ~U[i];
	//}
	//fptr.close();

	std::ifstream fptr(filename.c_str(), std::ios::in);
	fptr >> pathlen;
	std::cout << "Path len: " << pathlen << std::endl;
	for(int i = 0; i < X_DIM; ++i) { fptr >> x0[i]; }
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
	Xrange = 10.0; Yrange = 10.0; Zrange = 10.0;

	uMin[0] = 0.0; uMax[0] = 1.0;
	uMin[1] = 0.0; uMax[1] = 2.0*M_PI;
	uMin[2] = 0.0; uMax[2] = 0.25;
	
	goalr = 0.5;
	plane_l = 0.5;
	factor = 1.3;
	
	nsteps = 1;
	invnsteps = 1.0/(double)nsteps;

	for(int i = 0; i < X_DIM; ++i) { x0[i] = 0.0; }

	xGoal[0] = 0.5 * Xrange; xGoal[1] = 0.5 * Yrange; xGoal[2] = 1.0 * Zrange;
	xGoal[3] = xGoal[4] = xGoal[5] = 0.0; xGoal[6] = 1.0;
	
	Rint = identity<U_DIM>();
	Qint = identity<X_DIM>();
	QGoal = 100.0*identity<X_DIM>();
	QGoal(3,3) = QGoal(4,4) = QGoal(5,5) = QGoal(6,6) = 0.0;
	
	SqrtSigma0 = identity<X_DIM>() * sqrt(0.1);
	SqrtSigma0(3,3) = SqrtSigma0(4,4) = SqrtSigma0(5,5) = SqrtSigma0(6,6) = sqrt(0.25);
	
	CAL_Initialisation(true, true, true);

	// Set view params
	CAL_SetViewParams(0, 23.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 0.0f, 1.0f, 0.0f);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	
	CAL_CreateGroup (&cal_box, 0, false, "Box");
	CAL_SetGroupColor(cal_box, 0, 0, 0);
	float w = 0.05f;
	CAL_CreateBox(cal_box, (float)Xrange, w, w, 0.5*(float)Xrange, 0.0f, 0.0f);
	CAL_CreateBox(cal_box, (float)Xrange, w, w, 0.5*(float)Xrange, (float)Yrange, 0.0f);
	CAL_CreateBox(cal_box, (float)Xrange, w, w, 0.5*(float)Xrange, (float)Yrange, (float)Zrange);
	CAL_CreateBox(cal_box, (float)Xrange, w, w, 0.5*(float)Xrange, 0.0f, (float)Zrange);

	CAL_CreateBox(cal_box, w, (float)Yrange, w, 0.0f, 0.5*(float)Yrange, 0.0f);
	CAL_CreateBox(cal_box, w, (float)Yrange, w, (float)Xrange, 0.5*(float)Yrange, 0.0f);
	CAL_CreateBox(cal_box, w, (float)Yrange, w, (float)Xrange, 0.5*(float)Yrange, (float)Zrange);
	CAL_CreateBox(cal_box, w, (float)Yrange, w, 0.0f, 0.5*(float)Yrange, (float)Zrange);

	CAL_CreateBox(cal_box, w, w, (float)Zrange, 0.0f, 0.0f, 0.5*(float)Zrange);
	CAL_CreateBox(cal_box, w, w, (float)Zrange, (float)Xrange, 0.0f, 0.5*(float)Zrange);
	CAL_CreateBox(cal_box, w, w, (float)Zrange, (float)Xrange, (float)Yrange, 0.5*(float)Zrange);
	CAL_CreateBox(cal_box, w, w, (float)Zrange, 0.0f, (float)Yrange, 0.5*(float)Zrange);
	
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

	CAL_CreateCylinder(cal_obstacles, 0.7f, 6.0f, 5.0f, 5.0f, 7.0f);
	int obsid;
	CAL_CreateCylinder(cal_obstacles, 0.7f, 6.0f, 5.0f, 5.0f, 3.0f, &obsid);
	CAL_SetObjectOrientation(obsid, 0.0f, 0.0f, (float)(0.5*M_PI));

	// Visualize goal region
	CAL_CreateGroup(&cal_goal, 0, false, "Goal");
	CAL_CreateSphere(cal_goal, (float)goalr, (float)xGoal[0], (float)xGoal[1], (float)xGoal[2]);
	CAL_SetGroupColor(cal_goal, 0.0f, 1.0f, 0.0f, 0.5f);
	
	// Visualize sensing region
	int cal_sensor;
	CAL_CreateGroup(&cal_sensor, 0, false, "Sensor");

	// RRT planning routines
	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 0, 1, 0);

	CAL_CreateGroup(&cal_point, 0, true, "Point");
	CAL_CreateSphere(cal_point, 0, 0, 0, 0);

	//CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	//CAL_SetGroupColor(cal_ellipse, 0.2, 0.2, 0.2, 0.6f);

	// Initialize with collision-free trajectory
	initTrajectory();

	displayMeanPath(B);

	//x0 = B[0].subMatrix<X_DIM,1>(0,0);

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