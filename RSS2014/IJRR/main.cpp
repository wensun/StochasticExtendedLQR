#include "point.h"
#include "car.h"
#include "aircraft.h"

int main()
{
	// Floating point errors
	//unsigned int fp_control_word, new_fp_control_word;
	//_controlfp_s(&fp_control_word, 0, 0);
	////new_fp_control_word = fp_control_word & ~(_EM_INVALID | _EM_DENORMAL | _EM_ZERODIVIDE | _EM_OVERFLOW | _EM_UNDERFLOW /*| _EM_INEXACT*/);
	//new_fp_control_word = fp_control_word & ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	//_controlfp_s(&fp_control_word, new_fp_control_word, _MCW_EM);

	srand(time(0));

	/*
#if defined(CAR2D)
	Planner* ilopt = new CarPlanner();
#elif defined(POINT2D)
	Planner* ilopt = new PointPlanner();
#elif defined(AIRCRAFT3D)
	Planner* ilopt = new AircraftPlanner();
#endif

	ilopt->initPlanner();
	
	//ilopt->solveILQG();
	ilopt->solveDDP();

	//ilopt->simulate(10000, true);
	//ilopt->simulate(1, true);

	delete ilopt;	
	CAL_End();
	*/

	PointPlanner* popt = new PointPlanner();
	popt->initPlanner();

	popt->initProblem();
	//popt->solveILQG();
	popt->solveDDP();

	//for(int p = 0; p < 200; ++p) {
	//	std::cout << "Processing path: " << p+1 << std::endl;
	//	popt->initProblem();
	//	popt->solveDDP();
	//	popt->simulate(10000, false);
	//}

	//popt->simulate(10000, false);

	int num;
	std::cin >> num;

	delete popt;	
	CAL_End();

	return 0;
}
