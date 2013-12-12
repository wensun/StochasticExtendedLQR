#ifndef _OBSTACLES_
#define _OBSTACLES_


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
#include <time.h>
#include "callisto.h"
#include "matrix.h"
#include "utils.h"

class Obstacles
{
public:

	struct obs{
		Matrix<3> pos;
		double radius; 
		int dim;
	};

	std::vector<obs> obstacles;
	Matrix<3> Bottom;
	Matrix<3> Top;


	Obstacles()
	{
		obs obstacle;

		obstacle.pos[0] = -1.25; obstacle.pos[1] = 1.25; obstacle.pos[2] = 0;
		obstacle.radius = .4; obstacle.dim = 2;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = 1.25; obstacle.pos[1] = 1.25; obstacle.pos[2] = 0;
		obstacle.radius = .4; obstacle.dim = 2;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = 1.25; obstacle.pos[1] = -1.25; obstacle.pos[2] = 0;
		obstacle.radius = .4; obstacle.dim = 2;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = -1.25; obstacle.pos[1] = -1.25; obstacle.pos[2] = 0;
		obstacle.radius = .4; obstacle.dim = 2;
		obstacles.push_back(obstacle);


		obstacle.pos[0] = -1.25; obstacle.pos[1] = 0; obstacle.pos[2] = 1.25;
		obstacle.radius = .4; obstacle.dim = 1;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = 1.25; obstacle.pos[1] =  0; obstacle.pos[2] = 1.25;
		obstacle.radius = .4; obstacle.dim = 1;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = 1.25; obstacle.pos[1] = 0; obstacle.pos[2] = -1.25;
		obstacle.radius = .4; obstacle.dim = 1;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = -1.25; obstacle.pos[1] = 0; obstacle.pos[2] = -1.25;
		obstacle.radius = .4; obstacle.dim = 1;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = 0; obstacle.pos[1] = -1.25; obstacle.pos[2] = 1.25;
		obstacle.radius = .4; obstacle.dim = 0;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = 0; obstacle.pos[1] = 1.25;  obstacle.pos[2] = 1.25;
		obstacle.radius = .4; obstacle.dim = 0;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = 0; obstacle.pos[1] = 1.25; obstacle.pos[2] = -1.25;
		obstacle.radius = .4; obstacle.dim = 0;
		obstacles.push_back(obstacle);

		obstacle.pos[0] = 0; obstacle.pos[1] = -1.25; obstacle.pos[2] = -1.25;
		obstacle.radius = .4; obstacle.dim = 0;
		obstacles.push_back(obstacle);
	
		Top[0] = 3.0; Top[1] = 3.0; Top[2] = 3.0;
		Bottom[0] = -3.0; Bottom[1] = -3.0; Bottom[2] = -3.0;

	}

};


#endif