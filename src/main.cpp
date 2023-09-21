#include<iostream>
#include "include/SPH3D.h"
#include "include/SPH2D.h"
#include "include/utils.h"
using namespace std;
int main(int argc, char** argv)
{
#define GRID_OPTIMIZE
	SPH_2D(argc, argv,1,true);
	//checkFunc();
}