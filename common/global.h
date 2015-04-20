#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>


using namespace std;

#define pi 3.1415926

#include "..\common\random.h"

//******** Parameters in test instance *********************************************
int     nvar,      //  the number of variables
		nobj;      //  the number of objectives
fstream fout_temp;

vector<double> lowBound;
vector<double> uppBound;
vector<double> hv_ref;
double cubeHV;
enum TCompare { _Pareto_Dominating, _Pareto_Dominated, _Pareto_Nondominated, _Pareto_Equal };

//*************************************************
double acceptR = 1;
double T0 = 4.0;
double L = 2.0;
enum UpdateAlgorithmType { _CDP, _ADP, _PEN };
char *UpdateAlgorithmTypeStr[3] = { "CDP", "ADP","PEN" };
int updateType = _ADP;
bool isIGD;

char    strTestInstance[1024];
// *********************************************************************************

//******** Parameters in random number *********************************************
int     seed = 177;
long    rnd_uni_init;
//**********************************************************************************

//******** Common parameters in MOEAs **********************************************
int		max_gen = 500,    //  the maximal number of generations
max_run = 5,      //  the maximal number of runs
pops = 300,    //  the population size
nfes;             //  the number of function evluations
//**********************************************************************************


//*******  Parameters in MOEA/D ****************************************************
int	   niche = 20,
limit = 2,
unit = 33;

double scale[100];
vector <double> idealpoint;
char   strFunctionType[256],
strAlgorithmType[256];
//***************************************************************************************

//******** Parameters in SBX *******************************************************
int		etax = 20,
etam = 20;

double  realx, realm,    // probability in SBX crossover and polynomial mutation
realb = 0.9;     // probability of selecting mating parents from neighborhood

int     gID;
//**********************************************************************************

#endif