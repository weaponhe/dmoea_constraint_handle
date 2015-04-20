/*==========================================================================
//  Implementation of MOEA/D Based on Differential Evolution (DE) for Continuous Multiobjective
//  Optimization Problems with Complicate Pareto Sets (2007)
//
//  See the details of MOEA/D-DE and test problems in the following paper
//  H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  The component functions of each test instance can be found in "objective.h".
//
//  The source code of MOEA/D-DE and NSGA-II-DE were implemented by Hui Li and Qingfu Zhang
//
//  If you have any questions about the codes, please contact
//  Qingfu Zhang at qzhang@essex.ac.uk  or Hui Li at hzl@cs.nott.ac.uk
===========================================================================*/


#include "common\global.h"
#include "DMOEA\dmoeafunc.h"


void main()
{
	fstream configfile;
	configfile.open("TestInstance.txt", ios::in);
	int numInstance;
	configfile >> numInstance;
	char** ins = new char*[numInstance];
	int* var = new int[numInstance];
	int* obj = new int[numInstance];
	int* pop = new int[numInstance];
	int* gen = new int[numInstance];

	for (int i = 0; i < numInstance; i++)
	{
		ins[i] = new char[20];
		configfile >> ins[i];
		configfile >> var[i];
		configfile >> obj[i];
		configfile >> pop[i];
		configfile >> gen[i];
	}
	configfile.close();
	for (int i = 0; i < numInstance; i++)
	{
		isIGD = true;
		if (strcmp(ins[i], "CTP1") == 0 || strcmp(ins[i], "CTP2") == 0 || strcmp(ins[i], "CTP3") == 0 || strcmp(ins[i], "CTP4") == 0 || strcmp(ins[i], "CTP5") == 0 || strcmp(ins[i], "CTP6") == 0 || strcmp(ins[i], "CTP7") == 0 || strcmp(ins[i], "CTP8") == 0)
			isIGD = false;

		strcpy(strTestInstance, ins[i]);
		cout << strTestInstance << endl;
		nvar = var[i];
		nobj = obj[i];
		initTestInstance();

		// the parameter setting of MOEA/D-DE
		pops = pop[i];
		max_gen = gen[i];

		int H = 2;
		if (nobj == 2)
		{
			//pops
			H = pops - 1;
		}
		else
		{
			int mypop;
			while (1)
			{
				int numerator = 1;
				int denominator = 1;
				for (int i = H + nobj - 1; i > H; i--)
				{
					numerator *= i;
				}
				for (int i = nobj - 1; i > 0; i--)
				{
					denominator *= i;
				}
				mypop = numerator / denominator;
				if (mypop >= pops) {
					break;
				}
				H++;
			}
			pops = mypop;
		}

		niche = (int)(0.1 * pops);
		limit = (int)(0.01* pops);

		fstream fout_hv;
		char filename_hv[1024];
		sprintf(filename_hv, "LOG/%s/HV/HV_MOEAD-%s(%d)-p%d-g%d.dat", UpdateAlgorithmTypeStr[updateType], strTestInstance, nobj, pops, max_gen);
		fout_hv.open(filename_hv, ios::out);
		fstream fout_igd;
		if (isIGD){
			char filename_igd[1024];
			sprintf(filename_igd, "LOG/%s/IGD/IGD_MOEAD_%s(%d)-p%d-g%d.dat", UpdateAlgorithmTypeStr[updateType], strTestInstance, nobj, pops, max_gen);
			fout_igd.open(filename_igd, ios::out);
		}

		
		

		// Main loop
		for (int run = 1; run <= max_run; run++)
		{
			vector<double> hv;
			vector<double> igd;
			CMOEAD  MOEAD;

			char filename_temp[1024];
			sprintf(filename_temp, "LOG/%s/TEMP/%s.dat", UpdateAlgorithmTypeStr[updateType], strTestInstance);
			fout_temp.open(filename_temp, ios::out);

			MOEAD.execute(run, "_TCHE1", "_DE", hv, igd, H);

			fout_temp.close();

			for (int k = 0; k < hv.size();)
			{
				fout_hv << hv[k] << "\t" << hv[k + 1] << "\n";
				k = k + 2;
			}
			hv.clear();
			fout_hv << "\n";


			if (isIGD){
				for (int k = 0; k < igd.size();)
				{
					fout_igd << igd[k] << "\t" << igd[k + 1] << "\n";
					k = k + 2;
				}
				igd.clear();
				fout_igd << "\n";
				
			}
		}
		if (isIGD){
			fout_igd.close();
		}
		fout_hv.close();

	}
	cout << UpdateAlgorithmTypeStr[updateType] << endl;
}
