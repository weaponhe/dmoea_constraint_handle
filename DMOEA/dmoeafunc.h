#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include "..\common\global.h"
#include "..\common\recombination.h"
#include "..\common\common.h"
#include "dmoeaclass.h"
#include "..\common\wfg.h"


class CMOEAD
{

public:
	CMOEAD();
	virtual ~CMOEAD();


	void init_uniformweight();               // initialize the weights for subproblems
	void init_uniformweight(int sd);               // initialize the weights for subproblems
	void gen_uniformweight(int start_obj_index, int max_value_left, vector<int> coordinate, int sd);
	void init_neighbourhood();               // calculate the neighbourhood of each subproblem
	void init_population();                  // initialize the population


	void update_reference(CMOEADInd &ind);                 // update ideal point which is used in Tchebycheff or NBI method
	void update_problem(CMOEADInd &ind, int id, int type); // update current solutions in the neighbourhood
	void update_problem_PEN(CMOEADInd &ind, int id, int type);
	void update_problem_CDP(CMOEADInd &ind, int id, int type);
	void update_problem_ADP(CMOEADInd &indiv, int id, int type);
	bool isUpdate(CMOEADInd &indiv, CSUB &subP);

	void diffevolution();                                      // DE-based recombination
	void matingselection(vector<int> &list, int cid, int size, int type);  // select mating parents

	// execute MOEAD
	void execute(int run, char *strfunc, char *stralg, vector<double> &hv, vector<double> &igd, int sd);

	void read_front(char filename[1024]);
	double calc_IGD();
	double calc_HV();
	void population2front(vector <CSUB>  &mypopulation, vector <vector<double> >  &population_front);

	void save_front_pop(char savefilename[1024]);
	void save_front_pof(char savefilename[1024]);       // save the pareto front into files
	void save_front_archive(char savefilename[1024]);
	void save_ps(char savefilename[1024]);

	vector <CSUB>       population;
	vector <CMOEADInd>  ps;
	vector <CSUB>       offspring;
	vector <int>        array;
	CMOEADInd           *ind_arr;

	//double              distance;                   // generational distance from PF to solutions found
	int                 popsize;

	void operator=(const CMOEAD &moea);
};

CMOEAD::CMOEAD()
{

	ind_arr = new CMOEADInd[nobj];

	// initialize ideal point
	for (int n = 0; n < nobj; n++)
	{
		idealpoint.push_back(1.0e+30);
		ind_arr[n].rnd_init();
		ind_arr[n].obj_eval();
	}
}

CMOEAD::~CMOEAD()
{
	idealpoint.clear();
	delete[] ind_arr;
}


void CMOEAD::init_population()
{

	for (int i = 0; i < population.size(); i++)
	{
		population[i].indiv.rnd_init();
		population[i].indiv.obj_eval();
		update_reference(population[i].indiv);
		population[i].archive = population[i].indiv;
		nfes++;
	}
}

void CMOEAD::operator=(const CMOEAD &alg)
{

	population = alg.population;
	ps = alg.ps;
	ind_arr = alg.ind_arr;
	offspring = alg.offspring;
	popsize = alg.popsize;
}


// createt the weight vectors with uniform distribution
void CMOEAD::init_uniformweight()
{
	if (nobj == 2)
	{
		//vector<CMOEADInd> ws;
		//loadpfront("F6Weight500.dat",ws);
		//pops = 500;

		for (int n = 0; n < pops; n++)
		{
			CSUB sub;
			double a = 1.0*n / (pops - 1);
			sub.namda.push_back(a);
			sub.namda.push_back(1 - a);

			//load weight vectors from file
			//sub.namda.push_back(ws[n].y_obj[0]);
			//sub.namda.push_back(ws[n].y_obj[1]);

			population.push_back(sub);
		}
		popsize = pops;
	}
	else
	{
		for (int i = 0; i <= unit; i++)
		{
			for (int j = 0; j <= unit; j++)
			{
				if (i + j <= unit)
				{
					CSUB sub;
					sub.array.push_back(i);
					sub.array.push_back(j);
					sub.array.push_back(unit - i - j);
					for (int k = 0; k < sub.array.size(); k++)
						sub.namda.push_back(1.0*sub.array[k] / unit);
					population.push_back(sub);
				}
			}
		}

		popsize = population.size();
		pops = popsize;
	}
}
void CMOEAD::init_uniformweight(int sd)
{
	/*
	if(nobj==2)
	{
	//vector<CMOEADInd> ws;
	//loadpfront("F6Weight500.dat",ws);
	//pops = 500;

	for(int n=0; n<pops; n++)
	{
	CSUB sub;
	double a = 1.0*n/(pops - 1);
	sub.namda.push_back(a);
	sub.namda.push_back(1-a);

	//load weight vectors from file
	//sub.namda.push_back(ws[n].y_obj[0]);
	//sub.namda.push_back(ws[n].y_obj[1]);

	population.push_back(sub);
	}
	popsize = pops;
	}
	else
	{
	for(int i=0; i<=unit; i++)
	{
	for(int j=0; j<=unit; j++)
	{
	if(i+j<=unit)
	{
	CSUB sub;
	sub.array.push_back(i);
	sub.array.push_back(j);
	sub.array.push_back(unit-i-j);
	for(int k=0; k<sub.array.size(); k++)
	sub.namda.push_back(1.0*sub.array[k]/unit);
	population.push_back(sub);
	}
	}
	}

	popsize = population.size();
	pops    = popsize;
	}
	*/

	vector<int> coordinate(nobj, 0);
	gen_uniformweight(nobj - 1, sd, coordinate, sd);
	coordinate.clear();
	pops = population.size();
}

//分解子问题、并计算观察向量
void CMOEAD::gen_uniformweight(int start_obj_index, int max_value_left, vector<int> coordinate, int sd)
{
	//start_obj_index=nobj-1
	//max_value_left = sd(pops-1 when objs=2)
	if (0 == start_obj_index || 0 == max_value_left)
	{
		coordinate[start_obj_index] = max_value_left;
		CSUB sub;
		for (int k = 0; k < nobj; k++)
		{
			sub.array.push_back(coordinate[k]);
			sub.namda.push_back(1.0*sub.array[k] / sd);
		}
		population.push_back(sub);
		return;
	}

	for (int i = max_value_left; i >= 0; i--)
	{
		coordinate[start_obj_index] = i;
		gen_uniformweight(start_obj_index - 1, max_value_left - i, coordinate, sd);
	}
}

void CMOEAD::init_neighbourhood()
{
	double *x = new double[population.size()];
	int    *idx = new int[population.size()];
	for (int i = 0; i < population.size(); i++)
	{
		// calculate the distances based on weight vectors
		for (int j = 0; j < population.size(); j++)
		{
			x[j] = dist_vector(population[i].namda, population[j].namda);
			idx[j] = j;
		}

		// find 'niche' nearest neighboring subproblems
		minfastsort(x, idx, population.size(), niche);
		for (int k = 0; k < niche; k++)
		{
			population[i].table.push_back(idx[k]);
		}

	}
	delete[] x;
	delete[] idx;
}

void CMOEAD::update_problem(CMOEADInd &indiv, int id, int type)
{
	// indiv: child solution
	// id:   the id of current subproblem
	// type: update solutions in - neighborhood (1) or whole population (otherwise)
	int size, time = 0;
	if (type == 1)	size = population[id].table.size();
	else        size = population.size();
	int *perm = new int[size];
	random_permutation(perm, size);
	for (int i = 0; i < size; i++)
	{
		int k;
		if (type == 1) k = population[id].table[perm[i]];
		else        k = perm[i];

		// calculate the values of objective function regarding the current subproblem
		double f1, f2;
		f1 = fitnessfunction(population[k].indiv.y_obj, population[k].namda, ind_arr);
		f2 = fitnessfunction(indiv.y_obj, population[k].namda, ind_arr);
		if (f2 < f1)
		{
			population[k].indiv = indiv;
			time++;
		}
		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if (time >= limit)
		{
			return;
		}
	}
	delete[] perm;
}

void CMOEAD::update_problem_PEN(CMOEADInd &indiv, int id, int type)
{
	int size, time = 0;
	if (type == 1)	size = population[id].table.size();
	else        size = population.size();
	int *perm = new int[size];
	random_permutation(perm, size);
	for (int i = 0; i < size; i++)
	{
		int k;
		if (type == 1) k = population[id].table[perm[i]];
		else        k = perm[i];
		double f1, f2;
		f1 = fitnessfunction(population[k].indiv.y_obj, population[k].namda, ind_arr);
		f2 = fitnessfunction(indiv.y_obj, population[k].namda, ind_arr);
		f1 += abs(population[k].indiv.cv);
		f2 += abs(indiv.cv);
		if (f2 < f1)
		{
			population[k].indiv = indiv;
			time++;
		}
		if (time >= limit)
		{
			return;
		}
	}
	delete[] perm;
}
void CMOEAD::update_problem_CDP(CMOEADInd &indiv, int id, int type)
{
	// indiv: child solution
	// id:   the id of current subproblem
	// type: update solutions in - neighborhood (1) or whole population (otherwise)
	int size, time = 0;
	if (type == 1)	size = population[id].table.size();
	else        size = population.size();
	int *perm = new int[size];
	random_permutation(perm, size);
	for (int i = 0; i < size; i++)
	{
		int k;
		if (type == 1) k = population[id].table[perm[i]];
		else        k = perm[i];

		if (indiv.cv >= 0 && population[k].indiv.cv < 0)
		{
			population[k].indiv = indiv;
			time++;
		}
		else if (indiv.cv < 0 && population[k].indiv.cv < 0)
		{
			if (indiv.cv > population[k].indiv.cv)
			{
				population[k].indiv = indiv;
				time++;
			}
		}
		else if (indiv.cv >= 0 && population[k].indiv.cv >= 0)
		{
			// calculate the values of objective function regarding the current subproblem
			double f1, f2;
			f1 = fitnessfunction(population[k].indiv.y_obj, population[k].namda, ind_arr);
			f2 = fitnessfunction(indiv.y_obj, population[k].namda, ind_arr);
			if (f2 < f1)
			{
				population[k].indiv = indiv;
				time++;
			}
		}
		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if (time >= limit)
		{
			return;
		}
	}
	delete[] perm;
}
void CMOEAD::update_problem_ADP(CMOEADInd &indiv, int id, int type)
{
	// indiv: child solution
	// id:   the id of current subproblem
	// type: update solutions in - neighborhood (1) or whole population (otherwise)
	int size, time = 0;
	if (type == 1)	size = population[id].table.size();
	else        size = population.size();
	int *perm = new int[size];
	random_permutation(perm, size);

	for (int i = 0; i < size; i++)
	{
		int k;
		if (type == 1) k = population[id].table[perm[i]];
		else        k = perm[i];

		if (isUpdate(indiv, population[k]))
		{
			population[k].indiv = indiv;
			population[k].countLU = 0;
			time++;
		}
		else
			population[k].countLU++;

		//保存档案
		if (population[k].countLU == 0)
		{
			if (population[k].indiv.cv >= 0.0)
			{
				if (population[k].archive.cv < 0.0)
					population[k].archive = population[k].indiv;
				else//都是feasible
				{
					//为什么在这里要加入一个占优比较？
					TCompare result = ParetoCompare(population[k].archive.y_obj, population[k].indiv.y_obj);
					if (result == _Pareto_Dominated)
					{
						population[k].archive = population[k].indiv;
					}
					else if (result == _Pareto_Nondominated)
					{
						double f1, f2;
						f1 = fitnessfunction(population[k].indiv.y_obj, population[k].namda, ind_arr);
						f2 = fitnessfunction(population[k].archive.y_obj, population[k].namda, ind_arr);
						if (f2 > f1)
						{
							population[k].archive = population[k].indiv;
						}
					}
				}
			}
		}

		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if (time >= limit)
		{
			return;
		}
	}
	delete[] perm;
}

void CMOEAD::update_reference(CMOEADInd &ind)
{
	//ind: child solution
	for (int n = 0; n < nobj; n++)
	{
		if (ind.y_obj[n] < idealpoint[n])
		{
			idealpoint[n] = ind.y_obj[n];
			ind_arr[n] = ind;
		}
	}
}

void CMOEAD::matingselection(vector<int> &list, int cid, int size, int type){
	// list : the set of the indexes of selected mating parents
	// cid  : the id of current subproblem
	// size : the number of selected mating parents
	// type : 1 - neighborhood; otherwise - whole population
	int ss = population[cid].table.size(), r, p;
	while (list.size() < size)
	{
		if (type == 1){
			r = int(ss*rnd_uni(&rnd_uni_init));
			p = population[cid].table[r];
		}
		else
			p = int(population.size()*rnd_uni(&rnd_uni_init));

		bool flag = true;
		for (int i = 0; i < list.size(); i++)
		{
			if (list[i] == p) // p is in the list
			{
				flag = false;
				break;
			}
		}

		if (flag) list.push_back(p);
	}
}

void CMOEAD::diffevolution()
{
	pops = population.size();
	int *perm = new int[pops];
	random_permutation(perm, pops);

	for (int i = 0; i < pops; i++)
	{
		int n = perm[i];
		// or int n = i;
		int type;
		double rnd = rnd_uni(&rnd_uni_init);

		// mating selection based on probability
		if (rnd < realb)    type = 1;   // neighborhood
		else             type = 2;   // whole population

		// select the indexes of mating parents
		vector<int> p;
		matingselection(p, n, 2, type);  // neighborhood selection

		// produce a child solution
		CMOEADInd child;
		diff_evo_xover2(population[n].indiv, population[p[0]].indiv, population[p[1]].indiv, child);

		// apply polynomial mutation
		realmutation(child, 1.0 / nvar);

		// evaluate the child solution
		child.obj_eval();

		// update the reference points and other solutions in the neighborhood or the whole population
		update_reference(child);
		if (updateType == _ADP){
			update_problem_ADP(child, n, type);
		}
		else if (updateType == _CDP)
		{
			update_problem_CDP(child, n, type);
		}
		else if (updateType == _PEN)
		{
			update_problem_PEN(child, n, type);
		}
		p.clear(); 	nfes++;
	}

	delete[] perm;
}


void CMOEAD::execute(int run, char *strfunc, char *stralg, vector<double> &hv, vector<double> &igd, int sd)
{
	seed = (seed + 23) % 1377;
	rnd_uni_init = -(long)seed;

	strcpy(strFunctionType, strfunc);
	strcpy(strAlgorithmType, stralg);

	// initialization 
	int gen = 0;
	nfes = 0;
	init_uniformweight(sd);
	init_neighbourhood();
	init_population();

	char filename[1024];
	char filename_pf[1024];

	if (isIGD){
		sprintf(filename_pf, "ParetoFront/%s(%d).dat", strTestInstance, nobj);
		loadpfront(filename_pf, ps);
	}
	double hvvalue;
	double igdvalue;

	hvvalue = calc_HV();
	hv.push_back(gen);  hv.push_back(hvvalue);
	cout << "gen = " << gen << "  hv = " << hvvalue;

	if (isIGD){
		igdvalue = calc_IGD();
		igd.push_back(0);  igd.push_back(igdvalue);
		cout << "\t\t igd = " << igdvalue;
	}
	cout << endl;

	// evolution
	for (gen = 1; gen <= max_gen; gen++)
	{
		acceptR = 0.01 / (gen * gen);
		diffevolution();

		int dd = int(max_gen / 25.0);

		// calculate hv-values
		if (gen%dd == 0)
		{
			hvvalue = calc_HV();
			hv.push_back(gen);  hv.push_back(hvvalue);
			cout << "gen = " << gen << "  hv = " << hvvalue;

			if (isIGD)
			{
				igdvalue = calc_IGD();
				igd.push_back(gen); igd.push_back(igdvalue);
				cout << "\t\t  igd = " << igdvalue;
			}
			cout << endl;
		}

		// save the final population - F space
		if (gen%max_gen == 0)
		{
			sprintf(filename, "PF/%s/POP/MOEAD_%s(%d)_%d_%d_R%d.dat", UpdateAlgorithmTypeStr[updateType], strTestInstance, nobj, pops, max_gen, run);
			save_front_pop(filename);

			sprintf(filename, "PF/%s/POF/MOEAD_%s(%d)_%d_%d_R%d.dat", UpdateAlgorithmTypeStr[updateType], strTestInstance, nobj, pops, max_gen, run);
			save_front_pof(filename);

			sprintf(filename, "PF/%s/Archive/MOEAD_%s(%d)_%d_%d_R%d.dat", UpdateAlgorithmTypeStr[updateType], strTestInstance, nobj, pops, max_gen, run);
			save_front_archive(filename);
		}
	}

	population.clear();
	ps.clear();
	std::cout << " Outcome of the " << strTestInstance << " " << run << "th run:" << " nfes = " << nfes << endl;
}


void CMOEAD::save_front_pop(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n < population.size(); n++)
	{
		for (int k = 0; k < nobj; k++)
			fout << population[n].indiv.y_obj[k] << "  ";
		fout << "\n";
	}
	fout.close();
}

void CMOEAD::save_front_pof(char saveFilename[1024])
{

	vector<vector<double> > population_front;
	population2front(population, population_front);
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n < population_front.size(); n++)
	{
		for (int k = 0; k < population_front[n].size(); k++)
			fout << population_front[n][k] << "  ";
		fout << "\n";
	}
	fout.close();
}

void CMOEAD::save_front_archive(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n < population.size(); n++)
	{
		for (int k = 0; k < nobj; k++)
			fout << population[n].archive.y_obj[k] << "  ";
		fout << "\n";
	}
	fout.close();
}

void CMOEAD::save_ps(char saveFilename[1024])
{
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n < popsize; n++)
	{
		for (int k = 0; k < nvar; k++)
			fout << population[n].indiv.x_var[k] << "  ";
		fout << "\n";
	}
	fout.close();
}


double CMOEAD::calc_IGD()
{
	double distance = 0;
	for (int i = 0; i < ps.size(); i++)
	{
		double min_d = 1.0e+10;
		for (int j = 0; j < population.size(); j++)
		{
			double d = dist_vector(ps[i].y_obj, population[j].indiv.y_obj);
			if (d < min_d)  min_d = d;
		}
		distance += min_d;
	}
	distance /= ps.size();
	return distance;
}

double CMOEAD::calc_HV()
{
	wfg hv;
	vector<vector<double> > population_front;
	population2front(population, population_front);


	int i = 0;
	int j;
	int size = population_front.size();
	for (; i < size; i++)
	{
		for (j = 0; j<nobj; j++)
		{
			if (population_front[i][j] > hv_ref[j])
				population_front[i][j] = hv_ref[j];
		}
	}


	double result = hv.compute(population_front, hv_ref);
	population_front.clear();
	//cout<<result<<"\t"<<cubeHV<<endl;
	return result / cubeHV;
}

void CMOEAD::population2front(vector <CSUB>  &mypopulation, vector <vector<double> >  &population_front)
{
	vector<int> nDominated;
	for (int n = 0; n < mypopulation.size(); n++)
		nDominated.push_back(0);

	for (int k = 0; k < mypopulation.size(); k++)
	{
		if (mypopulation[k].indiv.cv < 0.0)
		{
			nDominated[k]++;
			continue;
		}
		for (int j = k + 1; j < mypopulation.size(); j++)
		{
			if (mypopulation[j].indiv.cv < 0.0)
			{
				nDominated[j]++;
				continue;
			}
			TCompare tresult = ParetoCompare(mypopulation[k].indiv.y_obj, mypopulation[j].indiv.y_obj);
			if (tresult == _Pareto_Dominated)
				nDominated[k]++;
			//对相等的解得处理会不会影响HV的值
			else if (tresult == _Pareto_Dominating || tresult == _Pareto_Equal)
				nDominated[j]++;
		}
	}
	for (int n = 0; n < mypopulation.size(); n++)
	if (nDominated[n] == 0)
		population_front.push_back(mypopulation[n].indiv.y_obj);

	nDominated.clear();
}

bool CMOEAD::isUpdate(CMOEADInd &indiv, CSUB &subP)
{
	if (subP.indiv.cv >= 0.0 && indiv.cv >= 0.0)
	{
		/*
		TCompare result = ParetoCompare(indiv.y_obj,subP.indiv.y_obj);
		if(result == _Pareto_Dominating)
		return true;
		else if(result == _Pareto_Nondominated)
		{
		double f1, f2;
		f1 = fitnessfunction(subP.indiv.y_obj, subP.namda, ind_arr);
		f2 = fitnessfunction(indiv.y_obj, subP.namda, ind_arr);
		if (f2 < f1)
		return true;
		else
		return false;
		}

		else
		return false;
		*/

		double f1, f2;
		f1 = fitnessfunction(subP.indiv.y_obj, subP.namda, ind_arr);
		f2 = fitnessfunction(indiv.y_obj, subP.namda, ind_arr);
		if (f2 < f1)
			return true;
		else
			return false;
	}
	else if (indiv.cv >= 0.0) //subP.indiv.cv < 0.0
		return true;
	else if (indiv.cv < 0.0 && subP.indiv.cv < 0.0)
	{
		if (indiv.cv > subP.indiv.cv)
			return true;
		else
		{

			double f1, f2;
			f1 = fitnessfunction(subP.indiv.y_obj, subP.namda, ind_arr);
			f2 = fitnessfunction(indiv.y_obj, subP.namda, ind_arr);

			//在这里只考虑到约束违反程度
			//要让约束违反程度大的，更有“潜力”的个体保留下来
			//T是根据进化代数的增加而降低
			//double T = T0 / (1.0 + subP.countLU / L);
			double T = T0 / (1.0 + subP.countLU / L);
			//replaceR随着温度的降低而降低
			//replaceR随着约束违反程度差距的增大而减小
			double replaceR = exp((-subP.indiv.cv + indiv.cv) / T);

			//记录约束违反程度:子cv/父cv/子f/父f
			fout_temp << subP.indiv.cv << "\t" << indiv.cv << "\t" << f2 << "\t" << f1 << "\t" << T << "\t" << replaceR << "\n";

			if (replaceR > acceptR)
			{
				return true;
			}
			else
				return false;
		}
	}
	else
		return false;
}
#endif