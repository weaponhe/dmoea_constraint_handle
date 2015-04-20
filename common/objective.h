#ifndef __OBJECTIVE_H_
#define __OBJECTIVE_H_

#include "..\common\global.h"
#include "CTP.h"
#include "cec09.h"

void initTestInstance(){
	if (!strcmp(strTestInstance, "CTP1") || !strcmp(strTestInstance, "CTP2") || !strcmp(strTestInstance, "CTP3") || !strcmp(strTestInstance, "CTP4") || !strcmp(strTestInstance, "CTP5") || !strcmp(strTestInstance, "CTP6") || !strcmp(strTestInstance, "CTP7") || !strcmp(strTestInstance, "CTP8"))
	{
		lowBound.resize(nvar);
		uppBound.resize(nvar);
		CTP::CTP_Vars_Allocate(&(*(lowBound.begin())), &(*(uppBound.begin())));

		hv_ref.resize(nobj, 11);
	}
	else if (!strcmp(strTestInstance, "CF1") || !strcmp(strTestInstance, "CF2") || !strcmp(strTestInstance, "CF3") || !strcmp(strTestInstance, "CF4") || !strcmp(strTestInstance, "CF5") || !strcmp(strTestInstance, "CF6") || !strcmp(strTestInstance, "CF7") || !strcmp(strTestInstance, "CF8") || !strcmp(strTestInstance, "CF9") || !strcmp(strTestInstance, "CF10"))
	{
		lowBound.resize(nvar);
		uppBound.resize(nvar);

		CEC09::CEC09_Vars_Allocate(&(*(lowBound.begin())), &(*(uppBound.begin())), nvar, CEC09::getCEC09Type(strTestInstance));

		hv_ref.resize(nobj, 11);
	}
	cubeHV = 1.0;
	for (int i = 0; i < nobj; i++)
	{
		cubeHV *= hv_ref[i];
	}
}
void objectives(vector<double> &x_var, vector <double> &y_obj, double &cv)
{

	if (!strcmp(strTestInstance, "CTP1"))
	{
		vector<double> c(2, 0.0);
		CTP::CTP1(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())));
		cv = 0;
		cv += (c[0] >= 0.0 ? 0.0 : c[0]);
		cv += (c[1] >= 0.0 ? 0.0 : c[1]);
		c.clear();
	}

	else if (!strcmp(strTestInstance, "CTP2"))
	{
		vector<double> c(1, 0.0);
		CTP::CTP2(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())));
		cv = (c[0] >= 0.0 ? 0.0 : c[0]);
		c.clear();
	}

	else if (!strcmp(strTestInstance, "CTP3"))
	{
		vector<double> c(1, 0.0);
		CTP::CTP3(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())));
		cv = (c[0] >= 0.0 ? 0.0 : c[0]);
		c.clear();
	}

	else if (!strcmp(strTestInstance, "CTP4"))
	{
		vector<double> c(1, 0.0);
		CTP::CTP4(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())));
		cv = (c[0] >= 0.0 ? 0.0 : c[0]);
		c.clear();
	}

	else if (!strcmp(strTestInstance, "CTP5"))
	{
		vector<double> c(1, 0.0);
		CTP::CTP5(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())));
		cv = (c[0] >= 0.0 ? 0.0 : c[0]);
		c.clear();
	}

	else if (!strcmp(strTestInstance, "CTP6"))
	{
		vector<double> c(1, 0.0);
		CTP::CTP6(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())));
		cv = (c[0] >= 0.0 ? 0.0 : c[0]);
		c.clear();
	}

	else if (!strcmp(strTestInstance, "CTP7"))
	{
		vector<double> c(1, 0.0);
		CTP::CTP7(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())));
		cv = (c[0] >= 0.0 ? 0.0 : c[0]);
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CTP8"))
	{
		vector<double> c(2, 0.0);
		CTP::CTP8(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())));
		cv = 0;
		cv += (c[0] >= 0.0 ? 0.0 : c[0]);
		cv += (c[1] >= 0.0 ? 0.0 : c[1]);
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CF1"))
	{
		vector<double> c(1, 0.0);
		CEC09::CF1(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);
		cv = c[0] >= 0.0 ? 0.0 : c[0];
		c.clear();
	}

	else if (!strcmp(strTestInstance, "CF2"))
	{
		vector<double> c(1, 0.0);
		CEC09::CF2(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = c[0] >= 0.0 ? 0.0 : c[0];
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CF3"))
	{
		vector<double> c(1, 0.0);
		CEC09::CF3(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = c[0] >= 0.0 ? 0.0 : c[0];
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CF4"))
	{
		vector<double> c(1, 0.0);
		CEC09::CF4(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = c[0] >= 0.0 ? 0.0 : c[0];
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CF5"))
	{
		vector<double> c(1, 0.0);
		CEC09::CF5(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = c[0] >= 0.0 ? 0.0 : c[0];
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CF6"))
	{
		vector<double> c(2, 0.0);
		CEC09::CF6(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = 0;
		cv += c[0] >= 0.0 ? 0.0 : c[0];
		cv += c[1] >= 0.0 ? 0.0 : c[1];
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CF7"))
	{
		vector<double> c(2, 0.0);
		CEC09::CF7(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = 0;
		cv += c[0] >= 0.0 ? 0.0 : c[0];
		cv += c[1] >= 0.0 ? 0.0 : c[1];
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CF8"))
	{
		vector<double> c(1, 0.0);
		CEC09::CF8(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = c[0] >= 0.0 ? 0.0 : c[0];
		c.clear();
	}
	else if (!strcmp(strTestInstance, "CF9"))
	{
		vector<double> c(1, 0.0);
		CEC09::CF9(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = c[0] >= 0.0 ? 0.0 : c[0];
		c.clear();
	}

	else if (!strcmp(strTestInstance, "CF10"))
	{
		vector<double> c(1, 0.0);
		CEC09::CF10(&(*(x_var.begin())), &(*(y_obj.begin())), &(*(c.begin())), nvar);

		cv = c[0] >= 0.0 ? 0.0 : c[0];
		c.clear();
	}
}
#endif