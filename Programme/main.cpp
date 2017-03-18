/*
 * main.cpp
 *
 *  Created on: 19 nov. 2015
 *      Author: haboud910e
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <cmath>

#include "Parser.h"



using namespace std;

bool vecotr(vector<int> o, vector<int> B){
		for (int k = 0; k < o.size(); ++k) {
			if (o[k]>B[k]) {
				return false;
			}
		}
		return true;
}

ILOLAZYCONSTRAINTCALLBACK4(Ctcoupe,IloArray<IloBoolVarArray>,x,IloIntVarArray,S,Parser,a,IloNum,eps){
	vector<int> L;
	for (int i = 0; i < a.jobs; ++i) {
		L.push_back(getValue(S[i])+a.duration[i]);
	}
	sort(L.begin(),L.end());
	for (int t = 0; t < a.jobs; ++t) {
		vector<int> o(a.renewable);
		vector<int> f;
		for (int j = 0; j < a.jobs; ++j) {
			if (getValue(S[j])<L[t] and getValue(S[j])+a.duration[j]>=L[t]){
				for (int k = 0; k < a.renewable; ++k) {
					o[k] += a.requests[j][k];
				}
				f.push_back(j);
				if (not vecotr(o,a.resource)) {
					IloRange cut;
					try{
						IloExpr coupe(getEnv());
						for (int i = 0; i < f.size(); ++i) {
							for (int ii = 0; ii < f.size(); ++ii) {
								if (i!=ii) {
									coupe += x[f[i]][f[ii]];
									coupe += x[f[ii]][f[i]];
								}
							}
						}
						//cout<<"he coucou"<<endl;
						cut = (coupe >=1);
						add(cut);

					}catch(...){
						throw;
					}

				}


			}

		}
	}

}


float mini(Parser a) {
	IloEnv env = IloEnv();
	IloModel m(env);
	float result = 0;
	try {
		IloInt M = a.horizon;
		IloArray<IloBoolVarArray> x(env, a.jobs);
		for (IloInt i = 0; i < a.jobs; ++i) {
			x[i] = IloBoolVarArray(env, a.jobs);
			for (IloInt j = 0; j < a.jobs; ++j) {
				if (i != j) {
					stringstream name1;
					name1 << "X_" << i << "_" << j;
					x[i][j] = IloBoolVar(env, name1.str().c_str());
				}
			}
		}
		IloIntVarArray S(env, a.jobs);
		for (IloInt i = 0; i < a.jobs; ++i) {
			stringstream name1;
			name1 << "S_" << i;
			S[i] = IloIntVar(env, name1.str().c_str());
		}
		IloArray<IloArray<IloIntVarArray> > f(env, a.jobs);
		for (IloInt i = 0; i < a.jobs; ++i) {
			f[i] = IloArray<IloIntVarArray>(env, a.jobs);
			for (IloInt j = 0; j < a.jobs; ++j) {
				f[i][j] = IloIntVarArray(env, a.renewable);
				for (IloInt k = 0; k < a.renewable; ++k) {
					if (i != j) {
						stringstream name2;
						name2 << "f_" << i << "_" << j << "_" << k;
						f[i][j][k] = IloIntVar(env, name2.str().c_str());
					}
				}
			}
		}

		IloArray<IloExpr> contr1(env, a.nbSuccesseur);
		IloInt cmp = 0;
		for (IloInt i = 0; i < a.successor.size() - 1; ++i) {
			for (IloInt j = 0; j < a.successor[i].size(); j++) {
				stringstream name;
				name << "Contr1_" << cmp;
				contr1[cmp] = IloExpr(env);
				contr1[cmp] += x[i][a.successor[i][j]];
				contr1[cmp].setName(name.str().c_str());
				m.add(contr1[cmp] == 1);
				cmp++;
			}
		}

		IloArray<IloExpr> contr2(env, (a.jobs * a.jobs - 1) / 2);
		cmp = 0;
		for (IloInt i = 0; i < a.jobs; ++i) {
			for (IloInt j = i + 1; j < a.jobs; j++) {
				contr2[cmp] = IloExpr(env);
				contr2[cmp] += x[i][j] + x[j][i];
				m.add(contr2[cmp] <= 1);
				cmp++;
			}
		}

		IloArray<IloExpr> contr3(env, pow(a.jobs, 3));

		for (IloInt i = 0; i < a.jobs; ++i) {
			for (IloInt j = 0; j < a.jobs; j++) {
				for (IloInt k = 0; k < a.jobs; k++) {
					if (i != j && k != j && k != j) {
						contr3[i * pow(a.jobs, 2) + j * a.jobs + k] = IloExpr(env);
						contr3[i * pow(a.jobs, 2) + j * a.jobs + k] += x[i][k]- x[i][j] - x[j][k];
						m.add(contr3[i * pow(a.jobs, 2) + j * a.jobs + k]>= -1);
					}
				}
			}
		}

		IloArray<IloExpr> contr4(env, pow(a.jobs, 2));
		for (IloInt i = 0; i < a.jobs; ++i) {
			for (IloInt j = 0; j < a.jobs; ++j) {
				if (i != j) {
					//M=abs(a.asap[j]-a.alap[i]);
					contr4[i * a.jobs + j] = IloExpr(env);
					contr4[i * a.jobs + j] += S[j] - S[i] + M- (a.duration[i] + M) * x[i][j];
					m.add(contr4[i * a.jobs + j] >= 0);
				}
			}
		}
		/*for (int i = 0; i < a.jobs; ++i) {
			m.add(S[i] >= a.asap[i]);
			m.add(S[i] <= a.alap[i]);
		}*/

		IloExpr obj(env);
		obj += S[a.jobs - 1];
		m.add(IloMinimize(env, obj));
		IloCplex cplex(m);
		cplex.use(Ctcoupe(env, x, S, a,cplex.getParam(IloCplex::EpRHS)));
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.exportModel("model.lp");
		cplex.setParam(cplex.TiLim,600);
		cplex.solve();
		cout<<"Temps : "<<cplex.getTime()<<" sec"<<endl;

		result = cplex.getObjValue();
		//m.end();

	} catch (IloException& e) {
		cerr << "Error" << e << endl;
		result = -1;
	}
	return result;
}
float conflits(Parser a){


IloEnv env=IloEnv();
	IloModel m(env);
	float result=0;
	try{
		IloInt T = a.horizon+1;
		IloArray<IloBoolVarArray> y(env, a.jobs);
		for (IloInt i = 0; i < a.jobs; ++i) {
			y[i] = IloBoolVarArray(env, T);
			for (IloInt j = 0; j < T; ++j) {
				stringstream name1;
				name1 << "Y_" << i << "_" << j;
				y[i][j] = IloBoolVar(env, name1.str().c_str());
				/*if (a.asap[i] > j or a.alap[i] < j) {
					m.add(y[i][j] == 0);
				}*/

			}
		}

		IloArray<IloExpr> contr1(env, a.jobs);
		for (IloInt i = 0; i < a.jobs; i++){
			contr1[i]= IloExpr(env);
			for (IloInt j = 0; j < T; j++)
			{
				contr1[i] += y[i][j];
			}
			m.add(contr1[i]==1);
		}

		IloArray<IloExpr> contr2(env, a.nbSuccesseur);
		IloInt cmp=0;
		for (IloInt i = 0; i < a.successor.size()-1; ++i) {
			for (IloInt j=0; j< a.successor[i].size();j++) {
				contr2[cmp] = IloExpr(env);
				for (IloInt t = 0; t < T; t++){
					contr2[cmp]+= t*(y[a.successor[i][j]][t] - y[i][t]);
				}
				m.add(contr2[cmp]>=a.duration[i]);
				cmp++;
			}
		}

		IloArray<IloExpr> contr3(env, a.renewable * T);
		for (IloInt k = 0; k < a.renewable; k++)
		{
			for (IloInt t = 0; t < T; t++)
			{
				contr3[k*T+t] = IloExpr(env);
				for (int i = 0; i < a.jobs; i++)
				{
					for (int r = t-a.duration[i] + 1; r <= t ; r++)
					{
						if (r>=0){//forcer r a etre positif
						contr3[k*T+t] += a.requests[i][k] * y[i][r];
						}
					}
				}
				m.add(contr3[k*T+t] <= a.resource[k]);
			}
		}
		IloExpr obj(env);
		for (int t = 0; t < T; t++)
		{
			obj += t * y[a.jobs-1][t];
		}
		
		m.add(IloMinimize(env, obj));
		IloCplex cplex(m);
		cplex.extract(m);
		cplex.setOut(env.getNullStream());
		cplex.setParam(cplex.TiLim,600);
		cplex.exportModel("model.lp");
		cplex.solve();
		cout<<"Temps : "<<cplex.getTime()<<" sec"<<endl;

		result = cplex.getObjValue();

		//m.end();

	}catch(IloException& e){
		cerr << "Error"<<e<<endl;
		result = -1;
	}
	return result;
}
float flot(Parser a){
	for (int i = 0; i < a.renewable; ++i) {
		a.requests[0][i] = a.resource[i];
		a.requests[a.jobs - 1][i] = a.resource[i];

	}
	IloEnv env=IloEnv();
	IloModel m(env);
	float result=0;
	try{
		IloInt M = a.horizon;
		IloArray<IloBoolVarArray> x(env, a.jobs);
		for (IloInt i = 0; i < a.jobs; ++i) {
			x[i] = IloBoolVarArray(env, a.jobs);
			for (IloInt j = 0; j < a.jobs; ++j) {
				if (i!=j){
					stringstream name1;
					name1 << "X_"<<i<<"_"<<j;
					x[i][j] = IloBoolVar(env,name1.str().c_str());
				}
			}
		}
		IloIntVarArray S (env, a.jobs);
		for (IloInt i = 0; i < a.jobs; ++i) {
			stringstream name1;
			name1 << "S_"<<i;
			S[i] = IloIntVar(env, name1.str().c_str());
		}
		IloArray<IloArray<IloIntVarArray> > f(env, a.jobs);
		for (IloInt i = 0; i < a.jobs; ++i) {
			f[i] = IloArray<IloIntVarArray>(env, a.jobs);
			for (IloInt j = 0; j < a.jobs; ++j) {
				f[i][j] = IloIntVarArray(env, a.renewable);
				for (IloInt k = 0; k < a.renewable; ++k) {
					if (i!=j){
						stringstream name2;
						name2 << "f_"<<i<<"_"<<j<<"_"<<k;
						f[i][j][k] = IloIntVar(env,name2.str().c_str());
					}
				}
			}
		}

		IloArray<IloExpr> contr1(env, a.nbSuccesseur);
		IloInt cmp=0;
		for (IloInt i = 0; i < a.successor.size()-1; ++i) {
			for (IloInt j=0; j< a.successor[i].size();j++) {
				stringstream name;
				name <<"Contr1_"<<cmp;
				contr1[cmp] = IloExpr(env);
				contr1[cmp] += x[i][a.successor[i][j]];
				contr1[cmp].setName(name.str().c_str());
				m.add(contr1[cmp]==1);
				cmp++;
			}
		}

		IloArray<IloExpr> contr2(env, (a.jobs*a.jobs-1)/2);
		cmp=0;
		for (IloInt i = 0; i < a.jobs; ++i) {
			for (IloInt j=i+1; j< a.jobs;j++) {
				contr2[cmp] = IloExpr(env);
				contr2[cmp] += x[i][j] + x[j][i];
				m.add(contr2[cmp] <= 1);
				cmp++;
			}
		}

		IloArray<IloExpr> contr3(env, pow(a.jobs,3));

		for (IloInt i = 0; i < a.jobs; ++i) {
			for (IloInt j=0; j< a.jobs;j++) {
				for (IloInt k=0; k< a.jobs;k++) {
					if (i!=j && k!=j && k!=j){
						contr3[i*pow(a.jobs,2)+j*a.jobs+k] = IloExpr(env);
						contr3[i*pow(a.jobs,2)+j*a.jobs+k] += x[i][k] - x[i][j] - x[j][k];
						m.add(contr3[i*pow(a.jobs,2)+j*a.jobs+k] >= -1);
					}
				}
			}
		}

		IloArray<IloExpr> contr4(env, pow(a.jobs,2));
		for (IloInt i = 0; i < a.jobs; ++i) {
			for (IloInt j = 0; j < a.jobs; ++j) {
				if (i!=j){
					//M=abs(a.asap[j]-a.alap[i]);
					contr4[i*a.jobs+j] = IloExpr(env);
					contr4[i*a.jobs+j] += S[j]-S[i]+M-(a.duration[i]+M)*x[i][j];
					m.add(contr4[i*a.jobs+j]>=0);
				}
			}
		}

		IloArray<IloExpr> contr5(env, a.jobs*a.renewable);
		IloArray<IloExpr> contr6(env, a.jobs*a.renewable);

		for (IloInt i = 0; i < a.jobs; ++i) {
			for (IloInt k = 0; k < a.renewable; ++k) {
				contr5[i*a.renewable+k] = IloExpr(env);
				contr6[i*a.renewable+k] = IloExpr(env);
				for (IloInt j = 0; j < a.jobs; ++j) {
					if (i!=j){
						if (j!=0 && i!=a.jobs-1){
							contr5[i*a.renewable+k] += f[i][j][k];
						}
						if (j!=a.jobs-1 && i!=0){
							contr6[i*a.renewable+k] += f[j][i][k];
						}
					}
				}
				if(i!=a.jobs-1){
					m.add(contr5[i*a.renewable+k]==a.requests[i][k]);
				}
				if(i!=0){
				m.add(contr6[i*a.renewable+k]==a.requests[i][k]);
				}
			}
		}

		IloArray<IloExpr> contr7(env, pow(a.jobs,2)*a.renewable);

		for (IloInt i = 0; i < a.jobs; ++i) {
			for (IloInt j = 0; j < a.jobs; ++j) {
				if (i!=j){
					for (IloInt k = 0; k < a.renewable; ++k) {
						contr7[i*a.jobs*a.renewable+j*a.renewable+k] = IloExpr(env);
						contr7[i*a.jobs*a.renewable+j*a.renewable+k] += f[i][j][k];
						if (a.requests[i][k]<a.requests[j][k]) {
							contr7[i*a.jobs*a.renewable+j*a.renewable+k] -= a.requests[i][k]*x[i][j];
						} else {
							contr7[i*a.jobs*a.renewable+j*a.renewable+k] -= a.requests[j][k]*x[i][j];
						}
						m.add(contr7[i*a.jobs*a.renewable+j*a.renewable+k]<=0);
					}
				}
			}
		}

		//test

		/*for (int i = 0; i < a.jobs; ++i) {
			m.add(S[i]>=a.asap[i]);
			m.add(S[i]<=a.alap[i]);
		}*/


		IloExpr obj(env);
		obj += S[a.jobs-1];
		m.add(IloMinimize(env, obj));
		IloCplex cplex(m);
		cplex.setOut(env.getNullStream());
		cplex.exportModel("model.lp");
		cplex.setParam(cplex.TiLim,600);
		cplex.solve();
		cout<<"Temps : "<<cplex.getTime()<<" sec"<<endl;
		result = cplex.getObjValue();
		//m.end();

	}catch(IloException& e){
		//cerr << "Error"<<e<<endl;
		result = -1;
	}
	return result;
}


int main(){


	/*for (int i = 0; i < a.jobs; ++i) {
		cout<<a.asap[i]<<"  "<<a.alap[i]<<endl;
	}*/

	cout<<"__________________________________"<<endl;
	for (int i = 1; i <= 1; ++i) {
		for (int j = 1; j <= 10; ++j) {

			stringstream nom;
			nom <<"/autofs/netapp/travail/haboud910e/workspace/GOPP/Instances/j30"<<i<<"_"<<j<<".sm";
			cout<<"Nom de l'instance : "<<nom.str()<<endl;
			Parser a(nom.str());
			cout<<"Formulation de flot"<<endl;
			cout<<"Valeur de l'objectif : "<<flot(a)<<endl;
			cout<<"Configuration interdites minimales"<<endl;
			cout<<"Valeur de l'objectif : "<<mini(a)<<endl;
			cout<<"Formulation indexé par le temps"<<endl;
			cout<<"Valeur de l'objectif : "<<conflits(a)<<endl;
			cout<<"__________________________________"<<endl;


		}
	}

	return 0;
}
