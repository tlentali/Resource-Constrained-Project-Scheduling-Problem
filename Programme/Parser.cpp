/*
 * Parser.cpp
 *
 *  Created on: 19 nov. 2015
 *      Author: haboud910e
 */

#include "Parser.h"

int Parser::heuristik(vector<vector<int> > pred){
	vector<vector<int> > cap(horizon,vector<int>(renewable));
	vector<int> s(jobs, 0);
	bool flag, flag2;
	for (int i = 0; i < horizon; ++i) {
		for (int k = 0; k < renewable; ++k) {
			cap[i][k] = resource[k];
		}
	}
	s[0] = 0;
	for (int i = 1; i < jobs; ++i) {
		int min = 0;
		for (int j = 0; j < pred[i].size(); ++j) {
			if (min<s[pred[i][j]] + duration[pred[i][j]]) {
				min = s[pred[i][j]] + duration[pred[i][j]];
			}
		}
		flag = true;
		flag2 = true;
		while(flag){
			int tmp = 0;
			flag2 = true;
			while(flag2 and tmp < duration[i]){
				for (int k = 0; k < renewable; ++k) {
					if (requests[i][k]>cap[min+tmp][k]) {
						flag2 = false;
					}
				}
				tmp++;
			}
			if (flag2) {
				flag = false;
				s[i] = min;
				for (int j = 0; j < duration[i]; ++j) {
					for (int k = 0; k < renewable; ++k) {
						cap[min+j][k]-=requests[i][k];
					}
				}
			}
			min++;
		}
	}
	/*for (int i = 0; i < jobs; ++i) {
		cout<<s[i]<<endl;
	}*/

	return s[jobs-1];
}
void Parser::pert(){
	vector<vector<int> > pred(jobs,vector<int>());
	vector<int> nbp(jobs);
	vector<int> nbs(jobs);
	for (int i = 0; i < jobs; ++i) {

		asap.push_back(-1);
		alap.push_back(horizon+1);
		for (int j = 0; j < successor[i].size(); ++j) {
			pred[successor[i][j]].push_back(i);
		}
	}
	horizon = heuristik(pred);
	for (int i = 0; i < jobs; ++i) {
		nbp[i] = pred[i].size();
		nbs[i] = successor[i].size();
	}

	queue<int> o;
	asap[0]=0;
	for (int i = 0; i < successor[0].size(); ++i) {
		nbp[successor[0][i]]--;
		if (nbp[successor[0][i]] == 0) {

			o.push(successor[0][i]);
		}
	}
	while(not o.empty()){
		int j = o.front();
		o.pop();
		for (int i = 0; i < successor[j].size(); ++i) {
			nbp[successor[j][i]]--;
			if (nbp[successor[j][i]] == 0) {
				o.push(successor[j][i]);
			}
		}
		int max = 0;
		for (int i = 0; i < pred[j].size(); ++i) {
			int tmp = pred[j][i];
			if (max<asap[tmp]+duration[tmp]){
				max = asap[tmp]+duration[tmp];
			}
			asap[j] = max;
		}
	}

	queue<int> oo;
	alap[jobs-1]=horizon;
	for (int i = 0; i < pred[jobs-1].size(); ++i) {
		nbs[pred[jobs-1][i]]--;
		if (nbs[pred[jobs-1][i]]==0) {
			oo.push(pred[jobs-1][i]);
		}
	}
	while(not oo.empty()){
		int j = oo.front();
		oo.pop();
		for (int i = 0; i < pred[j].size(); ++i) {
			nbs[pred[j][i]]--;
			if (nbs[pred[j][i]]==0){
				oo.push(pred[j][i]);
			}
		}
		int min = horizon;
		for (int i = 0; i < successor[j].size(); ++i) {
			int tmp = successor[j][i];
			if (min>alap[tmp]-duration[j]) {
				min = alap[tmp]-duration[j];
			}
			alap[j] = min;
		}
	}







}

Parser::Parser(string filename) {
	ifstream fichier(filename.c_str(), ios::in);
	if(fichier)
	{
		string ligne;
		for (int i = 0; i < 5; ++i) {
			getline(fichier,ligne);
		}
		fichier>>ligne>>ws>>ligne>>ws>>ligne>>ws>>ligne>>ws>>jobs>>ws;
		fichier>>ligne>>ws>>ligne>>ws>>horizon>>ws;
		getline(fichier,ligne);
		fichier>>ligne>>ws>>ligne>>ws>>ligne>>ws>>renewable;
		getline(fichier,ligne);
		fichier>>ligne>>ws>>ligne>>ws>>ligne>>ws>>nonrenewable;
		getline(fichier,ligne);
		fichier>>ligne>>ws>>ligne>>ws>>ligne>>ws>>ligne>>ws>>doubly_constrained;
		for (int i = 0; i < 4; ++i) {
			getline(fichier, ligne);
		}
		fichier>>ligne>>ws>>ligne>>ws>>rel_dat>>ws>>due_date>>ws>>tard_cost>>ws>>MPM_Time>>ws;
		for (int i = 0; i < 3; ++i) {
			getline(fichier, ligne);
		}
		int nbs;
		nbSuccesseur = 0;
		for (int i = 0; i < jobs; ++i) {
			successor.push_back(vector<int>());
			fichier>>ligne>>ws>>ligne>>ws>>nbs>>ws;
			nbSuccesseur+=nbs;
			for (int j = 0; j < nbs; ++j) {
				int tmp;
				fichier>>tmp>>ws;
				successor[i].push_back(tmp-1);
			}
		}
		for (int i = 0; i < 4; ++i) {
			getline(fichier, ligne);
		}
		for (int i = 0; i < jobs; ++i) {
			int d;
			fichier>>ligne>>ws>>ligne>>ws>>d>>ws;
			duration.push_back(d);
			requests.push_back(vector<int>());
			for (int j = 0; j < renewable; ++j) {
				int tmp;
				fichier >> tmp >> ws;
				requests[i].push_back(tmp);
			}
		}
		for (int i = 0; i < 3; ++i) {
			getline(fichier, ligne);
		}
		for (int i = 0; i < renewable; ++i) {
			int r;
			fichier>>r>>ws;
			resource.push_back(r);
		}
		fichier.close();

		//pert();

	}
	else{
		cerr << "Impossible d'ouvrir le fichier !"<<endl;
	}
}

Parser::~Parser() {}

