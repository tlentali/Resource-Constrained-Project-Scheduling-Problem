/*
 * Parser.h
 *
 *  Created on: 19 nov. 2015
 *      Author: haboud910e
 */
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <queue>
#ifndef PARSER_H_
#define PARSER_H_

using namespace std;


class Parser {
public:
	void pert();
	int heuristik(vector<vector<int> > pred);
	int jobs, horizon, renewable, nonrenewable, doubly_constrained, rel_dat, due_date, tard_cost, MPM_Time,nbSuccesseur;
	vector<vector<int> > successor,requests,gap;
	vector<int> duration,resource,asap,alap;
	Parser(string filename);
	virtual ~Parser();
	void initialize();

};

#endif /* PARSER_H_ */
