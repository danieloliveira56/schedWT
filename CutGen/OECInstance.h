#pragma once

#include<iostream>
using namespace std;

//Classe instancia
class Instance {
  public:
  string nome;	
  int n;
  int m;
  int T;
  vector<int> p; //Processing Time
  vector<int> w; //Weight
  vector<int> d; //Due date

  int rCount;
  int pTotal; //duracao de todas as tarefas somadas


	Instance(InstanceInfo instance_) 
	{		
		rCount=0;
		pTotal = 0;

		n = instance_.numNodes-1;
		m = instance_.nrootbranches;

		p.push_back(0);
		for (int i = 1; i <= n; i++)
		{
			p.push_back(instance_.demand[i]);
		}

		calculaT();
	}
		
	void calculaT(void) {
		int i;
		int pMax=0;
		T=0;
		for (i=0;i<(int)p.size();i++) {
			if(p[i]>pMax){	
				pMax=p[i];
			}
		}
		for (i=0;i<(int)p.size();i++) {
			T += p[i];
		}
		//T = int(floor((T - pMax)/m) + pMax);
		//T = int(floor(T/m + pMax*m)); //2850 para wt100-2m-1 
		T = int(floor((T - pMax)/m)+pMax); //2700 para wt100-2m-1
 	}
};