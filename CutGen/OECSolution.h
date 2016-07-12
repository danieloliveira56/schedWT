#pragma once

#define VAR_RELEVANTE 0
#include <iostream>

class ArcVariable {
	public:
		int i;
		int j;
		int t;
		double value;

	ArcVariable(int indice_i,int indice_j,int indice_t,double valor_var) {
		i=indice_i;
		j=indice_j;
		t=indice_t;
		value=valor_var;
	}

	string texto()
	{
		return "x" + to_string(i) + "_" + to_string(j) + "_" + to_string(t);
	}

	#ifndef DANIEL
	bool operator < (const ArcVariable& var) {
		return (t < var.t);    
    }
	#else
	bool operator < (ArcVariable& var) {
		return (t < var.t);       
    }
	#endif
  };

 class link_var_S {
 public:
	 int indice;
	 bool delta;

	 link_var_S(int indice_, bool delta_) {
		indice = indice_;
		delta = delta_;
	}

	 bool operator < (const link_var_S& rhs) const {
        return (indice < rhs.indice);
    }
 };
		

class Solution 
{
  public:
  int n;								//numero de tarefas para a solucao
  int T;								//tempo maximo da restricado para indexar
  string nome;							//nome da solucao
  vector <int> nFrag;					//guarda o numero de arcs x_i_j diferentes de 0 para a ArcVariable j
  //vector<vector<vector<double> > > x;	//guarda os valores de todas as arcs
  vector<ArcVariable> arcs;			//guarda todas as variaveis
  vector<vector<int> > indiceArcVariable;	//guarda os indices de todas as arcs indexadas a cada ArcVariable
  vector<int> valoresT;					//vetor guardando todos os valores de t presentes na solução, com repeticao
  vector<vector<link_var_S> > correspondencia_S;
  vector<vector<int> > entra;
  vector<vector<int> > sai;
  vector<vector<int> > adjMatrix;
  vector<vector<int> > adjMatrix_directed;
  vector<vector<double> > flowMatrix_directed;
  vector< vector<int> > adjVertex;
  vector< vector<double> > adjValue;

	Solution(Instance inst, LpSolution sol) 
	{	

		n = inst.n;
		T = inst.T;
		nFrag.resize(n+1);

		indiceArcVariable.resize(n+1);
		entra.resize(n+1);
	    sai.resize(n+1);

		for (int i = 0; i < sol.numArcsCap; i++)
		{
			leArcVariable(sol.arcsCap[i].j, sol.arcsCap[i].i, sol.arcsCap[i].d, sol.arcsCap[i].value);
		}

		ordena();
		gera_entra_sai();
		populaMatrizAdjacencias();
	}



  void popula_correspondencia_S(int tamanho_S) 
  {
	  //cout << "popula_correspondencia_S" << endl;
	 // cout << tamanho_S << endl;
	  correspondencia_S.resize(tamanho_S+1);


	  for (int i = 0; i < (int)arcs.size(); i++) {
		  correspondencia_S[arcs[i].i].push_back(link_var_S(i,true)); //true = S+		  
		  correspondencia_S[arcs[i].j].push_back(link_var_S(i,false)); //false = S-

//cout << arcs[i].i << ":" << i << "S+" << endl;
		 // cout << arcs[i].j << ":" << i << "S-" << endl;
	  }
  }
      
  void ordena() {
	  sort(arcs.begin(),arcs.end());
  }	
  
	void leArcVariable(int i,int j, int t, double valor) {
		//x[i][j][t]=valor;
		 nFrag[j]++;

		 //cria e insere variável no vetor de arcs
		 arcs.push_back(ArcVariable(i,j,t,valor));

		 //indexa a variável
		 indiceArcVariable[i].push_back(arcs.size()-1);
		 indiceArcVariable[j].push_back(arcs.size()-1);
		 		 
		 //insere o valor de t no vetor de t's
		 if (t>0)
			valoresT.push_back(t);
      }		

	void gera_entra_sai() {
		for (int t = 0; t < (int)arcs.size(); t++) {			
			//sai de i e entra em j
			sai[arcs[t].i].push_back(t);
			entra[arcs[t].j].push_back(t);
		}
	}

	void populaMatrizAdjacencias() 
	{
		double EccSetConnectEps = 0.01;

		vector<vector<int> > aux(n+1,vector<int>(n+1,0));
		adjMatrix = aux;		
		adjMatrix_directed = aux;		
		vector<vector<double> > aux2(n+1,vector<double>(n+1,0));		
		flowMatrix_directed = aux2;
		for (int t = 0; t < (int)arcs.size(); t++) 
		{
			if (arcs[t].value > VAR_RELEVANTE)
			{
				adjMatrix[arcs[t].i][arcs[t].j] = 1;
				adjMatrix[arcs[t].j][arcs[t].i] = 1;
				
				adjMatrix_directed[arcs[t].i][arcs[t].j] = 1;
				
				flowMatrix_directed[arcs[t].i][arcs[t].j]+= arcs[t].value;
			}
		}
		// Create the symmetric adjacency lists of vertices
	   adjVertex.resize(n+1);
	   adjValue.resize(n+1);
	   for (int t = 0; t < (int)arcs.size(); t++) 
	   {
		  if (arcs[t].value < EccSetConnectEps) continue;
		  int l;
		  int i = arcs[t].i;
		  int j = arcs[t].j;
		  if ((i == 0) || (j == 0)) continue;
		  for (l = 0; l < (int)adjVertex[i].size(); l++)
			 if (adjVertex[i][l] == j) break;
		  if (l == (int)adjVertex[i].size())
		  {
			 adjVertex[i].push_back(j);
			 adjValue[i].push_back(arcs[t].value);
		  }
		  else
		  {
			 adjValue[i][l] += arcs[t].value;
		  }
		  for (l = 0; l < (int)adjVertex[j].size(); l++)
			 if (adjVertex[j][l] == i) break;
		  if (l == (int)adjVertex[j].size())
		  {
			 adjVertex[j].push_back(i);
			 adjValue[j].push_back(arcs[t].value);
		  }
		  else
		  {
			 adjValue[j][l] += arcs[t].value;
		  } 
		}
	}

};
