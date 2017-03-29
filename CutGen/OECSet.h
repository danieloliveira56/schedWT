#pragma once

#include <math.h>
#include <cmath>
#include<map>
#include<vector>
#include<sstream>
#include "OECInstance.h"
#include "OECSolution.h"

using namespace std;

class RelevanciaT {
	public:
	int t;
	double relevancia;

	RelevanciaT(int t_) {
		t = t_;
		relevancia = 0;
	}

	RelevanciaT(int t_, double relevancia_) {
		t = t_;
		relevancia = relevancia_;
	}

	RelevanciaT() {
		t = 0;
		relevancia = 0;
	}

	bool operator < (const RelevanciaT& rhs) const {
		return (relevancia > rhs.relevancia); //operador invertido para ordenar relevancia decrescente
	}
};

class r {
public:
	int a;
	int b;

	r (int valor_a, int valor_b) {
		a = valor_a;
		b = valor_b;
	}

	double valor() {
		return a/double(b);
	}

	bool operator < (const r& ind) const
    {
        return (a/double(b) < ind.a/double(ind.b));
    }

};

class VarAgrup
{
public:
	VarAgrup();
	void adicionaArcVariable(int, double);
private:

};

VarAgrup::VarAgrup()
{
}

void VarAgrup::adicionaArcVariable(int t, double valor)
{

}

class conjuntoS {
public:
	vector<int> elementos;
	vector<int> anti_elementos; //vetor dos elementos que nao pertencem a S
	vector<int> combinacao;	
	vector<bool> pertenceS;

	vector<int> pontos_articulacao;
	vector<bool> ehPontoArticulacao;
	vector<vector<int> > comb;
	int tamanhoComb;
	int totalComb;
	vector<int> duracao;
	int pTotal;
	int tamanhoMax;
	bool cheio;
	int tamanho;
	int c;
	int Tcount;
	int TcountZ, TcountY;
	double somaarcsZ, somaarcsY;

	map<int, double> Y;
	map<int, double> Z;

	vector<double> Y_vec;
	vector<double> Z_vec;
	
	vector<RelevanciaT> relevanciaTs;
	vector<RelevanciaT> relevanciaTsY;
	vector<RelevanciaT> relevanciaTsZ;

	double somaarcs;
	int contVarNaoRaiz;

	double capIn;
	double capOut;
	double capAcum;
		int size() {

		return sizeof(c)*(
			10 + 
			elementos.capacity() +
			anti_elementos.capacity() +
			combinacao.capacity() +
			comb.capacity() +
			duracao.capacity() +
			pontos_articulacao.capacity() +
			Y.size() +
			Z.size() +
			relevanciaTs.capacity() +
			relevanciaTsY.capacity() +
			relevanciaTsZ.capacity()
			) +
			sizeof(cheio)*(
			1 + 
			ehPontoArticulacao.capacity() +
			pertenceS.capacity()
			) +
			sizeof(somaarcs)*(
			6 +
			Y.size() +
			Z.size() +
			relevanciaTs.capacity() +
			relevanciaTsY.capacity() +
			relevanciaTsZ.capacity() 
			);

	}

	string imprime() {	

		stringstream conjunto;

		std::sort(elementos.end(),elementos.end());

		for (int i=0;i<(int)elementos.size();i++){
			conjunto <<  "[" << elementos[i] << "]";
		}

		return conjunto.str();
	}

	conjuntoS (Instance& inst) 
	{
		pTotal=0;
		tamanho=0;
		cheio=false;
		Tcount = 0;

		tamanhoMax=inst.n;
		pertenceS.resize(tamanhoMax+1);
		duracao.resize(tamanhoMax+1);
		anti_elementos.resize(tamanhoMax);
		ehPontoArticulacao.resize(tamanhoMax+1);
		
		for (int i = 0; i < tamanhoMax+1; i++)
		{
			duracao[i] = inst.p[i];
		}
		//comentado em 22/01/2015 pois 0 deve ser um anti-elemento
		/*for (int i=0;i<tamanhoMax;i++) {
			anti_elementos[i] = i+1;
		}*/
		for (int i = 0; i < tamanhoMax; i++)
		{
			anti_elementos[i] = i+1;
		}
	}

	conjuntoS () {}; 
				
	conjuntoS (Instance& inst, vector<int> elementos_) 
	{
		pTotal=0;
		tamanho=0;
		cheio=false;
		Tcount = 0;

		tamanhoMax=inst.n;
		pertenceS.resize(tamanhoMax+1);
		duracao.resize(tamanhoMax+1);
		anti_elementos.resize(tamanhoMax);
		ehPontoArticulacao.resize(tamanhoMax+1);
		Y_vec.resize(tamanhoMax+1);
		Z_vec.resize(tamanhoMax+1);
		
		for (int i = 0; i < tamanhoMax+1; i++)
		{
			duracao[i] = inst.p[i];
		}

		for (int i = 0; i < tamanhoMax; i++)
		{
			anti_elementos[i] = i+1;
		}

		for (int i = 0; i < (int)elementos_.size(); i++)
		{
			incluiElemento(elementos_[i]);
		}
	}

	conjuntoS (const conjuntoS& clone) 
	{
		pTotal=clone.pTotal;
		tamanho=clone.tamanho;
		cheio=clone.cheio;
		Tcount = clone.Tcount;
		TcountY = clone.TcountY;
		TcountZ = clone.TcountZ;
		tamanhoMax=clone.tamanhoMax;
				
		elementos=clone.elementos;
		anti_elementos=clone.anti_elementos;
		combinacao=clone.combinacao;
		pertenceS=clone.pertenceS;
		ehPontoArticulacao=clone.ehPontoArticulacao;
		pontos_articulacao = clone.pontos_articulacao;
		comb=clone.comb;
		tamanhoComb=clone.tamanhoComb;
		totalComb=clone.totalComb;
		duracao=clone.duracao;
		pTotal=clone.pTotal;
		c=clone.c;
		Y=clone.Y;
		Z=clone.Z;
		relevanciaTsZ=clone.relevanciaTsZ;
		relevanciaTsY=clone.relevanciaTsY;
		contVarNaoRaiz=clone.contVarNaoRaiz;
		somaarcsY=clone.somaarcsY;
		somaarcsZ=clone.somaarcsZ;
		capAcum=clone.capAcum;
		capIn=clone.capIn;
		capOut=clone.capOut;
		
		Y_vec=clone.Y_vec;
		Z_vec=clone.Z_vec;
	}

	void esvazia(void) {
		cheio=false;
		tamanho=0;
		pTotal=0;
		elementos.resize(tamanho);
		fill(pertenceS.begin(), pertenceS.end(), false);
		anti_elementos.resize(tamanhoMax);

		for (int i = 0; i < tamanhoMax; i++)
		{
			anti_elementos[i] = i+1;
		}
	}
		
	void incluiElemento(int elemento) {
		elementos.push_back(elemento);
		pertenceS[elemento]=true;
		pTotal+=duracao[elemento];
		tamanho++;
		if (tamanho==tamanhoMax){cheio=true;}

		anti_elementos.erase(std::remove(anti_elementos.begin(), anti_elementos.end(), elemento), anti_elementos.end());
	}

	void excluiElemento(int elemento) {
		int i;
		pertenceS[elemento]=false;
		pTotal-=duracao[elemento];

		for (i=0;i<tamanho;i++) {
			if (elementos[i]==elemento) {
				elementos.erase(elementos.begin()+i);
				break;
			}
		}
		anti_elementos.push_back(elemento);
		
		tamanho--;
		cheio=false;
	}

	void incluiElementos(vector<int> elementos) {
		int length=elementos.size();
		 
		for(int i=0;i<length; i++) {
			incluiElemento(elementos[i]);
		}		
	}

	long factorial(int n) {
	  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
	}

	unsigned long long choose(unsigned long long n, unsigned long long k) {
    if (k > n) {
        return 0;
    }
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
	}

	void iniciaComb(int qtdElementos) {
		int i;

		esvazia();

		tamanhoComb=qtdElementos;
		//define o tamanho do vector
		totalComb = (int)choose(tamanhoMax,qtdElementos);
		comb.resize(totalComb);
		for (i=0;i<totalComb;i++) {
			comb[i].resize(qtdElementos);
			fill(comb[i].begin(), comb[i].end(), 0);
		}
		
		elementos.resize(0);
		//combinacao.resize(0);
		//for (int i = 0; i < qtdElementos; ++i) {combinacao.push_back(i+1);}
		c=0;
		go(0, qtdElementos,qtdElementos);
		
		c=0;
		for (i=0;i<qtdElementos;i++) {
			incluiElemento(comb[c][i]);
		}
	}

	void go(int offset, int k,int qtdElementos) {
		if (k == 0) {			
			adicionaComb(elementos,qtdElementos);
			c++;
			return;
		}
		for (int i = offset; i <= (int)tamanhoMax - k; ++i) {
			elementos.push_back(i+1);
			go(i+1, k-1,qtdElementos);
			elementos.pop_back();
		}
	}

	void adicionaComb(vector <int> elementosComb, int qtdElementos) {
		for (int i=0;i<qtdElementos;i++) {
			comb[c][i]=elementosComb[i];
		}
	}

	void stepComb(void) {
		c++;
		esvazia();

		if (c<totalComb)
		{
			for (int i=0;i<tamanhoComb;i++) {
				incluiElemento(comb[c][i]);
			}
		} else
		{
			cheio=true;
		}		
	}
	
	bool verificaPertinencia(int i){
		if(find(elementos.begin(), elementos.end(), i) != elementos.end()) {
			return true;
		} else {
			return false;
		}
	}

	void calculaYZ (Solution& sol)
	{	
		Y_vec.clear();
		Z_vec.clear();
		Y_vec.resize(sol.T+1);
		Z_vec.resize(sol.T+1);

		if ((int)elementos.size() <= tamanhoMax/2)
		{
			for (int s = 0; s < (int)elementos.size(); s++)
			{
				int i = elementos[s];
				for (int j = 0; j < (int)sol.entra[i].size(); j++)
				{
					ArcVariable &var = sol.arcs[sol.entra[i][j]];
					if (!pertenceS[var.i]) //Delta-
					{
						Z_vec[var.t] += var.value;
					}
				}
				for (int j = 0; j < (int)sol.sai[i].size(); j++)
				{
					ArcVariable &var = sol.arcs[sol.sai[i][j]]; 
					if (!pertenceS[var.j])  //Delta+
					{  
						Y_vec[var.t] += var.value;
					}
				}
			}
		} 
		else 
		{
			for (int s = 0; s < (int)anti_elementos.size()+1; s++) 
			{
				int i = s < (int)anti_elementos.size() ? anti_elementos[s] : 0;
				for (int j = 0; j < (int)sol.entra[i].size(); j++) 
				{
					ArcVariable &var = sol.arcs[sol.entra[i][j]]; 
					if (pertenceS[var.i]) //Delta+
					{ 
						Y_vec[var.t] += var.value;
					}
				}
				for (int j = 0; j < (int)sol.sai[i].size(); j++) 
				{
					ArcVariable &var = sol.arcs[sol.sai[i][j]];
					if (pertenceS[var.j]) //Delta-
					{  
						Z_vec[var.t] += var.value;
					}
				}
			}
		}
	}

	void calculaYZ_insertion(Solution& sol, int newJob)
	{					
		//for all variables representing flow going into newJob
		for (int j = 0; j < (int)sol.entra[newJob].size(); j++) 
		{
			ArcVariable &var = sol.arcs[sol.entra[newJob][j]]; 
			if (!pertenceS[var.i]) //Delta-
			{  
				Z_vec[var.t] += var.value;
			}
			else //was in Delta+
			{
				Y_vec[var.t] -= var.value;
			}
		}

		//for all variables representing flow going out of newJob
		for (int j = 0; j < (int)sol.sai[newJob].size(); j++)
		{
			ArcVariable &var = sol.arcs[sol.sai[newJob][j]]; 
			if (!pertenceS[var.j])  //Delta+
			{  
				Y_vec[var.t] += var.value;
			}
			else //was in Delta-
			{
				Z_vec[var.t] -= var.value;
			}
		}		
	}

	void calculaYZ_removal(Solution& sol, int oldJob)
	{					
		//for all variables representing flow going into newJob
		for (int j = 0; j < (int)sol.entra[oldJob].size(); j++) 
		{
			ArcVariable &var = sol.arcs[sol.entra[oldJob][j]]; 
			if (!pertenceS[var.i]) //Delta-
			{  
				Z_vec[var.t] -= var.value;
			}
			else //was in Delta+
			{
				Y_vec[var.t] += var.value;
			}
		}

		//for all variables representing flow going out of newJob
		for (int j = 0; j < (int)sol.sai[oldJob].size(); j++)
		{
			ArcVariable &var = sol.arcs[sol.sai[oldJob][j]]; 
			if (!pertenceS[var.j])  //Delta+
			{  
				Y_vec[var.t] -= var.value;
			}
			else //was in Delta-
			{
				Z_vec[var.t] += var.value;
			}
		}		
	}

};
