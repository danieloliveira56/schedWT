#pragma once

#include "OECSet.h"
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

const int igualdade = 0;
const int menor = 1;
const int maior = 2;

class coeficiente {
public:
	int i;
	int j;
	int t;
	double valor;
	double valor_var;

	coeficiente (int indice_i, int indice_j, int indice_t, double valor_coeficiente, double valor_var_) 
	{
		i = indice_i;
		j = indice_j;
		t = indice_t;
		valor = valor_coeficiente;
		valor_var = valor_var_;
	}

	int size () {
		return sizeof(i)*3 + sizeof(valor);
	}
};

class restricao {
public:
	vector<coeficiente> coeficientes;
	string nome;
	double rhs;
	double lhs;
	int tipo;
	int tamanho;
	string assinatura;

	restricao() {
		rhs = 0;
		lhs = 0;
		tamanho=0;
	}
	
	restricao(const restricao& clone) 
	{
		rhs = clone.rhs;
		lhs = clone.lhs;
		tamanho=clone.tamanho;
		tipo = clone.tipo;
		assinatura = clone.assinatura;
		nome = clone.nome;
		coeficientes = clone.coeficientes;
	}

	int size() {
//		return sizeof(rhs)*2 + coeficientes.capacity()*coeficientes[0].size() + sizeof(tipo)*2 + nome.size();
		return 0;
	}

	void resetCoeficientes()
	{
		coeficientes.clear();
		tamanho = 0;
	}

	void setCoeficiente (int &i, int &j, int &t, int valor, double &valor_var) {
		double dblvalor = valor;
		coeficientes.push_back(coeficiente(i,j,t,dblvalor,valor_var));
		tamanho++;
	}

	void setCoeficiente (int &i, int &j, int &t, double &valor, double &valor_var) {
		coeficientes.push_back(coeficiente(i,j,t,valor,valor_var));
		tamanho++;
	}

	void setCoeficiente (coeficiente coef) {
		coeficientes.push_back(coeficiente(coef.i,coef.j,coef.t,coef.valor,coef.valor_var));
		tamanho++;
	}

	void setTipo (int tipo_restricao) {
		tipo = tipo_restricao;
	}

	void defineNome(string nomeRestricao) {
	  nome = nomeRestricao;
	}

	double mediaCoeficientes() {
		double media=0;
		for (int i=0; i < (int)coeficientes.size(); i++) {
			media += coeficientes[i].valor;
		}
		return media/coeficientes.size();
	}

	double gap() {
		return lhs - rhs;
	}
	
	string imprime(bool imprimeValor = false){
		stringstream msg;
		stringstream aux;

		for (int i=0;i<tamanho;i++) {
			aux.str(std::string());
			aux << setprecision(5) << fabs(coeficientes[i].valor);
			if (coeficientes[i].valor > 0) {
				msg << " + " << aux.str() << " x" << coeficientes[i].i << "_" << coeficientes[i].j << "_" << coeficientes[i].t << " ";
				if (imprimeValor)
					msg << "(" << coeficientes[i].valor_var << ")";
			} else {
				msg << " - " << aux.str() << " x" << coeficientes[i].i << "_" << coeficientes[i].j << "_" << coeficientes[i].t << " ";
				if (imprimeValor)
					msg << "(" << coeficientes[i].valor_var << ")";
			}	
		}

		switch (tipo)
		{
		case igualdade:
			msg <<  " = ";
			break;
		case maior:
			msg <<  " > ";
			break;
		case menor:
			msg <<  " < ";
			break;
		}
		
		aux.str(std::string());
		aux << setprecision(5) << rhs;
		msg << aux.str();

		return msg.str();
	}
};


class InteDouble
{
public:
	int valor_int;
	double valor_double;

	InteDouble(int valor_int_, double valor_double_)
	{
		valor_int = valor_int_;
		valor_double = valor_double_;
	}

	InteDouble& operator += (double incrementa)
	{
	  valor_double += incrementa;
	  return (*this) ;
	}

	#ifdef _WIN32
		bool operator < (InteDouble& compara) {
			return (this->valor_double < compara.valor_double);     
		}
	#else
		bool operator < (const InteDouble& compara) {
			return (valor_double < compara.valor_double);       
		}
#endif

};
	
class CoefPathCut
{
public:
	int a, b, t;
	double valor;
	bool delta_mais;
	int distInteiro;
	double coefFrac;
	int coefFareyArred;

	
	CoefPathCut(){};

	int size () {
		return sizeof(a)*5 + sizeof(valor)*2 + sizeof(delta_mais);
	}

	CoefPathCut(int a_, int b_, int t_, double valor_, bool delta_mais_)
	{
		a = a_;
		b = b_;
		t = t_;
		valor = valor_;
		delta_mais = delta_mais_;
		coefFrac = (a*t)/double(b);

		if (a*t % b == 0)
		{
			if (delta_mais)
			{
				distInteiro = 0;
			}
			else
			{
				distInteiro = b;
			}
		}


		else
		{		
			///int aux = ((a*t)+b-1)/b; 
			distInteiro = b - (a*t % b);
		}
	}

	string imprime()
	{
		stringstream texto;

		if (delta_mais)
		{
			texto << "+Y(" << t << "*" << a << "/" << b << " = " << coefFrac << ")(distInteiro=" << distInteiro << ")(Y=" << valor << ")" << endl;
		}
		else
		{
			texto << "-Z(" << t << "*" << a << "/" << b << " = " << coefFrac << ")(distInteiro=" << distInteiro << ")(Z=" << valor << ")" << endl;
		}

		return texto.str();
	}
	

	#ifdef _WIN32
	bool operator < (CoefPathCut& coef) {
		if ( (this->distInteiro == coef.distInteiro) && (this->delta_mais != coef.delta_mais) )
		{
			//se o lado esquerdo do comparador < for do tipo delta+, ele perde depois de delta-
			return !this->delta_mais;
		}
		else
		{
			return (this->distInteiro < coef.distInteiro);
		}
    }
	#else
	bool operator < (const CoefPathCut& coef) {
		if ( (distInteiro == coef.distInteiro) && (this->delta_mais != coef.delta_mais) )
		{
			//se o lado esquerdo do comparador < for do tipo delta+, ele perde depois de delta-
			return !delta_mais;
		}
		else
		{
			return (distInteiro < coef.distInteiro);
		} 
    }
	#endif

};

class Individual 
{
public:
	string nome;
	restricao corte;
	conjuntoS S;
	int a;
	int b;
	int s;
	int t;
	string assinatura;
	string assinatura_completa;
	int Tcount;
	int i_Farey;
	double potencial;
	double custoAtual;
	double custoOtimo;
	vector<CoefPathCut> coeficientesPathCut;
	double ganho_s;
	int indice_s;

	//finds the best t for individual's current S
	void calculateOEC(Solution&  sol, Instance&  inst);
	double evaluateInsertionOEC(Solution& sol, Instance&  inst, int newJob);
	double evaluateRemovalOEC(Solution& sol, Instance&  inst, int oldJob);

	Individual(restricao & corte_indiv, conjuntoS & S_Individual, int valor_a, int valor_b, int valor_i_Farey) 
		:corte()
	{
		corte = restricao(corte_indiv);
		S = conjuntoS(S_Individual);
		a = valor_a;
		b = valor_b;
		s = 0;
		i_Farey = valor_i_Farey;
		Tcount = S_Individual.Tcount;
		buscaAssinatura();
		custoAtual = 0;
		custoOtimo = 0;
		potencial = 0;
		t = 1;
	}

	Individual( Instance inst, vector<int> elementos ) 
	{
		S = conjuntoS(inst, elementos);
		a = 1;
		b = 2;
		s = 0;
		i_Farey = 0;
		Tcount = 0;
		custoAtual = 0;
		custoOtimo = 0;
		potencial = 0;
		t = 1;
	}
	
	Individual(const Individual& clone) 
	{
		corte = restricao(clone.corte);
		S = conjuntoS(clone.S);
		a = clone.a;
		b = clone.b;
		s = clone.s;
		t = clone.t;
		i_Farey = clone.i_Farey;
		Tcount = clone.Tcount;
		buscaAssinatura();
		custoAtual = 0;
		custoOtimo = 0;
		potencial = 0;
		ganho_s = clone.ganho_s;
		indice_s = clone.indice_s; 
	}

	double calculaDistancia (Individual comparado) {
		double distancia = 0;
		double dif_r = a/(double)b- comparado.a/(double)comparado.b;
		distancia = fabs(dif_r)*2;
		for (int i = 0;i < (int)S.tamanhoMax;i++) 
		{
			if (S.pertenceS[i] != comparado.S.pertenceS[i]) 
			{
				distancia += 1;
			}
		}
		return distancia;
	}

	void buscaAssinatura() {
		std::sort(S.elementos.begin(),S.elementos.end());
		assinatura = "";
		for (int i = 0;i<S.tamanho;i++) {
			assinatura = assinatura + "[" + to_string(S.elementos[i]) + "]";
		}
		//comentado em 15-01-2015 retirando o r da assinatura
		//assinatura = assinatura + to_string(ceil(a/double(b)*1000000)/1000000);

		//inserido em 15-01-2015 inserindo o a/b na assinatura
		assinatura_completa = assinatura + to_string(a) + "/" + to_string(b);
	}

	
	//verifica se o conjunto S e igual iguais
	bool conjuntoSigual (Individual& ind) {
		for (int i = 0; i < (int)S.pertenceS.size(); i++)
		{
			if (S.pertenceS[i] != ind.S.pertenceS[i])	
				return false;			
		}		     

		return true;
	}

	#ifdef _WIN32
	//verifica se o conjunto S, a, b e s sao iguais
	bool operator == (Individual& ind) {
	#else
	bool operator == (const Individual& ind) {
	#endif
			if (t != ind.t)		
				return false;
		
		for (int i = 0; i < (int)S.pertenceS.size(); i++)
		{
			if (S.pertenceS[i] != ind.S.pertenceS[i])	
				return false;			
		}		     

		return true;
    }

	bool comparaS(Individual& ind) {
		for (int i = 0; i < (int)S.pertenceS.size(); i++)
		{
			if (S.pertenceS[i] != ind.S.pertenceS[i])	
				return false;			
		}		     

		return true;
    }



	#ifdef _WIN32
	bool operator < (Individual& ind) {
		if (abs(gap() - ind.gap()) > 0.001) {
			return (gap() < ind.gap());
		} else {
			return (S.Tcount < ind.S.Tcount);
			//return (this->gap2() < ind.gap2());
		}        
    }
	#else
	bool operator < (const Individual& ind) {
		if (abs(gap() - ind.gap()) > 0.001) {
			return (gap() < ind.gap());
		} else {
			return (S.Tcount < ind.S.Tcount);
			//return (gap2() < ind.gap2());
		}        
    }
	#endif
	

	int size () {
		return sizeof(a)*7 + sizeof(potencial)*4 + assinatura.size() + assinatura_completa.size() 
			+ nome.size() + corte.size() + S.size();// + coeficientesPathCut[0].size()*coeficientesPathCut.capacity();
	}

	double gap() {
		return corte.gap();

	}

	void nomeiaFilho (Individual pai, Individual mae) {
		nome = "(Filho de " + pai.nome + " e " + mae.nome + ")";
	}

	void preparaCalculo_s() 
	{
		coeficientesPathCut.resize(0);
		coeficientesPathCut.reserve(S.Y.size() + S.Z.size());

		for (pair<int, double> p : S.Y) //Delta+
		{
			coeficientesPathCut.push_back(CoefPathCut(a,b,p.first,p.second,true));
		}
		for (pair<int, double> p : S.Z) //Delta-
		{
			coeficientesPathCut.push_back(CoefPathCut(a,b,p.first,p.second,false));
		}

		sort(coeficientesPathCut.begin(),coeficientesPathCut.end());
	}

	string imprimeAnalise_s()
	{
		stringstream texto;
		preparaCalculo_s(); 

		texto << "Ganho com s: " << ganho_s << endl;
		texto << "Indice: " << indice_s << endl;
		texto << "s: " << s << endl;
		for (int i = 0; i < (int)coeficientesPathCut.size(); i++)
		{
			texto << coeficientesPathCut[i].imprime();
		}

		return texto.str();
	}

	double buscaMelhoria_s()
	{
		double soma=0;
		ganho_s = 0;
		indice_s = -1;
		s = 0;
		for (int i = 0; i < (int)coeficientesPathCut.size(); i++)
		{
			if (coeficientesPathCut[i].delta_mais)
			{
				soma += coeficientesPathCut[i].valor;
			}
			else
			{
				soma -= coeficientesPathCut[i].valor;
			}
			if (soma < ganho_s) //&& ( (i + 1 == coeficientesPathCut.size() ) || (coeficientesPathCut[i].distInteiro < coeficientesPathCut[i+1].distInteiro + 0.01) ) )
			{
				ganho_s = soma;
				indice_s = i;
			}
		}
		if (ganho_s < -0.05000001)
		{
			s = coeficientesPathCut[indice_s].distInteiro;
		}
		return ganho_s;
	}

	string imprimeCorte_s() {
		stringstream corte;
		for (int i = 0; i < (int)coeficientesPathCut.size(); i++)
		{
			if (i <= indice_s)
			{
				if (coeficientesPathCut[i].delta_mais)
				{
					corte << " + " <<  coeficientesPathCut[i].coefFareyArred +1 <<  "Y_" << coeficientesPathCut[i].t; 
				}
				else
				{
					corte << " - " <<  coeficientesPathCut[i].coefFareyArred + 1 <<  "Z_" << coeficientesPathCut[i].t; 
				}
			}
			else
			{
				if (coeficientesPathCut[i].delta_mais)
				{
					corte << " + " <<  coeficientesPathCut[i].coefFareyArred <<  "Y_" << coeficientesPathCut[i].t; 
				}
				else
				{
					corte << " - " <<  coeficientesPathCut[i].coefFareyArred <<  "Z_" << coeficientesPathCut[i].t; 
				}
			}
		}
		return corte.str();
	}

	string imprimeECCYZ(int a = 0, int b = 0)
	{
		stringstream ECC;
		for (pair<int, double> p : S.Y) //Delta+
		{
			
			ECC << " + Y_" << p.first << "*" << p.second;
			if (b != 0)
				ECC << "(" << (b - p.first*a % b)/double(b) << ")";
			ECC << endl;
		}
		for (pair<int, double> p : S.Z) //Delta-
		{
			ECC << " - Z_" << p.first << "*" << p.second;
			if (b != 0)
				ECC << "(" << (p.first*a % b)/double(b) << ")";
			ECC << endl;
		}
		return ECC.str();
	}

	double desvioMedio() {return 0;}

	bool equalS(const Individual& ind) 
	{
		if ( S.elementos.size() != ind.S.elementos.size() )	
			return false;

		for (int i = 0; i < (int)S.pertenceS.size(); i++)
		{
			if (S.pertenceS[i] != ind.S.pertenceS[i])	
				return false;			
		}
		return true;
	}

};

bool lexicographical_compareInd(const Individual& ind1, const Individual& ind2) {
   return lexicographical_compare( ind1.S.elementos.begin(), ind1.S.elementos.end(), ind2.S.elementos.begin(), ind2.S.elementos.end() );
}

bool comparaPotencial(const Individual& ind1, const Individual& ind2) {
   return ind1.potencial > ind2.potencial;
}