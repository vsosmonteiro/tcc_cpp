#include <iostream>
#include <stdlib.h>
#include <ilcplex/ilocplex.h>
#include <math.h>
#include <stdbool.h>
#include "PCInstanceReader.hpp"
// #include "../reader.hpp"
//#include <bits/stdc++.h>
//#include <chrono>




using namespace std;

void subset(bool *visited, int *lenneigb, int v)
{
	if (!visited[v]){
		visited[v]=true;
		for (int i=0; i<lenneigb[v]; i++){
			subset(visited, lenneigb, g[v][i]);
		}
	}
}

ILOLAZYCONSTRAINTCALLBACK2(CB_SubSet,
					IloArray<IloIntVarArray>, y, 
					int*, n_cut)
{			
	
	IloEnv env  = getEnv();
	
	//int solution[c][c-1];
	bool visited[c];
	int lenneigb[c];
	int cl1, cl2;
  
	for (int k=0; k<c; k++){
		visited[k]=false;
		lenneigb[k]=0;
	}
	
	for (int i=0; i<n-1; i++){
		for (int j=i+1; j<n; j++){
			if (graph[i][j]<0) continue;
			if (abs(getValue(y[i][j])-1)<0.1){
				cl1=cl[i];
				cl2=cl[j];
				g[cl1][lenneigb[cl1]]=cl2;
				g[cl2][lenneigb[cl2]]=cl1;
				lenneigb[cl1]++;
				lenneigb[cl2]++;
			}
		}
	}
	
	subset(visited, lenneigb, 0);
	
	for (int k=0; k<c; k++){
		if (!visited[k]){//Testa se algum vértice não foi visitado
			(*n_cut)++; //Aumenta o número de cortes0
			IloExpr corteSubSet(env);
			int cont=0;
			for (int e=0; e<m; e++){
				cl1=cl[v1[e]];
				cl2=cl[v2[e]];
				if(visited[cl1] != visited[cl2]){
					corteSubSet += y[v1[e]][v2[e]];
					cont++;
				}
			}
			add(corteSubSet >= IloNum(1)).end();
			corteSubSet.end();
			//cout << "\n"<<cont<<" arestas adicionadas no corte\n";
			break;
		}
	}	
  
	double ub     = getIncumbentObjValue();  // retorna a melhor solucao inteira (limite primal)
	double lb     = getBestObjValue();       // retorna o melhor limite dual
	double rlx    = getObjValue();           // quando chamada dentro do callback, retorna o valor da relaxacao do noh
	double nNodes = getNremainingNodes();    // retorna o numero restante de nos a serem analisados
	cout<<"***** LAZY CONSTRAINT: "<<"relax="<<rlx<<"\t bounds="<<lb<<"<->"<<ub<<"\t n_rest="<<nNodes<< "\t#lazy_con="<<*n_cut<<endl;

}

ILOLAZYCONSTRAINTCALLBACK2(CB_SubCiclo,
					IloArray<IloIntVarArray>, y, 
					int*, n_cut)
{			
	
	IloEnv env  = getEnv();
	
	//int solution[c][c-1];
	bool visited[c];
	int lenneigb[c];
	int cl1, cl2;
  
	for (int k=0; k<c; k++){
		visited[k]=false; //Todos os vértices são marcados como não visitados
		lenneigb[k]=0; //O número de vizinhos de cada vértices inicialmente é zero
	}
	
	//Monta a subsolução com c vértices
	for (int i=0; i<n-1; i++){
		for (int j=i+1; j<n; j++){
			if (graph[i][j]==-1) continue;
			if (abs(getValue(y[i][j])-1)<0.1){
				cl1=cl[i];
				cl2=cl[j];
				g[cl1][lenneigb[cl1]]=cl2;
				g[cl2][lenneigb[cl2]]=cl1;
				lenneigb[cl1]++;
				lenneigb[cl2]++;
			}
		}
	}
	
	subset(visited, lenneigb, 0);
	
	for (int k=0; k<c; k++){
		if (!visited[k]){//Testa se algum vértice não foi visitado
			(*n_cut)++; //Aumenta o número de cortes
			
			int nv=0;
			
			for (int i=0; i<c; i++){
				if (visited[i]) nv++;
			}
			
			IloExpr corteSubCiclo(env);
			IloExpr corteSubCiclo2(env);
			//IloExpr corteSet(env);
			int cont=0;
			for (int e=0; e<m; e++){
				cl1=cl[v1[e]];
				cl2=cl[v2[e]];
				if(visited[cl1] && visited[cl2]){
					corteSubCiclo += y[v1[e]][v2[e]];
					cont++;
				}
				if(!visited[cl1] && !visited[cl2]){
					corteSubCiclo2 += y[v1[e]][v2[e]];
				}
				//if(visited[cl1] != !visited[cl2]){
				//	corteSet += y[v1[e]][v2[e]];
				//}
			}
			
			
			add(corteSubCiclo <= IloNum(nv-1)).end();
			corteSubCiclo.end();
			add(corteSubCiclo2 <= IloNum(c-nv-1)).end();
			corteSubCiclo2.end();
			//add(corteSet >=1).end();
			break;
		}
	}	
  
	double ub     = getIncumbentObjValue();  // retorna a melhor solucao inteira (limite primal)
	double lb     = getBestObjValue();       // retorna o melhor limite dual
	double rlx    = getObjValue();           // quando chamada dentro do callback, retorna o valor da relaxacao do noh
	double nNodes = getNremainingNodes();    // retorna o numero restante de nos a serem analisados
	cout<<"***** LAZY CONSTRAINT: "<<"relax="<<rlx<<"\t bounds="<<lb<<"<->"<<ub<<"\t n_rest="<<nNodes<< "\t#lazy_con="<<*n_cut<<endl;

}


int main()
{
	cout << "Inicio\n";
	read();
	cout << "Leitura concluída.\n\n\n";
	PCInstance inst;
	
	//MODELO
	IloEnv env;
    IloModel model(env, "Problema da Árvore Geradora Mínima Generalizada");
    
    // Variáveis de decisão dos vértices: x[]
    IloIntVarArray x(env, n, 0, 1);
    
    //Variáveis de decisão das arestas: y[]
    //IloIntVarArray y(env, 2*m, 0, 1);//Uma variável por arco
    IloArray <IloIntVarArray> y(env, n);

    for(int i=0; i<n; i++)
	{
		IloIntVarArray vetor(env, n, 0, 1);
	    y[i] = vetor;
	}
	
	//Função Objetivo
	IloExpr custo(env);
	for (int i=0; i<n-1; i++){
		for(int j=i; j<n; j++){
			if (graph[i][j] > -1){
				custo += graph[i][j] * y[i][j];
			}
		}
	}
	IloExpr prize(env);
	for (int i=0; i<n; i++){
		prize += x[i] * pr[i];
	}
	
	model.add(IloMinimize(env, custo-prize));
	
	//RESTRIÇÕES
	
	IloExpr sumx(env);
	for (int i=0; i<n; i++){
		sumx += x[i];
	}
	model.add(sumx == IloNum(c));//Restrição do somatório de x[i]
	
	IloExpr sumy(env);
	for (int i=0; i<n-1; i++){
		for (int j=i; j<n; j++){
			if (graph[i][j] > -1){
				sumy += y[i][j];
			}
		}
	}
	model.add(sumy == IloNum(c-1)); //Restrição do somatório de y[i][j]
	
	for (int i=0; i<n-1; i++){
		for (int j=i; j<n; j++){
			if (graph[i][j]>-1){
				model.add(y[i][j] <= x[i]);//Restrições do uso de um arco apenas se os vértices forem usados
				model.add(y[i][j] <= x[j]);
			}
		}
	}
	
	for (int k=0; k<c; k++){
		IloExpr sumcluster(env);
		for (int i=0; i<n; i++){
			if (cl[i] == k){
				sumcluster += x[i];
			}
		}
		model.add(sumcluster >= inst); //Restrição do uso de um vértice por cluster
	}
	
	/*
	//Restrições de grau :: Rrestrição opicional
	for (int k=0; k<c; k++){
		IloExpr grau(env);
		for (int e=0; e<m; e++){
			if (cl[v1[e]] == k || cl[v2[e]] == k){
        			grau += y[v1[e]][v2[e]];
			}
		}
		model.add(grau >= IloNum(1));//Obriga que o vértice selecionado no cluster tenha pelo menos grau 1
	}*/
	
	
	//Máximo uma aresta entre dois clusters :: Restrição opicional
	for (int i=0; i<c; i++){
		for (int j=0; j<c; j++){
        		if (i==j) continue;
			IloExpr aux(env);
			for (int e=0; e<m; e++){
				if ((cl[v1[e]]==i && cl[v2[e]]==j) || (cl[v1[e]]==j && cl[v2[e]]==i)){
        				aux += y[v1[e]][v2[e]];
        			}
			}
			model.add(aux <=IloNum(1));//No máximo uma aresta entre dois clusters
		}
	}
	
	
    IloCplex cplex(model); //Alocação do solver
    
    //Criar a estrutura que será parâmetro do USERCALLBACK 
    int n_cut=0;
    cplex.use(CB_SubCiclo(env, y, &n_cut)); //CALLBACK que elimina os subciclos
    //cplex.use(CB_SubSet(env, y, &n_cut)); //CALLBACK que elimina os subciclos
    
    cplex.solve(); //Resolução
    
    double focost = cplex.getObjValue();
    
    cout << "Custo da FO: " << focost << endl << "Vértices: ";
    
    int sol[c], indx=0;
    for (int i=0; i<n; i++){
		if (abs(cplex.getValue(x[i])-1)<0.1){
			cout << i << "\t";
			sol[indx] = cl[i];
			indx++;
		}
	}
	cout << endl << "Clusters: ";
	for (int i=0; i<c; i++){
		cout << sol[i] << "\t";
	}
	
    cout << endl << "Arestas: ";
    for (int i=0; i<n-1; i++){
		for (int j=i; j<n; j++){
			if (graph[i][j]>-1 && abs(cplex.getValue(y[i][j])-1)<0.1) cout << "("<<i<<","<<j<<") ";
		}
	}
    cout <<endl;
	env.end();
	return 0;
}
