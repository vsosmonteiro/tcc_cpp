#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include "model.hpp"
#include "reader.hpp"
#include "preprocess.hpp" // Keep commented out as per previous context

ILOSTLBEGIN

// Globals
int n;
int c;
int m;
int **g;
std::vector<std::vector<int>> graph;
std::vector<int> cl;
std::vector<std::vector<int>> nodesPerCluster;
std::vector<int> v1;
std::vector<int> v2;
std::vector<int> prizes;


// Helper struct for adjacency list to keep track of edge indices
struct EdgeInfo
{
    int to;
    int index; // The index k in edgeVars
};

// ============================================================
// OPTIMIZED CALLBACK: GSEC
// ============================================================
void subset(bool *visited, int *lenneigb, int v)
{
	if (!visited[v]){
		visited[v]=true;
		for (int i=0; i<lenneigb[v]; i++){
			subset(visited, lenneigb, g[v][i]);
		}
	}
}
ILOLAZYCONSTRAINTCALLBACK3(CB_SubCiclo1,
                    IloIntVarArray, x,
					IloArray<IloIntVarArray>, y, 
					int*, n_cut)
{			
	
	IloEnv env  = getEnv();
	
	//int solution[c][c-1];
	bool visited[n];
	int lenneigb[n];
	int cl1, cl2;
    
    std::vector<bool> sol(n, false);
    int num_v=0;
    
    for (int i=0; i<n; i++){
        if (getValue(x[i]) > 0.5){
            num_v++;
            sol[i] = true;
        }
    }
    num_v=c;
  
	for (int k=0; k<n; k++){
		visited[k]=false; //Todos os vértices são marcados como não visitados
		lenneigb[k]=0; //O número de vizinhos de cada vértices inicialmente é zero
	}
	
    
    g = new int*[n]; //Grafo de suporte para a DFS
        for (int k=0; k<n; k++){
            g[k] = new int[n];
    }
    //g = std::vector<int> sol(n, -1);
    
	//Monta a subsolução com num_v vértices
	for (int i=0; i<n-1; i++){
		for (int j=i+1; j<n; j++){
			if (graph[i][j]==-1) continue;
			if (getValue(y[i][j])>0.5){
				//cl1=cl[i];
				//cl2=cl[j];
				g[i][lenneigb[i]]=j;
				g[j][lenneigb[j]]=i;
				lenneigb[i]++;
				lenneigb[j]++;
			}
		}
	}
    
    int start = 0;
    while (!sol[start]) start++;
	
	subset(visited, lenneigb, start);
    
    //cout << "\nVisitados:\n";
    std::vector<bool> visitedcluster(c, false);
	for (int i=0; i<n; i++){
        if (visited[i]){
            visitedcluster[cl[i]] = true;
    //        cout << i << " ";
        }
    }
   //cout << "\n";
    
	for (int k=0; k<n; k++){
		if (!visited[k] && sol[k]){//Testa se algum vértice não foi visitado
			(*n_cut)++; //Aumenta o número de cortes
			
			int nv=0;
			
			for (int i=0; i<n; i++){
				if (visited[i]) nv++;
			}
			
			IloExpr corteSubCiclo(env);
			IloExpr corteSubCiclo2(env);
			//IloExpr corteSet(env);
			//int cont=0;
			for (int e=0; e<m; e++){
				//cl1=cl[v1[e]];
				//cl2=cl[v2[e]];
				if(visited[v1[e]] && visited[v2[e]]){
					corteSubCiclo += y[v1[e]][v2[e]];
					//cont++;
				}
				if(!visited[v1[e]] && !visited[v2[e]]){
					corteSubCiclo2 += y[v1[e]][v2[e]];
				}
				//if(visited[cl1] != !visited[cl2]){
				//	corteSet += y[v1[e]][v2[e]];
				//}
			}
			
            IloExpr sumX1(env);
            IloExpr sumX2(env);
            
            for (int i=0; i<n; i++){
                //cl1 = cl[i];
                if (visited[i]) sumX1 += x[i];
                else sumX2 += x[i];
            }
			
            if (nv>2)
			add(corteSubCiclo <= sumX1-1).end();
			corteSubCiclo.end();
            sumX1.end();
			
            if (n-nv-1 > 2)
			add(corteSubCiclo2 <= sumX2-1).end();
			corteSubCiclo2.end();
            sumX2.end();
			
			//add(corteSet >=1).end();
			return;
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
    try
    {
        PCInstance inst = read_pcinstance();
        preprocess(inst); // Disabled to ensure correctness
        n = inst.N;
        c = inst.C;
        prizes = inst.prizes;
        graph.assign(n, vector<int>(n, -1));
        cl = inst.cluster_by_node;

        nodesPerCluster.assign(c, vector<int>());
        for (int i = 0; i < n; i++)
            nodesPerCluster[cl[i]].push_back(i);

        m = inst.edges.size();
        
        cout << "\nDetalhes da instância:\nN = " << n << "\nM = " << m << "\nc = " << c << "\n"; 
        
        v1.clear();
        v2.clear();

        for (Edge &e : inst.edges)
        {
            if (e.u < e.v){
                v1.push_back(e.u);
                v2.push_back(e.v);
            }else{
                v1.push_back(e.v);
                v2.push_back(e.u);
            }
            graph[e.u][e.v] = e.weight;
            graph[e.v][e.u] = e.weight;
        }

        IloEnv env;
        IloModel model(env, "Generalized MST / PCST");

        IloIntVarArray x(env, n, 0, 1);

        // Create edgeVars first (Linear)
        //IloIntVarArray edgeVars(env, m, 0, 1);

        // Map edgeVars to matrix y for convenience (Optional, but useful for constraints)
        // Note: We use IloNumVarArray in array for y to avoid double allocation issues if needed,
        // but here we just keep 'y' as a logical wrapper or build it manually if constraints need it.
        // For the model constraints, we can mostly use edgeVars.

        // Reconstructing y matrix just for model building consistency if you rely on y[i][j]
        // Warning: This creates aliasing. CPLEX handles it, but use edgeVars for global sums.
        IloArray<IloIntVarArray> y(env, n);
        for (int i = 0; i < n; i++)
            y[i] = IloIntVarArray(env, n, 0, 1);

        //for (int k = 0; k < m; k++)
        //{
        //    int u = v1[k];
        //    int v = v2[k];
        //    y[u][v] = edgeVars[k];
        //    y[v][u] = edgeVars[k];
        //}

        // Objective Function
        IloExpr custo(env);
        for (int k = 0; k < m; k++)
        {
            custo += graph[v1[k]][v2[k]] * y[v1[k]][v2[k]];
        }
        model.add(IloMinimize(env, custo));

        // Global Tree/Forest Constraint: SumEdges = SumNodes - 1
        IloExpr sumEdges(env), sumNodes(env);
        for (int k = 0; k < m; k++)
            sumEdges += y[v1[k]][v2[k]];
        for (int i = 0; i < n; i++)
            sumNodes += x[i];

        model.add(sumEdges == sumNodes - IloNum(1));

        sumEdges.end();
        sumNodes.end();

        // Edge-Node Coherence
        for (int k = 0; k < m; k++)
        {
            int u = v1[k];
            int v = v2[k];
            model.add(y[v1[k]][v2[k]] <= x[u]);
            model.add(y[v1[k]][v2[k]] <= x[v]);
        }

        // Minimum Prize per Cluster
        for (int k = 0; k < c; k++)
        {
            IloExpr expr(env);
            for (int v : nodesPerCluster[k])
                expr += inst.prizes[v] * x[v];
            model.add(expr >= inst.min_prize_by_cluster[k]);
            expr.end();
        }

        //for (int k = 0; k < c; k++)
        //{
            //IloExpr sumcluster(env);
            //for (int i = 0; i < n; i++)
            //{
                //if (cl[i] == k)
                //{
                    //sumcluster += x[i];
                //}
            //}
            //model.add(sumcluster == IloNum(1)); // Restrição do uso de um vértice por cluster
            //sumcluster.end();
        //}
        
        //Restrições de grau :: Restrição opcional
        for (int k=0; k<c; k++){
		    IloExpr grau(env);
            for (int e=0; e<m; e++){
                if (cl[v1[e]] == k || cl[v2[e]] == k){
        			grau += y[v1[e]][v2[e]];
                }
            }
            model.add(grau >= IloNum(1));//Obriga que o vértice selecionado no cluster tenha pelo menos grau 1
            grau.end();
        }
          for (int k=0; k<c; k++){
		    IloExpr grau(env);
            for (int e=0; e<m; e++){
                if (cl[v1[e]] == k || cl[v2[e]] == k){
        			grau += y[v1[e]][v2[e]];
                }
            }
            model.add(grau <= IloNum(5));//Obriga que o vértice selecionado no cluster tenha pelo menos grau 1
            grau.end();
        }
        
        //Restrição opcional, mas acelera a convergência
        for (int i = 0; i < c-1; i++)
        {
            IloExpr sumX1(env);
            for (int u: nodesPerCluster[i])
                sumX1 += x[u];
            for (int j = i+1; j < c; j++)
            {
                if (i == j)
                    continue;
                
                IloExpr sumX2(env);
                for (int u: nodesPerCluster[j])
                    sumX2 += x[u];
                
                IloExpr aux(env);
                for (int e = 0; e < m; e++)
                {
                    if ((cl[v1[e]] == i && cl[v2[e]] == j) || (cl[v1[e]] == j && cl[v2[e]] == i))
                    {
                        aux += y[v1[e]][v2[e]];
                    }
                }
                model.add(aux <= sumX1 + sumX2 -1); //O número de arestas entres dois clusters é no máximo o número de véritices menos 1
                aux.end();
                sumX2.end();
            }
            sumX1.end();
        }

        IloCplex cplex(model);

        // Parameters for Speed
        cplex.setParam(IloCplex::Param::Threads, 0); // Use all cores
        // cplex.setParam(IloCplex::Param::MIP::Strategy::Search, 1); // Traditional Branch & Cut (sometimes better with heavy callbacks)

        // CALLBACK
        int n_cut = 0;
        //cplex.use(CB_SubCiclo(env, x, y, &n_cut));
        cplex.use(CB_SubCiclo1(env, x, y, &n_cut));

        if (cplex.solve())
        {
            double focost = cplex.getObjValue();
            cout << "------------------------------------------------" << endl;
            cout << "Status: " << cplex.getStatus() << endl;
            cout << "Custo da FO: " << focost << endl;
            cout << "Cortes gerados: " << n_cut << endl;
            cout << "------------------------------------------------" << endl;

            cout << "Vertices Selecionados:" << endl;
            int vCount = 0;
            // Use solution retrieval again for printing
            IloNumArray final_x(env);
            cplex.getValues(final_x, x);
            for (int i = 0; i < n; i++)
            {
                if (final_x[i] > 0.9)
                {
                    cout << "["<< i << " c"<< cl[i] << "], ";
                    vCount++;
                }
            }
            cout << endl
                 << "Total Vertices: " << vCount << endl;

            cout << "Arestas Selecionadas:" << endl;
            int eCount = 0;
            //IloNumArray final_y(env);
            //cplex.getValues(final_y, edgeVars);
            for (int k = 0; k < m; k++)
            {
                float yvalue = cplex.getValue(y[v1[k]][v2[k]]);
                if (yvalue > 0.9)
                {
                    cout << "(" << v1[k] << ", " << v2[k] << ") Peso: " << graph[v1[k]][v2[k]] << endl;
                    eCount++;
                }
            }
            cout << "Total Arestas: " << eCount << endl;
        }
        else
        {
            cout << "Inviavel" << endl;
        }

        env.end();
    }
    catch (IloException &e)
    {
        std::cerr << "Erro CPLEX: " << e << "\n";
    }
    return 0;
}
