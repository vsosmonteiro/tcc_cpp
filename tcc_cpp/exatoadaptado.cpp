#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include "model.hpp"
#include "reader.hpp"
#include "preprocess.hpp"

ILOSTLBEGIN
int n;
vector<vector<int>> graph;
int c;
vector<int> cl;        // cluster de cada vértice
int m;                 // número de arestas
vector<vector<int>> g; // lista de adjacência dos clusters
vector<int> v1;
vector<int> v2;

void subset(bool *visited, int *lenneigb, int v)
{
    if (!visited[v])
    {
        visited[v] = true;
        for (int i = 0; i < lenneigb[v]; i++)
        {
            subset(visited, lenneigb, g[v][i]);
        }
    }
}


ILOLAZYCONSTRAINTCALLBACK2(CB_SubCiclo,
                           IloArray<IloIntVarArray>, y,
                           int *, n_cut)
{

    IloEnv env = getEnv();

    bool visited[c];
    int lenneigb[c];
    int cl1, cl2;

    for (int k = 0; k < c; k++)
    {
        visited[k] = false; // Todos os vértices são marcados como não visitados
        lenneigb[k] = 0;    // O número de vizinhos de cada vértices inicialmente é zero
        g[k].clear();
    }

    // Monta a subsolução com c vértices
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (graph[i][j] == -1)
                continue;
            if (abs(getValue(y[i][j]) - 1) < 0.1)
            {
                cl1 = cl[i];
                cl2 = cl[j];
                g[cl1].push_back(cl2);
                g[cl2].push_back(cl1);
                lenneigb[cl1]++;
                lenneigb[cl2]++;
            }
        }
    }

    subset(visited, lenneigb, 0);

    for (int k = 0; k < c; k++)
    {
        if (!visited[k])
        {               // Testa se algum vértice não foi visitado
            (*n_cut)++; // Aumenta o número de cortes

            int nv = 0;

            for (int i = 0; i < c; i++)
            {
                if (visited[i])
                    nv++;
            }

            IloExpr corteSubCiclo(env);
            IloExpr corteSubCiclo2(env);
            int cont = 0;
            for (int e = 0; e < m; e++)
            {
                cl1 = cl[v1[e]];
                cl2 = cl[v2[e]];
                if (visited[cl1] && visited[cl2])
                {
                    corteSubCiclo += y[v1[e]][v2[e]];
                    cont++;
                }
                if (!visited[cl1] && !visited[cl2])
                {
                    corteSubCiclo2 += y[v1[e]][v2[e]];
                }
            }
            /// callback correto

            add(corteSubCiclo <= IloNum(nv - 1)).end();
            corteSubCiclo.end();
            add(corteSubCiclo2 <= IloNum(c - nv - 1)).end();
            corteSubCiclo2.end();
            break;
        }
    }

    // double ub     = getIncumbentObjValue();  // retorna a melhor solucao inteira (limite primal)
    // double lb     = getBestObjValue();       // retorna o melhor limite dual
    // double rlx    = getObjValue();           // quando chamada dentro do callback, retorna o valor da relaxacao do noh
    // double nNodes = getNremainingNodes();    // retorna o numero restante de nos a serem analisados
    // cout<<"***** LAZY CONSTRAINT: "<<"relax="<<rlx<<"\t bounds="<<lb<<"<->"<<ub<<"\t n_rest="<<nNodes<< "\t#lazy_con="<<*n_cut<<endl;
}

int main()
{

    try
    {
        PCInstance inst = read_pcinstance();
        // preprocess(inst);
        n = inst.N;
        c = inst.C;
        graph.assign(n, vector<int>(n, -1));
        g.assign(c, vector<int>());
        cl = inst.cluster_by_node;
        m = inst.edges.size();
        for (Edge &e : inst.edges)
        {
            v1.push_back(e.u);
            v2.push_back(e.v);
            graph[e.u][e.v] = e.weight;
            graph[e.v][e.u] = e.weight;
        }

        IloEnv env;
        IloModel model(env, "Problema da Árvore Geradora Mínima Generalizada");

        IloIntVarArray x(env, n, 0, 1);
        IloArray<IloIntVarArray> y(env, n);

        for (int i = 0; i < n; i++)
        {
            IloIntVarArray vetor(env, n, 0, 1);
            y[i] = vetor;
        }

        // Função Objetivo
        IloExpr custo(env);
        for (int i = 0; i < n - 1; i++)
        {
            for (int j = i; j < n; j++)
            {
                if (graph[i][j] > -1)
                {
                    custo += graph[i][j] * y[i][j];
                }
            }
        }
        model.add(IloMinimize(env, custo));
        // restricao somatorio x[i]
        IloExpr sumx(env);
        for (int i = 0; i < n; i++)
        {
            sumx += x[i];
        }
        model.add(sumx == IloNum(c));
        // restricao somatorio y[i][j]
        IloExpr sumy(env);
        for (int i = 0; i < n - 1; i++)
        {
            for (int j = i; j < n; j++)
            {
                if (graph[i][j] > -1)
                {
                    sumy += y[i][j];
                }
            }
        }
        model.add(sumy == IloNum(c - 1));

        for (int i = 0; i < n - 1; i++)
        {
            for (int j = i; j < n; j++)
            {
                if (graph[i][j] > -1)
                {
                    model.add(y[i][j] <= x[i]); // Restrições do uso de um arco apenas se os vértices forem usados
                    model.add(y[i][j] <= x[j]);
                }
            }
        }

        for (int k = 0; k < c; k++)
        {
            IloExpr sumcluster(env);
            for (int i = 0; i < n; i++)
            {
                if (cl[i] == k)
                {
                    sumcluster += x[i];
                }
            }
            model.add(sumcluster == IloNum(1)); // Restrição do uso de um vértice por cluster
        }

        for (int i = 0; i < c; i++)
        {
            for (int j = 0; j < c; j++)
            {
                if (i == j)
                    continue;
                IloExpr aux(env);
                for (int e = 0; e < m; e++)
                {
                    if ((cl[v1[e]] == i && cl[v2[e]] == j) || (cl[v1[e]] == j && cl[v2[e]] == i))
                    {
                        aux += y[v1[e]][v2[e]];
                    }
                }
                model.add(aux <= IloNum(1)); // No máximo uma aresta entre dois clusters
            }
        }

        IloCplex cplex(model); // Alocação do solver

        // Criar a estrutura que será parâmetro do USERCALLBACK
        int n_cut = 0;
        cplex.use(CB_SubCiclo(env, y, &n_cut));          // CALLBACK que elimina os subciclos

        cplex.solve(); // Resolução
        double focost = cplex.getObjValue();

        cout << "Custo da FO: " << focost << endl
             << "Vértices: ";
        env.end();
    }
    catch (IloException &e)
    {
        std::cerr << e << "\n";
    }
    return 0;
}
