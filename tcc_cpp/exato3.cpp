#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stack>
#include "model.hpp"
#include "reader.hpp"
#include "preprocess.hpp"

ILOSTLBEGIN

// ================= GLOBALS =================
int n, c, m;
std::vector<std::vector<int>> graph;
std::vector<int> cl;
std::vector<int> v1, v2;
std::vector<int> prizes;
std::vector<std::vector<int>> nodesPerCluster;

// --- Callback buffers (REUSED, NO ALLOC INSIDE CALLBACK)
static std::vector<char> sol;
static std::vector<char> visited;
static std::vector<int> lenneigb;
static std::vector<std::vector<int>> g;

// ============================================================
// FAST DFS (ITERATIVE)
// ============================================================
inline void dfs(int start)
{
    std::stack<int> st;
    st.push(start);
    visited[start] = 1;

    while (!st.empty())
    {
        int u = st.top();
        st.pop();
        for (int i = 0; i < lenneigb[u]; i++)
        {
            int v = g[u][i];
            if (!visited[v])
            {
                visited[v] = 1;
                st.push(v);
            }
        }
    }
}

// ============================================================
// LAZY CALLBACK – SUBTOUR ELIMINATION
// ============================================================
ILOLAZYCONSTRAINTCALLBACK3(CB_SubCicloFast,
    IloIntVarArray, x,
    IloArray<IloIntVarArray>, y,
    int*, n_cut)
{
    IloEnv env = getEnv();

    // Reset buffers
    std::fill(sol.begin(), sol.end(), 0);
    std::fill(visited.begin(), visited.end(), 0);
    std::fill(lenneigb.begin(), lenneigb.end(), 0);

    int num_v = 0;
    for (int i = 0; i < n; i++)
    {
        if (getValue(x[i]) > 0.5)
        {
            sol[i] = 1;
            num_v++;
        }
    }

    // Build adjacency
    for (int k = 0; k < m; k++)
    {
        if (getValue(y[v1[k]][v2[k]]) > 0.5)
        {
            int u = v1[k], v = v2[k];
            g[u][lenneigb[u]++] = v;
            g[v][lenneigb[v]++] = u;
        }
    }

    // Find start
    int start = -1;
    for (int i = 0; i < n; i++)
        if (sol[i]) { start = i; break; }

    if (start == -1) return;

    dfs(start);

    // Check disconnected
    for (int k = 0; k < n; k++)
    {
        if (sol[k] && !visited[k])
        {
            (*n_cut)++;

            IloExpr cut1(env), cut2(env);
            IloExpr sumX1(env), sumX2(env);

            int nv = 0;
            for (int i = 0; i < n; i++)
            {
                if (visited[i])
                {
                    nv++;
                    sumX1 += x[i];
                }
                else
                    sumX2 += x[i];
            }

            for (int e = 0; e < m; e++)
            {
                int u = v1[e], v = v2[e];
                if (visited[u] && visited[v])
                    cut1 += y[u][v];
                if (!visited[u] && !visited[v])
                    cut2 += y[u][v];
            }

            if (nv > 2)
                add(cut1 <= sumX1 - 1).end();
            if (n - nv > 2)
                add(cut2 <= sumX2 - 1).end();

            cut1.end(); cut2.end();
            sumX1.end(); sumX2.end();

            double ub = getIncumbentObjValue();
            double lb = getBestObjValue();
            double rlx = getObjValue();
            double nodes = getNremainingNodes();

            // std::cout << "***** LAZY CUT "
            //           << " relax=" << rlx
            //           << " bounds=" << lb << "<->" << ub
            //           << " nodes=" << nodes
            //           << " #cuts=" << *n_cut << std::endl;
            return;
        }
    }
}

// ============================================================
// MAIN
// ============================================================
int main()
{
    try
    {
        PCInstance inst = read_pcinstance();
        preprocess(inst);

        n = inst.N;
        c = inst.C;
        prizes = inst.prizes;
        cl = inst.cluster_by_node;

        graph.assign(n, std::vector<int>(n, -1));
        nodesPerCluster.assign(c, {});

        for (int i = 0; i < n; i++)
            nodesPerCluster[cl[i]].push_back(i);

        for (auto &e : inst.edges)
        {
            int u = std::min(e.u, e.v);
            int v = std::max(e.u, e.v);
            v1.push_back(u);
            v2.push_back(v);
            graph[u][v] = graph[v][u] = e.weight;
        }

        m = v1.size();

        // Allocate callback buffers ONCE
        sol.resize(n);
        visited.resize(n);
        lenneigb.resize(n);
        g.assign(n, std::vector<int>(n));

        IloEnv env;
        IloModel model(env);

        IloIntVarArray x(env, n, 0, 1);
        IloArray<IloIntVarArray> y(env, n);
        for (int i = 0; i < n; i++)
            y[i] = IloIntVarArray(env, n, 0, 1);

        IloExpr obj(env);
        for (int k = 0; k < m; k++)
            obj += graph[v1[k]][v2[k]] * y[v1[k]][v2[k]];
        model.add(IloMinimize(env, obj));
        obj.end();

        IloExpr sumE(env), sumV(env);
        for (int k = 0; k < m; k++) sumE += y[v1[k]][v2[k]];
        for (int i = 0; i < n; i++) sumV += x[i];
        model.add(sumE == sumV - 1);
        sumE.end(); sumV.end();

        for (int k = 0; k < m; k++)
        {
            model.add(y[v1[k]][v2[k]] <= x[v1[k]]);
            model.add(y[v1[k]][v2[k]] <= x[v2[k]]);
        }

        for (int k = 0; k < c; k++)
        {
            IloExpr expr(env);
            for (int u : nodesPerCluster[k])
                expr += prizes[u] * x[u];
            model.add(expr >= inst.min_prize_by_cluster[k]);
            expr.end();
        }

        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::Threads, 0);
        cplex.setParam(IloCplex::Param::MIP::Strategy::Search, 1);

        int n_cut = 0;
        cplex.use(CB_SubCicloFast(env, x, y, &n_cut));

        if (cplex.solve())
        {
            std::cout << "\nSTATUS: " << cplex.getStatus() << std::endl;
            std::cout << "OBJ: " << cplex.getObjValue() << std::endl;
            std::cout << "CUTS: " << n_cut << std::endl;
        }
        else
            std::cout << "INVIAVEL\n";

        env.end();
    }
    catch (IloException &e)
    {
        std::cerr << "CPLEX ERROR: " << e << std::endl;
    }
}
