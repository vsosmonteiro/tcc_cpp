#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include <queue>
#include "model.hpp"
#include "reader.hpp"
#include "preprocess.hpp"

using namespace std;

/* =========================
   Lazy constraint: Connectivity
   ========================= */
ILOLAZYCONSTRAINTCALLBACK3(
    CB_Connectivity,
    IloArray<IloIntVarArray>, y,
    IloIntVarArray, x,
    const PCInstance*, inst)
{
    IloEnv env = getEnv();
    int n = inst->N;

    vector<int> active;
    for (int i = 0; i < n; i++)
        if (getValue(x[i]) > 0.5)
            active.push_back(i);

    if (active.size() <= 1) return;

    vector<bool> visited(n, false);
    queue<int> q;

    q.push(active[0]);
    visited[active[0]] = true;

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (const Edge& e : inst->edges) {
            int v = -1;
            if (e.u == u) v = e.v;
            else if (e.v == u) v = e.u;

            if (v != -1 &&
                getValue(y[e.u][e.v]) > 0.5 &&
                !visited[v]) {
                visited[v] = true;
                q.push(v);
            }
        }
    }

    for (int v : active) {
        if (!visited[v]) {
            IloExpr cut(env);
            for (const Edge& e : inst->edges) {
                if (visited[e.u] != visited[e.v])
                    cut += y[e.u][e.v];
            }
            add(cut >= 1).end();
            cut.end();
            return;
        }
    }
}

int main() {
    IloEnv env;
    try {
        /* =========================
           Instância de exemplo
           ========================= */
        PCInstance inst=read_pcinstance();
        preprocess(inst);
        /* =========================
           Modelo
           ========================= */
        IloModel model(env);

        // Decision variables
        IloIntVarArray x(env, inst.N, 0, 1);          // z_i
        IloArray<IloIntVarArray> y(env, inst.N);      // u_ij
        for (int i = 0; i < inst.N; i++)
            y[i] = IloIntVarArray(env, inst.N, 0, 1);

        /* Objective function */
        IloExpr obj(env);
        for (const Edge& e : inst.edges)
            obj += e.weight * y[e.u][e.v];
        model.add(IloMinimize(env, obj));
        obj.end();

        /* (1) At least one node per cluster */
        // for (int k = 0; k < inst.C; k++) {
        //     IloExpr expr(env);
        //     for (int v : inst.clusters[k])
        //         expr += x[v];
        //     model.add(expr >= 1);
        //     expr.end();
        // }

        /* (2) Minimum prize per cluster */
        for (int k = 0; k < inst.C; k++) {
            IloExpr expr(env);
            for (int v : inst.clusters[k])
                expr += inst.prizes[v] * x[v];
            model.add(expr >= inst.min_prize_by_cluster[k]);
            expr.end();
        }

        /* (3) Tree cardinality constraint */
        IloExpr sumEdges(env), sumNodes(env);
        for (const Edge& e : inst.edges)
            sumEdges += y[e.u][e.v];
        for (int i = 0; i < inst.N; i++)
            sumNodes += x[i];
        model.add(sumEdges == sumNodes - 1);
        sumEdges.end();
        sumNodes.end();

        /* (4) Edge–node linking */
        for (const Edge& e : inst.edges) {
            model.add(y[e.u][e.v] <= x[e.u]);
            model.add(y[e.u][e.v] <= x[e.v]);
        }

        /* =========================
           Solve
           ========================= */
        IloCplex cplex(model);
        cplex.use(CB_Connectivity(env, y, x, &inst));

        cplex.solve();

        /* =========================
           Output
           ========================= */
        cout << "FO = " << cplex.getObjValue() << endl;

        cout << "Vertices selecionados: ";
        for (int i = 0; i < inst.N; i++)
            if (cplex.getValue(x[i]) > 0.5)
                cout << i << " ";
        cout << endl;

        cout << "Arestas selecionadas: ";
        for (const Edge& e : inst.edges)
            if (cplex.getValue(y[e.u][e.v]) > 0.5)
                cout << "(" << e.u << "," << e.v << ") ";
        cout << endl;
    }
    catch (IloException& e) {
        cerr << e << endl;
    }
    env.end();
    return 0;
}
