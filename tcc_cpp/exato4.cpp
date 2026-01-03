#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <chrono>

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

// ================= UNION-FIND =================
struct UnionFind
{
    std::vector<int> parent, rank;

    UnionFind(int n = 0) { init(n); }

    void init(int n)
    {
        parent.resize(n);
        rank.assign(n, 0);
        for (int i = 0; i < n; i++)
            parent[i] = i;
    }

    int find(int x)
    {
        while (parent[x] != x)
        {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void unite(int a, int b)
    {
        a = find(a);
        b = find(b);
        if (a == b)
            return;
        if (rank[a] < rank[b])
            std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b])
            rank[a]++;
    }
};

// --- Callback buffers (reused)
static std::vector<char> sol;
static UnionFind uf;

// ============================================================
// LAZY CALLBACK – UNION-FIND VERSION
// ============================================================
ILOLAZYCONSTRAINTCALLBACK3(CB_SubCicloUF,
                           IloIntVarArray, x,
                           IloArray<IloIntVarArray>, y,
                           int *, n_cut)
{
    IloEnv env = getEnv();

    std::fill(sol.begin(), sol.end(), 0);
    uf.init(n);

    // Selected vertices
    int first = -1;
    for (int i = 0; i < n; i++)
    {
        if (getValue(x[i]) > 0.5)
        {
            sol[i] = 1;
            if (first == -1)
                first = i;
        }
    }
    if (first == -1)
        return;

    // Union edges
    for (int k = 0; k < m; k++)
    {
        if (getValue(y[v1[k]][v2[k]]) > 0.5)
        {
            uf.unite(v1[k], v2[k]);
        }
    }

    int root = uf.find(first);

    // Check disconnected vertex
    for (int i = 0; i < n; i++)
    {
        if (sol[i] && uf.find(i) != root)
        {
            (*n_cut)++;

            IloExpr cut1(env), cut2(env);
            IloExpr sumX1(env), sumX2(env);

            int nv = 0;
            for (int v = 0; v < n; v++)
            {
                if (sol[v] && uf.find(v) == root)
                {
                    sumX1 += x[v];
                    nv++;
                }
                else if (sol[v])
                {
                    sumX2 += x[v];
                }
            }

            for (int e = 0; e < m; e++)
            {
                int u = v1[e], v = v2[e];
                if (sol[u] && sol[v])
                {
                    if (uf.find(u) == root && uf.find(v) == root)
                        cut1 += y[u][v];
                    if (uf.find(u) != root && uf.find(v) != root)
                        cut2 += y[u][v];
                }
            }

            if (nv > 2)
                add(cut1 <= sumX1 - 1).end();
            if ((int)std::count(sol.begin(), sol.end(), 1) - nv > 2)
                add(cut2 <= sumX2 - 1).end();

            cut1.end();
            cut2.end();
            sumX1.end();
            sumX2.end();

            // double ub = getIncumbentObjValue();
            // double lb = getBestObjValue();
            // double rlx = getObjValue();
            // double nodes = getNremainingNodes();

            // std::cout << "***** LAZY CUT (UF) "
            //           << " relax=" << rlx
            //           << " bounds=" << lb << "<->" << ub
            //           << " nodes=" << nodes
            //           << " #cuts=" << *n_cut << std::endl;
            return;
        }
    }
}
namespace fs = std::filesystem;

// Check if result file already exists
bool resultExists(const std::string &instanceName)
{
    fs::path dir("exato_results");
    fs::path file = dir / (instanceName + ".txt");
    return fs::exists(file);
}

// Save solution
void saveResult(const std::string &instanceName,
                double obj,
                IloCplex &cplex,
                IloIntVarArray &x,
                IloArray<IloIntVarArray> &y)
{
    fs::create_directory("exato_results");

    fs::path file = fs::path("exato_results") / (instanceName + ".txt");
    std::ofstream out(file);

    out << "Objective: " << obj << "\n\n";

    out << "Selected vertices:\n";
    for (int i = 0; i < n; i++)
    {
        if (cplex.getValue(x[i]) > 0.5)
            out << i << "\n";
    }

    out << "\nSelected edges:\n";
    for (int k = 0; k < m; k++)
    {
        int u = v1[k], v = v2[k];
        if (cplex.getValue(y[u][v]) > 0.5)
            out << u << " " << v << "\n";
    }

    out.close();
}

void saveResultCSV(const std::string &instanceName,
                   double obj,
                   double runtime_sec,
                   IloCplex &cplex,
                   IloIntVarArray &x,
                   IloArray<IloIntVarArray> &y)
{
    namespace fs = std::filesystem;

    fs::create_directory("exato_results");
    fs::path file = fs::path("exato_results") / (instanceName + ".csv");

    std::ofstream out(file);

    // CSV header
    out << "type,id_u,id_v,value\n";

    // Objective
    out << "objective,,,";
    out << obj << "\n";

    // Runtime (seconds)
    out << "runtime_seconds,,," << runtime_sec << "\n";

    // Selected vertices
    for (int i = 0; i < n; i++)
    {
        if (cplex.getValue(x[i]) > 0.5)
        {
            out << "vertex," << i << ",,1\n";
        }
    }

    // Selected edges
    for (int k = 0; k < m; k++)
    {
        int u = v1[k], v = v2[k];
        if (cplex.getValue(y[u][v]) > 0.5)
        {
            out << "edge," << u << "," << v << "," << graph[u][v] << "\n";
        }
    }

    out.close();
}

// ============================================================
// MAIN
// ============================================================
int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: ./exato4 <instance_file>\n";
        return 1;
    }

    namespace fs = std::filesystem;

    fs::path input_path(argv[1]);

    if (!fs::exists(input_path) || !fs::is_regular_file(input_path))
    {
        std::cerr << "File does not exist: " << input_path << "\n";
        return 1;
    }

    // ===== instance name from file =====
    std::string instanceName = input_path.stem().string();

    fs::create_directory("exato_results");
    fs::path out_path = fs::path("exato_results") / (instanceName + ".csv");

    if (fs::exists(out_path))
    {
        std::cout << "Result already exists: " << out_path << "\n";
        return 0;
    }

    // ===== redirect file to cin =====
    std::ifstream in(input_path);
    if (!in.is_open())
    {
        std::cerr << "Cannot open file: " << input_path << "\n";
        return 1;
    }

    std::streambuf *cin_backup = std::cin.rdbuf();
    std::cin.rdbuf(in.rdbuf());

    try
    {
        // ===== normal solver flow (UNCHANGED) =====
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

        v1.clear();
        v2.clear();
        for (auto &e : inst.edges)
        {
            int u = std::min(e.u, e.v);
            int v = std::max(e.u, e.v);
            v1.push_back(u);
            v2.push_back(v);
            graph[u][v] = graph[v][u] = e.weight;
        }
        m = v1.size();

        sol.assign(n, 0);
        uf.init(n);

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
        for (int k = 0; k < m; k++)
            sumE += y[v1[k]][v2[k]];
        for (int i = 0; i < n; i++)
            sumV += x[i];
        model.add(sumE == sumV - 1);
        sumE.end();
        sumV.end();

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

        for (int k = 0; k < c; k++)
        {
            IloExpr sumcluster(env);
            for (int i = 0; i < n; i++)
                if (cl[i] == k)
                    sumcluster += x[i];
            model.add(sumcluster >= 1);
            sumcluster.end();
        }

        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::Threads, 0);
        cplex.setParam(IloCplex::Param::TimeLimit, 3600.0);

        int n_cut = 0;
        cplex.use(CB_SubCicloUF(env, x, y, &n_cut));
        auto t_start = std::chrono::steady_clock::now();

        bool solved = cplex.solve();

        auto t_end = std::chrono::steady_clock::now();
        double runtime_sec =
            std::chrono::duration<double>(t_end - t_start).count();

        if (solved)
        {

            std::cout << "OBJ = " << cplex.getObjValue()
                      << " | CUTS = " << n_cut << "\n";

            saveResultCSV(instanceName,
                          cplex.getObjValue(),
                            runtime_sec,
                          cplex,
                          x,
                          y);
        }
        else
        {
            std::cout << "INFEASIBLE\n";
        }

        env.end();
    }
    catch (IloException &e)
    {
        std::cerr << "CPLEX ERROR: " << e << "\n";
    }

    // ===== restore cin =====
    std::cin.rdbuf(cin_backup);
    in.close();
    return 0;
}
