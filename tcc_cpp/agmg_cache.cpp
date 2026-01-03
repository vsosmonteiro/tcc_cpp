#include <random>
#include <iostream>
#include <set>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include "model.hpp"
#include "reader.hpp"
#include "preprocess.hpp"

static const double TIME_LIMIT = 60.0; // seconds
std::chrono::steady_clock::time_point START_TIME;
std::vector<std::vector<Edge>> adj; // Global adjacency list
std::vector<Edge> sorted_edges;
std::vector<Edge> out_edges;
int n=0;
int best_cost_g = INT_MAX;

inline bool time_limit_reached()
{
    auto now = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = now - START_TIME;
    return elapsed.count() >= TIME_LIMIT;
}

struct dsu
{
    std::vector<int> id, sz;

    dsu(int n) : id(n), sz(n, 1) { iota(id.begin(), id.end(), 0); }

    int find(int a) { return a == id[a] ? a : id[a] = find(id[a]); }

    void uni(int a, int b)
    {
        a = find(a), b = find(b);
        if (a == b)
            return;
        if (sz[a] < sz[b])
            std::swap(a, b);
        sz[a] += sz[b], id[b] = a;
    }
    bool unite(int a, int b)
    {
        a = find(a);
        b = find(b);

        if (a == b)
            return false; // cycle detected

        if (sz[a] < sz[b])
            std::swap(a, b);

        id[b] = a;
        sz[a] += sz[b];
        return true;
    }
};


// Helper: Calculate MST cost from a pre-sorted edge list
long long fast_kruskal(const std::vector<Edge>& edges, int n_nodes) {
    if (n_nodes <= 1) return 0;
    
    dsu global_dsu = dsu(n_nodes);
    long long cost = 0;
    int edges_count = 0;
    
    for (const auto& e : edges) {
        if (global_dsu.unite(e.u, e.v)) {
            cost += e.weight;
            edges_count++;
        }
    }
    if (edges_count < n_nodes - 1) return INT_MAX; // disconnected
    return cost; 
}

// Logic to evaluate "ADD NODE" by merging new edges
long long eval_add(int new_node,std::vector<bool> isActive) {
    out_edges.clear();
    out_edges.reserve(sorted_edges.size() + adj[new_node].size());

    // 1. Identify new edges connecting new_node to current solution
    // Since adj is pre-sorted, these are naturally sorted!
    std::vector<Edge> new_connections;
    new_connections.reserve(adj[new_node].size());
    
    for (const auto& e : adj[new_node]) {
        if ((isActive[e.v] || e.v==new_node) && (isActive[e.u] || e.u==new_node)) {
            new_connections.push_back(e);
        }
    }

    // 2. Merge existing sorted edges with new sorted edges (Linear O(N))
    // This replaces std::sort
    std::merge(sorted_edges.begin(), sorted_edges.end(),
               new_connections.begin(), new_connections.end(),
               std::back_inserter(out_edges),
               [](const Edge& a, const Edge& b) { return a.weight < b.weight; });

    return fast_kruskal(out_edges, n);
}

// // Logic to evaluate "DROP NODE" by filtering edges
// long long eval_drop(const SolutionState& current, int drop_node, std::vector<Edge>& out_edges) {
//     out_edges.clear();
//     out_edges.reserve(current.sorted_edges.size());

//     // Filter out edges touching the dropped node (Linear O(N))
//     for (const auto& e : current.sorted_edges) {
//         if (e.u != drop_node && e.v != drop_node) {
//             out_edges.push_back(e);
//         }
//     }

//     return fast_kruskal(out_edges, current.nodes.size() - 1);
// }

std::pair<int, std::vector<Edge>>
mst(PCInstance& inst, const std::vector<int>& subgraph)
{
    dsu uf(inst.N);
    std::vector<Edge> edges;
    std::vector<bool> in_subgraph(inst.N, false);

    for (int u : subgraph)
        in_subgraph[u] = true;

    for (int u : subgraph) {
        for (const Edge& e : adj[u]) {

            // ensure both endpoints are selected
            if (!in_subgraph[e.u] || !in_subgraph[e.v])
                continue;

            // avoid duplicates (only u < v)
            if (e.u != u)
                continue;

            edges.push_back(e);
        }
    }

    std::sort(edges.begin(), edges.end(),
              [](const Edge& a, const Edge& b) {
                  return a.weight < b.weight;
              });

    int cost = 0;
    std::vector<Edge> mst_edges;

    for (const Edge& e : edges) {
        if (uf.unite(e.u, e.v)) {
            mst_edges.push_back(e);
            cost += e.weight;
            if (cost > best_cost_g) {
                return {INT_MAX, mst_edges}; // Early exit if cost exceeds current best
            }
        }
    }

    return {cost, mst_edges};
}


std::vector<int> build_random_subgraph(PCInstance &inst)
{
    std::vector<int> subgraph;

    for (int i = 0; i < inst.C; i++)
    {
        auto nodes = inst.clusters[i];
        std::shuffle(nodes.begin(), nodes.end(), rng);
        int prize = 0;
        int j = 0;
        while (prize < inst.min_prize_by_cluster[i])
        {
            prize += inst.prizes[nodes[j]];
            subgraph.push_back(nodes[j]);
            j++;
        }
    }
    return subgraph;
}

int calc_solution_cost(PCInstance &inst, std::vector<int> &subgraph)
{
    auto [mst_weight, _] = mst(inst, subgraph);
    return mst_weight;
}

bool has_cycle(int num_vertices, const std::vector<Edge> &selected_edges)
{

    dsu uf(num_vertices);

    for (const Edge &e : selected_edges)
    {
        if (!uf.unite(e.u, e.v))
        {
            std::cout << "❌ Cycle detected when adding edge: "
                      << e.u << " - " << e.v
                      << " (weight = " << e.weight << ")\n";
            return true;
        }
    }

    return false;
}

int local_search(PCInstance &inst, std::vector<int> &current_subgraph, int current_cost)
{
    bool improved = true;
    std::vector<int> p(inst.C);
    std::iota(p.begin(), p.end(), 0);

    std::vector<bool> is_selected(inst.N, false);
    for (int u : current_subgraph)
        is_selected[u] = true;

    while (improved)
    {
        improved = false;
        std::shuffle(p.begin(), p.end(), rng);

        for (int c : p)
        {
            std::vector<int> nodes_in, nodes_out;
            long long cluster_prize_sum = 0;

            for (int u : inst.clusters[c])
            {
                if (is_selected[u])
                {
                    nodes_in.push_back(u);
                    cluster_prize_sum += inst.prizes[u];
                }
                else
                {
                    nodes_out.push_back(u);
                }
            }

            // 1. DROP (Remover)
            if (!nodes_in.empty())
            {
                for (int u : nodes_in)
                {
                    if (cluster_prize_sum - inst.prizes[u] >= inst.min_prize_by_cluster[c])
                    {
                        std::vector<int> candidate;
                        candidate.reserve(current_subgraph.size());
                        for (int x : current_subgraph)
                            if (x != u)
                                candidate.push_back(x);

                        long long new_cost = calc_solution_cost(inst, candidate);
                        if (new_cost < current_cost)
                        {
                            current_subgraph = candidate;
                            current_cost = new_cost;
                            is_selected[u] = false;
                            improved = true;
                            break;
                        }
                    }
                }
                if (improved)
                    break;
            }

            // 2. ADD (Adicionar)
            for (int v : nodes_out)
            {
                std::vector<int> candidate = current_subgraph;
                candidate.push_back(v);
                

                long long new_cost = eval_add(v, is_selected);
                if (new_cost < current_cost)
                {
                    current_subgraph = candidate;
                    current_cost = new_cost;
                    is_selected[v] = true;
                    improved = true;
                    break; // <--- EARLY RETURN
                }
            }
            if (improved)
                break;

            // 3. SWAP (Trocar)
            for (int u : nodes_in)
            {
                for (int v : nodes_out)
                {
                    if (cluster_prize_sum - inst.prizes[u] + inst.prizes[v] >= inst.min_prize_by_cluster[c])
                    {
                        std::vector<int> candidate;
                        candidate.reserve(current_subgraph.size());
                        for (int x : current_subgraph)
                            if (x != u)
                                candidate.push_back(x);
                        candidate.push_back(v);

                        long long new_cost = calc_solution_cost(inst, candidate);
                        if (new_cost < current_cost)
                        {
                            current_subgraph = candidate;
                            current_cost = new_cost;
                            is_selected[u] = false;
                            is_selected[v] = true;
                            improved = true;
                            break;
                        }
                    }
                }
                if (improved)
                    break;
            }
            if (improved)
                break;
        }
    }
    return current_cost;
}

void ils(PCInstance &inst)
{
    auto current_subgraph = build_random_subgraph(inst);
    auto best_subgraph = current_subgraph;

    int max_iterations = 50;

    int lim_p = 0;
    for (int i = 0; i < inst.C; i++)
    {
        if (inst.clusters[i].size() > 1)
            lim_p++;
    }

    int perturb = std::min(lim_p, 4);
    int current_cost = calc_solution_cost(inst, current_subgraph);
    int best_cost = current_cost;
    best_cost_g = best_cost;
    int count = 0;
    while (count < max_iterations)
    {
        count += 1;
        current_cost = local_search(inst, current_subgraph, current_cost);

        if (current_cost < best_cost)
        {
            best_cost = current_cost;
            best_cost_g = best_cost;
            best_subgraph = current_subgraph;
            // std::cout << "Nova melhor solucao na iter " << iteration << ": " << best_cost << "\n";
        }

        // 3. Perturbação (Kick) para escapar de mínimos locais
        // Sorteia alguns clusters e re-gera seus vértices aleatoriamente
        std::vector<int> next_subgraph = current_subgraph; // Copia o estado atual (não o best)

        for (int k = 0; k < perturb; k++)
        {
            int c = uniform(0, inst.C - 1); // Escolhe cluster aleatório

            // Remove vértices desse cluster da solução atual
            for (int u : inst.clusters[c])
            {
                next_subgraph.erase(std::remove(next_subgraph.begin(), next_subgraph.end(), u), next_subgraph.end());
            }

            // Adiciona novos vértices aleatórios até bater a meta
            std::vector<int> cluster_nodes = inst.clusters[c];
            std::shuffle(cluster_nodes.begin(), cluster_nodes.end(), rng);

            int p_sum = 0;
            for (int node : cluster_nodes)
            {
                if (p_sum >= inst.min_prize_by_cluster[c])
                    break;
                next_subgraph.push_back(node);
                p_sum += inst.prizes[node];
            }
        }

        // Aceita a perturbação incondicionalmente para explorar
        current_subgraph = next_subgraph;
        current_cost = calc_solution_cost(inst, current_subgraph);
    }
    for (int u : best_subgraph)
    {
        std::cout << "Selected vertice " << u << " \n";
    }
    std::vector<Edge> edges_in_mst = mst(inst, best_subgraph).second;
    int cost_mst = 0;
    for (auto e : edges_in_mst)
        cost_mst += e.weight;
    for (Edge e : edges_in_mst)
    {
        std::cout << "Selected edge " << e.u << " " << e.v << " " << e.weight << "\n";
    }
    std::cout << edges_in_mst.size() << " edges in MST\n";
    std::cout << best_subgraph.size() << " vertices selected\n";
    if (has_cycle(inst.N, edges_in_mst))
    {
        std::cout << "❌ Cycle detected in the MST!\n";
    }
    else
    {
        std::cout << "✅ No cycles in the MST.\n";
    }
    std::vector<int> prizes_mst(inst.C, 0);
    for (int u : best_subgraph)
    {
        int cluster = inst.cluster_by_node[u];
        prizes_mst[cluster] += inst.prizes[u];
    }
    for (int c = 0; c < inst.C; c++)
    {
        if (prizes_mst[c] < inst.min_prize_by_cluster[c])
        {
            std::cout << "❌ Prize constraint violated in cluster " << c
                      << ": required " << inst.min_prize_by_cluster[c]
                      << ", got " << prizes_mst[c] << "\n";
        }
        else
        {
            std::cout << "✅ Cluster " << c
                      << " meets prize requirement: required " << inst.min_prize_by_cluster[c]
                      << ", got " << prizes_mst[c] << "\n";
        }
    }

    std::cout << "Custo da MST na melhor solução: " << cost_mst << "\n";
    std::cout << "Melhor custo encontrado: " << best_cost << "\n";
}
void create_adjacency_list(PCInstance &inst)
{
    adj.resize(inst.N);
    for (const auto &edge : inst.edges)
    {
        adj[edge.u].push_back(edge);
        adj[edge.v].push_back(edge); // Aresta reversa
    }
}
void sort_adjacency_list(int N) {
    for (int i = 0; i < N; i++) {
        std::sort(adj[i].begin(), adj[i].end(), [](const Edge& a, const Edge& b) {
            return a.weight < b.weight;
        });
    }
}
int main()
{
    auto inst = read_pcinstance();
    START_TIME = std::chrono::steady_clock::now();
    std::cout << inst.edges.size() << "\n";
    preprocess(inst);
    create_adjacency_list(inst);
    sort_adjacency_list(inst.N);
    out_edges.reserve(inst.edges.size());
    n= inst.N;
    ils(inst);
    auto now = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = now - START_TIME;
    std::cout << "Finished execution within time of " << elapsed.count() << " seconds.\n";
    return 0;
}