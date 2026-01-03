#include <random>
#include <iostream>
#include <set>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include "model.hpp"
#include "reader.hpp"
#include "preprocess.hpp"

static const double TIME_LIMIT = 60.0; 
std::chrono::steady_clock::time_point START_TIME;
std::vector<std::vector<Edge>> adj; 

// --- STATE MANAGEMENT ---
std::vector<Edge> sorted_edges; // MUST stay in sync with current_subgraph
std::vector<Edge> out_edges;    // Temporary buffer
int best_cost_g = INT_MAX;

inline bool time_limit_reached() {
    auto now = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = now - START_TIME;
    return elapsed.count() >= TIME_LIMIT;
}

struct dsu {
    std::vector<int> id, sz;
    dsu(int n) : id(n), sz(n, 1) { std::iota(id.begin(), id.end(), 0); }
    int find(int a) { return a == id[a] ? a : id[a] = find(id[a]); }
    bool unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return false;
        if (sz[a] < sz[b]) std::swap(a, b);
        sz[a] += sz[b]; id[b] = a;
        return true;
    }
};

// --- HELPER TO KEEP STATE SYNCED ---
// Call this whenever you change the subgraph via "Slow" methods (Drop/Swap/Init)
void rebuild_sorted_edges(const std::vector<int>& subgraph, int N) {
    sorted_edges.clear();
    std::vector<bool> active(N, false);
    for(int u : subgraph) active[u] = true;

    for (int u : subgraph) {
        for (const auto& e : adj[u]) {
            // Only add edge if u < v (avoid duplicates) and both are active
            if (active[e.v] && u < e.v) {
                sorted_edges.push_back(e);
            }
        }
    }
    // We must sort because we just rebuilt from scratch
    std::sort(sorted_edges.begin(), sorted_edges.end(), [](const Edge& a, const Edge& b){
        return a.weight < b.weight;
    });
}

// Fixed Kruskal: Takes 'active_nodes_count', NOT total 'inst.N'
long long fast_kruskal(const std::vector<Edge>& edges, int active_nodes_count, int total_nodes_N) {
    if (active_nodes_count <= 1) return 0;
    
    dsu global_dsu(total_nodes_N);
    long long cost = 0;
    int edges_count = 0;
    
    for (const auto& e : edges) {
        if (global_dsu.unite(e.u, e.v)) {
            cost += e.weight;
            edges_count++;
            if (cost > best_cost_g) {
                return INT_MAX; // Early exit if cost exceeds current best
            }
        }
    }
    
    // FIX 1: Check against active subgraph size, not total graph
    if (edges_count < active_nodes_count - 1) return INT_MAX; 
    return cost; 
}

// Eval Add
long long eval_add(int new_node, const std::vector<bool>& isActive, int current_subgraph_size, int total_N) {
    out_edges.clear();
    // Heuristic reserve
    out_edges.reserve(sorted_edges.size() + adj[new_node].size());

    static std::vector<Edge> new_connections; // static to reduce allocs
    new_connections.clear();
    
    for (const auto& e : adj[new_node]) {
        int neighbor = (e.u == new_node) ? e.v : e.u;
        if (isActive[neighbor]) {
            new_connections.push_back(e);
        }
    }

    // Merge: sorted_edges (Global State) + new_connections
    std::merge(sorted_edges.begin(), sorted_edges.end(),
               new_connections.begin(), new_connections.end(),
               std::back_inserter(out_edges),
               [](const Edge& a, const Edge& b) { return a.weight < b.weight; });

    // FIX 2: Pass correct size (current + 1)
    return fast_kruskal(out_edges, current_subgraph_size + 1, total_N);
}

// "Slow" MST for verification/initialization
std::pair<int, std::vector<Edge>> mst(PCInstance& inst, const std::vector<int>& subgraph) {
    dsu uf(inst.N);
    std::vector<Edge> edges;
    std::vector<bool> in_subgraph(inst.N, false);
    for (int u : subgraph) in_subgraph[u] = true;

    for (int u : subgraph) {
        for (const Edge& e : adj[u]) {
            if (in_subgraph[e.v] && u < e.v) edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b){ return a.weight < b.weight; });

    int cost = 0;
    std::vector<Edge> mst_edges;
    for (const Edge& e : edges) {
        if (uf.unite(e.u, e.v)) {
            mst_edges.push_back(e);
            cost += e.weight;
        }
    }
    return {cost, mst_edges};
}

// --- STANDARD HELPER FUNCTIONS ---
std::vector<int> build_random_subgraph(PCInstance &inst) {
    std::vector<int> subgraph;
    std::vector<bool> added(inst.N, false);
    // Simple random pick (likely disconnected, corrected by ILS usually)
    for (int i = 0; i < inst.C; i++) {
        auto nodes = inst.clusters[i];
        std::shuffle(nodes.begin(), nodes.end(), rng);
        int prize = 0; int j = 0;
        while (prize < inst.min_prize_by_cluster[i] && j < nodes.size()) {
            if(!added[nodes[j]]) {
                prize += inst.prizes[nodes[j]];
                subgraph.push_back(nodes[j]);
                added[nodes[j]] = true;
            }
            j++;
        }
    }
    return subgraph;
}

int calc_solution_cost(PCInstance &inst, std::vector<int> &subgraph) {
    auto [mst_weight, _] = mst(inst, subgraph);
    return mst_weight;
}

int local_search(PCInstance &inst, std::vector<int> &current_subgraph, int current_cost) {
    bool improved = true;
    std::vector<int> p(inst.C);
    std::iota(p.begin(), p.end(), 0);

    std::vector<bool> is_selected(inst.N, false);
    for (int u : current_subgraph) is_selected[u] = true;

    while (improved) {
        improved = false;
        std::shuffle(p.begin(), p.end(), rng);

        for (int c : p) {
            std::vector<int> nodes_in, nodes_out;
            long long cluster_prize_sum = 0;
            for (int u : inst.clusters[c]) {
                if (is_selected[u]) { nodes_in.push_back(u); cluster_prize_sum += inst.prizes[u]; }
                else { nodes_out.push_back(u); }
            }

            // 1. DROP (Use slow calc, rebuild state)
            if (!nodes_in.empty()) {
                for (int u : nodes_in) {
                    if (cluster_prize_sum - inst.prizes[u] >= inst.min_prize_by_cluster[c]) {
                        std::vector<int> candidate;
                        candidate.reserve(current_subgraph.size());
                        for (int x : current_subgraph) if (x != u) candidate.push_back(x);

                        // Slow check (Safe)
                        long long new_cost = calc_solution_cost(inst, candidate);
                        if (new_cost < current_cost) {
                            current_subgraph = candidate;
                            current_cost = new_cost;
                            is_selected[u] = false;
                            
                            // SYNC FIX: Rebuild sorted edges because we did a "Slow" Drop
                            rebuild_sorted_edges(current_subgraph, inst.N);
                            
                            improved = true;
                            break;
                        }
                    }
                }
                if (improved) break;
            }

            // 2. ADD (Use Fast Eval)
            for (int v : nodes_out) {
                // Pass subgraph size, not inst.N
                long long new_cost = eval_add(v, is_selected, current_subgraph.size(), inst.N);
                
                if (new_cost < current_cost) {
                    current_subgraph.push_back(v);
                    current_cost = new_cost;
                    is_selected[v] = true;
                    
                    // SYNC FIX: Update sorted_edges with the result we just calculated
                    sorted_edges = out_edges; 
                    
                    improved = true;
                    break; 
                }
            }
            if (improved) break;

            // 3. SWAP (Slow calc)
            for (int u : nodes_in) {
                for (int v : nodes_out) {
                    if (cluster_prize_sum - inst.prizes[u] + inst.prizes[v] >= inst.min_prize_by_cluster[c]) {
                        std::vector<int> candidate;
                        candidate.reserve(current_subgraph.size());
                        for (int x : current_subgraph) if (x != u) candidate.push_back(x);
                        candidate.push_back(v);

                        long long new_cost = calc_solution_cost(inst, candidate);
                        if (new_cost < current_cost) {
                            current_subgraph = candidate;
                            current_cost = new_cost;
                            is_selected[u] = false; is_selected[v] = true;
                            
                            // SYNC FIX: Rebuild
                            rebuild_sorted_edges(current_subgraph, inst.N);
                            
                            improved = true;
                            break;
                        }
                    }
                }
                if (improved) break;
            }
            if (improved) break;
        }
    }
    return current_cost;
}

void ils(PCInstance &inst) {
    auto current_subgraph = build_random_subgraph(inst);
    
    // FIX: Initialize sorted_edges before starting search
    rebuild_sorted_edges(current_subgraph, inst.N);
    
    int current_cost = calc_solution_cost(inst, current_subgraph);
    
    auto best_subgraph = current_subgraph;
    int best_cost = current_cost;
    best_cost_g = best_cost;

    int max_iterations = 1000;
    int lim_p = 0;
    for (int i = 0; i < inst.C; i++) if (inst.clusters[i].size() > 1) lim_p++;
    int perturb = std::min(lim_p, 4);
    
    int count = 0;
    while (count < max_iterations && !time_limit_reached()) {
        count++;
        
        current_cost = local_search(inst, current_subgraph, current_cost);

        if (current_cost < best_cost) {
            best_cost = current_cost;
            best_cost_g = best_cost;
            best_subgraph = current_subgraph;
            std::cout << "New Best: " << best_cost << "\n";
        }

        // Perturbation
        std::vector<int> next_subgraph = current_subgraph;
        for (int k = 0; k < perturb; k++) {
            int c = uniform(0, inst.C - 1);
            // Remove cluster
            for (int u : inst.clusters[c]) {
                next_subgraph.erase(std::remove(next_subgraph.begin(), next_subgraph.end(), u), next_subgraph.end());
            }
            // Add random
            std::vector<int> cluster_nodes = inst.clusters[c];
            std::shuffle(cluster_nodes.begin(), cluster_nodes.end(), rng);
            int p_sum = 0;
            for (int node : cluster_nodes) {
                if (p_sum >= inst.min_prize_by_cluster[c]) break;
                next_subgraph.push_back(node);
                p_sum += inst.prizes[node];
            }
        }

        current_subgraph = next_subgraph;
        // FIX: Rebuild state after perturbation chaos
        rebuild_sorted_edges(current_subgraph, inst.N);
        current_cost = calc_solution_cost(inst, current_subgraph);
    }
    
    // Final Report
    std::cout << "--- Final Result ---\n";
    std::cout << "Best Cost: " << best_cost << "\n";
}

void create_adjacency_list(PCInstance &inst) {
    adj.resize(inst.N);
    for (const auto &edge : inst.edges) {
        adj[edge.u].push_back(edge);
        adj[edge.v].push_back(edge);
    }
}

void sort_adjacency_list(int N) {
    for (int i = 0; i < N; i++) {
        std::sort(adj[i].begin(), adj[i].end(), [](const Edge& a, const Edge& b) {
            return a.weight < b.weight;
        });
    }
}

int main() {
    auto inst = read_pcinstance();
    START_TIME = std::chrono::steady_clock::now();
    
    preprocess(inst);
    create_adjacency_list(inst);
    sort_adjacency_list(inst.N);
    out_edges.reserve(inst.edges.size());
    
    ils(inst);
    std::cout << "Time" << std::chrono::duration<double>(std::chrono::steady_clock::now() - START_TIME).count() << " seconds.\n";
    
    return 0;
}