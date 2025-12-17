#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>

// --- Struct definitions ---
struct Edge {
    int u, v;
    int weight;

    bool operator<(const Edge& other) const {
        if (weight != other.weight) return weight < other.weight;
        if (u != other.u) return u < other.u;
        return v < other.v;
    }
};

struct PCInstance {
    int N; // Number of nodes
    int C; // Number of clusters
    
    std::vector<std::pair<int, int>> coords; 
    std::vector<std::vector<int>> clusters;  
    std::vector<int> cluster_by_node;        
    std::vector<int> centers;                
    std::vector<int> prizes;                 
    std::vector<int> min_prize_by_cluster;
    std::vector<Edge> edges;                 
};

// --- Solver Function ---
void solvePCInstance(const PCInstance& instance) {
    IloEnv env;
    try {
        IloModel model(env);
        
        int n_nodes = instance.N;
        int n_edges = instance.edges.size();

        // 1. Variables
        IloNumVarArray y(env, n_nodes, 0, 1, ILOBOOL); // Node active
        IloNumVarArray x(env, n_edges, 0, 1, ILOBOOL); // Edge active
        
        // Flow variables for connectivity
        IloNumVarArray f_uv(env, n_edges, 0, IloInfinity, ILOFLOAT); 
        IloNumVarArray f_vu(env, n_edges, 0, IloInfinity, ILOFLOAT); 
        IloNumVarArray f_root(env, n_nodes, 0, IloInfinity, ILOFLOAT); // Flow from virtual root
        IloNumVarArray is_root(env, n_nodes, 0, 1, ILOBOOL); // Which node connects to virtual root

        // 2. Objective: Minimize Edge Costs
        IloExpr objExpr(env);
        for (int i = 0; i < n_edges; ++i) {
            objExpr += x[i] * instance.edges[i].weight;
        }
        model.add(IloMinimize(env, objExpr));
        objExpr.end();

        // 3. Constraints

        // A. Cluster Prize Constraints & At least one vertex per cluster
        for (int c = 0; c < instance.C; ++c) {
            IloExpr clusterPrizeSum(env);
            IloExpr clusterPresence(env);
            
            for (int node_id : instance.clusters[c]) {
                clusterPrizeSum += y[node_id] * instance.prizes[node_id];
                clusterPresence += y[node_id];
            }
            
            model.add(clusterPrizeSum >= instance.min_prize_by_cluster[c]);
            model.add(clusterPresence >= 1); 
            
            clusterPrizeSum.end();
            clusterPresence.end();
        }

        // B. Edge Validity (Edge needs both nodes active)
        for (int i = 0; i < n_edges; ++i) {
            int u = instance.edges[i].u;
            int v = instance.edges[i].v;
            model.add(x[i] <= y[u]);
            model.add(x[i] <= y[v]);
        }

        // C. Tree Property: Edges = Nodes - 1 (prevents cycles if connected)
        IloExpr sumEdges(env);
        IloExpr sumNodes(env);
        for(int i=0; i<n_edges; ++i) sumEdges += x[i];
        for(int i=0; i<n_nodes; ++i) sumNodes += y[i];
        model.add(sumEdges == sumNodes - 1); 
        sumEdges.end();
        sumNodes.end();

        // D. Connectivity (Flow formulation)
        // Only one node connects to virtual root
        IloExpr sumIsRoot(env);
        for(int i=0; i<n_nodes; ++i) sumIsRoot += is_root[i];
        model.add(sumIsRoot == 1);
        sumIsRoot.end();

        for(int i=0; i<n_nodes; ++i) {
            model.add(is_root[i] <= y[i]);
        }

        // Flow capacities
        double M = n_nodes; 
        for (int i = 0; i < n_edges; ++i) {
            model.add(f_uv[i] <= M * x[i]);
            model.add(f_vu[i] <= M * x[i]);
        }
        for(int i=0; i<n_nodes; ++i) {
            model.add(f_root[i] <= M * is_root[i]);
        }

        // Flow Conservation
        std::vector<std::vector<int>> adj(n_nodes);
        for(int i=0; i<n_edges; ++i) {
            adj[instance.edges[i].u].push_back(i);
            adj[instance.edges[i].v].push_back(i);
        }

        for (int i = 0; i < n_nodes; ++i) {
            IloExpr flowIn(env);
            IloExpr flowOut(env);

            flowIn += f_root[i];

            for (int e_idx : adj[i]) {
                int u = instance.edges[e_idx].u;
                int v = instance.edges[e_idx].v;

                if (u == i) { // Edge is i -- v
                    flowIn += f_vu[e_idx];
                    flowOut += f_uv[e_idx];
                } else { // Edge is u -- i
                    flowIn += f_uv[e_idx];
                    flowOut += f_vu[e_idx];
                }
            }

            model.add(flowIn - flowOut == y[i]);
            flowIn.end();
            flowOut.end();
        }

        // 4. Solve
        IloCplex cplex(model);
        // cplex.setOut(env.getNullStream()); // Uncomment to silence CPLEX output
        
        if (cplex.solve()) {
            std::cout << "\n-----------------------------------" << std::endl;
            std::cout << "Solution Status: " << cplex.getStatus() << std::endl;
            std::cout << "Objective Cost : " << cplex.getObjValue() << std::endl;
            std::cout << "-----------------------------------" << std::endl;

            std::cout << "Active Nodes:" << std::endl;
            for (int i = 0; i < n_nodes; ++i) {
                if (cplex.getValue(y[i]) > 0.5) {
                    std::cout << "  Node " << i << " (Cluster " << instance.cluster_by_node[i] 
                              << ", Prize " << instance.prizes[i] << ")" << std::endl;
                }
            }
            std::cout << "Active Edges:" << std::endl;
            for (int i = 0; i < n_edges; ++i) {
                if (cplex.getValue(x[i]) > 0.5) {
                    std::cout << "  (" << instance.edges[i].u << ", " << instance.edges[i].v 
                              << ") Cost: " << instance.edges[i].weight << std::endl;
                }
            }
        } else {
            std::cout << "No solution found." << std::endl;
        }

    } catch (IloException& e) {
        std::cerr << "Concert exception caught: " << e << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception caught" << std::endl;
    }
    env.end();
}
PCInstance read_pcinstance() {
    PCInstance inst;

    std::cin >> inst.N;

    inst.coords.resize(inst.N);
    for (int i = 0; i < inst.N; i++) {
        std::cin >> inst.coords[i].first >> inst.coords[i].second;
    }

    std::cin >> inst.C;

    inst.clusters.resize(inst.C);
    inst.min_prize_by_cluster.resize(inst.C);
    inst.cluster_by_node.resize(inst.N);

    // 4. Read Clusters structure
    // Format: [size] [node1] [node2] ...
    for (int c = 0; c < inst.C; c++) {
        int size; 
        std::cin >> size;
        
        inst.clusters[c].reserve(size);
        for (int k = 0; k < size; k++) {
            int node_idx;
            std::cin >> node_idx;
            
            // Adjust 1-based input to 0-based index
            node_idx--; 
            
            inst.clusters[c].push_back(node_idx);
            inst.cluster_by_node[node_idx] = c;
        }
    }

    // 5. Read Centers
    inst.centers.resize(inst.C);
    for (int c = 0; c < inst.C; c++) {
        std::cin >> inst.centers[c];
        inst.centers[c]--; // Adjust to 0-based
    }

    // 6. Validation -999
    int validation;
    std::cin >> validation;
    if (validation != -999) {
        throw std::runtime_error("Inconsistent input: Missing first -999 marker");
    }

    // 7. Read Prizes
    inst.prizes.resize(inst.N);
    for (int i = 0; i < inst.N; i++) {
        std::cin >> inst.prizes[i];
        inst.prizes[i]=1; // Increment prize by 1
    }

    // 8. Validation -999
    std::cin >> validation;
    if (validation != -999) {
        throw std::runtime_error("Inconsistent input: Missing second -999 marker");
    }

    std::vector<int> cluster_sum(inst.C), cluster_min(inst.C, INF);
    // 9. Generate Edges (Euclidean distance)
    // Connect all nodes i, j UNLESS they are in the same cluster
    for (int i = 0; i < inst.N; i++) {
        cluster_sum[inst.cluster_by_node[i]] += inst.prizes[i];
        cluster_min[inst.cluster_by_node[i]] = std::min(cluster_min[inst.cluster_by_node[i]], inst.prizes[i]);
        
        for (int j = i + 1; j < inst.N; j++) {
            // Skip if in the same cluster
            if (inst.cluster_by_node[i] == inst.cluster_by_node[j] && !allowIntra) continue;

            double dx = inst.coords[i].first - inst.coords[j].first;
            double dy = inst.coords[i].second - inst.coords[j].second;
            
            // Round to nearest integer per Reinelt (1991) logic
            int weight = (int)std::round(std::sqrt(dx*dx + dy*dy));

            inst.edges.push_back({i, j, weight});
        }
    }
    std::sort(inst.edges.begin(), inst.edges.end());

    for (int i = 0; i < inst.C; i++) {
        inst.min_prize_by_cluster[i] =1;
    }

    return inst;
}


// --- Main Entry Point ---
int main() {
    // create a dummy instance to test compilation and logic
    PCInstance inst;
    inst.N = 5;
    inst.C = 2;
    
    // Nodes 0, 1 are Cluster 0. Nodes 2, 3, 4 are Cluster 1.
    inst.clusters = {{0, 1}, {2, 3, 4}};
    inst.cluster_by_node = {0, 0, 1, 1, 1};
    inst.prizes = {10, 10, 5, 5, 20}; // Node prizes
    
    // Cluster 0 needs 10 prize (1 node). Cluster 1 needs 20 prize (Node 4 or Nodes 2+3+more).
    inst.min_prize_by_cluster = {10, 20}; 
    
    // Simple line graph 0-1-2-3-4
    inst.edges = {
        {0, 1, 5},
        {1, 2, 10},
        {2, 3, 5},
        {3, 4, 5},
        {0, 2, 100} // expensive shortcut
    };
    PCInstance inst2= read_pcinstance();

    std::cout << "Solving PCInstance..." << std::endl;
    solvePCInstance(inst);

    return 0;
}