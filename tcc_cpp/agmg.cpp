#include <random>
#include <iostream>
#include <set>
#include <iomanip>

const int SEED = 1;
const int INF = 0x3f3f3f3f;
std::mt19937_64 rng(SEED);
bool allowIntra = false;
int uniform(int l, int r){
    std::uniform_int_distribution<int> uid(l, r);
	return uid(rng);
}

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
    
    std::vector<std::pair<int, int>> coords; // Coordinates (x, y)
    std::vector<std::vector<int>> clusters;  // clusters[c] = list of nodes in cluster c
    std::vector<int> cluster_by_node;        // Map: node_id -> cluster_id
    std::vector<int> centers;                // Center node of each cluster
    std::vector<int> prizes;                 // Prize for each node
    std::vector<int> min_prize_by_cluster;
    std::vector<Edge> edges;                 // Adjacency list/Edge list
};

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
        inst.min_prize_by_cluster[i] = uniform(cluster_min[i], cluster_sum[i]);
    }

    return inst;
}

struct dsu {
	std::vector<int> id, sz;

	dsu(int n) : id(n), sz(n, 1) { iota(id.begin(), id.end(), 0); }

	int find(int a) { return a == id[a] ? a : id[a] = find(id[a]); }

	void uni(int a, int b) {
		a = find(a), b = find(b);
		if (a == b) return;
		if (sz[a] < sz[b]) std::swap(a, b);
		sz[a] += sz[b], id[b] = a;
	}
};

std::pair<int, std::vector<Edge>> mst(PCInstance& inst, std::vector<int> subgraph) {
    dsu union_find(inst.N);
	
    std::set<int> vertices;
    for (auto x : subgraph) vertices.insert(x);

	int cost = 0;
	std::vector<Edge> mst;
	for (auto edge : inst.edges) {
        if (!vertices.count(edge.u) || !vertices.count(edge.v)) continue;
        
        if (union_find.find(edge.u) != union_find.find(edge.v)) {
            mst.push_back(edge);
            cost += edge.weight;
            union_find.uni(edge.u, edge.v);
        }
    }
	return {cost, mst};
}

std::vector<int> build_random_subgraph(PCInstance& inst) {
    std::vector<int> subgraph;

    for (int i = 0; i < inst.C; i++) {
        auto nodes = inst.clusters[i];
        std::shuffle(nodes.begin(), nodes.end(), rng);
        int prize = 0;
        int j = 0;
        while (prize < inst.min_prize_by_cluster[i]) {
            prize += inst.prizes[nodes[j]];
            subgraph.push_back(nodes[j]);
            j++;
        }
    }
    return subgraph;
} 

int calc_solution_cost(PCInstance& inst, std::vector<int>& subgraph) {
    auto [mst_weight, _] = mst(inst, subgraph);
    int total_prize = 0;
    for (int node: subgraph) {
        total_prize += inst.prizes[node];
    }
    return mst_weight - total_prize;
}

int local_search(PCInstance& inst, std::vector<int>& current_subgraph, int current_cost) {
    bool improved = true;
    std::vector<int> p(inst.C);
    std::iota(p.begin(), p.end(), 0); 

    std::vector<bool> is_selected(inst.N, false);
    for(int u : current_subgraph) is_selected[u] = true;

    while (improved) {
        improved = false;
        std::shuffle(p.begin(), p.end(), rng);

        for (int c : p) {
            std::vector<int> nodes_in, nodes_out;
            long long cluster_prize_sum = 0;

            for (int u : inst.clusters[c]) {
                if (is_selected[u]) {
                    nodes_in.push_back(u);
                    cluster_prize_sum += inst.prizes[u];
                } else {
                    nodes_out.push_back(u);
                }
            }

            // 1. DROP (Remover)
            if (!nodes_in.empty()) {
                for (int u : nodes_in) {
                    if (cluster_prize_sum - inst.prizes[u] >= inst.min_prize_by_cluster[c]) {
                        std::vector<int> candidate;
                        candidate.reserve(current_subgraph.size());
                        for(int x : current_subgraph) if(x != u) candidate.push_back(x);

                        long long new_cost = calc_solution_cost(inst, candidate);
                        if (new_cost < current_cost) {
                            current_subgraph = candidate;
                            current_cost = new_cost;
                            is_selected[u] = false;
                            improved = true;
                            break;
                        }
                    }
                }
                if (improved) break;
            }

            // 2. ADD (Adicionar)
            for (int v : nodes_out) {
                std::vector<int> candidate = current_subgraph;
                candidate.push_back(v);

                long long new_cost = calc_solution_cost(inst, candidate);
                if (new_cost < current_cost) {
                    current_subgraph = candidate;
                    current_cost = new_cost;
                    is_selected[v] = true;
                    improved = true;
                    break; // <--- EARLY RETURN
                }
            }
            if (improved) break; 

            // 3. SWAP (Trocar)
            for (int u : nodes_in) {
                for (int v : nodes_out) {
                    if (cluster_prize_sum - inst.prizes[u] + inst.prizes[v] >= inst.min_prize_by_cluster[c]) {
                        std::vector<int> candidate;
                        candidate.reserve(current_subgraph.size());
                        for(int x : current_subgraph) if(x != u) candidate.push_back(x);
                        candidate.push_back(v);

                        long long new_cost = calc_solution_cost(inst, candidate);
                        if (new_cost < current_cost) {
                            current_subgraph = candidate;
                            current_cost = new_cost;
                            is_selected[u] = false;
                            is_selected[v] = true;
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


void ils(PCInstance& inst) {
    auto current_subgraph = build_random_subgraph(inst);
    auto best_subgraph = current_subgraph;

    int max_iterations = inst.N * inst.N;

    int lim_p = 0;
    for (int i = 0; i < inst.C; i++) {
        if (inst.clusters[i].size() > 1) lim_p++;
    }

    int perturb = std::min(lim_p, 4);
    int current_cost = calc_solution_cost(inst, current_subgraph);
    int best_cost = current_cost;

    for (int iteration = 0; iteration < max_iterations; iteration++) {

        local_search(inst, current_subgraph, current_cost);

        if (current_cost < best_cost) {
            best_cost = current_cost;
            best_subgraph = current_subgraph;
            std::cout << "Nova melhor solucao na iter " << iteration << ": " << best_cost << "\n";
        }

        // 3. Perturbação (Kick) para escapar de mínimos locais
        // Sorteia alguns clusters e re-gera seus vértices aleatoriamente
        std::vector<int> next_subgraph = current_subgraph; // Copia o estado atual (não o best)
        
        for (int k = 0; k < perturb; k++) {
            int c = uniform(0, inst.C - 1); // Escolhe cluster aleatório
            
            // Remove vértices desse cluster da solução atual
            for (int u : inst.clusters[c]) {
                next_subgraph.erase(std::remove(next_subgraph.begin(), next_subgraph.end(), u), next_subgraph.end());
            }

            // Adiciona novos vértices aleatórios até bater a meta
            std::vector<int> cluster_nodes = inst.clusters[c];
            std::shuffle(cluster_nodes.begin(), cluster_nodes.end(), rng);
            
            int p_sum = 0;
            for (int node : cluster_nodes) {
                if (p_sum >= inst.min_prize_by_cluster[c]) break;
                next_subgraph.push_back(node);
                p_sum += inst.prizes[node];
            }
        }
        
        // Aceita a perturbação incondicionalmente para explorar
        current_subgraph = next_subgraph;
        current_cost = calc_solution_cost(inst, current_subgraph);
    }
    


}
std::vector<std::vector<int>> build_weight_matrix(const PCInstance& inst) {
    std::vector<std::vector<int>> weight_matrix(inst.N, std::vector<int>(inst.N, INF));

    // Inicializa a diagonal com 0 (distância de um nó a si mesmo)
    for (int i = 0; i < inst.N; ++i) {
        weight_matrix[i][i] = 0;
    }

    // Preenche com os pesos das arestas existentes
    for (const auto& edge : inst.edges) {
        weight_matrix[edge.u][edge.v] = edge.weight;
        weight_matrix[edge.v][edge.u] = edge.weight;
    }

    // Marca distâncias intra-cluster com INF
    // (Opcional, pois inst.edges já não as inclui, mas é bom para clareza)
    if(!allowIntra) {
    for (int c = 0; c < inst.C; ++c) {
        for (size_t i = 0; i < inst.clusters[c].size(); ++i) {
            for (size_t j = i + 1; j < inst.clusters[c].size(); ++j) {
                int u = inst.clusters[c][i];
                int v = inst.clusters[c][j];
                weight_matrix[u][v] = INF;
                weight_matrix[v][u] = INF;
            }
        }
    }
}
    return weight_matrix;
}
int eliminate_bottleneck_edges(PCInstance& inst, const std::vector<std::vector<int>>& weights) {
    int N = inst.N;
    int C = inst.C;
    int eliminated_count = 0;

    // A. Identificar arestas a serem removidas
    // Não podemos modificar inst.edges enquanto iteramos sobre ela.
    std::set<std::pair<int, int>> edges_to_keep;
    
    // A complexidade do teste é O(N*C). Executando para M arestas: O(M*N*C).
    for (const auto& edge : inst.edges) {
        int u = edge.u;
        int v = edge.v;
        int w_uv = edge.weight;
        bool is_bottleneck = false; // Flag para indicar se a aresta DEVE ser eliminada
        
        // Se a aresta foi eliminada em uma iteração anterior, continue
        // Esta verificação é opcional, mas evita re-processamento em implementações mais complexas.
        
        // Itera sobre todos os clusters 'c'
        for (int c = 0; c < C; ++c) {
            // Se u ou v pertencem ao cluster c, pula o cluster
            if ((inst.cluster_by_node[u] == c || inst.cluster_by_node[v] == c) && !allowIntra) {
                continue;
            }

            bool cluster_test_passed = true;
            
            // Testa se w_uv é maior que o gargalo de TODOS os nós 'x' no cluster 'c'
            for (int x : inst.clusters[c]) {
                int w_ux = weights[u][x];
                int w_vx = weights[v][x];

                // Condição de gargalo: w_ux ou w_vx deve ser menor que w_uv
                // E w_ux e w_vx DEVEM ser distâncias finitas (não intra-cluster)
                if (w_ux == INF || w_vx == INF) {
                    // A aresta não pode ser eliminada se o caminho por x for "infinito"
                    cluster_test_passed = false;
                    break; 
                }
                
                // Se w_uv NÃO for maior que AMBOS w_ux E w_vx, o cluster 'c' não elimina (bottleneck) a aresta.
                // A condição para ELIMINAR é: PARA TODO x in c, w_uv > max(w_ux, w_vx)
               if (w_ux > w_uv || w_vx > w_uv) { 
                    cluster_test_passed = false;
                    break;
                }
            }

            if (cluster_test_passed) {
                // Se o teste passou para PELO MENOS UM cluster 'c', a aresta é um gargalo e deve ser removida.
                is_bottleneck = true;
                break; // Não precisa testar outros clusters
            }
        }

        if (!is_bottleneck) {
            // Se não for gargalo, mantemos a aresta
            if (u < v) {
                edges_to_keep.insert({u, v});
            } else {
                edges_to_keep.insert({v, u});
            }
        } else {
            eliminated_count++;
        }
    }

    // B. Reconstituir inst.edges
    std::vector<Edge> new_edges;
    for (const auto& edge : inst.edges) {
        int u = edge.u;
        int v = edge.v;
        if (u > v) std::swap(u, v); // Normaliza (u, v)
        
        if (edges_to_keep.count({u, v})) {
            new_edges.push_back(edge);
        }
    }
    
    inst.edges = new_edges;
    return eliminated_count;
}

void preprocess(PCInstance& inst) {
    std::cout << "-> Iniciando pre-processamento de eliminacao de arestas...\n";
    
    // 1. Constrói a matriz de pesos para acesso rápido
    // Nota: O cálculo da distância euclidiana é feito apenas uma vez em read_pcinstance.
    // Esta matriz armazena esses pesos.
    std::vector<std::vector<int>> weights = build_weight_matrix(inst);
    
    // 2. Executa a eliminação de arestas
    int initial_M = inst.edges.size();
    int eliminated = eliminate_bottleneck_edges(inst, weights);

    if (eliminated > 0) {
        std::cout << "-> Arestas eliminadas: " << eliminated << " de " << initial_M
                  << " (" << std::fixed << std::setprecision(2) << ((double)eliminated / initial_M) * 100.0 << "%)\n";
        std::cout << "-> Numero de arestas restantes: " << inst.edges.size() << "\n";
    } else {
        std::cout << "-> Nenhuma aresta foi eliminada. Total: " << initial_M << "\n";
    }
}


int main() {
    allowIntra = false;
    auto inst = read_pcinstance();
    std::cout << inst.edges.size() << "\n";
    preprocess(inst);
    ils(inst);
    return 0;
    
}