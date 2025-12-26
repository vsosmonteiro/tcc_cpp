#include "preprocess.hpp"
#include "model.hpp"
#include <iostream>
#include <iomanip>
#include <set>
const int INF = 0x3f3f3f3f;
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
//     if(!allowIntra) {
//     for (int c = 0; c < inst.C; ++c) {
//         for (size_t i = 0; i < inst.clusters[c].size(); ++i) {
//             for (size_t j = i + 1; j < inst.clusters[c].size(); ++j) {
//                 int u = inst.clusters[c][i];
//                 int v = inst.clusters[c][j];
//                 weight_matrix[u][v] = INF;
//                 weight_matrix[v][u] = INF;
//             }
//         }
//     }
// }
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
            if ((inst.cluster_by_node[u] == c || inst.cluster_by_node[v] == c) && false) {
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

