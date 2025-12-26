#ifndef MODEL_HPP
#define MODEL_HPP

#include <vector>
#include <utility>

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
    
    std::vector<std::pair<int,int>> coords;
    std::vector<std::vector<int>> clusters;
    std::vector<int> cluster_by_node;
    std::vector<int> centers;
    std::vector<int> prizes;
    std::vector<int> min_prize_by_cluster;
    std::vector<Edge> edges;
};

#endif
