#include "reader.hpp"

#include <random>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace std;

const int SEED = 1;

mt19937_64 rng(SEED);

int uniform(int l, int r)
{
    uniform_int_distribution<int> uid(l, r);
    return uid(rng);
}

PCInstance read_pcinstance()
{
    PCInstance inst;
    cerr << "Reading N..." << endl;
    
    cin >> inst.N;
    cerr << "Reading N..." << endl;

    inst.coords.resize(inst.N);
    for (int i = 0; i < inst.N; i++)
        cin >> inst.coords[i].first >> inst.coords[i].second;

    cin >> inst.C;

    inst.clusters.resize(inst.C);
    inst.min_prize_by_cluster.resize(inst.C);
    inst.cluster_by_node.resize(inst.N);

    for (int c = 0; c < inst.C; c++)
    {
        int size;
        cin >> size;
        inst.clusters[c].resize(size);
        for (int i = 0; i < size; i++)
        {
            int v;
            cin >> v;
            v--;
            inst.clusters[c][i] = v;
            inst.cluster_by_node[v] = c;
        }
    }
    cerr << "Reading clusters..." << endl;
    inst.centers.resize(inst.C);
    for (int c = 0; c < inst.C; c++)
    {
        cin >> inst.centers[c];
        inst.centers[c]--;
    }

    int validation;
    cin >> validation;
    if (validation != -999)
        throw runtime_error("Missing -999 marker");

    inst.prizes.resize(inst.N);
    for (int i = 0; i < inst.N; i++)
    {
        cin >> inst.prizes[i];
    }

    cin >> validation;
    if (validation != -999)
        throw runtime_error("Missing second -999 marker");

    inst.edges.clear();

    for (int i = 0; i < inst.N; i++)
    {
        for (int j = i + 1; j < inst.N; j++)
        {
            // if (inst.cluster_by_node[i] == inst.cluster_by_node[j] )
            //     continue;

            double dx = inst.coords[i].first - inst.coords[j].first;
            double dy = inst.coords[i].second - inst.coords[j].second;
            int w = (int)round(sqrt(dx * dx + dy * dy));
            inst.edges.push_back({i, j, w});
        }
    }

    sort(inst.edges.begin(), inst.edges.end());

    for (int c = 0; c < inst.C; c++)
    {
        cin >> inst.min_prize_by_cluster[c];
    }

    return inst;
}
