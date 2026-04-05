#ifndef EDQC_H
#define EDQC_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <string>
#include <queue>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <limits>
#include <functional>
#include <set>
using namespace std;

struct Graph {
    int n, m; // v_num and e_num
    vector<vector<int>> adj; // adjacency list
    vector<int> degree; // degree list
    unordered_map<uint64_t, bool> edge_exists; // edge existence map
    Graph(int _n = 0);
    void add_edge(int u, int v);
    bool has_edge(int u, int v) const;
};

inline uint64_t edge_key(int u, int v) {
    if (u > v) std::swap(u, v);
    return ((uint64_t)u << 32) | (uint64_t)v;
}

Graph::Graph(int _n) : n(_n), m(0) {
    adj.resize(n);
    degree.resize(n, 0);
}

void Graph::add_edge(int u, int v) {
    if (u == v) return; 
    uint64_t key = edge_key(u, v);
    if (edge_exists.count(key)) return; 
    
    adj[u].push_back(v);
    adj[v].push_back(u);
    degree[u]++;
    degree[v]++;
    edge_exists[key] = true;
    m++;
}

bool Graph::has_edge(int u, int v) const {
    return edge_exists.count(edge_key(u, v)) > 0;
}

Graph read_graph(const string& filename);
double compute_local_clustering(const Graph& G, int node);
vector<int> edqc(const Graph& G, double density_threshold, int verbose = 1);
vector<int> extract_QC(const vector<double>& energy, const vector<int>& active_nodes, 
                                      const Graph& G, double density_threshold);
int find_energy_breakpoint(const vector<double>& energy_spectrum);
double compute_density(const Graph& G, const vector<int>& nodes);
bool reduction_rule(const Graph& G, int v, double gamma);
vector<int> compute_kcore(const Graph& G);

Graph read_graph(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    string line, type;
    int n = 0, m = 0;
    
    while (getline(file, line)) {
        istringstream iss(line);
        iss >> type;
        if (type == "p") {
            string format;
            iss >> format >> n >> m;
            break;
        }
    }
    Graph G(n);
    file.clear();
    file.seekg(0);
    while (getline(file, line)) {
        istringstream iss(line);
        iss >> type;
        if (type == "e") {
            int u, v;
            iss >> u >> v;
            G.add_edge(u-1, v-1);
       }
    }
    return G;
}

#endif // EDQC