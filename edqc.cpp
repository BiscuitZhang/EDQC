#include "edqc.h"

double compute_density(const Graph &G, const vector<int> &nodes)
{
    if (nodes.size() <= 1)
        return 0.0;

    static vector<bool> in_set;
    if (in_set.size() < G.n)
    {
        in_set.resize(G.n);
    }
    fill(in_set.begin(), in_set.begin() + G.n, false);

    for (int v : nodes)
        in_set[v] = true;

    int edges = 0;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        int u = nodes[i];
        for (int v : G.adj[u])
        {
            if (in_set[v] && u < v)
            {
                edges++;
            }
        }
    }

    double possible_edges = (nodes.size() * (nodes.size() - 1)) / 2.0;
    return edges / possible_edges;
}

namespace
{
    // Perform one round of the adaptive energy diffusion process.
    // The release rate of each vertex is adjusted by the energy in its neighborhood
    // to concentrate scores within locally cohesive regions.
    static void OneRoundDiffuse(const Graph &G,
                                const std::vector<double> &Ein,
                                std::vector<double> &Eout,
                                const std::vector<int> &actives)
    {
        Eout = Ein; // E_new = E_cur
        for (int node : actives)
        {
            if (G.adj[node].empty())
                continue;

            std::vector<double> w;
            w.reserve(G.adj[node].size());
            double sumNbrE = 0.0, sumw = 0.0;
            for (int nbr : G.adj[node])
            {
                sumNbrE += Ein[nbr];
                double U = (std::rand() + 1.0) / (RAND_MAX + 1.0); // U∈(0,1]
                double g = -std::log(U);                           // Exp(1)
                w.emplace_back(g);
                sumw += g;
            }
            double denom = sumNbrE + 2.0 * Ein[node];
            double eta = (denom > 0.0) ? ((sumNbrE + Ein[node]) / denom) : 0.0;
            double emit = Ein[node] * eta;

            double remaining = emit;
            const int deg = (int)G.adj[node].size();
            for (int i = 0; i < deg; ++i)
            {
                int nbr = G.adj[node][i];
                double share = 0.0;
                if (i + 1 < deg)
                {
                    share = (sumw > 0.0) ? (emit * (w[i] / sumw)) : 0.0;
                    remaining -= share;
                }
                else
                {
                    share = remaining;
                }
                Eout[nbr] += share;
            }
            Eout[node] -= emit;
        }
    }

    static void BuildActive(const std::vector<double> &e,
                            double theta,
                            std::vector<int> &out,
                            std::vector<char> &mask)
    {
        out.clear();
        if ((int)mask.size() < (int)e.size())
        {
            mask.assign(e.size(), 0);
        }
        else
        {
            std::fill(mask.begin(), mask.end(), 0);
        }
        for (int v = 0; v < (int)e.size(); ++v)
        {
            if (e[v] > theta)
            {
                out.push_back(v);
                mask[v] = 1;
            }
        }
    }

  
    inline double DensityBy(int s, long long e)
    {
        return (s <= 1) ? 0.0 : (2.0 * (double)e) / (1.0 * s * (s - 1));
    }
    inline bool DensityOK(int s, long long e, double gamma)
    {
        return (s > 1) && (DensityBy(s, e) + 1e-15 >= gamma);
    }
    inline int NeedEdges(int s, long long e, double gamma)
    {
        // γ*C(s+1,2)
        long double target = (long double)gamma * ((long double)(s + 1) * s / 2.0L);
        long double need = target - (long double)e;
        if (need <= 0.0L)
            return 0;
        return (int)std::ceil(need);
    }

    // remove S[pos] and maintain eS / deg_in / inS
    static void RemoveAt(std::vector<int> &S, int pos,
                         std::vector<char> &inS,
                         std::vector<int> &deg_in,
                         long long &eS,
                         const Graph &G)
    {
        const int v = S[pos];
        inS[v] = 0;
        for (int u : G.adj[v])
            if (inS[u])
            {
                --deg_in[u];
                --eS;
            }
        S[pos] = S.back();
        S.pop_back();
    }

    // Single-round expansion: Try adding vertices from R to S by sorting by (connections, energy).
    static bool ExpandOnce(const Graph &G,
                           const std::vector<double> &energy,
                           double gamma,
                           std::vector<int> &S,
                           std::vector<int> &R,
                           std::vector<char> &inS,
                           std::vector<int> &deg_in,
                           long long &eS)
    {
        struct Cand
        {
            int conn;
            double en;
            int v;
        };
        std::vector<Cand> cand;
        cand.reserve(R.size());

        int s = (int)S.size();
        long long e = eS;
        int Tneed = NeedEdges(s, e, gamma);

        for (int v : R)
        {
            if (s == 0)
                continue;
            if (Tneed > 0 && std::min(G.degree[v], s) < Tneed)
                continue; 
            int conn = 0;
            for (int u : G.adj[v])
                if (inS[u])
                    ++conn;
            cand.push_back({conn, energy[v], v});
        }
        if (cand.empty())
            return false;

        std::sort(cand.begin(), cand.end(),
                  [](const Cand &a, const Cand &b)
                  {
                      if (a.conn != b.conn)
                          return a.conn > b.conn;
                      return a.en > b.en;
                  });

        std::vector<char> inR(G.n, 0);
        for (int v : R)
            inR[v] = 1;

        bool any = false;
        for (const auto &c : cand)
        {
            const int v = c.v;
            if (!inR[v])
                continue;

            const int conn = c.conn;
            const int s2 = s + 1;
            const long long e2 = e + conn;

            if (DensityOK(s2, e2, gamma))
            {
                inS[v] = 1;
                S.push_back(v);
                deg_in[v] = 0;
                for (int u : G.adj[v])
                    if (inS[u] && u != v)
                    {
                        ++deg_in[u];
                        ++deg_in[v];
                    }
                eS = e = e2;
                s = (int)S.size();
                Tneed = NeedEdges(s, e, gamma);
                inR[v] = 0;
                any = true;
            }
        }

        if (any)
        {
            std::vector<int> R2;
            R2.reserve(R.size());
            for (int v : R)
                if (!inS[v])
                    R2.push_back(v);
            R.swap(R2);
        }
        return any;
    }
}

// Converts the energy function ranking into a feasible gamma-quasi-clique.
// This executes the candidate extraction and density-constrained refinement steps:
// Stage 1: Conductance-based sweep cut to extract an explicit candidate region
// Stage 2: Pruning and greedy expansion to strictly satisfy the explicit density constraint
std::vector<int> extract_QC(const std::vector<double> &energy, const std::vector<int> &active_nodes, const Graph &G, double density_threshold)
{
    // 1.1 Rank vertices by their remaining energy to localize promising candidate regions
    std::vector<std::pair<double, int>> P;
    P.reserve(active_nodes.size());
    for (int v : active_nodes)
        if (energy[v] > 0.0)
            P.emplace_back(energy[v], v);

    if (P.size() < 2)
        return {};
    std::sort(P.begin(), P.end(), std::greater<std::pair<double, int>>());

    // 1.2 Prepare for conductance-based sweep cut extraction
    std::vector<int> S_initial;
    double min_conductance = std::numeric_limits<double>::max();

    std::vector<int> current_S;
    current_S.reserve(P.size());
    std::vector<char> in_current_S(G.n, 0);

    long long current_cut = 0;
    long long current_vol = 0;
    long long total_vol = 0; 
    for (const auto &pe : P)
        total_vol += G.degree[pe.second];
    if (total_vol == 0)
        return {};

    // 1.3 Conductance sweep cut to turn the ranking into an explicit candidate region
    for (const auto &pe : P)
    {
        int new_node = pe.second;
        long long new_internal_edges = 0; 

        for (int neighbor : G.adj[new_node])
        {
            if (in_current_S[neighbor])
            {
                new_internal_edges++;
            }
        }

        current_S.push_back(new_node);
        in_current_S[new_node] = 1;

        current_cut = current_cut + G.degree[new_node] - 2 * new_internal_edges;
        current_vol += G.degree[new_node];

        if (current_vol > 0 && total_vol - current_vol > 0)
        {
            double conductance = (double)current_cut / std::min(current_vol, total_vol - current_vol);

            if (current_S.size() >= 2 && conductance < min_conductance)
            {
                min_conductance = conductance;
                S_initial = current_S;
            }
        }
    }

    if (S_initial.empty())
    {
        if (P.size() >= 2)
        {
            for (size_t i = 0; i < std::min((size_t)P.size(), (size_t)2); ++i)
                S_initial.push_back(P[i].second);
        }
        else
        {
            return {};
        }
    }

    // 2.1 Initialize explicit candidate region extracted from the sweep cut prefix
    std::vector<int> S = S_initial;
    std::vector<char> inS(G.n, 0);
    std::vector<int> deg_in(G.n, 0);
    for (int v : S)
        inS[v] = 1;

    long long eS = 0;
    for (int v : S)
    {
        for (int u : G.adj[v])
        {
            if (inS[u] && u < v)
            {
                ++eS;
                ++deg_in[v];
                ++deg_in[u];
            }
        }
    }

    // 2.2 Pruning phase: Enforce feasibility under gamma explicitly.
    while ((int)S.size() > 1 && !DensityOK((int)S.size(), eS, density_threshold))
    {
        int min_pos = 0;
        int min_val = deg_in[S[0]];
        double min_en_tiebreak = energy[S[0]];
        for (int i = 1; i < (int)S.size(); ++i)
        {
            int v = S[i];
            if (deg_in[v] < min_val || (deg_in[v] == min_val && energy[v] < min_en_tiebreak))
            {
                min_val = deg_in[v];
                min_en_tiebreak = energy[v];
                min_pos = i;
            }
        }
        RemoveAt(S, min_pos, inS, deg_in, eS, G);
    }

    if ((int)S.size() < 2)
        return {};
    if (!DensityOK((int)S.size(), eS, density_threshold))
        return {};

    // 2.3 Greedy expansion: Enforce feasibility while maximizing size.
    std::vector<int> R;
    for (const auto &pe : P)
        if (!inS[pe.second])
            R.push_back(pe.second);

    while (true)
    {
        bool has_expanded = ExpandOnce(G, energy, density_threshold, S, R, inS, deg_in, eS);
        if (!has_expanded)
            break;
    }

    return S;
}


std::vector<int> edqc(const Graph &G, double density_threshold, int max_iterations, double activation_threshold, int verbose, double &best_find_time)
{
    auto print_message = [verbose](const std::string &msg, int level)
    {
        if (verbose >= level)
            std::cout << "[EDQC] " << msg << std::endl;
    };

    std::clock_t start_time = std::clock();

    // Process each source vertex v in non-increasing order of degree
    int max_degree = 0;
    for (int i = 0; i < G.n; ++i)
        max_degree = std::max(max_degree, G.degree[i]);

    std::vector<std::vector<int>> buckets(max_degree + 1);
    for (int i = 0; i < G.n; ++i)
        if (G.degree[i] >= 2)
            buckets[G.degree[i]].push_back(i);

    std::vector<int> sorted_nodes;
    sorted_nodes.reserve(G.n);
    for (int d = max_degree; d >= 0; --d)
    {
        if (!buckets[d].empty())
        {
            sorted_nodes.insert(sorted_nodes.end(), buckets[d].begin(), buckets[d].end());
        }
    }
    double sorting_time = double(std::clock() - start_time) / CLOCKS_PER_SEC;
    print_message("Node bucket-sort completed in " + std::to_string(sorting_time) +
                      "s (" + std::to_string(sorted_nodes.size()) + " nodes sorted)",
                  1);
    // End sort

    std::vector<int> best_quasi_clique;
    print_message("Diffusion with T (max_iterations)=" + std::to_string(max_iterations) +
                      ", theta=" + std::to_string(activation_threshold),
                  1);

    // Traversing nodes
    for (int seed : sorted_nodes)
    {
        if (double(std::clock() - start_time) / CLOCKS_PER_SEC > 60.0)
            break;

        std::vector<double> energy(G.n, 0.0);
        energy[seed] = 1.0;

        std::vector<int> active_nodes(1, seed);
        std::vector<char> mask_cur(G.n, 0), mask_next(G.n, 0);
        mask_cur[seed] = 1;

        std::vector<double> e_cur = energy;
        std::vector<int> A_cur = active_nodes;

        // Run the adaptive energy diffusion process for T synchronous rounds
        int round_id = 0;
        while (round_id < max_iterations)
        {
            ++round_id;

            // Run one round diffusion and obtain e_new
            std::vector<double> e_new;
            OneRoundDiffuse(G, e_cur, e_new, A_cur);

            // Build new active set A_new after one round diffusion
            std::vector<int> A_new;
            A_new.reserve(A_cur.size() * 2);
            BuildActive(e_new, activation_threshold, A_new, mask_next);

            // Whether has new active nodes
            bool has_new = (A_new.size() > A_cur.size());
            if (!has_new)
            {
                for (int v : A_new)
                    if (!mask_cur[v])
                    {
                        has_new = true;
                        break;
                    }
            }

            e_cur.swap(e_new);
            A_cur.swap(A_new);
            mask_cur.swap(mask_next);

            // shutdown
            // if (!has_new) break;
            // if (double(std::clock() - start_time) / CLOCKS_PER_SEC > 60.0) break;
        }

        energy.swap(e_cur);
        active_nodes.swap(A_cur);

        if ((int)active_nodes.size() <= (int)best_quasi_clique.size())
            continue;

        std::vector<int> quasi_clique = extract_QC(energy, active_nodes, G, density_threshold);

        if (quasi_clique.size() > best_quasi_clique.size())
        {
            best_quasi_clique.swap(quasi_clique);
            best_find_time = double(std::clock() - start_time) / CLOCKS_PER_SEC;
            if (verbose > 0)
            {
                double cur_den = compute_density(G, best_quasi_clique);
                print_message("New best quasi-clique: " + std::to_string(best_quasi_clique.size()) +
                                  " nodes, density: " + std::to_string(cur_den) +
                                  ", seed: " + std::to_string(seed) + ", " +
                                  std::to_string(best_find_time) + " s",
                              1);
            }
        }
    }

    return best_quasi_clique;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cout << "Usage: " << argv[0] << " <graph_file> <density_threshold> [max_iterations] [activation_threshold] [seed]" << endl;
        return 1;
    }

    // Input pararemeters
    string filename = argv[1];
    double density_threshold = atof(argv[2]);
    int max_iterations = (argc > 3) ? atoi(argv[3]) : 2;              // Default to 2 iterations
    double activation_threshold = (argc > 4) ? atof(argv[4]) : 0.0001; // Default to 0.0001
    int seed = (argc > 5) ? atoi(argv[5]) : 1;

    srand(seed);
    Graph G = read_graph(filename);

    double best_time = 0.0;

    // EDQC: 0/1-> print/no print
    vector<int> result = edqc(G, density_threshold, max_iterations, activation_threshold, 0, best_time);

    double density = compute_density(G, result);
    printf("%s, %d, %.6f, %.6f\n", filename.c_str(), (int)result.size(), density, best_time);

    return 0;
}
