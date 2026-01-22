//
//  FastNBSim.cpp
//
//  Created by on 2023/4/5.
//

#include "FastNBSim.h"
#include <string>
#include <ctime>

void read_graph(string filename) {
    clock_t begin = clock();
    ifstream ifile;
    ifile.open(filename);
    if(ifile.fail())
    {
        cout << "### Error Open, File Name: " << filename << endl;
        return;
    }
    
    string type, keyword;
    int u, v;
    
    // 读取第一行 "p edge 节点数 边数"
    ifile >> type;
    if (type != "p") {
        cout << "### Error: Not a valid DIMACS format, expected 'p' at the beginning" << endl;
        return;
    }
    
    ifile >> keyword >> n >> m;
    if (keyword != "edge") {
        cout << "### Error: Not a valid DIMACS edge format" << endl;
        return;
    }
    
    initial_variable();
    
    // 读取每条边 "e 节点1 节点2"
    for (int i = 0; i < m; ++i) {
        ifile >> type >> u >> v;
        if (type != "e") {
            cout << "### Error: Expected edge definition 'e'" << endl;
            return;
        }
        
        // 如果节点编号是从1开始的，可能需要减1来适应从0开始的数组索引
        // 根据您的需求决定是否需要调整节点编号
        
        degree[u-1]++;
        degree[v-1]++;
        neighbor[u-1].emplace_back(v-1);
        neighbor[v-1].emplace_back(u-1);
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
}

void initial_variable() {
    degree.resize(n);
    neighbor.resize(n);
    kcore_size.resize(n);
    vis.resize(n);

    if (use_signature) {
        vec_a.resize(k);
        vec_b.resize(k);
        signature.resize(n);
    }
    for (int i = 0; i < k; i++) {
        vec_a[i] = rand() % (p - 1) + 1;
        vec_b[i] = rand();
    }
}

void get_k_signatures(int node) {
    for (int index = 0; index < k; ++index) {
        int minhash = p;
        for (int j : neighbor[node]) {
            minhash = min(minhash, int(((long long)(vec_a[index]) * (long long)(j) + (long long)(vec_b[index])) % p));
        }
        minhash = min(minhash, int(((long long)(vec_a[index]) * (long long)(node) + (long long)(vec_b[index])) % p));
        signature[node].emplace_back(minhash);
    }
}

void sort_vertex() {
    for (int i = 0; i < n; ++i) {
        int num = core_num(i);
        if (num > 0) sorted_vertex.emplace_back(pair<int,int>(num, i));
        
    }
    sort(sorted_vertex.begin(), sorted_vertex.end(), greater<pair<int,int>>());
}

void iterate_vertex() {
    for (int i = 0; i < sorted_vertex.size(); ++i) {
        auto p = sorted_vertex[i];
        if (p.first > cur_sol_size) {
            operations++;
            int num = sift_num(p.second);
            if (num > cur_sol_size) {
                cur_sol_size = num;
                cur_sol_vertex = p.second;
            }
        }
    }
}

int core_num(int node) {
    int num = 0;
    for (int v : neighbor[node]) {
        if ((degree[v]+1) >= (degree[node]+1) * ts) num++;
    }
    return num;
}

double ct_score(int v1, int v2) {
    int res = 0;
    unordered_set<int> st;
    st.insert(v1);
    for (int num : neighbor[v1]) st.insert(num);
    for (int num : neighbor[v2]) if (st.count(num)) res++;
    if (st.count(v2)) res++;
    return (double)res / (double)(degree[v1] + 1);
}

double ct_score_minhash(int v1, int v2) {
    int res = 0;
    for (int i = 0; i < k; ++i) {
        if (signature[v1][i] == signature[v2][i]) res++;
    }
    double j_sim = (double)res / (double)k;
    return (double)(degree[v1] + degree[v2] + 2) * j_sim / (j_sim + 1) / (double)(degree[v1] + 1);
}

int sift_num(int node) {
    int res = 0;
    if (use_signature){
       if (vis[node] == 0) {
           get_k_signatures(node);
           vis[node] = 1;
       }
       for (int v : neighbor[node]) {
           if (vis[v] == 0) {
               get_k_signatures(v);
               vis[v] = 1;
           }
           if (ct_score_minhash(node, v) > ts) res++;
       }
    }
    else {
        for (int v : neighbor[node]) {
          if (ct_score(node, v) >= ts) res++;
        }
    }
    double tmp = (double)res / (double)(degree[node]+1);
    if (tmp < b) return 0;
    return res + 1;
}

// void get_res() {
//     sort_vertex();
//     iterate_vertex();
//     cur_sol.emplace_back(cur_sol_vertex);
//     for (int v : neighbor[cur_sol_vertex]) {
//         if (use_signature) {
//             if (ct_score_minhash(cur_sol_vertex, v) > ts) cur_sol.emplace_back(v);
//         }
//         else {
//             if (ct_score(cur_sol_vertex, v) > ts) cur_sol.emplace_back(v);
//         }
//     }
    
//     for (int i : cur_sol) {
//         for (int i2 : neighbor[i]) {
//             if(find(cur_sol.begin(), cur_sol.end(), i2) != cur_sol.end()) cur_sol_edges++;
//         }
//     }
    
//     cur_sol_quasi = double(cur_sol_edges) / double(cur_sol_size * (cur_sol_size-1));
// }

void get_res() {
    clock_t start_time = clock();  // 记录开始时间
    
    sort_vertex();
    
    // 修改 iterate_vertex 函数来包含时间检查
    for (int i = 0; i < sorted_vertex.size(); ++i) {
        // 检查是否超时
        if ((double)(clock() - start_time) / CLOCKS_PER_SEC > 60) {
            // cout << "Time limit of 60 seconds reached. Returning best solution found so far." << endl;
            break;  // 终止循环，使用当前找到的最佳解
        }
        
        auto p = sorted_vertex[i];
        if (p.first > cur_sol_size) {
            operations++;
            int num = sift_num(p.second);
            if (num > cur_sol_size) {
                cur_sol_size = num;
                cur_sol_vertex = p.second;
            }
        }
    }
    
    cur_sol.emplace_back(cur_sol_vertex);
    for (int v : neighbor[cur_sol_vertex]) {
        if (use_signature) {
            if (ct_score_minhash(cur_sol_vertex, v) > ts) cur_sol.emplace_back(v);
        }
        else {
            if (ct_score(cur_sol_vertex, v) > ts) cur_sol.emplace_back(v);
        }
    }
    
    for (int i : cur_sol) {
        for (int i2 : neighbor[i]) {
            if(find(cur_sol.begin(), cur_sol.end(), i2) != cur_sol.end()) cur_sol_edges++;
        }
    }
    
    cur_sol_quasi = double(cur_sol_edges) / double(cur_sol_size * (cur_sol_size-1));
}

int main(int argc, const char * argv[]) {
    
    if (argc < 3 || argc > 6){
        cout << "wrong input format!!" << endl;
        return 0;
    }
    if (argc == 5) {
        k = atoi(argv[4]);
        use_signature = true;
    }
    string filename;
    filename = argv[1];
    ts = atof(argv[2]);
    b = atof(argv[3]);
    
    if (argc == 6) {
        int seed = atoi(argv[5]);
        srand(seed);
        k = atoi(argv[4]);
        use_signature = true;
    }
    
    read_graph(filename);

    clock_t begin = clock();
    get_res();
    clock_t end = clock();
    
    double run_secs = double(end - begin) / CLOCKS_PER_SEC;
    if(cur_sol_size){
        cout << filename << " " << ts << " " << b << " " << k << " " << cur_sol_size << " " << cur_sol_quasi << " ";
        cout << setprecision(4) << run_secs << endl;
    }
    // cout << "********************************" << endl;
    return 0;
}
