#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;

const int min_L = 64;           // 短辺の最小サイズ
const int n_samples = 10000;      // 各アスペクト比でのサンプル数
const double beta_c = 0.44068679350977147533;     // (1/2) * log(1 + sqrt(2))

// 調査するアスペクト比 r = Lx / Ly のリスト
vector<double> r_list = {0.25, 0.5, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0};


// 高速な連結判定のためのUnion-Find構造体
struct DSU {
    vector<int> parent;
    DSU(int n) : parent(n) {
        iota(parent.begin(), parent.end(), 0);
    }
    int find(int i) {
        if (parent[i] == i) return i;
        return parent[i] = find(parent[i]);
    }
    void unite(int i, int j) {
        int root_i = find(i);
        int root_j = find(j);
        if (root_i != root_j) parent[root_i] = root_j;
    }
    bool connected(int i, int j) {
        return find(i) == find(j);
    }
};

class IsingSimulation {
    int Lx, Ly;
    vector<int> spins;
    double p_bond;
    mt19937 gen;
    uniform_real_distribution<double> dist;

    int get_idx(int x, int y) const { return y * Lx + x; }

public:
    IsingSimulation(int lx, int ly, double beta, int seed)
        : Lx(lx), Ly(ly), gen(seed), dist(0.0, 1.0) {
        p_bond = 1.0 - exp(-2.0 * beta);
        spins.assign(Lx * Ly, 0);
        
        // 境界条件: 左右(x方向)を+1, 上下(y方向)を-1に固定
        for (int y = 0; y < Ly; ++y) {
            for (int x = 0; x < Lx; ++x) {
                if (x == 0 || x == Lx - 1) spins[get_idx(x, y)] = 1;
                else if (y == 0 || y == Ly - 1) spins[get_idx(x, y)] = -1;
                else spins[get_idx(x, y)] = (dist(gen) < 0.5) ? 1 : -1;
            }
        }
    }

    void update_wolff() {
        uniform_int_distribution<int> distX(1, Lx - 2);
        uniform_int_distribution<int> distY(1, Ly - 2);
        int sx = distX(gen), sy = distY(gen);
        int start_idx = get_idx(sx, sy);
        int old_s = spins[start_idx];
        int new_s = -old_s;

        vector<int> q = {start_idx};
        spins[start_idx] = new_s;
        int head = 0;

        while(head < (int)q.size()){
            int curr = q[head++];
            int cx = curr % Lx, cy = curr / Lx;
            int dx[] = {1, -1, 0, 0}, dy[] = {0, 0, 1, -1};

            for(int i=0; i<4; ++i){
                int nx = cx + dx[i], ny = cy + dy[i];
                if(nx >= 0 && nx < Lx && ny >= 0 && ny < Ly){
                    int next_idx = get_idx(nx, ny);
                    bool is_boundary = (nx == 0 || nx == Lx - 1 || ny == 0 || ny == Ly - 1);
                    if(spins[next_idx] == old_s && dist(gen) < p_bond){
                        if(!is_boundary){
                            spins[next_idx] = new_s;
                            q.push_back(next_idx);
                        }
                    }
                }
            }
        }
    }

    pair<bool, bool> check_crossings() {
        DSU dsu_spin(Lx * Ly + 2);
        DSU dsu_fk(Lx * Ly + 2);
        int LEFT = Lx * Ly, RIGHT = Lx * Ly + 1;

        for (int y = 0; y < Ly; ++y) {
            dsu_spin.unite(get_idx(0, y), LEFT);
            dsu_spin.unite(get_idx(Lx - 1, y), RIGHT);
            dsu_fk.unite(get_idx(0, y), LEFT);
            dsu_fk.unite(get_idx(Lx - 1, y), RIGHT);
        }

        for (int y = 0; y < Ly; ++y) {
            for (int x = 0; x < Lx; ++x) {
                int curr = get_idx(x, y);
                int dx_list[] = {1, 0}, dy_list[] = {0, 1};
                for(int i=0; i<2; ++i){
                    int nx = x + dx_list[i], ny = y + dy_list[i];
                    if(nx < Lx && ny < Ly){
                        int next = get_idx(nx, ny);
                        if(spins[curr] == spins[next]){
                            if(spins[curr] == 1) dsu_spin.unite(curr, next);
                            if(dist(gen) < p_bond) dsu_fk.unite(curr, next);
                        }
                    }
                }
            }
        }
        return {dsu_spin.connected(LEFT, RIGHT), dsu_fk.connected(LEFT, RIGHT)};
    }
};

int main() {
    string filename = "crossing_probs_L" + to_string(min_L) + ".csv";
    ofstream ofs(filename);
    ofs << "L_min,r,Lx,Ly,P_spin,P_fk" << endl;

    auto start = chrono::high_resolution_clock::now();
    
    cout << fixed << setprecision(3);
    for (double r : r_list) {
        int Lx, Ly;
        if (r >= 1.0) {
            Ly = min_L;
            Lx = (int)round(min_L * r);
        } else {
            Lx = min_L;
            Ly = (int)round(min_L / r);
        }

        IsingSimulation sim(Lx, Ly, beta_c, 42);

        // 熱浴
        for(int i=0; i<1000; ++i) sim.update_wolff();

        int spin_cnt = 0, fk_cnt = 0;
        for(int i=0; i<n_samples; ++i){
            for(int j=0; j<10; ++j) sim.update_wolff();
            auto [s_cross, fk_cross] = sim.check_crossings();
            if(s_cross) spin_cnt++;
            if(fk_cross) fk_cnt++;
        }

        double p_s = (double)spin_cnt / n_samples;
        double p_fk = (double)fk_cnt / n_samples;

        cout << "r=" << r << " (" << Lx << "x" << Ly << ") -> P_spin: " << p_s << ", P_fk: " << p_fk << endl;
        ofs << min_L << "," << r << "," << Lx << "," << Ly << "," << p_s << "," << p_fk << endl;
    }

    ofs.close();
    
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    // ログファイルへの追記
    ofstream log_file("log.txt", ios::app);
    if (log_file) {
        // 現在時刻の取得（オプション）
        auto now = chrono::system_clock::to_time_t(chrono::system_clock::now());
        
        log_file << "--- Simulation Log ---" << endl;
        log_file << "Date: " << ctime(&now); // 実行日時
        log_file << "L_min: " << min_L << ", n_samples: " << n_samples << endl;
        log_file << "Num_r: " << r_list.size() << " (" << *std::min_element(r_list.begin(), r_list.end()) << " to " << *std::max_element(r_list.begin(), r_list.end()) << ")" << endl;
        log_file << "Elapsed time: " << fixed << setprecision(2) << elapsed.count() << " seconds" << endl;
        log_file << "-----------------------" << endl << endl;
        log_file.close();
    }

    cout << "Simulation completed in " << elapsed.count() << " seconds." << endl;
    
    cout << "\nResults saved to crossing_probs.csv" << endl;

    return 0;
}
