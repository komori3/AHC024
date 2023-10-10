#define _CRT_NONSTDC_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
#include <bits/stdc++.h>
#include <random>
#include <unordered_set>
#include <array>
#include <optional>
#ifdef _MSC_VER
#include <opencv2/core.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <conio.h>
#include <ppl.h>
#include <omp.h>
#include <filesystem>
#include <intrin.h>
/* g++ functions */
int __builtin_clz(unsigned int n) { unsigned long index; _BitScanReverse(&index, n); return 31 - index; }
int __builtin_ctz(unsigned int n) { unsigned long index; _BitScanForward(&index, n); return index; }
namespace std { inline int __lg(int __n) { return sizeof(int) * 8 - 1 - __builtin_clz(__n); } }
/* enable __uint128_t in MSVC */
//#include <boost/multiprecision/cpp_int.hpp>
//using __uint128_t = boost::multiprecision::uint128_t;
#else
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif

/** compro io **/
namespace aux {
    template<typename T, unsigned N, unsigned L> struct tp { static void output(std::ostream& os, const T& v) { os << std::get<N>(v) << ", "; tp<T, N + 1, L>::output(os, v); } };
    template<typename T, unsigned N> struct tp<T, N, N> { static void output(std::ostream& os, const T& v) { os << std::get<N>(v); } };
}
template<typename... Ts> std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) { os << '{'; aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t); return os << '}'; } // tuple out
template<class Ch, class Tr, class Container> std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x); // container out (fwd decl)
template<class S, class T> std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) { return os << '{' << p.first << ", " << p.second << '}'; } // pair out
template<class S, class T> std::istream& operator>>(std::istream& is, std::pair<S, T>& p) { return is >> p.first >> p.second; } // pair in
std::ostream& operator<<(std::ostream& os, const std::vector<bool>::reference& v) { os << (v ? '1' : '0'); return os; } // bool (vector) out
std::ostream& operator<<(std::ostream& os, const std::vector<bool>& v) { bool f = true; os << '{'; for (const auto& x : v) { os << (f ? "" : ", ") << x; f = false; } os << '}'; return os; } // vector<bool> out
template<class Ch, class Tr, class Container> std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) { bool f = true; os << '{'; for (auto& y : x) { os << (f ? "" : ", ") << y; f = false; } return os << '}'; } // container out
template<class T, class = decltype(std::begin(std::declval<T&>())), class = typename std::enable_if<!std::is_same<T, std::string>::value>::type> std::istream& operator>>(std::istream& is, T& a) { for (auto& x : a) is >> x; return is; } // container in
template<typename T> auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) { out << t.stringify(); return out; } // struct (has stringify() func) out
/** io setup **/
struct IOSetup { IOSetup(bool f) { if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); } std::cout << std::fixed << std::setprecision(15); } }
iosetup(true); // set false when solving interective problems
/** string formatter **/
template<typename... Ts> std::string format(const std::string& f, Ts... t) { size_t l = std::snprintf(nullptr, 0, f.c_str(), t...); std::vector<char> b(l + 1); std::snprintf(&b[0], l + 1, f.c_str(), t...); return std::string(&b[0], &b[0] + l); }
/** dump **/
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<']'<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }
/** timer **/
class Timer {
    double t = 0, paused = 0, tmp;
public:
    Timer() { reset(); }
    static double time() {
#ifdef _MSC_VER
        return __rdtsc() / 2.9e9;
#else
        unsigned long long a, d;
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
        return (d << 32 | a) / 2.9e9;
#endif
    }
    void reset() { t = time(); }
    void pause() { tmp = time(); }
    void restart() { paused += time() - tmp; }
    double elapsed_ms() const { return (time() - t - paused) * 1000.0; }
};
/** rand **/
struct Xorshift {
    static constexpr uint64_t M = INT_MAX;
    static constexpr double e = 1.0 / M;
    uint64_t x = 88172645463325252LL;
    Xorshift() {}
    Xorshift(uint64_t seed) { reseed(seed); }
    inline void reseed(uint64_t seed) { x = 0x498b3bc5 ^ seed; for (int i = 0; i < 20; i++) next(); }
    inline uint64_t next() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
    inline int next_int() { return next() & M; }
    inline int next_int(int mod) { return next() % mod; }
    inline int next_int(int l, int r) { return l + next_int(r - l + 1); }
    inline double next_double() { return next_int() * e; }
};
/** shuffle **/
template<typename T> void shuffle_vector(std::vector<T>& v, Xorshift& rnd) { int n = v.size(); for (int i = n - 1; i >= 1; i--) { int r = rnd.next_int(i); std::swap(v[i], v[r]); } }
/** split **/
std::vector<std::string> split(const std::string& str, const std::string& delim) {
    std::vector<std::string> res;
    std::string buf;
    for (const auto& c : str) {
        if (delim.find(c) != std::string::npos) {
            if (!buf.empty()) res.push_back(buf);
            buf.clear();
        }
        else buf += c;
    }
    if (!buf.empty()) res.push_back(buf);
    return res;
}
/** misc **/
template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) { std::fill((T*)array, (T*)(array + N), val); } // fill array
template<typename T, typename ...Args> auto make_vector(T x, int arg, Args ...args) { if constexpr (sizeof...(args) == 0)return std::vector<T>(arg, x); else return std::vector(arg, make_vector<T>(x, args...)); }
template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }

/* fast queue */
class FastQueue {
    int front = 0;
    int back = 0;
    int v[4096];
public:
    inline bool empty() { return front == back; }
    inline void push(int x) { v[front++] = x; }
    inline int pop() { return v[back++]; }
    inline void reset() { front = back = 0; }
    inline int size() { return front - back; }
};

class RandomQueue {
    int sz = 0;
    int v[4096];
public:
    inline bool empty() const { return !sz; }
    inline int size() const { return sz; }
    inline void push(int x) { v[sz++] = x; }
    inline void reset() { sz = 0; }
    inline int pop(int i) {
        std::swap(v[i], v[sz - 1]);
        return v[--sz];
    }
    inline int pop(Xorshift& rnd) {
        return pop(rnd.next_int(sz));
    }
};

#if 1
inline double get_temp(double stemp, double etemp, double t, double T) {
    return etemp + (stemp - etemp) * (T - t) / T;
};
#else
inline double get_temp(double stemp, double etemp, double t, double T) {
    return stemp * pow(etemp / stemp, t / T);
};
#endif

constexpr int N = 50;
constexpr int M = 100;
constexpr int L = 52;
constexpr int LL = L * L;

constexpr int d4[] = { 1, -L, -1, L };

// 境界を動かす遷移
// 隣接する二点 p1, p2 について、c1 = grid[p1], c2 = grid[p2] とする
// c1 != c2 のとき、境界である
// p2 の色を c1 に変更する
// c2 の領域が分断されてはいけない
// c1 境界の隣接関係は変化しない
// c2 境界の隣接関係は変化しない

// 空白を増やす遷移
// 隣接する二点 p1, p2 について、c1 = grid[p1] = 0, c2 = grid[p2] とする
// p2 の色を 0 に変更する
// c2 の領域が分断されてはいけない
// c1 境界の隣接関係は変化しない
// c2 境界の隣接関係は変化しない

struct State2 {

    std::array<int, LL> grid{}; // char にする？

    std::array<bool, LL> wall{};
    std::array<std::array<bool, M + 1>, M + 1> adj{};

    State2() {}

    State2(std::istream& in) {
        {
            int buf;
            in >> buf >> buf; // ignore N, M
        }
        wall.fill(true);
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                in >> grid[i * L + j];
                wall[i * L + j] = false;
            }
        }
        create_adj(adj);
    }

    bool is_valid(int pos, int c, int nc) const {
        int cc = 0;
        for (int d = 0; d < 4; d++) {
            cc += (grid[pos + d4[d]] == c);
        }
        if (cc == 0) return false; // 孤立点は消せない
        FastQueue fqu;
        std::array<bool, M + 1> a{};
        std::array<bool, LL> seen{};
        a[c] = a[nc] = true;
        seen[pos] = true;
        bool searched = false;
        for (int d = 0; d < 4; d++) {
            int spos = pos + d4[d];
            if (grid[spos] != c) continue;
            if (!searched) {
                // bfs
                fqu.reset();
                fqu.push(spos);
                seen[spos] = true;
                while (!fqu.empty()) {
                    int u = fqu.pop();
                    for (int d = 0; d < 4; d++) {
                        int v = u + d4[d];
                        if (seen[v] || grid[v] != c) {
                            a[grid[v]] = true;
                            continue;
                        }
                        seen[v] = true;
                        fqu.push(v);
                    }
                }
                searched = true;
            }
            else if (!seen[spos]) {
                return false;
            }
        }
        return a == adj[c];
    }

    bool is_valid_wall(int pos, int c, int nc) const {
        assert(c == 0 && nc != 0);
        if (wall[pos]) return false; // 簡単のため壁は除外する
        int cc = 0;
        for (int d = 0; d < 4; d++) {
            cc += (grid[pos + d4[d]] == c);
        }
        if (cc == 0) return false; // 孤立点は消せない
        FastQueue fqu;
        std::array<bool, M + 1> a{};
        std::array<bool, LL> seen{};
        a[c] = a[nc] = true;
        seen[pos] = true;
        bool searched = false;
        for (int d = 0; d < 4; d++) {
            int spos = pos + d4[d];
            if (grid[spos] != c) continue;
            if (wall[spos]) continue;
            if (!searched) {
                // bfs
                fqu.reset();
                fqu.push(spos);
                seen[spos] = true;
                while (!fqu.empty()) {
                    int u = fqu.pop();
                    for (int d = 0; d < 4; d++) {
                        int v = u + d4[d];
                        if (seen[v] || wall[v] || grid[v] != c) {
                            a[grid[v]] = true;
                            continue;
                        }
                        seen[v] = true;
                        fqu.push(v);
                    }
                }
                searched = true;
            }
            else if (!seen[spos]) {
                return false;
            }
        }
        return a == adj[c];
    }

    bool can_erode(int p1, int p2) const {
        if (wall[p2]) return false;
        if (grid[p1] == grid[p2]) return false;
        int c1 = grid[p1], c2 = grid[p2];
        // 浸食する側のチェックはすぐできる
        for (int d = 0; d < 4; d++) {
            int nc2 = grid[p2 + d4[d]];
            if (!adj[c1][nc2]) return false;
        }
        // 連結性と隣接関係を同時にチェック
        if (grid[p2] == 0) {
            // wall erode
            if (!is_valid_wall(p2, c2, c1)) return false;
        }
        else {
            if (!is_valid(p2, c2, c1)) return false;
        }
        return true;
    }

    void erode(int p1, int p2) {
        grid[p2] = grid[p1];
    }

    //void find_contour_o4(int spos) {
    //    int c = grid[spos];
    //    while (grid[spos] == c) spos -= L;
    //    // spos は境界上にある
    //    dump(spos / L, spos % L);
    //    int stop_dir = -1;
    //    int prev_dir = -1;
    //    int cpos = -1;
    //    for (int d = 0; d < 4; d++) {
    //        cpos = spos + d4[d];
    //        if (grid[cpos] == c) continue;
    //        prev_dir = stop_dir = (d + 2) & 3;
    //        break;
    //    }
    //    dump(cpos / L, cpos % L);
    //    do {
    //        int start_dir = (prev_dir + 1) & 3;
    //        for (int dir = start_dir; dir < start_dir + 4; dir++) {
    //            int d = dir & 3;
    //            int npos = cpos + d4[d], nc = grid[npos];
    //            if (nc == c) continue;
    //            prev_dir = (d + 2) & 3;
    //            cpos = npos;
    //            break;
    //        }
    //        dump(cpos / L, cpos % L);
    //    } while ((cpos + d4[prev_dir]) != spos);
    //}

    void create_adj(std::array<std::array<bool, M + 1>, M + 1>& adj) const {
        for (auto& v : adj) v.fill(0);
        for (int i = 0; i <= M; i++) {
            adj[i][i] = true;
        }
        for (int i = 0; i < L; i++) {
            for (int j = 1; j < L; j++) {
                int u = i * L + j, v = u - 1;
                int cu = grid[u], cv = grid[v];
                if (!adj[cu][cv]) {
                    adj[cu][cv] = adj[cv][cu] = true;
                }
            }
        }
        for (int i = 1; i < L; i++) {
            for (int j = 0; j < L; j++) {
                int u = i * L + j, v = u - L;
                int cu = grid[u], cv = grid[v];
                if (!adj[cu][cv]) {
                    adj[cu][cv] = adj[cv][cu] = true;
                }
            }
        }
    }

    void print() const {
        std::string str;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                str += format("%3d ", (int)grid[i * L + j]);
            }
            str += '\n';
        }
        str += '\n';
        for (int i = 0; i < M + 1; i++) {
            for (int j = 0; j < M + 1; j++) {
                str += adj[i][j] ? '1' : '0';
            }
            str += '\n';
        }
        str += '\n';
        std::cerr << str;
    }

    void output(std::ostream& out) const {
        std::string str;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                str += std::to_string(grid[i * L + j]) + ' ';
            }
            str += '\n';
        }
        out << str;
    }

};

void batch_execution() {

    //    int num_seeds = 50;
    //
    //    // grid_search
    //    int intervals[] = { 1,2,4,8,16 };
    //    int num_trials[] = { 1,2,4,8,16,32,64,128,256 };
    //    for (int interval : intervals) {
    //        for (int num_trial : num_trials) {
    //            std::vector<Metrics> metrics_list(num_seeds);
    //
    //            int progress = 0;
    //            int64_t score_sum = 0;
    //            int64_t min_score = INT64_MAX, max_score = INT64_MIN;
    //
    //#pragma omp parallel for num_threads(6)
    //            for (int seed = 0; seed < num_seeds; seed++) {
    //                //std::string input_file(format("../../tools_win/in/%04d.txt", seed));
    //                //std::string output_file(format("../../tools_win/out/%04d.txt", seed));
    //                //auto judge = std::make_shared<FileJudge>(input_file, output_file);
    //                auto judge = std::make_shared<LocalJudge>(seed, -1, -1, 4);
    //                Solver solver(judge);
    //                //solver.set_params(interval, num_trial);
    //                solver.set_params_opt();
    //                solver.solve();
    //#pragma omp critical(crit_sct)
    //                {
    //                    auto metrics = *solver.get_metrics();
    //                    progress++;
    //                    score_sum += metrics.score;
    //                    chmin(min_score, metrics.score);
    //                    chmax(max_score, metrics.score);
    //                    std::cerr << format("\rprogress=%3d/%3d, avg=%13.2f, min=%11lld, max=%11lld, interval=%3d, num_trial=%3d", progress, num_seeds, (double)score_sum / progress, min_score, max_score, interval, num_trial);
    //                    metrics_list[seed] = metrics;
    //                }
    //            }
    //            std::cerr << '\n';
    //            //std::cerr << format("\ninterval=%3d, num_trial=%3d, avg=%13.2f, min=%11lld, max=%11lld\n", interval, num_trial, (double)score_sum / progress, min_score, max_score);
    //            //std::cerr << '\n';
    //            //for (int seed = 0; seed < num_seeds; seed++) {
    //            //    std::cerr << format("seed=%3d, ", seed) << metrics_list[seed] << '\n';
    //            //}
    //        }
    //    }

}

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

    Timer timer;

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

#if 0
    std::ifstream ifs("../../tools_win/in/0000.txt");
    std::istream& in = ifs;
    std::ofstream ofs("../../tools_win/out/0000.txt");
    std::ostream& out = ofs;
#else
    std::istream& in = std::cin;
    std::ostream& out = std::cout;
#endif

    Xorshift rnd;

    State2 state(in);

    std::vector<std::pair<int, int>> cands;
    for (int i = 0; i < L; i++) {
        for (int j = 1; j < L; j++) {
            int u = i * L + j, v = u - 1;
            cands.emplace_back(v, u);
            cands.emplace_back(u, v);
        }
    }
    for (int i = 1; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int u = i * L + j, v = u - L;
            cands.emplace_back(v, u);
            cands.emplace_back(u, v);
        }
    }

    int best_score = -1;
    State2 best_state;

    int score = 0;
    int loop = 0;
    double start_time = timer.elapsed_ms(), now_time, end_time = 1980;
    while ((now_time = timer.elapsed_ms()) < end_time) {
        auto [u, v] = cands[rnd.next_int(cands.size())];
        if (state.grid[u] == state.grid[v]) continue;
        //if (state.grid[v] == 0 && rnd.next_int(5)) continue;
        loop++;
        int diff = int(state.grid[u] == 0) - int(state.grid[v] == 0);
        double temp = get_temp(2.0, 0.01, now_time - start_time, end_time - start_time);
        double prob = exp(diff / temp);
        if (rnd.next_double() > prob) continue;
        if (state.can_erode(u, v)) {
            state.erode(u, v);
            score += diff;
            if (chmax(best_score, score)) {
                best_state = state;
            }
        }
        if (!(loop & 0xFFFF)) dump(loop, score);
    }
    dump(loop, score, best_score);

    state = best_state;

    state.output(out);

    return 0;
}