#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <climits>
#include <stdexcept>
#include <utility>
#include <fstream>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

/* =========================
   Random helper
   ========================= */
int uniform_int(int l, int r, mt19937_64& rng) {
    uniform_int_distribution<int> dist(l, r);
    return dist(rng);
}

/* =========================
   Process ONE file safely
   ========================= */
void process_file(const fs::path& input_path,
                  const fs::path& output_path)
{
    ifstream in(input_path);
    if (!in.is_open()) {
        cerr << "Skipping (cannot open): " << input_path << "\n";
        return;
    }

    ofstream out(output_path);
    if (!out.is_open()) {
        cerr << "Cannot write to: " << output_path << "\n";
        return;
    }

    mt19937_64 rng(42);

    try {
        /* ===== 1. Read N ===== */
        int N;
        if (!(in >> N) || N <= 0 || N > 100000)
            throw runtime_error("Invalid N");

        out << N << "\n";

        vector<pair<int,int>> coords(N);
        for (int i = 0; i < N; i++) {
            if (!(in >> coords[i].first >> coords[i].second))
                throw runtime_error("Invalid coordinates");
            out << coords[i].first << " " << coords[i].second << "\n";
        }

        /* ===== 2. Read clusters ===== */
        int C;
        if (!(in >> C) || C <= 0 || C > N)
            throw runtime_error("Invalid C");

        out << C << "\n";

        vector<vector<int>> clusters(C);
        vector<int> cluster_of(N, -1);

        for (int c = 0; c < C; c++) {
            int sz;
            if (!(in >> sz) || sz <= 0)
                throw runtime_error("Invalid cluster size");

            out << sz;
            clusters[c].resize(sz);

            for (int i = 0; i < sz; i++) {
                int v;
                if (!(in >> v))
                    throw runtime_error("Invalid cluster node");
                out << " " << v;

                v--; // to 0-based
                if (v < 0 || v >= N)
                    throw runtime_error("Cluster node out of range");

                clusters[c][i] = v;
                cluster_of[v] = c;
            }
            out << "\n";
        }

        /* ===== 3. Centers ===== */
        vector<int> centers(C);
        for (int c = 0; c < C; c++) {
            if (!(in >> centers[c]))
                throw runtime_error("Invalid center");
            out << centers[c] << "\n";
        }

        /* ===== 4. First -999 ===== */
        int marker;
        if (!(in >> marker) || marker != -999)
            throw runtime_error("Missing first -999");

        out << "-999\n";

        /* ===== 5. Prizes ===== */
        vector<int> prizes(N);
        for (int i = 0; i < N; i++) {
            if (!(in >> prizes[i]))
                throw runtime_error("Invalid prize");
            out << prizes[i] << "\n";
        }

        /* ===== 6. Second -999 ===== */
        if (!(in >> marker) || marker != -999)
            throw runtime_error("Missing second -999");

        out << "-999\n";

        /* ===== 7. Compute min prizes ===== */
        vector<int> min_prize(C);
        for (int c = 0; c < C; c++) {
            int sumP = 0;
            int minP = INT_MAX;

            for (int v : clusters[c]) {
                sumP += prizes[v];
                minP = min(minP, prizes[v]);
            }

            min_prize[c] = uniform_int(minP, sumP, rng);
        }

        /* ===== 8. Append min prizes ===== */
        for (int c = 0; c < C; c++) {
            out << min_prize[c] << "\n";
        }
    }
    catch (const exception& e) {
        cerr << "Skipping invalid file: "
             << input_path.filename()
             << " (" << e.what() << ")\n";
    }
}

/* =========================
   MAIN
   ========================= */
int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "Usage: " << argv[0]
             << " <input_folder> <output_folder>\n";
        return 1;
    }

    fs::path input_dir(argv[1]);
    fs::path output_dir(argv[2]);

    if (!fs::exists(input_dir) || !fs::is_directory(input_dir)) {
        cerr << "Input folder does not exist.\n";
        return 1;
    }

    if (!fs::exists(output_dir)) {
        fs::create_directories(output_dir);
    }

    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (!entry.is_regular_file())
            continue;

        // Skip hidden files (.DS_Store, etc.)
        string name = entry.path().filename().string();
        if (!name.empty() && name[0] == '.')
            continue;

        fs::path out_path = output_dir / entry.path().filename();
    

        cout << "Processing: " << entry.path().filename() << "\n";
        process_file(entry.path(), out_path);
    }

    cout << "Done.\n";
    return 0;
}
