
#include <limits>
#include <random>
#include <algorithm>
#include <map>
#include <fstream>

#include <satgirgs/Generator.h>

using namespace std;
using namespace satgirgs;

void makePDF(
        const std::vector<std::vector<double>> &clausePositions,
        const std::vector<std::vector<double>> &variablePositions,
        const std::vector<std::vector<int>> &clauses, const std::string &file) {

    auto n = variablePositions.size();
    auto m = clausePositions.size();
    auto D = variablePositions.front().size();
    std::ofstream f{file + ".dot"};
    if(!f.is_open())
        throw std::runtime_error{"Error: failed to open file \"" + file + '\"'};

    f << "graph girg {\n\tscale=2000\n\n";
    f << "\tviewport=\"2000,2000\"\n\n";
    f << std::fixed;

    // duplicate nodes 4 times
    auto duplicate4 = [](const vector<vector<double>>& pos) {
        vector<vector<double>> res;
        for(auto& coord : pos) {
            for(int i=0;i<4;++i) {
                res.push_back(vector{coord[0] + (i&1), coord[1] + ((i>>1)&1)});
                // TODO possibly copy color, label and other meta info
            }
        }
        return res;
    };

    auto clsPosTorus = duplicate4(clausePositions);
    auto varPosTorus = duplicate4(variablePositions);

    // expand edges and filter to dist<0.5
    auto dist = [](auto a, auto b) { return max(abs(a[0]-b[0]), abs(a[1]-b[1])); };
    vector<vector<int>> clsTorus(4*m);
    for(int c=0; c<4*m; ++c)
        for(auto v : clauses[c/4])
            for(int i=0; i<4; ++i)
                if(dist(clsPosTorus[c], varPosTorus[v*4+i]) < 0.5)
                    clsTorus[c].push_back(v*4+i);

    n *= 4;
    m *= 4;

    // write variables
    for (int i = 0; i < n; ++i) {
        f << '\t' << i << " [shape=circle pos=\"";
        for (auto d = 0u; d < D; ++d)
            f << (d == 0 ? "" : ",") << varPosTorus[i][d];
        f << '\"';
        //f << " label=\"" << i/4 << "v\"";
        f << " label=\"\"";
        f << " style =\"filled\" fillcolor=\"lightblue\"";
        f << "];\n";
    }
    // write clauses
    for (int i = 0; i < m; ++i) {
        f << '\t' << n+i << " [shape=point, pos=\"";
        for (auto d = 0u; d < D; ++d)
            f << (d == 0 ? "" : ",") << clsPosTorus[i][d];
        f << '\"';
        //f << " label=\"" << i/4 << "c\"";
        f << " label=\"\"";
        f << " style =\"filled\" fillcolor=\"black\"";
        f << "];\n";
    }
    f << '\n';
    // write edges
    for (int i = 0; i < m; ++i) {
        if(clsTorus[i].empty()) continue;
        f << '\t' << n+i << "\t-- {";
        for(auto var : clsTorus[i])
            f << var << ' ';
        f << "};\n";
    }
    f << "}\n";

    f.close();

    string command = "neato -n -Tpdf " + file + ".dot -o " + file + ".pdf";
    system(command.c_str());
}


std::vector<std::vector<int>> generate(int n, int m, int deg, int seed, double T, bool univar, bool unicls) {
    mt19937 gen(seed);
    uniform_int_distribution dist(0,numeric_limits<int>::max());
    int wseed1 = dist(gen);
    int wseed2 = dist(gen);
    int pseed1 = dist(gen);
    int pseed2 = dist(gen);
    int sseed = dist(gen);

    double ple = 2.8;
    int d = 2;
    double alpha = T>0 ? 1.0/T : numeric_limits<double>::infinity();

    auto varW = univar ? vector(n,1.0) : generateWeights(n,ple,wseed1);
    auto varP = generatePositions(n,d,pseed1);

    auto clsW = unicls ? vector(m,1.0) : generateWeights(m,ple,wseed2);
    auto clsP = generatePositions(m,d,pseed2);

    auto sum_deg = [&](double scaling) {
        auto scaled_clsW = clsW;
        for(auto& e : scaled_clsW) e *= scaling;
        auto res = generateClauses(scaled_clsW, clsP, varW, varP, alpha, sseed);
        long long avg = 0;
        for(auto& cls : res) avg += cls.size();
        //cerr << scaling << " :" << 1.0*avg/m << endl;
        return avg;
    };

    double low = 1.0, high = low;
    while(sum_deg(low) > deg*m) low /= 2;
    while(sum_deg(high)< deg*m) high *= 2;
    while(high-low > 0.01) {
        auto mid = (high+low)/2;
        if(sum_deg(mid) <= deg*m) low = mid;
        else high = mid;
    }


    for(auto& e : clsW) e *= low;
    auto res = generateClauses(clsW, clsP, varW, varP, alpha, sseed);

    long long avg = 0;
    for(auto& cls : res) avg += cls.size();
    printf("%f", 1.0 * avg / m);

    makePDF(clsP,varP,res, "graph_torus");
    saveDot(clsP,varP,res, "graph.dot");
    string command = "neato -n -Tpdf graph.dot -o graph.pdf";
    system(command.c_str());

    return res;
}



int main(int argc, char* argv[]) {

    int n = 200;
    int m = 200;
    int deg = 5;
    int seed = 12;
    double T = 0.5;
    bool univar = true;
    bool unicls = true;

    auto graph = generate(n, m, deg, seed, T, univar, unicls);

    return 0;
}
