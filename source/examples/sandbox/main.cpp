
#include <iostream>
#include <cassert>
#include <limits>

#include <satgirgs/SpatialTree.h>
#include <satgirgs/Generator.h>
#include <fstream>

using namespace std;
using namespace satgirgs;



int main(int argc, char* argv[]) {

    // params
    auto n = 30000;
    auto m = 130000;
    const auto d = 2;
    // auto alpha = std::numeric_limits<double>::infinity();
    auto alpha = 1.3;

    auto edgeSeed = 1338;
    auto posSeedC = 378;
    auto posSeedV = 2049;
    auto weiSeedV = 378;

    auto claLengths = vector<int>(m,3);
    auto varWeights = generateWeights(n, 3.0, weiSeedV);
    auto claPos = generatePositions(m, d, posSeedC);
    auto varPos = generatePositions(n, d, posSeedV);

    // sample stuff
    auto sat = generateCVIG(
            claLengths, claPos, // clauses
            varWeights, varPos, // variables
            alpha, edgeSeed); // sampling parameter

    // debug variable distribution
    auto varOccs = vector<pair<int,double>>(n);
    for (int i = 0; i < n; ++i) varOccs[i] = {0,varWeights[i]};
    for(auto& clause : sat)
        for(int var : clause)
            varOccs[var].first++;

    sort(varOccs.begin(), varOccs.end(), greater<>());
    for(auto [deg, weight] : varOccs)
        cout << deg << '\t' << weight << endl;
    cout << endl;

    // print degree distribution
    auto current = 0;
    auto occ = 0;
    auto min_w = accumulate(varWeights.begin(), varWeights.end(), 0.0);
    auto max_w = 0.0;
    auto sum_w = 0.0;
    for(int i = varOccs.size() -1; i>=0; --i) {
        auto [deg, weight] = varOccs[i];
        if(deg == current) {
            occ++;
            min_w = min(min_w, weight);
            max_w = max(max_w, weight);
            sum_w += weight;
        } else {
            cout << current << '\t' << occ << '\t' << min_w << '\t' << max_w << '\t' << sum_w / occ << endl;
            current = deg;
            occ = 1;
            sum_w = max_w = min_w = weight;
        }
    }
    cout << current << '\t' << occ << '\t' << min_w << '\t' << max_w << '\t' << sum_w / occ << endl;

    // saveDot(claPos, varPos, sat, "graph.dot");
    // system("neato -n -Tpdf graph.dot -o graph.pdf");

    // write dimacs
    ofstream f{"graph.sat"};
    f << "p cnf " << n << ' ' << m << endl;
    for(auto& clause : sat) {
        for(auto v : clause)
            f << (rand() < RAND_MAX/2 ? v+1 : -(v+1)) << ' ';
        f << "0\n";
    }

    return 0;
}
