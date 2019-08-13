
#include <iostream>
#include <cassert>
#include <limits>

#include <satgirgs/SpatialTree.h>
#include <satgirgs/Generator.h>

using namespace std;
using namespace satgirgs;


int main(int argc, char* argv[]) {

    const auto d = 2;

    for (double alpha : {4.0, numeric_limits<double>::infinity()}) {
        for (int n = 128; n <= 1u<<13u; n*=2) {

            // params
            auto m = 4*n;
            auto seed = 1337+n;
            auto posSeedC = seed+378;
            auto posSeedV = seed*seed+4;

            auto pc = generatePositions(m, d, posSeedC);
            auto wc = vector<double>(m, 3.0);

            auto pv = generatePositions(n, d, posSeedV);
            auto wv = generateWeights(n, 2.5, seed);

            auto sat = vector<vector<int>>(m);
            auto addVar = [&sat](int c, int v) { sat[c].push_back(v); };
            {
                ScopedTimer timer("   " + to_string(n)+'\t');
                makeSpatialTree<d>(wc, pc, wv, pv, alpha, addVar).generateEdges(seed);
            }
        }
    }

    return 0;
}
