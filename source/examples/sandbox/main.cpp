
#include <iostream>
#include <cassert>
#include <limits>

#include <satgirgs/SpatialTree.h>
#include <satgirgs/Generator.h>

using namespace std;
using namespace satgirgs;

double distance(const vector<double>& a, const vector<double>& b) {
    auto D = a.size();
    auto result = 0.0;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(a[d] - b[d]);
        dist = std::min(dist, 1.0-dist);
        result = std::max(result, dist);
    }
    return result;
}



int main(int argc, char* argv[]) {

    // params
    auto n = 30;
    auto m = 50;
    const auto d = 2;
    //auto alpha = std::numeric_limits<double>::infinity();
    auto alpha = 7.0;

    // seeds
    auto seed = 1337;
    auto posSeedC = 378;
    auto posSeedV = seed*seed+4;

    auto clauseLengths = vector<int>(m,3);
    auto pc = generatePositions(m, d, posSeedC);
    auto wc = vector<double>(m);
    for (int i = 0; i < m; ++i)
        wc[i] = clauseLengths[i];

    auto pv = generatePositions(n, d, posSeedV);
    auto wv = generateWeights(n, 2.5, seed);
    auto W = accumulate(wv.begin(), wv.end(), 0.0);

    auto sat = vector<vector<int>>(m);

    // increase clause weights until each clause has > desired length
    auto scaled = 1.0;
    while(true) {

        for(auto& clause : sat)
            clause.clear();
        auto addVar = [&sat](int c, int v) { sat[c].push_back(v); };
        makeSpatialTree<d>(wc, pc, wv, pv, alpha, addVar).generateEdges(seed);

        auto works = true;
        for (int i = 0; i < m; ++i)
            works = works && sat[i].size() >= clauseLengths[i];
        if(works) break;

        scaled *= 2;
        auto scaling = 2.0;
        if(alpha == numeric_limits<double>::infinity())
            scaling = pow(2.0, 1.0/alpha);
        for(auto& w : wc)
            w *= scaling;
    }

    // clean formulas
    if(alpha==numeric_limits<double>::infinity()) {
        for (int i = 0; i < m; ++i) {
            auto& clause = sat[i];
            auto& clausePos = pc[i];
            assert(clause.size() >= clauseLengths[i]);
            // can be made linear with n-th element, and partition
            sort(clause.begin(), clause.end(), [&clausePos, &pv,&wv](int var1, int var2) {
                auto var1Goodness = wv[var1]/distance(clausePos, pv[var1]);
                auto var2Goodness = wv[var2]/distance(clausePos, pv[var2]);
                return var1Goodness > var2Goodness;
            });
            clause.resize(clauseLengths[i]);
        }
    } else {
        satgirgs::default_random_engine gen(seed);
        for (int i = 0; i < m; ++i) {
            auto& clause = sat[i];
            auto& clausePos = pc[i];
            assert(clause.size() >= clauseLengths[i]);

            // prepare vars
            vector<pair<int,double>> vars; // id, selection weight
            vars.reserve(clause.size());
            for(auto v : clause) {
                auto actual = pow((wc[i]*wv[v]/W)/distance(clausePos, pv[v]), alpha);
                auto girgDid = min(actual*scaled, 1.0);
                auto goodness = actual/girgDid;
                // if actual*scaled < 1 then goodness is 1/scaled (i.e. uniform dist)
                vars.emplace_back(v, goodness);
            }

            // take some vars
            vector<int> selection;
            while(selection.size() < clauseLengths[i]) {

                auto totalGoodness = 0.0;
                for(auto [_,g] : vars) totalGoodness += g;
                auto dist = uniform_real_distribution<>(0.0, totalGoodness);
                auto choice = dist(gen);
                // find choice, remove from vars and add to selection

                double acc = 0;
                auto toSelect = vars.end();
                for(auto it=vars.begin(); it!=vars.end(); ++it) {
                    acc += it->second;
                    if(acc >= choice) {
                        toSelect = it;
                        break;
                    }
                }
                assert(toSelect != vars.end());
                selection.push_back(toSelect->first);
                vars.erase(toSelect);
            }
            // clause i is done here
            clause = selection;
        }
    }

    // debug
    auto varOccs = vector<pair<int,double>>(n);
    for (int i = 0; i < n; ++i) varOccs[i] = {0,wv[i]};
    for(auto& clause : sat)
        for(int var : clause)
            varOccs[var].first++;

    sort(varOccs.begin(), varOccs.end(), greater<>());
    for(auto [occ, weight] : varOccs)
        cout << occ << '\t' << weight << endl;

    saveDot(wc, pc, wv, pv, sat, "graph.dot");
    system("neato -n -Tpdf graph.dot -o graph.pdf");
    return 0;
}
