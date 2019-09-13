
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <numeric>
#include <limits>

#include <satgirgs/Generator.h>
#include <satgirgs/SpatialTree.h>

namespace satgirgs {

double distance(const std::vector<double>& a, const std::vector<double>& b) {
    auto D = a.size();
    auto result = 0.0;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(a[d] - b[d]);
        dist = std::min(dist, 1.0-dist);
        result = std::max(result, dist);
    }
    return result;
}

std::vector<double> generateWeights(int n, double ple, int weightSeed) {
    auto result = std::vector<double>(n);
    auto gen = default_random_engine{weightSeed >= 0 ? weightSeed : std::random_device()()};
    auto dist = std::uniform_real_distribution<>{};
    for (int i = 0; i < n; ++i)
        result[i] = std::pow((std::pow(n, -ple + 1) - 1) * dist(gen) + 1, 1 / (-ple + 1));
    return result;
}

std::vector<std::vector<double>> generatePositions(int n, int dimension, int positionSeed) {
    auto result = std::vector<std::vector<double>>(n, std::vector<double>(dimension));
    auto gen = default_random_engine{positionSeed >= 0 ? positionSeed : std::random_device()()};
    auto dist = std::uniform_real_distribution<>{};
    for(int i=0; i<n; ++i)
        for (int d=0; d<dimension; ++d)
            result[i][d] = dist(gen);
    return result;
}

std::vector<std::vector<int>> generateClauses(
        const std::vector<double> &clauseWeights, const std::vector<std::vector<double>> &clausePositions,
        const std::vector<double> &variableWeights, const std::vector<std::vector<double>> &variablePositions,
        double alpha, int seed) {
    auto n = variableWeights.size();
    auto m = clauseWeights.size();
    auto d = clausePositions.front().size();
    auto clauses = std::vector<std::vector<int>>(m);
    auto addVarToClause = [&clauses](int c, int v) {
        clauses[c].push_back(v);
    };

    assert(variablePositions.size()==n);
    assert(clausePositions.size()==m);
    assert(variablePositions.front().size() == d);

    if(d==1) makeSpatialTree<1>(clauseWeights, clausePositions, variableWeights, variablePositions, alpha, addVarToClause).generateEdges(seed);
    if(d==2) makeSpatialTree<2>(clauseWeights, clausePositions, variableWeights, variablePositions, alpha, addVarToClause).generateEdges(seed);
    if(d==3) makeSpatialTree<3>(clauseWeights, clausePositions, variableWeights, variablePositions, alpha, addVarToClause).generateEdges(seed);
    if(d==4) makeSpatialTree<4>(clauseWeights, clausePositions, variableWeights, variablePositions, alpha, addVarToClause).generateEdges(seed);
    if(d==5) makeSpatialTree<5>(clauseWeights, clausePositions, variableWeights, variablePositions, alpha, addVarToClause).generateEdges(seed);

    return clauses;
}

std::vector<std::vector<int>> generateCVIG(
        const std::vector<int> &clauseLengths, const std::vector<std::vector<double>> &clausePositions,
        const std::vector<double> &variableWeights, const std::vector<std::vector<double>> &variablePositions,
        double alpha, int seed) {

    auto m = clausePositions.size();
    auto W = accumulate(variableWeights.begin(), variableWeights.end(), 0.0);
    auto d = clausePositions.front().size();

    std::vector<std::vector<int>> sat;

    // increase clause weights until each clause has > desired length
    auto scaled = 1.0;
    auto clauseWeights = std::vector<double>(m);
    for (int i = 0; i < m; ++i) clauseWeights[i] = clauseLengths[i] / 3.0; // this 3.0 is a magic tuning parameter

    while(true) {

        sat = generateClauses(clauseWeights, clausePositions, variableWeights, variablePositions, alpha, seed);

        auto works = true;
        for (int i = 0; i < m; ++i)
            works = works && sat[i].size() >= clauseLengths[i];
        if(works) break;

        scaled *= 2;
        auto scaling = 2.0; // correct scaling such that all edge probs actually double
        // this correctly doubles all probs, but do we want this?
        // for near-threshold it does nothing. we actually want to double expected avg degrees
        //if(alpha != std::numeric_limits<double>::infinity())
        //    scaling = pow(2.0, 1.0/alpha);
        for(auto& w : clauseWeights)
            w *= scaling;
    }

    // trim clauses to desired length
    if(alpha==std::numeric_limits<double>::infinity()) {
        for (int i = 0; i < m; ++i) {
            auto &clause = sat[i];
            auto &clausePos = clausePositions[i];
            assert(clause.size() >= clauseLengths[i]);
            // can be made linear with n-th element, and partition
            sort(clause.begin(), clause.end(), [&](int var1, int var2) {
                auto var1Goodness = variableWeights[var1] / pow(distance(clausePos, variablePositions[var1]), d);
                auto var2Goodness = variableWeights[var2] / pow(distance(clausePos, variablePositions[var2]), d);
                return var1Goodness > var2Goodness;
            });
            clause.resize(clauseLengths[i]);
        }
        return sat;
    }

    // binomial model trimming is the tricky part
    satgirgs::default_random_engine gen(seed);
    for (int i = 0; i < m; ++i) {
        auto& clause = sat[i];
        auto& clausePos = clausePositions[i];
        assert(clause.size() >= clauseLengths[i]);

        // prepare vars
        std::vector<std::pair<int,double>> vars; // id, selection weight
        vars.reserve(clause.size());
        for(auto v : clause) {
            auto actual = pow(
                    (clauseWeights[i]*variableWeights[v]/W) / pow(distance(clausePos, variablePositions[v]), d),
                    alpha);
            assert(actual != 0.0);
            auto girgDid = std::min(actual*scaled, 1.0);
            auto goodness = actual/girgDid;
            if(goodness == std::numeric_limits<double>::infinity()) {
                std::cout << "Accuracy insufficient to sample this close to threshold model\n"
                    << "Use T=0 instead\n";
                exit(0);
            }

            // if actual*scaled < 1 then goodness is 1/scaled (i.e. uniform dist)
            vars.emplace_back(v, goodness);
        }

        // take some vars
        std::vector<int> selection;
        while(selection.size() < clauseLengths[i]) {

            auto totalGoodness = 0.0;
            for(auto [_,g] : vars) totalGoodness += g;
            auto dist = std::uniform_real_distribution<>(0.0, totalGoodness);
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
    return sat;
}

std::vector<std::vector<int>> generateCVIG_for_real(
        const std::vector<int> &clauseLengths, const std::vector<std::vector<double>> &clausePositions,
        const std::vector<double> &variableWeights, const std::vector<std::vector<double>> &variablePositions,
        double alpha, int seed) {
    using namespace std;

    auto n = variableWeights.size();
    auto m = clauseLengths.size();
    auto d = clausePositions.front().size();
    auto W = accumulate(variableWeights.begin(), variableWeights.end(), 0.0);

    auto cvig = std::vector(m, std::vector<int>());

    mt19937 gen(seed);
    for (int c = 0; c < m; ++c) {
        auto length = clauseLengths[c];
        vector<int>& clause = cvig[c];

        // threshold
        if(alpha > 20.0) {
            vector<pair<int,double>> dists;
            for (int v = 0; v < n; ++v) {
                auto w_term = variableWeights[v] / W;
                auto d_term = pow(distance(clausePositions[c], variablePositions[v]), d);
                dists.emplace_back(v, w_term / d_term);
            }
            nth_element(dists.begin(), dists.begin()+(length-1), dists.end() , [](auto a, auto b) {
                return a.second > b.second;
            });
            for (int i = 0; i < length; ++i) {
                clause.push_back(dists[i].first);
            }
            continue;
        }

        // binomial
        vector<pair<int,double>> choices; // (var, weight)
        for (int v = 0; v < n; ++v) {
            auto w_term = variableWeights[v] / W;
            auto d_term = pow(distance(clausePositions[c], variablePositions[v]), d);
            auto goodness = pow(w_term/d_term, alpha);
            choices.emplace_back(v, goodness);
        }

        // repeat size of clause times
        for (int rep = 0; rep < length; ++rep) {

            // compute total weight
            auto total = accumulate(choices.begin(), choices.end(), 0.0, [](double acc, pair<int,double> choice){
                return acc + choice.second;
            });
            if(total == numeric_limits<double>::infinity()) {
                cout << "Double Overflow ... ABORT" << endl;
                return cvig;
            }

            // throw coin
            auto dist = std::uniform_real_distribution<>(0.0, total);
            auto choice = dist(gen);

            // find choice
            auto acc = 0.0;
            auto toSelect = choices.end();
            for(auto it=choices.begin(); it!=choices.end(); ++it) {
                acc += it->second;
                if(acc >= choice) {
                    toSelect = it;
                    break;
                }
            }

            // remove from list and add to clause
            assert(toSelect != choices.end());
            clause.push_back(toSelect->first);
            choices.erase(toSelect);
        }
    }

    return cvig;
}


void saveDot(
        const std::vector<std::vector<double>> &clausePositions,
        const std::vector<std::vector<double>> &variablePositions,
        const std::vector<std::vector<int>> &clauses, const std::string &file) {

    auto n = variablePositions.size();
    auto m = clausePositions.size();
    auto D = variablePositions.front().size();
    std::ofstream f{file};
    if(!f.is_open())
        throw std::runtime_error{"Error: failed to open file \"" + file + '\"'};

    f << "digraph sat {\n\toverlap=scale;\n\n" << std::fixed;
    // write variables
    for (int i = 0; i < n; ++i) {
        f << '\t' << i << " [pos=\"";
        for (auto d = 0u; d < D; ++d)
            f << (d == 0 ? "" : ",") << variablePositions[i][d];
        f << "\"];\n";
    }
    // write clauses
    for (int i = 0; i < m; ++i) {
        f << '\t' << n+i << " [shape=box, pos=\"";
        for (auto d = 0u; d < D; ++d)
            f << (d == 0 ? "" : ",") << clausePositions[i][d];
        f << "\"];\n";
    }
    f << '\n';
    // write edges
    for (int i = 0; i < m; ++i) {
        if(clauses[i].empty()) continue;
        f << '\t' << n+i << "\t-> {";
        for(auto var : clauses[i])
            f << var << ' ';
        f << "};\n";
    }
    f << "}\n";
}

} // namespace satgirgs
