#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <mutex>
#include <ios>

#include <omp.h>

#include <satgirgs/Generator.h>
#include <satgirgs/SpatialTree.h>

namespace satgirgs {

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

void saveDot(
        const std::vector<double> &clauseWeights, const std::vector<std::vector<double>> &clausePositions,
        const std::vector<double> &variableWeights, const std::vector<std::vector<double>> &variablePositions,
        const std::vector<std::vector<int>> &clauses, const std::string &file) {

    auto n = variableWeights.size();
    auto m = clauseWeights.size();
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
