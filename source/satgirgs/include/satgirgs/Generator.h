#pragma once

#include <vector>
#include <string>

#include <satgirgs/satgirgs_api.h>


namespace satgirgs {

SATGIRGS_API double distance(const std::vector<double>& a, const std::vector<double>& b);

SATGIRGS_API std::vector<double> generateWeights(int n, double ple, int weightSeed);

SATGIRGS_API std::vector<std::vector<double>> generatePositions(int n, int dimension, int positionSeed);

// returns clauses where the prob that var v is in clause c is:
// min(1, (wc*wv/Wv/dist(c,v))^alpha)
SATGIRGS_API std::vector<std::vector<int>> generateClauses(
        const std::vector<double>& clauseWeights, const std::vector<std::vector<double>>& clausePositions,
        const std::vector<double>& variableWeights, const std::vector<std::vector<double>>& variablePositions,
        double alpha, int seed);

SATGIRGS_API std::vector<std::vector<int>> generateCVIG(
        const std::vector<int>& clauseLengths, const std::vector<std::vector<double>>& clausePositions,
        const std::vector<double>& variableWeights, const std::vector<std::vector<double>>& variablePositions,
        double alpha, int seed);

SATGIRGS_API void saveDot(
        const std::vector<std::vector<double>> &clausePositions,
        const std::vector<std::vector<double>> &variablePositions,
        const std::vector<std::vector<int>> &clauses, const std::string &file);

} // namespace satgirgs
