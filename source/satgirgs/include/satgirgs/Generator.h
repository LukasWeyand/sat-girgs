#pragma once

#include <vector>
#include <string>

#include <satgirgs/satgirgs_api.h>


namespace satgirgs {

SATGIRGS_API std::vector<double> generateWeights(int n, double ple, int weightSeed);

SATGIRGS_API std::vector<std::vector<double>> generatePositions(int n, int dimension, int positionSeed);

SATGIRGS_API void saveDot(
        const std::vector<double> &clauseWeights, const std::vector<std::vector<double>> &clausePositions,
        const std::vector<double> &variableWeights, const std::vector<std::vector<double>> &variablePositions,
        const std::vector<std::vector<int>> &clauses, const std::string &file);

} // namespace satgirgs
