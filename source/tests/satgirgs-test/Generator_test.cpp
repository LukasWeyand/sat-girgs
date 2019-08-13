
#include <algorithm>
#include <cmath>
#include <numeric>

#include <gmock/gmock.h>

#include <satgirgs/Generator.h>

using namespace std;

// FWD for distance function. Declared in main.
double distance(const std::vector<double>& a, const std::vector<double>& b);


class Generator_test: public testing::Test
{
protected:
    int seed = 1337;
};


bool connected(int c, int v, const vector<vector<int>>& graph) {
    return find(graph[c].begin(), graph[c].end(), v) != graph[c].end();
}


TEST_F(Generator_test, testThresholdModel)
{
    const auto n = 100;
    const auto m = 100;
    const auto alpha = numeric_limits<double>::infinity();
    const auto ple = 2.8;

    // weights
    auto variabW = satgirgs::generateWeights(n, ple, seed);
    auto clauseW = vector<double>(m, 1.0);
    auto W = accumulate(variabW.begin(), variabW.end(), 0.0);

    for(auto d=1u; d<5; ++d){

        auto variabP = satgirgs::generatePositions(n, d, seed+d);
        auto clauseP = satgirgs::generatePositions(m, d, 2*seed+d);
        auto clauses = satgirgs::generateClauses(clauseW, clauseP, variabW, variabP, alpha, 0);

        // dist^d < wv/W
        for(int v=0; v<n; ++v){
            for(int c=0; c<m; ++c){

                const auto dist = distance(clauseP[c], variabP[v]);
                const auto d_term = pow(dist, d);
                const auto w_term = clauseW[c] * variabW[v] / W;

                if(d_term < w_term) {
                    EXPECT_TRUE(connected(c,v, clauses)) << "edge should be present";
                } else {
					EXPECT_FALSE(connected(c,v, clauses)) << "edge should be absent";
                }
            }
        }
    }
}

TEST_F(Generator_test, testGeneralModel)
{
    const auto n = 300;
    const auto m = 600;
    const auto alpha = 2.5;
    const auto ple = 2.5;

    auto variabW = satgirgs::generateWeights(n, ple, seed);
    auto clauseW = vector<double>(m, 1.0);

    for(auto d=1u; d<5; ++d){

        // generate and hope for no failing assertions
        auto variabP = satgirgs::generatePositions(n, d, seed+d);
        auto clauseP = satgirgs::generatePositions(m, d, 2*seed+d);
        auto clauses = satgirgs::generateClauses(clauseW, clauseP, variabW, variabP, alpha, 3*seed+d);
    }

    // if no assertion failed we are fine here :)
}


TEST_F(Generator_test, testCompleteGraph)
{
    const auto n = 100;
    const auto m = 200;
    const auto alpha = 0.0; // each edge prob will be 100% now
    const auto ple = 2.5;

    auto variabW = satgirgs::generateWeights(n, ple, seed);
    auto clauseW = vector<double>(m, 1.0);

    for(auto d=1u; d<5; ++d) {

        auto variabP = satgirgs::generatePositions(n, d, seed+d);
        auto clauseP = satgirgs::generatePositions(m, d, 2*seed+d);
        auto clauses = satgirgs::generateClauses(clauseW, clauseP, variabW, variabP, alpha, 3*seed+d);

        auto edges = accumulate(clauses.begin(), clauses.end(), 0ll, [](auto acc, auto& c){return acc + c.size();});

		// check for the correct number of edges
		EXPECT_EQ(edges, n*m) << "expect a complete graph without self loops";

        // check that each node is connected to all other nodes
        for (int v = 0; v < n; ++v) {
            for (int c = 0; c < m; ++c) {
                EXPECT_TRUE(connected(c,v,clauses));
            }
        }
    }
}


TEST_F(Generator_test, testWeightSampling)
{
    auto n = 10000;
    auto ple = 2.1;
    int runs = 10;

    for(int i=0; i<runs; ++i){

        auto weights = satgirgs::generateWeights(n, ple, seed+i);
        for(auto each : weights) {
            EXPECT_GE(each, 1.0);
            EXPECT_LT(each, n);
        }
        auto max_weight = *max_element(weights.begin(), weights.end());
        EXPECT_GT(max_weight * max_weight, n) << "max weight should be large";
    }
}


TEST_F(Generator_test, testReproducible)
{
    auto n = 1000;
    auto m = 1000;
    auto ple = 2.4;
    auto weight_seed    = 1337;
    auto position_seed  = 42;

    auto alphas = { 1.5, std::numeric_limits<double>::infinity() };
    auto dimensions = { 1, 2 };

    for (auto alpha : alphas) {
        for (auto d : dimensions) {

            auto clauses1 = satgirgs::generateClauses(
                    vector<double>(m,1.0),
                    satgirgs::generatePositions(m, d, position_seed),
                    satgirgs::generateWeights(m, ple, weight_seed),
                    satgirgs::generatePositions(n, d, position_seed*2+1),
                    alpha, weight_seed+position_seed+d);
            auto clauses2 = satgirgs::generateClauses(
                    vector<double>(m,1.0),
                    satgirgs::generatePositions(m, d, position_seed),
                    satgirgs::generateWeights(m, ple, weight_seed),
                    satgirgs::generatePositions(n, d, position_seed*2+1),
                    alpha, weight_seed+position_seed+d);

            // same edges
            EXPECT_EQ(clauses1, clauses2);
        }
    }
}
