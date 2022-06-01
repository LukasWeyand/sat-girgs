
#include <iostream>
#include <limits>
#include <fstream>

#include <satgirgs/SpatialTree.h>
#include <satgirgs/Generator.h>

using namespace std;
using namespace satgirgs;


std::vector<std::vector<int>> generate(
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
        auto res = generate(scaled_clsW, clsP, varW, varP, alpha, sseed);
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
    return generate(clsW, clsP, varW, varP, alpha, sseed);
}


int main(int argc, char* argv[]) {

    int n = 100;
    int m = 1000;
    int deg = 5;
    int seed = 12;
    double T = 0;
    bool univar = false;
    bool unicls = true;
    if(argc==1) cerr << "USAGE: n m deg seed T univar unicls" << endl, exit(0);

    if(argc>1) n = stoi(argv[1]);
    if(argc>2) m = stoi(argv[2]);
    if(argc>3) deg = stoi(argv[3]);
    if(argc>4) seed = stoi(argv[4]);
    if(argc>5) T = stod(argv[5]);
    if(argc>6) univar = argv[6] != "0"s;
    if(argc>7) unicls = argv[7] != "0"s;

    auto graph = generate(n, m, deg, seed, T, univar, unicls);

    int zero_cnt = 0;
    for(auto& cls : graph) zero_cnt += empty(cls);
    if(zero_cnt) cerr << "WARNING: " << zero_cnt << " empty clauses" << endl;

    cout << n << ' ' << m-zero_cnt << endl;
    for(auto& cls : graph) {
        if(empty(cls)) continue;
        cout << size(cls);
        for(auto var : cls) cout << ' ' << var;
        cout << endl;
    }

    return 0;
}
