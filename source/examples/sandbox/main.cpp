
#include <iostream>
#include <cassert>
#include <limits>
#include <fstream>
#include <future>

#include <minisat/simp/SimpSolver.h>

#include <satgirgs/SpatialTree.h>
#include <satgirgs/Generator.h>

using namespace std;
using namespace satgirgs;

mt19937 gen;
bool randomBit() {
    return gen() < gen.max() / 2;
}
void seedBit(int seed) {
    gen.seed(seed);
}

struct RunData {
    int restarts;
    int conflicts;
    int decisions;
    int propagations;
    Minisat::lbool res;
    double time; // in ms
};

RunData runMinisat(const vector<vector<int>>& cvig) {
    using namespace Minisat;
    auto s = make_shared<SimpSolver>();
    seedBit(378);
    for(auto& clause : cvig) {
        vec<Lit> v(clause.size());
        for (int i = 0; i < clause.size(); ++i) {
            v[i] = mkLit(clause[i], !randomBit());
            while(clause[i] >= s->nVars()) s->newVar(); // irgs
        }
        s->addClause(v);
    }
    auto fut = std::async(launch::async, [s](){
        ScopedTimer timer;
        s->eliminate(true);
        auto res = s->solveLimited({});
        RunData run;
        run.restarts = s->starts;
        run.conflicts = s->conflicts;
        run.decisions = s->decisions;
        run.propagations = s->propagations;
        run.res = res;
        run.time = timer.elapsed();
        return run;
    });
    if(fut.wait_for(chrono::seconds(5)) == future_status::timeout) // let it work 5 seconds
        s->interrupt();
    return fut.get(); // wait for the result
}


void printDistributionDebug(const vector<vector<int>>& sat, const vector<double>& varWeights) {
    auto n = varWeights.size();

    // debug variable distribution
    auto varOccs = vector<pair<int,double>>(n);
    for (int i = 0; i < n; ++i) varOccs[i] = {0,varWeights[i]};
    for(auto& clause : sat)
        for(int var : clause)
            varOccs[var].first++;

    sort(varOccs.begin(), varOccs.end(), greater<>());
    // for(auto [deg, weight] : varOccs) cout << deg << '\t' << weight << endl;

    // print degree distribution
    auto current = 0;
    auto occ = 0;
    auto min_w = accumulate(varWeights.begin(), varWeights.end(), 0.0);
    auto max_w = 0.0;
    auto sum_w = 0.0;
    cout << "deg\tocc\tmin\tmax\tavg\n";
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
}

vector<double> concretePowerLaw(int n, double ple) {
    auto res = vector(n,0.0);
    for (int i = 1; i <= n; ++i)
        res[n-i] = pow((double(n)/i),1/(ple-1));
    return res;
}

double searchSatisfiabilityThreshold(int n, int k, int d, double T, double ple) {
    auto alpha = (T==0 ? std::numeric_limits<double>::infinity() : 1/T);

    auto edgeSeed = 1338;
    auto posSeedC = 0; // todo change
    auto posSeedV = 50; // todo change

    auto varPoitions = generatePositions(n, d, posSeedV);
    auto varWeights = concretePowerLaw(n, ple);
    vector<int> clauseLengths;

    auto low = 0.1;
    auto high = 10.0;
    while(high-low >= 0.01 || int(n*high) == int(n*low)) {
        auto mid = (low+high)/2.0;
        // sample mid 50 times
        clauseLengths.resize(mid*n, k);

        auto sats = 0;
        auto unsats = 0;
        auto maxRuntime = 0.0;
        for (int i = 0; i < 10; ++i) {
            auto clausePositions = generatePositions(n*mid, d, posSeedC+i);
            auto sat = generateCVIG_for_real(
                    clauseLengths, clausePositions, // clauses
                    varWeights, varPoitions, // variables
                    alpha, edgeSeed); // sampling parameter
            cout << "generated " << i << flush;
            auto run = runMinisat(sat);
            if(run.res == Minisat::lbool(true))
                sats++;
            if(run.res == Minisat::lbool(false))
                unsats++;
            maxRuntime = max(maxRuntime, run.time);
            cout << " solved " << i << '\r' << flush;
        }

        // debug
        cout << "                                           \r" << flush;
        cout << mid << '\t' << sats << '\t' << unsats << '\t' << maxRuntime << endl;

        if(abs(sats-unsats) <= 1) return mid;

        if(sats>unsats)
            low = mid;
        else
            high = mid;
    }

    return low;
}



int main(int argc, char* argv[]) {
    cout.precision(4);

    // params
    auto n = 1000;
    const auto d = 2;
    auto k = 3;
    auto T = atof(argv[1]);
    // auto T = 0.01;
    auto alpha = (T==0 ? std::numeric_limits<double>::infinity() : 1/T);

    auto threshold = searchSatisfiabilityThreshold(n, k, d, T, 3.0);
    auto m = static_cast<int>(n*threshold);

    auto edgeSeed = 1338;
    auto posSeedC = 0;
    auto posSeedV = 50;

    auto claLengths = vector<int>(m,k);
    auto varWeights = concretePowerLaw(n, 3.0);
    auto claPos = generatePositions(m, d, posSeedC);
    auto varPos = generatePositions(n, d, posSeedV);

    // sample stuff
    auto sat = generateCVIG_for_real(
            claLengths, claPos, // clauses
            varWeights, varPos, // variables
            alpha, edgeSeed); // sampling parameter


    saveDot(claPos, varPos, sat, "graph.dot");
    system("neato -n -Tpdf graph.dot -o graph.pdf");
    printDistributionDebug(sat, varWeights);

    // write dimacs
    /*
    seedBit(378);
    ofstream f{"graph.sat"};
    f << "p cnf " << n << ' ' << m << endl;
    for(auto& clause : sat) {
        for(auto v : clause)
            f << (randomBit() ? v+1 : -(v+1)) << ' ';
        f << "0\n";
    }
    */

    auto run = runMinisat(sat);

    cout << run.restarts << endl;
    cout << run.conflicts << endl;
    cout << run.decisions << endl;
    cout << run.propagations << endl;
    cout << toInt(run.res) << endl;
    cout << run.time << endl;

    return 0;
}
