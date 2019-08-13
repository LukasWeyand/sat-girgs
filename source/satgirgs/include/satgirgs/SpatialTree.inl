
#include <satgirgs/ScopedTimer.h>

namespace satgirgs {

template<size_t exp>
double pow_to_the(double x) {
    double res = 1.0;
    for (int i = 0; i < exp; ++i)
        res *= x;
    return res;
}


template<unsigned int D, typename EdgeCallback>
SpatialTree<D, EdgeCallback>::SpatialTree(
        const std::vector<double>& weightsC, const std::vector<std::vector<double>>& positionsC,
        const std::vector<double>& weightsV, const std::vector<std::vector<double>>& positionsV,
        double alpha, EdgeCallback& edgeCallback, bool profile)
: m_EdgeCallback(edgeCallback)
, m_profile(profile)
, m_alpha(alpha)
{
    ScopedTimer timer("Preprocessing", profile);
    assert(weightsC.size() == positionsC.size());
    assert(weightsV.size() == positionsV.size());
    assert(!positionsC.empty() && positionsC.front().size() == D);
    assert(!positionsV.empty() && positionsV.front().size() == D);

    m_m = weightsC.size();
    m_n = weightsV.size();
    auto w0C = *std::min_element(weightsC.begin(), weightsC.end());
    auto wnC = *std::max_element(weightsC.begin(), weightsC.end());
    auto WC = std::accumulate(weightsC.begin(), weightsC.end(), 0.0);
    auto w0V = *std::min_element(weightsV.begin(), weightsV.end());
    auto wnV = *std::max_element(weightsV.begin(), weightsV.end());
    auto WV = std::accumulate(weightsV.begin(), weightsV.end(), 0.0);
    m_w0 = std::min(w0C, w0V);
    m_wn = std::max(wnC, wnV);
    m_W = WV;

    m_baseLevelConstant = static_cast<int>(std::log2(m_W/m_w0/m_w0)); // log2(W/w0^2)
    m_layers = static_cast<unsigned int>(floor(std::log2(m_wn/m_w0))) + 1;
    m_levels = partitioningBaseLevel(0,0) + 1; // (log2(W/w0^2) - 2) / d

    // determine which layer pairs to sample in which level
    m_layer_pairs.resize(m_levels);
    for (auto i = 0u; i < m_layers; ++i)
        for (auto j = 0u; j < m_layers; ++j)
            m_layer_pairs[partitioningBaseLevel(i, j)].emplace_back(i,j);

    // sort clause and variable weights into exponentially growing layers
    {
        ScopedTimer timer("Build DS", profile);
        m_weight_layersC = buildPartition(weightsC, positionsC, true);
        m_weight_layersV = buildPartition(weightsV, positionsV, false);
    }
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::generateEdges(int seed) {

    // one random generator
    m_gen.seed(seed >= 0 ? seed : std::random_device()());

#ifndef NDEBUG
    // ensure that all node pairs are compared either type 1 or type 2
    m_type1_checks = 0;
    m_type2_checks = 0;
#endif // NDEBUG

    // sample all edges
    visitCellPair(0, 0, 0);
    assert(m_type1_checks + m_type2_checks == m_m*m_n);
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level) {

    if(!CoordinateHelper::touching(cellA, cellB, level)) { // not touching
        // sample all type 2 occurrences with this cell pair
        #ifdef NDEBUG
		if (m_alpha == std::numeric_limits<double>::infinity()) return; // dont trust compilter optimization
        #endif // NDEBUG
        for(auto l=level; l<m_levels; ++l)
            for(auto& [i,j] : m_layer_pairs[l]) {
                sampleTypeII(cellA, cellB, level, i, j);
                sampleTypeII(cellB, cellA, level, i, j);
            }
        return;
    }

    // sample all type 1 occurrences with this cell pair
    for(auto& [i,j] : m_layer_pairs[level]) {
        sampleTypeI(cellA, cellB, level, i, j);
        if(cellA != cellB)
            sampleTypeI(cellB, cellA, level, i, j);
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    for(auto a = CoordinateHelper::firstChild(cellA); a<=CoordinateHelper::lastChild(cellA); ++a)
        for(auto b = cellA == cellB ? a : CoordinateHelper::firstChild(cellB); b<=CoordinateHelper::lastChild(cellB); ++b)
            visitCellPair(a, b, level+1);
}




template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::sampleTypeI(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{
    assert(partitioningBaseLevel(i, j) == level);

    // in this method we sample CiA x VjB
    auto [beginA, endA] = m_weight_layersC[i].cellIterators(cellA, level); // clauses in cell A and layer i
    auto [beginB, endB] = m_weight_layersV[j].cellIterators(cellB, level); // variables in cell B and layer j

    if (beginA == endA || beginB == endB)
        return;

#ifndef NDEBUG
    {
        const auto sizeV_i_A = std::distance(beginA, endA);
        const auto sizeV_j_B = std::distance(beginB, endB);
        m_type1_checks += sizeV_i_A * sizeV_j_B;
    }
#endif // NDEBUG

    std::uniform_real_distribution<> dist;
    const auto inThresholdMode = m_alpha == std::numeric_limits<double>::infinity();

    for(auto pointerA = beginA; pointerA != endA; ++pointerA) {
        for (auto pointerB = beginB; pointerB != endB; ++pointerB) {

            const Node<D>& clause = *pointerA;
            const auto& variable = *pointerB;

            // pointer magic gives same results
            assert(clause.index == m_weight_layersC[i].kthPoint(cellA, level, std::distance(beginA, pointerA)).index);
            assert(variable.index == m_weight_layersV[j].kthPoint(cellB, level, std::distance(beginB, pointerB)).index);

            // points are in correct cells
            assert(cellA - CoordinateHelper::firstCellOfLevel(level) == CoordinateHelper::cellForPoint(clause.coord, level));
            assert(cellB - CoordinateHelper::firstCellOfLevel(level) == CoordinateHelper::cellForPoint(variable.coord, level));

            // points are in correct weight layer
            // TODO change later
            assert(i == static_cast<unsigned int>(std::log2(clause.weight/m_w0)));
            assert(j == static_cast<unsigned int>(std::log2(variable.weight/m_w0)));

            const auto distance = clause.distance(variable);
            const auto w_term = clause.weight*variable.weight/m_W;
            const auto d_term = pow_to_the<D>(distance);

            if(inThresholdMode) {
                if(d_term < w_term)
                    m_EdgeCallback(clause.index, variable.index);
            } else {
                auto edge_prob = std::pow(w_term/d_term, m_alpha); // we don't need min with 1.0 here
                if(dist(m_gen) < edge_prob)
                    m_EdgeCallback(clause.index, variable.index);
            }
        }
    }
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::sampleTypeII(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{
    assert(partitioningBaseLevel(i, j) >= level);

    // sample CiA x VjB
    auto [beginA, endA] = m_weight_layersC[i].cellIterators(cellA, level); // clauses in cell A and layer i
    auto [beginB, endB] = m_weight_layersV[j].cellIterators(cellB, level); // variables in cell B and layer j

    if (beginA == endA || beginB == endB)
        return;

    const auto sizeV_i_A = std::distance(beginA, endA);
    const auto sizeV_j_B = std::distance(beginB, endB);
#ifndef NDEBUG
    m_type2_checks += sizeV_i_A * sizeV_j_B;
#endif // NDEBUG

    // get upper bound for probability
    const auto w_upper_bound = m_w0*(1<<(i+1)) * m_w0*(1<<(j+1)) / m_W;
    const auto cell_distance = CoordinateHelper::dist(cellA, cellB, level);
    const auto dist_lower_bound = pow_to_the<D>(cell_distance);
    const auto max_connection_prob = std::min(std::pow(w_upper_bound/dist_lower_bound, m_alpha), 1.0);
    assert(dist_lower_bound > w_upper_bound); // in threshold model we would not sample anything
    const auto num_pairs = sizeV_i_A * sizeV_j_B;
    const auto expected_samples = num_pairs * max_connection_prob;

    if(expected_samples < 1e-6)
        return;

    // init geometric distribution
    auto geo = std::geometric_distribution<unsigned long long>(max_connection_prob);
    auto dist = std::uniform_real_distribution<>(0, max_connection_prob);

    for (auto r = geo(m_gen); r < num_pairs; r += 1 + geo(m_gen)) {
        // determine the r-th pair
        const Node<D>& clause = beginA[r%sizeV_i_A];
        const Node<D>& variable = beginB[r/sizeV_i_A];

        // points are in correct weight layer
        assert(i == static_cast<unsigned int>(std::log2(clause.weight/m_w0)));
        assert(j == static_cast<unsigned int>(std::log2(variable.weight/m_w0)));

        // points are in correct cells
        assert(cellA - CoordinateHelper::firstCellOfLevel(level) == CoordinateHelper::cellForPoint(clause.coord, level));
        assert(cellB - CoordinateHelper::firstCellOfLevel(level) == CoordinateHelper::cellForPoint(variable.coord, level));

        const auto rnd = dist(m_gen);

        // get actual connection probability
        const auto distance = clause.distance(variable);
        const auto w_term = clause.weight*variable.weight/m_W;
        const auto d_term = pow_to_the<D>(distance);
        const auto connection_prob = std::pow(w_term/d_term, m_alpha); // we don't need min with 1.0 here
        assert(w_term < w_upper_bound);
        assert(d_term >= dist_lower_bound);

        if(rnd < connection_prob) {
            m_EdgeCallback(clause.index, variable.index);
        }
    }
}


template<unsigned int D, typename EdgeCallback>
unsigned int SpatialTree<D, EdgeCallback>::weightLayerTargetLevel(int layer) const {
    // -1 coz w0 is the upper bound for layer 0 in paper and our layers are shifted by -1
    auto result = std::max((m_baseLevelConstant - layer - 1) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct insertion level
        assert(0 <= layer && layer < m_layers);
        assert(0 <= result && result <= m_levels); // note the result may be one larger than the deepest level (hence the <= m_levels)
        auto volume_requested  = m_w0*m_w0*std::pow(2,layer+1)/m_W; // v(i) = w0*wi/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current || volume_requested >= 1.0); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}


template<unsigned int D, typename EdgeCallback>
unsigned int SpatialTree<D, EdgeCallback>::partitioningBaseLevel(int layer1, int layer2) const {

    // we do the computation on signed ints but cast back after the max with 0
    // m_baseLevelConstant is just log(W/w0^2)
    auto result = std::max((m_baseLevelConstant - layer1 - layer2 - 2) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct comparison level
        assert(0 <= layer1 && layer1 < m_layers);
        assert(0 <= layer2 && layer2 < m_layers);
        auto volume_requested  = m_w0*std::pow(2,layer1+1) * m_w0*std::pow(2,layer2+1) / m_W; // v(i,j) = wi*wj/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current || volume_requested >= 1.0); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}

template<unsigned int D, typename EdgeCallback>
std::vector<WeightLayer<D>> SpatialTree<D, EdgeCallback>::buildPartition(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions, bool clauses) {

    auto& m_nodes = (clauses ? m_clauses : m_variables);
    auto& m_first_in_cell = (clauses ? m_firstC_in_cell : m_firstV_in_cell);

    const auto n = weights.size();
    assert(positions.size() == n);

    auto weight_to_layer = [=] (double weight) {
        return std::log2(weight / m_w0);
    };

    const auto first_cell_of_layer = [&] {
        std::vector<unsigned int> first_cell_of_layer(m_layers + 1);
        unsigned int sum = 0;
        for (auto l = 0; l < m_layers; ++l) {
            first_cell_of_layer[l] = sum;
            sum += CoordinateHelper::numCellsInLevel(weightLayerTargetLevel(l));
        }
        first_cell_of_layer.back() = sum;
        return first_cell_of_layer;
    }();
    const auto max_cell_id = first_cell_of_layer.back();

    // Node<D> should incur no init overhead; checked on godbolt
    m_nodes = std::vector<Node<D>>(n); 
    // compute the cell a point belongs to
    {
        ScopedTimer timer("Classify points & precompute coordinates", m_profile);

        for (int i = 0; i < n; ++i) {
            const auto layer = weight_to_layer(weights[i]);
            const auto level = weightLayerTargetLevel(layer);
            m_nodes[i] = Node<D>(positions[i], weights[i], i);
            m_nodes[i].cell_id = first_cell_of_layer[layer] + CoordinateHelper::cellForPoint(m_nodes[i].coord, level);
            assert(m_nodes[i].cell_id < max_cell_id);
        }
    }

    // Sort points by cell-ids
    {
        ScopedTimer timer("Sort points", m_profile);

        auto compare = [](const Node<D> &a, const Node<D> &b) { return a.cell_id < b.cell_id; };

        std::sort(m_nodes.begin(), m_nodes.end(), compare);

        assert(std::is_sorted(m_nodes.begin(), m_nodes.end(), compare));
    }


    // compute pointers into points
    constexpr auto gap_cell_indicator = std::numeric_limits<unsigned int>::max();
    m_first_in_cell = std::vector<unsigned int>(max_cell_id + 1, gap_cell_indicator);
    {
        ScopedTimer timer("Find first point in cell", m_profile);

        m_first_in_cell[max_cell_id] = n;

        // First, we mark the begin of cells that actually contain points
        // and repair the gaps (i.e., empty cells) later. In the mean time,
        // the values of those gaps will remain at gap_cell_indicator.
        m_first_in_cell[m_nodes[0].cell_id] = 0;
        for (int i = 1; i < n; ++i) {
            if (m_nodes[i - 1].cell_id != m_nodes[i].cell_id) {
                m_first_in_cell[m_nodes[i].cell_id] = i;
            }
        }

        // Now repair gaps: since first_point_in_cell shell contain
        // a prefix sum, we simply replace any "gap_cell_indicator"
        // with its nearest non-gap successor. In the main loop,
        // this is always the direct successors since we're iterating
        // from right to left.
        for(int i=max_cell_id-1; i>=0; --i) {
            m_first_in_cell[i] = std::min(m_first_in_cell[i], m_first_in_cell[i+1]);
        }

        #ifndef NDEBUG
        {
            assert(m_nodes[n-1].cell_id < max_cell_id);

            // assert that we have a prefix sum starting at 0 and ending in n
            assert(m_first_in_cell[0] == 0);
            assert(m_first_in_cell[max_cell_id] == n);
            assert(std::is_sorted(m_first_in_cell.begin(), m_first_in_cell.end()));

            // check that each point is in its right cell (and that the cell boundaries are correct)
            for (auto cid = 0u; cid != max_cell_id; ++cid) {
                const auto begin = m_first_in_cell[cid];
                const auto end = m_first_in_cell[cid + 1];
                for (auto idx = begin; idx != end; ++idx)
                    assert(m_nodes[idx].cell_id == cid);
            }
        }
        #endif
    }

    // build spatial structure and find insertion level for each layer based on lower bound on radius for current and smallest layer
    std::vector<WeightLayer<D>> weight_layers;
    weight_layers.reserve(m_layers);
    {
        ScopedTimer timer("Build data structure", m_profile);
        for (auto layer = 0u; layer < m_layers; ++layer) {
            weight_layers.emplace_back(weightLayerTargetLevel(layer), m_nodes.data(), m_first_in_cell.data() + first_cell_of_layer[layer]);
        }
    }

    return weight_layers;
}


} // namespace satgirgs
