
#include <cassert>

namespace satgirgs {

template<size_t D>
std::array<uint32_t, D> extract(uint32_t cell) {
    std::array<uint32_t, D> result;
    result.fill(0);
    auto bit = 0u;
    for(auto l = 0u; (l+1)*D <= 32; ++l) {
        for(auto d = 0u; d<D; ++d) {
            result[d] |= ((cell >> bit++) & 1u) << l;
        }
    }
    return result;
}

template<size_t D>
uint32_t deposit(std::array<uint32_t, D> coords) {
    if (D == 1) return coords.front();

    uint32_t result = 0u;
    uint32_t bit = 0;
    // stop if next level would be (partly) unrepresentable
    // thus we encode the same number of bits for each coordinate
    for(auto l = 0u; (l+1)*D <= 32; l++) {
        for(auto d = 0u; d < D; d++) {
            result |= ((coords[d] >> l) & 1u) << bit++;
        }
    }
    return result;
}


template<unsigned int D>
unsigned int SpatialTreeCoordinateHelper<D>::cellOfLevel(unsigned cell) noexcept {
    // sets all bits below the most significant bit set in x
    auto assertLower = [] (uint32_t x) {
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return x;
    };

    // cool c++17 constexpr lambdas
    constexpr auto mask = [](){
        auto res = 0u;
        for(auto bit = 0u; bit < 32; bit += D)
            res |= 1u << bit;
        return res;
    }();

    auto firstCellInLayer = mask & assertLower(cell);

    if (cell < firstCellInLayer)
        firstCellInLayer >>= D;

    return cell - firstCellInLayer;
}

template<unsigned int D>
std::array<std::pair<double, double>, D> SpatialTreeCoordinateHelper<D>::bounds(unsigned int cell, unsigned int level) noexcept {
    const auto diameter = 1.0 / (1<<level);
    const auto coord = extract<D>(cellOfLevel(cell));

    auto result = std::array<std::pair<double, double>, D>();
    for(auto d=0u; d<D; ++d)
        result[d]= {coord[d]*diameter, (coord[d]+1)*diameter };

    return result;
}

template<unsigned int D>
unsigned int SpatialTreeCoordinateHelper<D>::cellForPoint(const std::array<double, D>& position, unsigned int targetLevel) noexcept {
    const auto diameter = static_cast<double>(1 << targetLevel);

    std::array<uint32_t, D> coords;
    for (auto d = 0u; d < D; ++d)
        coords[d] = static_cast<uint32_t>(position[d] * diameter);

    return deposit(coords);
}

template<unsigned int D>
bool SpatialTreeCoordinateHelper<D>::touching(unsigned int cellA, unsigned int cellB, unsigned int level) noexcept  {
    const auto coordA = extract<D>(cellOfLevel(cellA));
    const auto coordB = extract<D>(cellOfLevel(cellB));

    auto touching = true;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(static_cast<int>(coordA[d]) - static_cast<int>(coordB[d]));
        dist = std::min(dist, (1<<level) - dist);
        touching &= (dist <= 1);
    }

    return touching;
}

template<unsigned int D>
double SpatialTreeCoordinateHelper<D>::dist(unsigned int cellA, unsigned int cellB, unsigned int level) noexcept  {
    // first work with integer d dimensional index
    const auto coordA = extract<D>(cellOfLevel(cellA));
    const auto coordB = extract<D>(cellOfLevel(cellB));

    auto result = 0;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(static_cast<int>(coordA[d]) - static_cast<int>(coordB[d]));
        dist = std::min(dist, (1<<level) - dist);
        result = std::max(result, dist);
    }

    // then apply the diameter
    auto diameter = 1.0 / (1<<level);
    return std::max(0.0, (result-1) * diameter); // TODO if cellA and cellB are not touching, this max is irrelevant
}
} // namespace satgirgs
