#ifndef VELO_PHYSICS_ENERGY_H
#define VELO_PHYSICS_ENERGY_H

#include <vector>
#include <array>

namespace velo::physics::constraints
{
template <typename ConstraintT>
struct Energy : public ConstraintT {
    std::vector<std::array<int, ConstraintT::NbParticles>> indices;
};
}  // namespace velo::physics::constraints

#endif  //  VELO_PHYSICS_ENERGY_H
