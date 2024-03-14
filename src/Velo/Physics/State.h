#ifndef VELO_PHYSICS_STATE_H
#define VELO_PHYSICS_STATE_H

#include <Eigen/Dense>

namespace velo
{
template <typename TMaterialSpace, typename TWorldSpace>
struct MechanicalState {
    // The mechanical state of an object is composed of the rest position of the object
    // the world positions and the velocity
    using MaterialScalar = typename TMaterialSpace::Scalar;
    using WorldScalar = typename TWorldSpace::Scalar;

    using RestPoseT = typename Eigen::Matrix<MaterialScalar, Eigen::Dynamic, TMaterialSpace::Dim>;
    using PositionT = typename Eigen::Matrix<WorldScalar, Eigen::Dynamic, TWorldSpace::Dim>;
    using VelocitiesT = typename Eigen::Matrix<WorldScalar, Eigen::Dynamic, TWorldSpace::Dim>;

    // Positions in the rest configuration for this set of particles
    RestPoseT x0;

    // Positions and velocities in the world for this set of particles
    PositionT x;
    VelocitiesT v;

    // TODO move this to its own components
    // Mass for each particle
    Eigen::Matrix<WorldScalar, Eigen::Dynamic, 1> mass;

    // Is the ith particle fixed?
    std::vector<int> fixed;
};
}  // namespace velo

#endif  //  VELO_PHYSICS_STATE_H