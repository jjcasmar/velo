#ifndef VELO_PHYSICS_STVKCONSTRAINT_H
#define VELO_PHYSICS_STVKCONSTRAINT_H

#include <Eigen/Dense>

#include <Velo/Core/Types.h>
#include <Velo/Core/Space.h>

#include <Velo/Physics/Energy.h>

#include <Velo/Physics/velophysics_export.h>

namespace velo::physics
{

namespace constraints
{

template <typename WorldSpaceT, typename MaterialSpaceT>
struct StVK {
    // Number of particles this constraint type affects
    // For example, for StVK in a tetrahedron, this should be 4

    // For now this is assuming we are working with a tetrahedron, which assumes Material and World spaces
    // TODO Make it generic
    // Number of particles should depend on the simplex of the material space
    static constexpr int NbParticles = 4;

    // The rest pose for this energy
    using RestPoseT = Eigen::Matrix<typename MaterialSpaceT::Scalar, Eigen::Dynamic, 3>;
    using PoseT = Eigen::Matrix<typename WorldSpaceT::Scalar, Eigen::Dynamic, 3>;

    // Number of actual constraint terms this constraint has
    // For example, in StVK this would be 2, the hidrostatic and the deviatoric
    static constexpr int CSize = 2;

    // These are helpers, they should only be modified by the initialization
    std::vector<Eigen::Matrix<typename MaterialSpaceT::Scalar, MaterialSpaceT::Dim, NbParticles - 1>> invDm;
    std::vector<typename WorldSpaceT::Scalar> volumes;
    std::vector<Eigen::Matrix<float, WorldSpaceT::Dim *(NbParticles - 1), NbParticles * WorldSpaceT::Dim>> dFdx;

    // TODO Add a material manager and specify which material each element uses
    std::vector<float> lambda;
    std::vector<float> mu;
};

void VELOPHYSICS_EXPORT initialize(Energy<StVK<Space3D, Space3D>> &constraints,
                                   const StVK<Space3D, Space3D>::RestPoseT &x0);

std::pair<Eigen::Vector<Space3D::Scalar, Energy<StVK<Space3D, Space3D>>::CSize>,
          Eigen::Matrix<Space3D::Scalar,
                        Energy<StVK<Space3D, Space3D>>::NbParticles * Space3D::Dim,
                        Energy<StVK<Space3D, Space3D>>::CSize>>
    VELOPHYSICS_EXPORT gradients(const Energy<StVK<Space3D, Space3D>> &constraints,
                                 const StVK<Space3D, Space3D>::PoseT &x,
                                 int cId);

Eigen::Matrix<Space3D::Scalar, 2, 2> VELOPHYSICS_EXPORT compliance(const Energy<StVK<Space3D, Space3D>> &constraints,
                                                                   int cId);

}  // namespace constraints
}  // namespace velo::physics

#endif  // VELO_PHYSICS_STVKCONSTRAINT_H