#ifndef VELO_PHYSICS_XPBD_H
#define VELO_PHYSICS_XPBD_H

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Velo/Core/Space.h>

#include <Velo/Physics/State.h>

#include <Velo/Physics/velophysics_export.h>

namespace velo
{

namespace constraints
{
struct StVK;
}

namespace xpbd
{

void VELOPHYSICS_EXPORT step(velo::Real dt,
                             int substeps,
                             velo::MechanicalState<velo::Space3D, velo::Space3D> &state,
                             const velo::constraints::StVK &constraints);

}  // namespace xpbd
}  // namespace velo

#endif  // VELO_PHYSICS_XPBD_H