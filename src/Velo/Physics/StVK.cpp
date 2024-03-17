#include <Velo/Physics/StVK.h>

#include <Velo/Core/Space.h>

#include <Velo/Physics/Energy.h>
#include <Velo/Physics/MechanicalState.h>

#include <Eigen/Core>

#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>

namespace velo::physics::constraints
{

void initialize(Energy<StVK<Space3D, Space3D>> &constraints, const Energy<StVK<Space3D, Space3D>>::RestPoseT &x0)
{
    constraints.invDm.resize(constraints.indices.size());
    constraints.dFdx.resize(constraints.indices.size());
    constraints.volumes.resize(constraints.indices.size());
    constraints.lambda.resize(constraints.indices.size());
    constraints.mu.resize(constraints.indices.size());
    for (int cId{0UL}; cId < constraints.indices.size(); ++cId) {
        const auto &indices = constraints.indices[cId];

        auto dm = (x0({indices[1], indices[2], indices[3]}, Eigen::indexing::all).rowwise() -
                   x0(indices[0], Eigen::indexing::all))
                      .transpose()
                      .eval();

        constraints.volumes[cId] = 1.0F / 6.0F * dm.determinant();

        constraints.invDm[cId] = dm.inverse();

        constraints.lambda[cId] = 210000 * 0.33 / ((1 + 0.33) * (1 - 2 * 0.33));
        constraints.mu[cId] = 210000 / (2 * (1 + 0.33));

        // Apendix E from T.Kim Siggraph course
        const Space3D::Scalar m = constraints.invDm[cId](0, 0);
        const Space3D::Scalar n = constraints.invDm[cId](0, 1);
        const Space3D::Scalar o = constraints.invDm[cId](0, 2);
        const Space3D::Scalar p = constraints.invDm[cId](1, 0);
        const Space3D::Scalar q = constraints.invDm[cId](1, 1);
        const Space3D::Scalar r = constraints.invDm[cId](1, 2);
        const Space3D::Scalar s = constraints.invDm[cId](2, 0);
        const Space3D::Scalar t = constraints.invDm[cId](2, 1);
        const Space3D::Scalar u = constraints.invDm[cId](2, 2);

        const Space3D::Scalar t1 = -m - p - s;
        const Space3D::Scalar t2 = -n - q - t;
        const Space3D::Scalar t3 = -o - r - u;

        auto &PFPx{constraints.dFdx[cId]};
        PFPx.setZero();
        PFPx(0, 0) = t1;
        PFPx(0, 3) = m;
        PFPx(0, 6) = p;
        PFPx(0, 9) = s;
        PFPx(1, 1) = t1;
        PFPx(1, 4) = m;
        PFPx(1, 7) = p;
        PFPx(1, 10) = s;
        PFPx(2, 2) = t1;
        PFPx(2, 5) = m;
        PFPx(2, 8) = p;
        PFPx(2, 11) = s;
        PFPx(3, 0) = t2;
        PFPx(3, 3) = n;
        PFPx(3, 6) = q;
        PFPx(3, 9) = t;
        PFPx(4, 1) = t2;
        PFPx(4, 4) = n;
        PFPx(4, 7) = q;
        PFPx(4, 10) = t;
        PFPx(5, 2) = t2;
        PFPx(5, 5) = n;
        PFPx(5, 8) = q;
        PFPx(5, 11) = t;
        PFPx(6, 0) = t3;
        PFPx(6, 3) = o;
        PFPx(6, 6) = r;
        PFPx(6, 9) = u;
        PFPx(7, 1) = t3;
        PFPx(7, 4) = o;
        PFPx(7, 7) = r;
        PFPx(7, 10) = u;
        PFPx(8, 2) = t3;
        PFPx(8, 5) = o;
        PFPx(8, 8) = r;
        PFPx(8, 11) = u;
    }
}

// Computes the constraint values and the gradient of the constraints wrt particle positions
std::pair<Eigen::Vector<Space3D::Scalar, Energy<StVK<Space3D, Space3D>>::CSize>,
          Eigen::Matrix<Space3D::Scalar,
                        Energy<StVK<Space3D, Space3D>>::NbParticles * Space3D::Dim,
                        Energy<StVK<Space3D, Space3D>>::CSize>>
gradients(const Energy<StVK<Space3D, Space3D>> &constraints, const StVK<Space3D, Space3D>::PoseT &x, int cId)
{
    // For StVK we can write two constraints
    // Ch = det(F) - i
    // Cd = sqrt(tr(F.T * F))
    // Neo Hookean energy is
    // E = V * (Ch * lambda/2 * Ch + Cd * mu * Cd)

    const auto &indices = constraints.indices[cId];

    auto poseMatrix =
        (x({indices[1], indices[2], indices[3]}, Eigen::all).rowwise() - x(indices[0], Eigen::all)).transpose().eval();

    const auto F = (poseMatrix * constraints.invDm[cId]).eval();

    const auto cH = 0.5 * ((F.transpose() * F).trace() - 3);
    const auto &dcHdF = F;
    const auto dcHdx = (constraints.dFdx[cId].transpose() * dcHdF.reshaped()).eval();

    const auto E = (0.5 * (F.transpose() * F - Eigen::Matrix<Space3D::Scalar, 3, 3>::Identity())).eval();
    const auto cD = std::sqrt((E * E).trace());
    const auto &dcDdF = ((1.0 / (2 * (cD + 1e-8))) * (F * F.transpose() * F - F)).eval();
    const auto dcDdx = (constraints.dFdx[cId].transpose() * dcDdF.reshaped()).eval();

    auto c = Eigen::Matrix<Space3D::Scalar, 2, 1>{cH, cD};

    auto dC = Eigen::Matrix<Space3D::Scalar, 12, 2>::Zero().eval();
    dC.block<12, 1>(0, 0) = dcHdx;
    dC.block<12, 1>(0, 1) = dcDdx;

    return {c, dC};
}

Eigen::Matrix<Space3D::Scalar, 2, 2> compliance(const Energy<StVK<Space3D, Space3D>> &constraints, int cId)
{
    const float alphaH = 1.0F / (constraints.volumes[cId] * constraints.lambda[cId]);
    const float alphaD = 1.0F / (constraints.volumes[cId] * 2 * constraints.mu[cId]);

    return Eigen::Matrix<Space3D::Scalar, 2, 2>{{alphaH, 0}, {0, alphaD}};
}

}  // namespace velo::constraints