#include "gtsam_optimizer.h"

namespace gtsam_wrapper
{
//------------------------------------------------------------------------------
// sugar for derivative blocks
#define D_R_R(H) (H)->block<3, 3>(0, 0)
#define D_R_t(H) (H)->block<3, 3>(0, 3)
#define D_R_v(H) (H)->block<3, 3>(0, 6)
#define D_t_R(H) (H)->block<3, 3>(3, 0)
#define D_t_t(H) (H)->block<3, 3>(3, 3)
#define D_t_v(H) (H)->block<3, 3>(3, 6)
#define D_v_R(H) (H)->block<3, 3>(6, 0)
#define D_v_t(H) (H)->block<3, 3>(6, 3)
#define D_v_v(H) (H)->block<3, 3>(6, 6)

    NavState PreintegratedImuMeasurementsGS::predict(const NavStateGS &state_i,
                                                     const imuBias::ConstantBias &bias_i, const Vector3 &gdir, OptionalJacobian<9, 9> H1,
                                                     OptionalJacobian<9, 6> H2, OptionalJacobian<9, 3> H3) const
    {
        // TODO(frank): make sure this stuff is still correct
        Matrix96 D_biasCorrected_bias;
        Vector9 biasCorrected = biasCorrectedDelta(bias_i,
                                                   H2 ? &D_biasCorrected_bias : 0);

        // Correct for initial velocity and gravity
        Matrix9 D_delta_state, D_delta_biasCorrected;
        Matrix93 D_delta_g;
        Vector9 xi = state_i.correctPIM(biasCorrected, deltaTij_, gdir,
                                        p().omegaCoriolis, p().use2ndOrderCoriolis, H1 ? &D_delta_state : 0,
                                        H2 ? &D_delta_biasCorrected : 0, H3 ? &D_delta_g : 0);

        // Use retract to get back to NavStateGS manifold
        Matrix9 D_predict_state, D_predict_delta;
        NavState state_j = state_i.retract(xi, D_predict_state, D_predict_delta);
        if (H1)
            *H1 = D_predict_state + D_predict_delta * D_delta_state;
        if (H2)
            *H2 = D_predict_delta * D_delta_biasCorrected * D_biasCorrected_bias;
        if (H3)
            *H3 = D_predict_delta * D_delta_g;
        return state_j;
    }

    gtsam::Matrix9 PreintegratedImuMeasurementsGS::PreintMeasCov() const
    {
        return preintMeasCov_;
    }

    gtsam::Vector9 PreintegratedImuMeasurementsGS::computeError(const NavStateGS &state_i,
                                                                const NavState &state_j,
                                                                const imuBias::ConstantBias &bias_i,
                                                                const Vector3 &gdir,
                                                                OptionalJacobian<9, 9> H1,
                                                                OptionalJacobian<9, 9> H2,
                                                                OptionalJacobian<9, 6> H3,
                                                                OptionalJacobian<9, 3> H4) const
    {
        // Predict state at time j
        Matrix9 D_predict_state_i;
        Matrix96 D_predict_bias_i;
        Matrix93 D_predict_g;
        NavState predictedState_j = predict(
            state_i, bias_i, gdir, H1 ? &D_predict_state_i : 0, H3 ? &D_predict_bias_i : 0, H4 ? &D_predict_g : 0);

        // Calculate error
        Matrix9 D_error_state_j, D_error_predict;
        Vector9 error =
            state_j.localCoordinates(predictedState_j, H2 ? &D_error_state_j : 0,
                                     H1 || H3 ? &D_error_predict : 0);

        if (H1)
            *H1 << D_error_predict * D_predict_state_i;
        if (H2)
            *H2 << D_error_state_j;
        if (H3)
            *H3 << D_error_predict * D_predict_bias_i;
        if (H4)
            *H4 << D_error_predict * D_predict_g;
        return error;
    }

    //------------------------------------------------------------------------------
    gtsam::Vector9 PreintegratedImuMeasurementsGS::computeErrorAndJacobians(const Pose3 &pose_i, const Vector3 &vel_i, const Pose3 &pose_j, const Vector3 &vel_j,
                                                                            const imuBias::ConstantBias &bias_i, const gtsam::Vector3 gdir, OptionalJacobian<9, 6> H1,
                                                                            OptionalJacobian<9, 3> H2, OptionalJacobian<9, 6> H3,
                                                                            OptionalJacobian<9, 3> H4, OptionalJacobian<9, 6> H5, OptionalJacobian<9, 3> H6) const
    {
        // Note that derivative of constructors below is not identity for velocity, but
        // a 9*3 matrix == Z_3x3, Z_3x3, state.R().transpose()
        NavStateGS state_i(pose_i, vel_i);
        NavState state_j(pose_j, vel_j);

        // Predict state at time j
        Matrix9 D_error_state_i, D_error_state_j;
        Vector9 error = computeError(state_i, state_j, bias_i, gdir,
                                     H1 || H2 ? &D_error_state_i : 0, H3 || H4 ? &D_error_state_j : 0, H5, H6);

        // Separate out derivatives in terms of 5 arguments
        // Note that doing so requires special treatment of velocities, as when treated as
        // separate variables the retract applied will not be the semi-direct product in NavStateGS
        // Instead, the velocities in nav are updated using a straight addition
        // This is difference is accounted for by the R().transpose calls below
        if (H1)
            *H1 << D_error_state_i.leftCols<6>();
        if (H2)
            *H2 << D_error_state_i.rightCols<3>() * state_i.R().transpose();
        if (H3)
            *H3 << D_error_state_j.leftCols<6>();
        if (H4)
            *H4 << D_error_state_j.rightCols<3>() * state_j.R().transpose();

        return error;
    }

    Vector9 NavStateGS::correctPIM(const Vector9 &pim, double dt,
                                   const Vector3 &n_gravity, const boost::optional<Vector3> &omegaCoriolis,
                                   bool use2ndOrderCoriolis, OptionalJacobian<9, 9> H1,
                                   OptionalJacobian<9, 9> H2, OptionalJacobian<9, 3> H3) const
    {
        const Rot3 &nRb = gtsam::Rot3(R());
        const Velocity3 &n_v = v(); // derivative is Ri !
        const double dt22 = 0.5 * dt * dt;

        Vector9 xi;
        Matrix3 D_dP_Ri1, D_dP_Ri2, D_dP_nv, D_dV_Ri;
        Matrix3 D_dP_g, D_dV_g;

        dR(xi) = dR(pim);
        dP(xi) = dP(pim) + dt * nRb.unrotate(n_v, H1 ? &D_dP_Ri1 : 0, H2 ? &D_dP_nv : 0) + dt22 * nRb.unrotate(n_gravity, H1 ? &D_dP_Ri2 : 0, H3 ? &D_dP_g : 0);
        dV(xi) = dV(pim) + dt * nRb.unrotate(n_gravity, H1 ? &D_dV_Ri : 0, H3 ? &D_dV_g : 0);

        if (omegaCoriolis)
        {
            xi += coriolis(dt, *omegaCoriolis, use2ndOrderCoriolis, H1);
        }

        if (H1 || H2)
        {
            Matrix3 Ri = nRb.matrix();

            if (H1)
            {
                if (!omegaCoriolis)
                    H1->setZero(); // if coriolis H1 is already initialized
                D_t_R(H1) += dt * D_dP_Ri1 + dt22 * D_dP_Ri2;
                D_t_v(H1) += dt * D_dP_nv * Ri;
                D_v_R(H1) += dt * D_dV_Ri;
            }
            if (H2)
            {
                H2->setIdentity();
            }
        }

        gtsam::Matrix93 Jac_xi_g;
        Jac_xi_g.setZero();
        if (H3)
        {
            Jac_xi_g.block<3, 3>(3, 0) = dt22 * D_dP_g;
            Jac_xi_g.block<3, 3>(6, 0) = dt * D_dP_nv;
            *H3 = Jac_xi_g;
        }

        return xi;
    }
} // namespace gtsam_wrapper

namespace gtsam_optimizer
{
    int GtsamOptimizer::InertialOptimization(ORB_SLAM3::Map *pMap, Eigen::Matrix3d &Rwg, double &scale, Eigen::Vector3d &bg,
                                             Eigen::Vector3d &ba, bool bMono, Eigen::MatrixXd &covInertial, bool bFixedVel, bool bGauss, float priorG, float priorA)
    {
    }

    gtsam::Vector InertialEgdeGS::evaluateError(gtsam::Pose3 &pose_i, gtsam::Vector3 &vel_i, gtsam::Pose3 &pose_j, gtsam::Vector3 &vel_j,
                                                const Rot3 &rwg, const double &scale, imuBias::ConstantBias &bias_i,
                                                boost::optional<Matrix &> H_pi,
                                                boost::optional<Matrix &> H_vi,
                                                boost::optional<Matrix &> H_pj,
                                                boost::optional<Matrix &> H_vj,
                                                boost::optional<Matrix &> H_rwg,
                                                boost::optional<Matrix &> H_scale,
                                                boost::optional<Matrix &> H_b) const
    {
        gtsam::Matrix96 _H1, _H3, _H5;
        gtsam::Matrix93 _H2, _H4, _H6;
        gtsam::Pose3 pose_i_scale(pose_i.rotation(), pose_i.translation() * scale);
        gtsam::Pose3 pose_j_scale(pose_j.rotation(), pose_j.translation() * scale);

        gtsam::Matrix3 H_1_rwg;
        gtsam::Vector3 gravity_rot = rwg.rotate(gI_, H_1_rwg);
        gtsam::Vector9 imu_pim_error = _PIM_.computeErrorAndJacobians(pose_i_scale, scale * vel_i, pose_j_scale, scale * vel_j, bias_i,
                                                                      gravity_rot, _H1, _H2, _H3, _H4, _H5, _H6);

        if (H_pi)
        {
            //*H_pi = ;
        }
        if (H_vi)
        {
            *H_vi = _H2 * scale;
        }
        if (H_pj)
        {
            //*H_pj = ;
        }
        if (H_vj)
        {
            *H_vj = _H4 * scale;
        }
        if (H_rwg)
        {
            *H_rwg = _H6 * H_1_rwg;
        }
        if (H_b)
        {
            *H_b = _H5;
        }
        if (H_scale)
        {
        }

        return imu_pim_error;
    }

}; // namespace gtsam_optimizer