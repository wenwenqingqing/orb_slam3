#pragma once
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/navigation/ImuFactor.h>
#include <gtsam/navigation/NavState.h>
#include <gtsam/navigation/PreintegrationBase.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/NonlinearOptimizer.h>

#include "Map.h"
#include "MapPoint.h"
#include "KeyFrame.h"
#include "LoopClosing.h"
#include "Frame.h"

using namespace gtsam;

namespace gtsam_wrapper
{
    /** A convenient base class for creating your own NoiseModelFactor with 7
    * variables.  To derive from this class, implement evaluateError(). */
    template <class VALUE1, class VALUE2, class VALUE3, class VALUE4, class VALUE5, class VALUE6, class VALUE7>
    class NoiseModelFactor7 : public NoiseModelFactor
    {

    public:
        // typedefs for value types pulled from keys
        typedef VALUE1 X1;
        typedef VALUE2 X2;
        typedef VALUE3 X3;
        typedef VALUE4 X4;
        typedef VALUE5 X5;
        typedef VALUE6 X6;
        typedef VALUE7 X7;

    protected:
        typedef NoiseModelFactor Base;
        typedef NoiseModelFactor7<VALUE1, VALUE2, VALUE3, VALUE4, VALUE5, VALUE6, VALUE7> This;

    public:
        /**
         * Default Constructor for I/O
         */
        NoiseModelFactor7() {}

        /**
         * Constructor
         * @param noiseModel shared pointer to noise model
         * @param j1 Key of the first variable
         * @param j2 Key of the second variable
         * @param j3 Key of the third variable
         * @param j4 Key of the fourth variable
         * @param j5 Key of the fifth variable
         * @param j6 Key of the sixth variable
         * @param j7 Key of the seventh variable
         */
        NoiseModelFactor7(const SharedNoiseModel &noiseModel, Key j1, Key j2, Key j3, Key j4, Key j5, Key j6, Key j7) : Base(noiseModel, cref_list_of<7>(j1)(j2)(j3)(j4)(j5)(j6)(j7)) {}

        virtual ~NoiseModelFactor7() {}

        /** methods to retrieve keys */
        inline Key key1() const { return keys_[0]; }
        inline Key key2() const { return keys_[1]; }
        inline Key key3() const { return keys_[2]; }
        inline Key key4() const { return keys_[3]; }
        inline Key key5() const { return keys_[4]; }
        inline Key key6() const { return keys_[5]; }
        inline Key key7() const { return keys_[6]; }

        /** 
         * Calls the 7-Key specific version of evaluateError, which is pure virtual
         * so must be implemented in the derived class. 
         */
        virtual Vector unwhitenedError(const Values &x, boost::optional<std::vector<Matrix> &> H = boost::none) const
        {
            if (this->active(x))
            {
                if (H)
                    return evaluateError(x.at<X1>(keys_[0]), x.at<X2>(keys_[1]), x.at<X3>(keys_[2]), x.at<X4>(keys_[3]), x.at<X5>(keys_[4]), x.at<X6>(keys_[5]), x.at<X7>(keys_[6]), (*H)[0], (*H)[1], (*H)[2], (*H)[3], (*H)[4], (*H)[5], (*H)[6]);
                else
                    return evaluateError(x.at<X1>(keys_[0]), x.at<X2>(keys_[1]), x.at<X3>(keys_[2]), x.at<X4>(keys_[3]), x.at<X5>(keys_[4]), x.at<X6>(keys_[5]), x.at<X7>(keys_[6]));
            }
            else
            {
                return Vector::Zero(this->dim());
            }
        }

        /**
         *  Override this method to finish implementing a 7-way factor.
         *  If any of the optional Matrix reference arguments are specified, it should compute
         *  both the function evaluation and its derivative(s) in X1 (and/or X2, X3).
         */
        virtual Vector
        evaluateError(const X1 &, const X2 &, const X3 &, const X4 &, const X5 &, const X6 &, const X7 &,
                      boost::optional<Matrix &> H1 = boost::none,
                      boost::optional<Matrix &> H2 = boost::none,
                      boost::optional<Matrix &> H3 = boost::none,
                      boost::optional<Matrix &> H4 = boost::none,
                      boost::optional<Matrix &> H5 = boost::none,
                      boost::optional<Matrix &> H6 = boost::none,
                      boost::optional<Matrix &> H7 = boost::none) const = 0;

    private:
        /** Serialization function */
        friend class boost::serialization::access;
        template <class ARCHIVE>
        void serialize(ARCHIVE &ar, const unsigned int /*version*/)
        {
            ar &boost::serialization::make_nvp("NoiseModelFactor",
                                               boost::serialization::base_object<Base>(*this));
        }
    }; // \class NoiseModelFactor7

    class GTSAM_EXPORT NavStateGS : public gtsam::NavState
    {
    public:
        NavStateGS(const gtsam::Pose3 p, const gtsam::Vector v)
            : NavState(p, v)
        {
        }

        Vector9 correctPIM(const Vector9 &pim, double dt,
                           const Vector3 &n_gravity, const boost::optional<Vector3> &omegaCoriolis,
                           bool use2ndOrderCoriolis, OptionalJacobian<9, 9> H1,
                           OptionalJacobian<9, 9> H2, OptionalJacobian<9, 3> H3) const;
    };

    class GTSAM_EXPORT PreintegratedImuMeasurementsGS : public PreintegratedImuMeasurements
    {
    public:
        NavState predict(const NavStateGS &state_i,
                         const imuBias::ConstantBias &bias_i,
                         const Vector3 &gdir, OptionalJacobian<9, 9> H1,
                         OptionalJacobian<9, 6> H2, OptionalJacobian<9, 3> H3) const;

        gtsam::Matrix9 PreintMeasCov() const;
        gtsam::Vector9 computeError(const NavStateGS &state_i,
                                    const NavState &state_j,
                                    const imuBias::ConstantBias &bias_i,
                                    const Vector3 &gdir,
                                    OptionalJacobian<9, 9> H1,
                                    OptionalJacobian<9, 9> H2,
                                    OptionalJacobian<9, 6> H3,
                                    OptionalJacobian<9, 3> H4) const;

        gtsam::Vector9 computeErrorAndJacobians(const Pose3 &pose_i,
                                                const Vector3 &vel_i, const Pose3 &pose_j, const Vector3 &vel_j,
                                                const imuBias::ConstantBias &bias_i, const gtsam::Vector3 gdir, OptionalJacobian<9, 6> H1,
                                                OptionalJacobian<9, 3> H2, OptionalJacobian<9, 6> H3,
                                                OptionalJacobian<9, 3> H4, OptionalJacobian<9, 6> H5, OptionalJacobian<9, 3> H6) const;
    };
}; // namespace gtsam_wrapper

namespace gtsam_optimizer
{
    class GTSAM_EXPORT InertialEgdeGS : public gtsam_wrapper::NoiseModelFactor7<gtsam::Pose3, gtsam::Vector3, gtsam::Pose3, gtsam::Vector3,
                                                                                gtsam::Rot3, double, imuBias::ConstantBias>
    {
        typedef InertialEgdeGS This;
        typedef gtsam_wrapper::NoiseModelFactor7<gtsam::Pose3, gtsam::Vector3, gtsam::Pose3, gtsam::Vector3,
                                                 gtsam::Rot3, double, imuBias::ConstantBias>
            Base;
        typedef gtsam_wrapper::PreintegratedImuMeasurementsGS PimGS;

        InertialEgdeGS() {}

    public:
        InertialEgdeGS(gtsam::Key posei_key, gtsam::Key veli_key, gtsam::Key posej_key, gtsam::Key velj_key,
                       gtsam::Key rwg_key, gtsam::Key scale_key, gtsam::Key b_key, const PimGS &p)
            : Base(noiseModel::Gaussian::Covariance(p.PreintMeasCov()), posei_key, veli_key, posej_key, velj_key, rwg_key, scale_key, b_key), _PIM_(p), dt_(p.deltaTij())
        {
        }

        gtsam::Vector evaluateError(gtsam::Pose3 &pose_i, gtsam::Vector3 &vel_i, gtsam::Pose3 &pose_j, gtsam::Vector3 &vel_j,
                                    const Rot3 &rwg, const double &scale, imuBias::ConstantBias &bias_i,
                                    boost::optional<Matrix &> H_pi,
                                    boost::optional<Matrix &> H_vi,
                                    boost::optional<Matrix &> H_pj,
                                    boost::optional<Matrix &> H_vj,
                                    boost::optional<Matrix &> H_rwg,
                                    boost::optional<Matrix &> H_scale,
                                    boost::optional<Matrix &> H_b) const;

    private:
        PimGS _PIM_;
        double dt_;
        gtsam::Vector3 gI_;
    };

    class GtsamOptimizer
    {
        int InertialOptimization(ORB_SLAM3::Map *pMap, Eigen::Matrix3d &Rwg, double &scale, Eigen::Vector3d &bg,
                                 Eigen::Vector3d &ba, bool bMono, Eigen::MatrixXd &covInertial, bool bFixedVel, bool bGauss, float priorG, float priorA);
    };

}; // namespace gtsam_optimizer