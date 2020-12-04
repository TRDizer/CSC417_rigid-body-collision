#include <inverse_rigid_body.h>

void inverse_rigid_body(Eigen::Vector3d &X, Eigen::Ref<const Eigen::Vector3d> x_world, 
                        Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p) {

    X = R.lu().solve(x_world - p);

}