#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> X) {

    Eigen::Matrix3d cross_X;
    cross_X <<     0,     -X.z(),   X.y(),
                  X.z(),    0,     -X.x(),
                 -X.y(),   X.x(),    0;

    // temp1 = | [X]^T   I |
    Eigen::Matrix36d temp1;
    temp1 << cross_X.transpose(), Eigen::Matrix3d::Identity();

    //         | R^T  0 |
    // temp2 = |        |
    //         |  0  R^T|
    Eigen::Matrix66d temp2;
    temp2.setZero();
    temp2.block<3,3>(0,0) = R.transpose();
    temp2.block<3,3>(3,3) = R.transpose();

    J = R * temp1 * temp2;
}

