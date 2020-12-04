#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {
    // Assumption: omega is pre-multiplied with dt
    Eigen::Vector3d axis = omega.normalized();
    Eigen::Matrix3d cross_axis;
    cross_axis <<    0,     -axis.z(),   axis.y(),
                  axis.z(),    0,       -axis.x(),
                 -axis.y(),  axis.x(),    0;
    double theta = omega.norm();
    
    R = Eigen::Matrix3d::Identity() + std::sin(theta) * cross_axis + (1 - std::cos(theta)) * cross_axis * cross_axis;
}