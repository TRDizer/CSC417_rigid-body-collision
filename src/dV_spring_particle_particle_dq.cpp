#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {

   f.setZero();
    double dV_constant = stiffness * (1 - l0 / (q1-q0).norm());
    // f = dV(q)/dq = dV(q)/d(q0_x, q0_y, q0_z, q1_x, q1_y, q1_z)
    // hence the the last term of q0_? would yield a -1 to flip the difference term
    // f << dV_constant * (q0 - q1)(0),
    //      dV_constant * (q0 - q1)(1),
    //      dV_constant * (q0 - q1)(2),
    //      dV_constant * (q1 - q0)(0),
    //      dV_constant * (q1 - q0)(1),
    //      dV_constant * (q1 - q0)(2);
    f << dV_constant * (q0 - q1), 
         dV_constant * (q1 - q0);
    
}