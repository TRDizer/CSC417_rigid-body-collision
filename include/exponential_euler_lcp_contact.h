#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>
#include <rigid_body_jacobian.h>
#include <inverse_rigid_body.h>
#include <unordered_map>
#include <cmath>
#include <iostream>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//  n - list of collision normals
//  x - list of world space collision points
//  obj - list of collision object ids 
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 
static const int floor_id = -1;

inline void exponential_euler_lcp_contact(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces,
                            std::vector<Eigen::Vector3d> &n, std::vector<Eigen::Vector3d> &x, std::vector<std::pair<int,int> > &obj) {
    
    Eigen::VectorXd qdot_unc;
    qdot_unc.resizeLike(qdot);
    qdot_unc.setZero();

    // calculate unconstrained qdot
    for (int i = 0; i < qdot_unc.size() / 6; i++) {
        // Get all the current parameters
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12*i).data());
        // Eigen::Vector3d p = q.segment<3>(12*i+9);
        Eigen::Vector3d omega = qdot.segment<3>(6*i);
        Eigen::Vector3d pdot = qdot.segment<3>(6*i+3);
        Eigen::Matrix3d I = masses[i].block<3,3>(0,0);
        double mass = masses[i](3,3);
        Eigen::Vector3d t_torq = forces.segment<3>(6*i);
        Eigen::Vector3d f_ext = forces.segment<3>(6*i+3);

        // do not update here, need to use collision adjusted qdot
        // Eigen::Matrix3d update_rotation, R_next_t;
        // rodrigues(update_rotation, omega * dt);
        // R_next_t = update_rotation * R;
        // q.segment<9>(12*i) = Eigen::Map<const Eigen::Vector9d>(R_next_t.data());
        // q.segment<3>(12*i+9) = p + dt * pdot;

        qdot_unc.segment<3>(6*i) = (R*I*R.transpose()).inverse() * ((R*I*R.transpose())*omega - dt * omega.cross((R*I*R.transpose())*omega) + dt * t_torq);
        qdot_unc.segment<3>(6*i+3) = (mass * pdot + dt * f_ext) / mass;
    }

    // At this point, the unconstrained qdot for the system is computed
    unsigned int iteration = 0;
    Eigen::VectorXd alpha(n.size());
    alpha.setZero();
    std::unordered_map<int, Eigen::Matrix66d> mass_inverses;
    mass_inverses[-1] = Eigen::Matrix66d::Zero(); // just add it in for completeness, even though computation of velocity for floor will just be skipped
    std::unordered_map<int, Eigen::Matrix36d> jacobians;
    jacobians[-1] = Eigen::Matrix36d::Zero(); // for completeness 

    auto get_mass_inverse = [&](int obj_id) {
        if (mass_inverses.find(obj_id) == mass_inverses.end()) {
            // If inverse not stored, calculate and store
            Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12*obj_id).data());
            //                    I     0
            // masses[obj_id] =              // hence need to fill in R terms to compute M^-1
            //                    0    mI
            Eigen::Matrix66d M = masses[obj_id];
            M.block<3,3>(0,0) = R * M.block<3,3>(0,0) * R.transpose();
            mass_inverses[obj_id] = M.inverse();
        }

        return mass_inverses[obj_id];
    };

    auto get_jacobian = [&](int obj_id, Eigen::Ref<const Eigen::Vector3d> x_world) {
        if (jacobians.find(obj_id) == jacobians.end()) {
            // If jacobian not found, calculate and store
            Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12*obj_id).data());
            Eigen::Vector3d p = q.segment<3>(12*obj_id+9);
            Eigen::Vector3d X; // X bar over here
            inverse_rigid_body(X, x_world, R, p);

            Eigen::Matrix36d J;
            rigid_body_jacobian(J, R, p, X);
            jacobians[obj_id] = J;

            return J;
        }
        else {
            return jacobians[obj_id];
        }
    };

    // can use any negative current_collision to collect the net force from all contacts
    auto collect_other_contact = [&](int current_obj, int current_collision) {
        Eigen::Vector6d return_qdot;
        return_qdot.setZero();

        for (int i = 0; i < n.size(); i++) {
            if (i == current_collision) {
                continue;
            }

            Eigen::Vector3d signed_normal;
            if ((obj[i]).first == current_obj){
                // current obj is A in the collision 
                signed_normal = -n[i];
            }
            else if ((obj[i]).second == current_obj) {
                // current obj is B in the collision 
                signed_normal = n[i];
            }
            else {
                // not a relevant collision 
                continue;
            }

            Eigen::Matrix36d J = get_jacobian(current_obj, x[i]);
            return_qdot += alpha(i) * J.transpose() * signed_normal;
        }

        return return_qdot;
    };

    // std::cout << "num collision points: " << n.size() << std::endl;

    // address collisions
    if (n.size() > 0) {
        while(iteration < 10) {
            // alpha update
            for (int i = 0; i < n.size(); i++) {
                // Get all parameters for the current contact
                int obj_a = obj[i].first;
                int obj_b = obj[i].second;

                Eigen::Vector3d c_normal = n[i]; // collision normal, pointing towards B
                Eigen::Vector3d c_world = x[i]; // collision point in world coordinate

                // they are (g_i)^T * M^-1 * g_i, which are needed for delta_i
                double current_contact_velocity_a, current_contact_velocity_b;
                // they are (g_i)^T * (qdot_unc + dt * M^-1 * f_i), which are needed for gamma_i
                double other_velocity_a, other_velocity_b;

                // encode rule: floor does not move
                if (obj_a == floor_id) {
                    current_contact_velocity_a = 0;
                    other_velocity_a = 0;
                }
                // else we are looking at an object
                else {
                    Eigen::Matrix66d M_a_inverse = get_mass_inverse(obj_a);
                    Eigen::Vector3d signed_normal = -c_normal;
                    Eigen::Matrix36d J_a = get_jacobian(obj_a, c_world);
                    Eigen::Vector6d qdot_unc_a = qdot_unc.segment<6>(6*obj_a);

                    // g_a is the current contact force and f_a is the net contact force for other collisions
                    Eigen::Vector6d g_a = J_a.transpose() * signed_normal;
                    Eigen::Vector6d f_a = collect_other_contact(obj_a, i);
                    current_contact_velocity_a = g_a.transpose() * M_a_inverse * g_a;
                    other_velocity_a = g_a.transpose() * (qdot_unc_a + dt * M_a_inverse * f_a);
                }

                // encode rule: floor does not move
                if (obj_b == floor_id) {
                    current_contact_velocity_b = 0;
                    other_velocity_b = 0;
                }
                else {
                    Eigen::Matrix66d M_b_inverse = get_mass_inverse(obj_b);
                    Eigen::Vector3d signed_normal = c_normal; // normal points from A to B, hence the raw normal repels B
                    Eigen::Matrix36d J_b = get_jacobian(obj_b, c_world);
                    Eigen::Vector6d qdot_unc_b = qdot_unc.segment<6>(6*obj_b);

                    // g_b is the current contact force and f_b is the net contact force for other collisions
                    Eigen::Vector6d g_b = J_b.transpose() * signed_normal;
                    Eigen::Vector6d f_b = collect_other_contact(obj_b, i);
                    current_contact_velocity_b = g_b.transpose() * M_b_inverse * g_b;
                    other_velocity_b = g_b.transpose() * (qdot_unc_b + dt * M_b_inverse * f_b);
                }

                double delta = dt * (current_contact_velocity_a + current_contact_velocity_b);
                double gamma = other_velocity_a + other_velocity_b;
                alpha(i) = std::max(-gamma / delta, 0.0);
            }

            iteration++;

            // early termination check
            // bool violated = false;
            // for (int i = 0; i < n.size(); i++) {
            //     int obj_a = obj[i].first;
            //     int obj_b = obj[i].second;
            //     Eigen::Vector3d c_normal = n[i]; // collision normal, pointing towards B
            //     Eigen::Vector3d c_world = x[i]; // collision point in world coordinate

            //     // point_velocity = J * qdot (aka. general velocity)
            //     Eigen::Vector3d point_velocity_a, point_velocity_b;

            //     if (obj_a == floor_id) {
            //         point_velocity_a.setZero();
            //     }
            //     else {
            //         Eigen::Vector6d qdot_unc_a = qdot_unc.segment<6>(6*obj_a);
            //         point_velocity_a = get_jacobian(obj_a, c_world) * (qdot_unc_a + dt * get_mass_inverse(obj_a) * collect_other_contact(obj_a, -1));
            //     }

            //     if (obj_b == floor_id) {
            //         point_velocity_b.setZero();
            //     }
            //     else {
            //         Eigen::Vector6d qdot_unc_b = qdot_unc.segment<6>(6*obj_b);
            //         point_velocity_b = get_jacobian(obj_b, c_world) * (qdot_unc_b + dt * get_mass_inverse(obj_b) * collect_other_contact(obj_b, -1));
            //     }

            //     // alpha >= 0 is satisfied by design
            //     if (c_normal.dot(point_velocity_b - point_velocity_a) < 0) {
            //         // there are still objects kissing
            //         violated = true;
            //         break;
            //     }
            //     if (alpha(i) * c_normal.dot(point_velocity_b - point_velocity_a) != 0) {
            //         // there are still objects kissing
            //         violated = true;
            //         break;
            //     }
            // }

            // if (violated == false) {
            //     break;
            // }
        }
    }
    // std::cout << "contraint satisfies after " << iteration << "iterations" << std::endl;

    // compute qdot at t+1
    qdot = qdot_unc;
    
    for (int obj_i = 0; obj_i < qdot.size() / 6; obj_i++) {
        qdot.segment<6>(6*obj_i) += dt * get_mass_inverse(obj_i) * collect_other_contact(obj_i, -1);
    }

    // std::cout << (qdot - qdot_unc).norm() << std::endl;

    // for (int i = 0; i < n.size(); i++) {
    //     int obj_a = obj[i].first;
    //     int obj_b = obj[i].second;
    //     Eigen::Vector3d c_normal = n[i]; // collision normal, pointing towards B
    //     Eigen::Vector3d c_world = x[i]; // collision point in world coordinate

    //     if (obj_a != floor_id) {
    //         Eigen::Vector3d signed_normal = -c_normal;
    //         qdot.segment<6>(6*obj_a) += dt * alpha(i) * get_mass_inverse(obj_a) * get_jacobian(obj_a, c_world).transpose() * signed_normal;
    //     }

    //     if (obj_b != floor_id) {
    //         Eigen::Vector3d signed_normal = c_normal;
    //         qdot.segment<6>(6*obj_b) += dt * alpha(i) * get_mass_inverse(obj_b) * get_jacobian(obj_b, c_world).transpose() * signed_normal;
    //     }
    // }

    // update q
    for (int i = 0; i < qdot.size() / 6; i++) {
        // Get all the current parameters
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12*i).data());
        Eigen::Vector3d p = q.segment<3>(12*i+9);
        Eigen::Vector3d omega = qdot.segment<3>(6*i);
        Eigen::Vector3d pdot = qdot.segment<3>(6*i+3);

        Eigen::Matrix3d update_rotation, R_next_t;
        rodrigues(update_rotation, omega * dt);
        R_next_t = update_rotation * R;
        q.segment<9>(12*i) = Eigen::Map<const Eigen::Vector9d>(R_next_t.data());
        q.segment<3>(12*i+9) = p + dt * pdot;
    }
}