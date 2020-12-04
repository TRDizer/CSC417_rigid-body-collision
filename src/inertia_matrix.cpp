#include <inertia_matrix.h>
#include <cassert>
#include <igl/doublearea.h>
#include <cmath>
#include <iostream>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {

    I.setZero();
    center.setZero();
    mass = 0;

    Eigen::VectorXd A;
    igl::doublearea(V, F, A);

    // interation 1 for getting mass and center of mass
    Eigen::Vector3d X0, X1, X2, N;
    Eigen::Array3d X0_a, X1_a, X2_a;
    for (int i = 0; i < F.rows(); i++) {
        X0 = V.row(F(i,0)).transpose();
        X1 = V.row(F(i,1)).transpose();
        X2 = V.row(F(i,2)).transpose();

        // for element wise computation
        X0_a = X0.array();
        X1_a = X1.array();
        X2_a = X2.array();

        N = ((X1 - X0).cross(X2 - X0)).normalized();
        // Barycentric integration over surface for X_x * N_x can be rewritten to the following and taking out constant N_x
        // N_x * 2 * Area * [sigma{0 to 1} sigma{0 to 1-phi1} X0_x * phi0 + X1_x * phi1 + X2_x * phi2 dphi2 dphi1]
        // which evaluates to (x0 + x1 + x2)/3 * A * N_x (note that A is just the area here)
        mass += (X0.x() + X1.x() + X2.x()) * N.x() * A(i); // A here is double area, and density and quotient will be multipled at the end
        center += ((X0_a.square() + X0_a * X1_a + X0_a * X2_a + X1_a.square() + X1_a * X2_a + X2_a.square()) * N.array()).matrix() * A(i); // Add density and quotient at the end
    }

    mass *= density / 6.0;
    center *= density / 24.0 / mass;

    // std::cout << "mass: " << mass << std::endl;
    // std::cout << "center:\n" << center << std::endl;

    // Shift over to center of mass coordinate
    Eigen::MatrixXd V_bar = V.rowwise() - center.transpose();
    // iteration 2 for [X_bar] * [X_bar]^T
    // will reuse X0, X1, X2, N again but now they all have a bar above them
    Eigen::Vector3d Xxyz2, Xxyz2Xyzx; // translate to vector for [Xx^2 | Xy^2 | Xz^2] and [Xx^2 * Xy | Xy^2 * Xz |  Xz^2 * Xx]
    Eigen::Array3d X0csu, X1csu, X2csu; // They are the circular shifted up by 1 (|y z x| instead of |x y z|) equivalence of X0 X1 and X2
    for (int i = 0; i < F.rows(); i++) {
        X0 = V_bar.row(F(i,0)).transpose();
        X1 = V_bar.row(F(i,1)).transpose();
        X2 = V_bar.row(F(i,2)).transpose();

        N = ((X1 - X0).cross(X2 - X0)).normalized();

        // for element wise computation
        X0_a = X0.array();
        X1_a = X1.array();
        X2_a = X2.array();
        X0csu << X0.y(), X0.z(), X0.x();
        X1csu << X1.y(), X1.z(), X1.x();
        X2csu << X2.y(), X2.z(), X2.x();

        // Will integrate several terms seperately and combine to form the inertia matrix
        // (X0^3 + X0^2*x1 + X0^2*X2 + X0*X1^2 + X0*X1*X2 + X0*X2^2 + X1^3 + X1^2*X2 + X1*X2^2 + X2^3) / 20 * 1/3 * 2 * A * N * rho [again, quotient and density will be added last]
        Xxyz2 = ((X0_a.pow(3) + X0_a.square()*X1_a + X0_a.square()*X2_a + X0_a*X1_a.square() + X0_a*X1_a*X2_a + X0_a*X2_a.square() + X1_a.pow(3) + X1_a.square()*X2_a + X1_a*X2_a.square() + X2_a.pow(3)) * N.array()).matrix() * A(i);
        
        // (3*x0_x^2*x0_y + x0_y*x1_x^2 + x0_x^2*x1_y + x0_y*x2_x^2 + x0_x^2*x2_y + 3*x1_x^2*x1_y + x1_y*x2_x^2 + x1_x^2*x2_y + 3*x2_x^2*x2_y + 2*x0_x*x0_y*x1_x + 2*x0_x*x0_y*x2_x + 2*x0_x*x1_x*x1_y + x0_x*x1_x*x2_y + x0_x*x1_y*x2_x + x0_y*x1_x*x2_x + 2*x0_x*x2_x*x2_y + 2*x1_x*x1_y*x2_x + 2*x1_x*x2_x*x2_y) / 60 / 2 * 2 * A * N * rho
        // just like terms above, quotient and density will be added last
        Xxyz2Xyzx = ((3*X0_a.square()*X0csu + X0csu*X1_a.square() + X0_a.square()*X1csu + X0csu*X2_a.square() + X0_a.square()*X2csu + 3*X1_a.square()*X1csu + X1csu*X2_a.square() + X1_a.square()*X2csu + 3*X2_a.square()*X2csu + 2*X0_a*X0csu*X1_a + 2*X0_a*X0csu*X2_a + 2*X0_a*X1_a*X1csu + X0_a*X1_a*X2csu + X0_a*X1csu*X2_a + X0csu*X1_a*X2_a + 2*X0_a*X2_a*X2csu + 2*X1_a*X1csu*X2_a + 2*X1_a*X2_a*X2csu) * N.array()).matrix() * A(i);
    
        // Now we have all the needed terms. Time to combine!
        // Xy^2+Xz^2    -XxXy      -XxXz
        //  -XxXy     Xx^2+Xz^2    -XyXz
        //  -XxXz       -XyXz    Xx^2+Xy^2
        I(0,0) += Xxyz2(1) + Xxyz2(2);
        I(1,1) += Xxyz2(0) + Xxyz2(2);
        I(2,2) += Xxyz2(0) + Xxyz2(1);

        I(0,1) -= Xxyz2Xyzx(0);
        I(1,0) -= Xxyz2Xyzx(0);

        I(0,2) -= Xxyz2Xyzx(2);
        I(2,0) -= Xxyz2Xyzx(2);

        I(1,2) -= Xxyz2Xyzx(1);
        I(2,1) -= Xxyz2Xyzx(1);
    }

    // quotient for diagonal entry
    I.diagonal() *= density / 60.0;
    // quotient for off-diagonal entry
    I(0,1) *= density / 120.0;
    I(1,0) *= density / 120.0;
    I(0,2) *= density / 120.0;
    I(2,0) *= density / 120.0;
    I(1,2) *= density / 120.0;
    I(2,1) *= density / 120.0;

    // std::cout << "inertia:\n" << I << std::endl;
}