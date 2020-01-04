#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in; // object state
    P_ = P_in; // object covariance matrix
    F_ = F_in; // state transition matrix
    H_ = H_in; // measurment matrix
    R_ = R_in; // measurement covariance matrix
    Q_ = Q_in; // process covariance matrix
}

void KalmanFilter::Predict() {

    MatrixXd Ft = F_.transpose();

    x_ = F_ * x_;
    P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {

    // innovation
    VectorXd y = z - H_ * x_;

    CommonUpdate(y);

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

    VectorXd h = VectorXd(3);

    double ro     = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
    double theta  = atan2(x_(1), x_(0));
    double ro_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / ro;

    h << ro, theta, ro_dot;

    // innovation
    VectorXd y = z - h;

    CommonUpdate(y);
}

void KalmanFilter::CommonUpdate(const VectorXd & y) {

    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    MatrixXd Ht  = H_.transpose();
    MatrixXd S   = H_ * P_ * Ht + R_;   // innovation covariance matrix
    MatrixXd Si  = S.inverse();
    MatrixXd K   = P_ * Ht * Si;        // kalman gain

    // calculate Mahalanobis's Distance
    MatrixXd Xp = x_;
    MatrixXd Xz = x_ + (K * y);
    MatrixXd Pi = P_.inverse();
    MatrixXd chi = (Xp - Xz).transpose() * Pi * (Xp - Xz);

    cout << "chi: " << chi << endl;
    if (chi(0, 0) > 4) {
        is_update = false;
    }else {
        is_update = true;
    }

    // new estimated state
    if (is_update) {
        x_ = x_ + (K * y);
        P_ = (I - K * H_) * P_;
    }

}