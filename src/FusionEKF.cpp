#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0,      0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0,      0,
                0,    0.0009, 0,
                0,    0,      0.09;

    // process noises - laser
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    // initial transition matrix
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    // initial covariance Matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1,  0,  0,  0,
               0,  1,  0,  0,
               0,  0,  1,  0,
               0,  0,  0,  1;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /**
     * Initialization
     */
    if (!is_initialized_) {

        cout << "Initialization" << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

            cout << "First mesurment is RADAR" << endl;

            double ro     = measurement_pack.raw_measurements_[0];
            double theta  = measurement_pack.raw_measurements_[1];
            double ro_dot = measurement_pack.raw_measurements_[2];
            double phi = theta;

            // Compute cordinate
            double x  = ro * cos(phi);     // px
            double y  = ro * sin(phi);     // py
            double vx = ro_dot * cos(phi); // vx
            double vy = ro_dot * sin(phi); // vy

            ekf_.x_ << x, y, vx, vy;

        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

            cout << "First mesurment is LASER" << endl;

            double px = measurement_pack.raw_measurements_[0];
            double py = measurement_pack.raw_measurements_[1];

            ekf_.x_ << px, py, 0, 0;

        }

        // Save initial timestamp
        previous_timestamp_ = measurement_pack.timestamp_;

        // Change done flag
        is_initialized_ = true;

        cout << "Done initializing: " << ekf_.x_ << endl;
        return;
    }

    /**
     * Prediction
     */

    cout << "Prediction step" << endl;

    // Compute dt [s] from difference timestamp [ns]
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    // Modify state transition matrix
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    // Compute the Noise Covariance Matrix
    double noise_ax = 9.0;
    double noise_ay = 9.0;
    double dt_2     = dt   * dt;
    double dt_3     = dt_2 * dt;
    double dt_4     = dt_3 * dt;
    double dt_4d4   = dt_4 / 4;
    double dt_3d2   = dt_3 / 2;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4d4*noise_ax, 0,               dt_3d2*noise_ax, 0,
               0,               dt_4d4*noise_ay, 0,               dt_3d2*noise_ay,
               dt_3d2*noise_ax, 0,               dt_2*noise_ax,   0,
               0,               dt_3d2*noise_ay, 0,               dt_2*noise_ay;
    
    ekf_.Predict();

    cout << "Done prediction:" << ekf_.x_ << endl;

    /**
     * Update
     */

    cout << "Update step" << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        cout << "RADAR measurment update" << endl;
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        cout << "LASER measurment update" << endl;
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ =" << ekf_.x_ << endl;
    cout << "P_ =" << ekf_.P_ << endl;
}
