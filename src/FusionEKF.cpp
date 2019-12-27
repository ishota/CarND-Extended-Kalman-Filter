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

    // initial covariance Matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 10,  0,  0,   0,
                0,  10, 0,   0,
                0,  0,  100, 0,
                0,  0,  0,   100;

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
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

        cout << "Initialization" << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

            cout << "First mesurment is RADAR" << endl;

            float ro     = measurement_pack.raw_measurements_[0];
            float theta  = measurement_pack.raw_measurements_[1];
            float ro_dot = measurement_pack.raw_measurements_[2];

            // Normalize theta to phi in range [-pi, pi]
            float phi = theta;
            while (phi >  M_PI) phi -= 2.0 * M_PI;
            while (phi < -M_PI) phi += 2.0 * M_PI;

            // Compute cordinate
            float x  = ro * cos(phi);
            float y  = ro * sin(phi);
            float vx = ro_dot * cos(phi);
            float vy = ro_dot * sin(phi);

            ekf_.x_ << x, y, vx, vy;

        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

            cout << "First mesurment is LASER" << endl;

            float px = measurement_pack.raw_measurements_[0];
            float py = measurement_pack.raw_measurements_[1];

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
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    // Update state transition matrix
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1,  0,  dt, 0,
               0,  1,  0,  dt,
               0,  0,  1,  0,
               0,  0,  0,  1;
    
    // Compute the Noise Covariance Matrix
    ekf_.Q_ = MatrixXd(4, 4);
    float noise_ax = 9.0;
    float noise_ay = 9.0;
    float dt_2     = dt * dt;
    float dt_3     = dt * dt_2;
    float dt_4     = dt * dt_3;
    float dt_4d4   = dt_4 / 4;
    float dt_3d2   = dt_3 / 2;

    ekf_.Q_ << dt_4d4*noise_ax, 0,               dt_3d2*noise_ax, 0,
               0,               dt_4d4*noise_ay, 0,               dt_3d2*noise_ay,
               dt_3d2*noise_ax, 0,               dt_2*noise_ax,   0,
               0,               dt_3d2*noise_ay, 0,               dt_2*noise_ay;
    
    ekf_.Predict();

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
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
