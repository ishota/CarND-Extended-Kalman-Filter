# Extended Kalman Filter Project Starter Code

[flow]: ./Pics/sensorflow.png "flow"
[result1]: ./Pics/result_dataset1.png "result1"
[result2]: ./Pics/result_dataset2.png "result2"

Self-Driving Car Engineer Nanodegree Program

In this project, A kalman filter estimates the state of a moving object of interest with noisy lidar and radar measurments.
I improved the robustness against large deviations of observations by using Bayesian filters.

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases).

This repository includes two files that can be used to set up and install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. Please see the uWebSocketIO Starter Guide page in the classroom within the EKF Project lesson for the required version and installation scripts.

Once the install for uWebSocketIO is complete, the main program can be built and run by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ExtendedKF

**INPUT**: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurement that the simulator observed (either lidar or radar)

**OUTPUT**: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x

["estimate_y"] <= kalman filter estimated position y

["rmse_x"] <= root mean squared errora of position x

["rmse_y"] <= root mean squared errora of position y 

["rmse_vx"] <= root mean squared errora of Velocity in the direction of travel

["rmse_vy"] <= root mean squared errora of Velocity in the direction of lateral

---

## Other Important Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make` 
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./ExtendedKF `

## Program flow

![alt text][flow]

### initialize EKF matrices
Initialize metrices related to the Extended Kalman filter.
The metrices to be initialized are as flollows.

* R_laser: measurement covariance matrix - laser
* R_radar: measurement covariance matrix - radar
* H_laser: process noises - laser
* P: initial covariance matrix

### FIRST MEASUREMENT? -> YES
If input is the first observation, initialize the following values.

* x: state of a car. [position x, position y, velocity x, velocity y]
* save initial update timestamp

If the first observation was a laser, velocity is initialized to 0.

### FIRST MEASUREMENT? -> No
Repeat the "Predict" and "Update" after the first observation.

#### Predict

* update the elapsed time since the last update to the present
* compute F: transition matrix
* compute Q: process covariance matrix
* predict x and P

#### Update

**LASER**
Update with the Kalman filter for linear state transition
* H: measurment matrix for linear state transition

**RADAR**
Update with the EKF:Extended Kalman filter for non-linear state transition
* h: measurment matrix for non-linear state transition usint by Jacobian matrix

After computing the observation residual y, update the state estimate based on the Kalman gain.
* K: kalman gain

## **Proposed Method**: Bayesian filter

The Kalman filter can deal with some steay-observation errors.
However, an impossibly large error may occur in the RADAR and LASER measurment observation value as compared with the steady errors.
If this error is directly incorporated into the Kalman filter, the filter estimation becomes unstable, and in the worst case, the filter may diverge.

Implemented bayesian filter was introduced to prevent such outliers in the observed values.
It calculates the statistical distance (Maharanobis distance) between the predicted position and updated position and thresholds the distance according to the Chi-square probability distribution.

```cpp
    // calculate Mahalanobis's Distance
    MatrixXd Xp = x_;
    MatrixXd Xz = x_ + (K * y);
    MatrixXd Pi = P_.inverse();
    MatrixXd chi = (Xp - Xz).transpose() * Pi * (Xp - Xz);

    // degrees of freedom is 4
    // significance level is 5%
    if (chi(0, 0) > 9.49) {
        is_update = false;
    }else {
        is_update = true;
    }

    // new estimated state
    if (is_update) {
        x_ = x_ + (K * y);
        P_ = (I - K * H_) * P_;
    }
```

In this project case, a degrees of freedom is 4.
If it is determined that the distribution of the estimated value is out of the distribution at the significance level of 5%, the update with the observation value at that time is not perfomed.

## Simulator Details And Outcome

The simulator is a birds eye view of a car.
The following two figures show the results for datasets 1 and 2, respectively.
The green triangles are the estimated position from the kalman filter with bayesian filter, the red circles represent laser measurements, and blue circles are radar.

As can be seen below both RMSE results, proposed method's result underpassed threshold [0.11, 0.11, 0.52, 0.52] for this project rubric.

![alt text][result1]
![alt text][result2]

## Conclustion
In this project, the position of a moving car was estimated using the Kalman filter.
When the input observation result was an outlier for the estimation result, the estimation result sometimes become abnormal.
On the other hand, by using the bayesian filter using the distribution of the estimation results, it was detected as an abnormal value, and the estimation result was stabilized.
As a result, the RMSE was able to meet the criteria for the two datasets.