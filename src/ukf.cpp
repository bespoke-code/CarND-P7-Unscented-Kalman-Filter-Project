#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_x_;

    // initial state vector
    x_ = VectorXd(n_x_);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.3; // 30? Really?

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.3; // 30? Really?

    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

    /**
    TODO:

    Complete the initialization. See ukf.h for other member properties.

    Hint: one or more values initialized above might be wildly off...
    */

    weights_ = VectorXd(2*n_aug_+1);
    weights_.fill(0.0);

    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2*n_aug_+1; ++i) {
        weights_(i) = lambda_ / (2 * (lambda_ + n_aug_));
    }

    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_.fill(0.0);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
        // LASER
    }
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // LASER
    }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
    TODO:

    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */

    // augment X
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(n_x_) = x_;
    x_aug(n_x_) = 0;
    x_aug(n_x_+1) = 0;

    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    // adding the Q matrix to the augmented P
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;

    // Augmented sigma points matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    // Predicted sigma points according to the model

    MatrixXd A = P_aug.llt().matrixL();

    // Calculating sigma points
    Xsig_aug.col(0) = x_;

    for(int i=0; i < n_x_; i++) {
        Xsig_aug.col(i+1) = x_ + std::sqrt(lambda_ + n_aug_) * A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_ - std::sqrt(lambda_ + n_aug_) * A.col(i);
    }

    // Predict the sigma points after the process model
    double px, py, v, psi, yaw_rate, ni_acc, ni_yaw;
    for(int i=0; i < 2*n_aug_+1; ++i) {
        px = Xsig_aug.col(i)[0];
        py = Xsig_aug.col(i)[1];
        v = Xsig_aug.col(i)[2];
        psi = Xsig_aug.col(i)[3];
        yaw_rate = Xsig_aug.col(i)[4];
        ni_acc = Xsig_aug.col(i)[5];
        ni_yaw = Xsig_aug.col(i)[6];

        VectorXd sigma_pt = VectorXd(n_x_);

        //avoid division by zero
        // if yaw rate is zero
        if(Xsig_aug.col(i)(4) == 0) {
            // calculate px, py when yaw rate is zero
            sigma_pt(0) = px + v * cos(psi)*delta_t + pow(delta_t, 2)/2 * cos(psi)* ni_acc;
            sigma_pt(1) = py + v * sin(psi)*delta_t + pow(delta_t, 2)/2 * sin(psi)* ni_acc;
        }
        else {
            // calculate px, py otherwise
            sigma_pt(0) = px + v/yaw_rate * (sin(psi + yaw_rate*delta_t) - sin(psi)) + pow(delta_t, 2)/2 * cos(psi)* ni_acc;
            sigma_pt(1) = py + v/yaw_rate * (-cos(psi + yaw_rate*delta_t) + cos(psi)) + pow(delta_t, 2)/2 * sin(psi)* ni_acc;
        }

        // calculate v, psi, yaw_rate
        sigma_pt(2) = v + delta_t * ni_acc;
        sigma_pt(3) = psi + yaw_rate * delta_t + pow(delta_t, 2)/2 * ni_yaw;
        sigma_pt(4) = yaw_rate + delta_t * ni_yaw;

        //write predicted sigma points into the corresponding column
        Xsig_pred_.col(i) = sigma_pt;
    }

    // Use the predicted sigma points to calculate the predicted state's mean and covariance
    //create vector for predicted state
    VectorXd x = VectorXd(n_x_);
    x.fill(0.0);

    //create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x_, n_x_);
    P.fill(0.0);

    //predict the state mean
    for(int i=0; i < 2*n_aug_+1; ++i) {
        x += weights_(i) * Xsig_pred_.col(i);
    }

    // predict state covariance matrix
    for (int i = 0; i < 2*n_aug_+1; ++i) {

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        //angle normalization check
        if (x_diff(3) > M_PI || x_diff(3) < -M_PI)
            x_diff(3) = std::fmod(x_diff(3), M_PI);
        P += weights_(i) * x_diff * x_diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */
}
