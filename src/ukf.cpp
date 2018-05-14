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
    lambda_ = 3 - n_aug_; // TODO: check with n_aug?

    // initial state vector
    x_ = VectorXd(n_x_);
    x_.fill(0.0);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_.fill(0.0);

    // TODO: tune the process noise parameters std_a_ and std_yawdd_
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1; // 30? Really?

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.5; // 30? Really?

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
    TODO: Complete the initialization. See ukf.h for other member properties.

    Hint: one or more values initialized above might be wildly off...
    */

    is_initialized_ = false;
    time_us_ = 0;

    weights_ = VectorXd(2*n_aug_+1);
    weights_.fill(0.0);

    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2*n_aug_+1; ++i) {
        weights_(i) = 1.0 / (2.0 * (lambda_ + n_aug_)); // was lambda_ / (2 * (lambda_ + n_aug_))
    }

    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_.fill(0.0);

    R_laser_ = MatrixXd(2,2);
    R_radar_ = MatrixXd(3,3);

    R_laser_ << std_laspx_*std_laspx_, 0,
                0, std_laspy_*std_laspy_;

    R_radar_ << std_radr_*std_radr_, 0, 0,
                0, std_radphi_*std_radphi_, 0,
                0, 0, std_radrd_*std_radrd_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    TODO: Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    if(!is_initialized_) {
        // initialize the state vector x and state covariance matrix P with appropriate values.
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            // LASER
            double px = meas_package.raw_measurements_(0);
            double py = meas_package.raw_measurements_(1);
            x_ << px, py, 4.16, 0, 0; // 4.16m/s = 15km/h predicted speed
            P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
                    0, std_laspy_*std_laspy_, 0, 0, 0,
                    0, 0, 1, 0, 0,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1;

        }
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            // RADAR
            double px = meas_package.raw_measurements_(0) * std::cos(meas_package.raw_measurements_(1));
            double py = meas_package.raw_measurements_(0) * std::sin(meas_package.raw_measurements_(1));

            double rho, phi, rho_dot;
            rho = meas_package.raw_measurements_(0);
            phi = meas_package.raw_measurements_(1);
            rho_dot = meas_package.raw_measurements_(2);
            double v = std::sqrt(std::pow(rho_dot*std::cos(phi), 2) +
                                 std::pow(rho_dot*std::sin(phi), 2));
            x_ << px, py, v, 0, 0;
            P_ << std_radr_*std_radr_, 0, 0, 0, 0,
                    0, std_radr_*std_radr_, 0, 0, 0,
                    0, 0, 1, 0, 0,
                    0, 0, 0, std_radphi_, 0,
                    0, 0, 0, 0, std_radphi_;
        }
        // Initialization is done only once
        is_initialized_ = true;
    }
    else {
        double delta_t = (double)(meas_package.timestamp_ - time_us_) / 1000000.0; // time in seconds
        std::cout << delta_t << std::endl;
        this->Prediction(delta_t);
        // update last timestamp for step k+1
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            this->UpdateLidar(meas_package);
        }
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            this->UpdateRadar(meas_package);
        }
    }
    // set time of last measurement
    time_us_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
    TODO: Complete this function! Estimate the object's location. Modify the state
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
    Xsig_aug.fill(0.0);

    MatrixXd A = P_aug.llt().matrixL();

    // Calculating sigma points
    Xsig_aug.col(0) = x_aug; //

    for(int i=0; i < n_aug_; ++i) {
        Xsig_aug.col(i+1) = x_aug + std::sqrt(lambda_ + n_aug_) * A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - std::sqrt(lambda_ + n_aug_) * A.col(i);
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
        if(Xsig_aug.col(i)(4) <= 0.001) {
            // calculate px, py when yaw rate is zero
            sigma_pt(0) = px + v * std::cos(psi)*delta_t + std::pow(delta_t, 2)/2 * std::cos(psi)* ni_acc;
            sigma_pt(1) = py + v * std::sin(psi)*delta_t + std::pow(delta_t, 2)/2 * std::sin(psi)* ni_acc;
        }
        else {
            // calculate px, py otherwise
            sigma_pt(0) = px + v/yaw_rate * (std::sin(psi + yaw_rate*delta_t) - std::sin(psi)) + std::pow(delta_t, 2)/2 * std::cos(psi)* ni_acc;
            sigma_pt(1) = py + v/yaw_rate * (-std::cos(psi + yaw_rate*delta_t) + std::cos(psi)) + std::pow(delta_t, 2)/2 * std::sin(psi)* ni_acc;
        }

        // calculate v, psi, yaw_rate
        sigma_pt(2) = v + delta_t * ni_acc;
        sigma_pt(3) = psi + yaw_rate * delta_t + std::pow(delta_t, 2)/2 * ni_yaw;
        normalizeAngle(sigma_pt, 3);
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
    normalizeAngle(x, 3);

    // predict state covariance matrix
    for (int i = 0; i < 2*n_aug_+1; ++i) {

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x;
        normalizeAngle(x_diff, 3);

        P += weights_(i) * x_diff * x_diff.transpose();
    }

    // Just a possibility for the future: save x, P in other variables to use in a more complex scenario
    // where the predicted position should be accessible anytime from the main program (for example, in SLAM scenarios)
    x_ = x;
    P_ = P;
    std::cout << "v       = " << x_(2) << std::endl;
    std::cout << "phi_dot = " << x_(4) << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    TODO: Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */

    // TODO: CHECK FOR ERRORS!
    VectorXd z = meas_package.raw_measurements_;
    int n_z = 2;

    double px, py; // v, yaw_rate, vx, vy

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
    for (int i=0; i < 2*n_aug_+1; ++i) {  //2n+1 sigma points
        Zsig(0,i) = Xsig_pred_(0,i);  // px
        Zsig(1,i) = Xsig_pred_(1,i);  // py
    }

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; ++i) {
        z_pred += weights_(i) * Zsig.col(i);
    }

    //innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);

    for (int i = 0; i < 2*n_aug_+1; ++i) {  //2n+1 sigma points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        S += weights_(i) * z_diff * z_diff.transpose();
    }
    // add measurement noise covariance matrix
    S += R_laser_;


    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  //2n+1 sigma points

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        normalizeAngle(x_diff, 3);

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();

    // calculate NIS
    double nis_l = z_diff.transpose() * S.inverse() * z_diff;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    TODO: Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */

    VectorXd z = meas_package.raw_measurements_;
    int n_z = 3; // number of dimensions for the Radar Sensor (rho, phi, r_dot)

    double px, py, v, yaw_angle, vx, vy;

    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
    Zsig.fill(0.0);

    // transform sigma points into measurement space
    for (int i=0; i < 2*n_aug_+1; ++i) {  //2n+1 sigma points

        // extract values for better readability
        px = Xsig_pred_(0,i);
        py = Xsig_pred_(1,i);
        v  = Xsig_pred_(2,i);
        yaw_angle = Xsig_pred_(3,i);

        vx = cos(yaw_angle)*v;
        vy = sin(yaw_angle)*v;

        // measurement model
        Zsig(0,i) = sqrt(px*px + py*py);              //r
        Zsig(1,i) = atan2(py,px);                     //phi
        Zsig(2,i) = (px*vx + py*vy ) / Zsig(0,i);     //r_dot
    }

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; ++i) {
        z_pred += weights_(i) * Zsig.col(i);
    }

    // innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);

    for (int i = 0; i < 2*n_aug_+1; ++i) {  // 2n+1 sigma points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        normalizeAngle(z_diff, 1);

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    S += R_radar_;

    // create cross correlation matrix Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2*n_aug_+1; ++i) {  //2n+1 sigma points

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        normalizeAngle(z_diff, 1);

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        normalizeAngle(x_diff, 3);
        Tc += weights_(i) * x_diff * z_diff.transpose();

    }
    std::cout << "Bingo" << std::endl;

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    normalizeAngle(z_diff, 1);

    //update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();

    // calculate NIS
    double nis_r = z_diff.transpose() * S.inverse() * z_diff;
}

void UKF::normalizeAngle(VectorXd &z_diff, int index) {
    //angle normalization check
    if (z_diff(index) > M_PI || z_diff(index) < -M_PI)
        z_diff(index) = std::fmod(z_diff(index), M_PI);
}
