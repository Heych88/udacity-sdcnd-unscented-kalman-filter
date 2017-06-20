#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;
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

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    previous_timestamp_ = meas_package.timestamp_;
    
    // State dimension
    n_x_ = 5;
    // Augmented state dimension
    n_aug_ = 7;
    n_aug_size_ = 2 * n_aug_ + 1;
    // Sigma point spreading parameter
    lambda_ = 3 - n_aug_;
    
    // initial state vector
    x_ = VectorXd(n_x_);
    x_.fill(1.0);
    x_aug_ = VectorXd(n_aug_);
    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    // Init covariance matrix to be large so as to rely on the measured data to 
    // start with
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1; 
    P_aug_ = MatrixXd(n_aug_, n_aug_);
    sigma_pts_ = MatrixXd(n_aug_, n_aug_size_);
    //create matrix with predicted sigma points as columns
    predicted_sigma_pts_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    // Q prediction noise matrix
    Q_ = MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
    Q_ << std_a_*std_a_, 0.0,
          0.0, std_yawdd_*std_yawdd_;  
    
    ///* Weights of sigma points
    weights_ = VectorXd(n_aug_size_);
    for(int i=0; i < n_aug_size_; i++){
      if(i < 1){
        weights_(i) = lambda_ / (lambda_ + n_aug_);
      } else {
        weights_(i) = 1 / (2 * (lambda_ + n_aug_));
      }
    }
        
    // sigma point matrix for radar measurement
    n_lidar_ = 2;
    // Measurment noise
    R_lidar_ = MatrixXd(n_lidar_, n_lidar_);
    R_lidar_ << std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      //Initialize state.  
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  // dt - expressed in seconds
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; 
  previous_timestamp_ = meas_package.timestamp_;
  
  UKF::Prediction(dt);
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // Update the predicted position with the measured sensor data
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    // Use normal kalman filter on measured data of the same coordinates system
    UKF::UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UKF::UpdateRadar(meas_package);
  }

  // print the output
  std::cout << "x_ = " << x_ << std::endl;
  std::cout << "P_ = " << P_ << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // update the state augmentation vector
  x_aug_ << x_, 0, 0;
  // Calculate the augmentation matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_.bottomRightCorner(n_aug_-n_x_, n_aug_-n_x_) = Q_;
  
  sigma_pts_.fill(0.0);
  sigma_pts_.col(0) = x_aug_;
  //create square root matrix
  MatrixXd sqrt_P_ = P_aug_.llt().matrixL();
  for(int i=0; i < n_aug_; i++) {
    // columns 1 -> n_x_ = x + sqrt((lambda + n_x_) * P_) 
    sigma_pts_.col(i+1) = x_aug_ + sqrt(abs(lambda_ + n_aug_)) * sqrt_P_.col(i);
    // columns n_x_+1 -> 2*n_x_+1 = x - sqrt((lambda + n_x_) * P_)
    sigma_pts_.col(i+1+n_aug_) = x_aug_ - sqrt(abs(lambda_ + n_aug_)) * sqrt_P_.col(i);
  }
  
  //predict sigma points
  for(int i=0; i < n_aug_size_; i++){
      double pos_x = sigma_pts_.col(i)[0];
      double pos_y = sigma_pts_.col(i)[1];
      double vel = sigma_pts_.col(i)[2];
      double thi = sigma_pts_.col(i)[3];
      double thi_dot = sigma_pts_.col(i)[4];
      double a_x = sigma_pts_.col(i)[5];
      double a_thi = sigma_pts_.col(i)[6];
      
      double pred_pos_x, pred_pos_y;
    // Check for angle rate of change and prevent a divide by zero
      if(thi_dot > 0.0001){
          pred_pos_x = (vel/thi_dot) * (sin(thi + thi_dot * delta_t) - sin(thi));
          pred_pos_y = (vel/thi_dot) * (-cos(thi + thi_dot * delta_t) + cos(thi));
      } else {
          pred_pos_x = vel * (sin(thi) - sin(thi));
          pred_pos_y = vel * (-cos(thi) + cos(thi));
      }
        
      double noise_px = 0.5 * (delta_t * delta_t) * cos(thi) *a_x;
      predicted_sigma_pts_.col(i)[0] = pos_x + pred_pos_x + noise_px;

      double noise_py = 0.5 * (delta_t * delta_t) * sin(thi) *a_x;
      predicted_sigma_pts_.col(i)[1] = pos_y + pred_pos_y + noise_py;

      predicted_sigma_pts_.col(i)[2] = vel + delta_t * a_x;
      predicted_sigma_pts_.col(i)[3] = thi + thi_dot * delta_t + 0.5 * (delta_t * delta_t) * a_thi;
      predicted_sigma_pts_.col(i)[4] = thi_dot + delta_t * a_thi;
  }
  
  //predict state mean
  x_.fill(0.0);
  for(int i=0; i < 2 * n_aug_ + 1; i++){
      x_ += predicted_sigma_pts_.col(i) * weights_(i);
  }
  
  //predict state covariance matrix
  P_.fill(0.0);
  for(int i=0; i < n_aug_size_; i++) {
    
    VectorXd diff = predicted_sigma_pts_.col(i) - x_;
    //angle normalization
    while (diff(3)> M_PI) diff(3) -= 2.*M_PI;
    while (diff(3)<-M_PI) diff(3) += 2.*M_PI;

    P_ += weights_(i) * diff * diff.transpose();
  }  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  std::cout << "UKF::UpdateLidar Start" << std::endl;
  
    
 
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  std::cout << "UKF::UpdateRadar Start" << std::endl;
}
