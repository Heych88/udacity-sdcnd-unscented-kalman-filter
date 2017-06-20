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
  use_radar_ = true;
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;
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
    x_.fill(0.1);
    
    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    // Init covariance matrix to be large so as to rely on the measured data to 
    // start with
    P_.fill(0.0);
    /*P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1; */
    //P_aug_ = MatrixXd(n_aug_, n_aug_);
    //sigma_pts_ = MatrixXd(n_aug_, n_aug_size_);
    //create matrix with predicted sigma points as columns
    predicted_sigma_pts_ = MatrixXd(n_x_, n_aug_size_);
    
    // Q prediction noise matrix
    //Q_ = MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
    //Q_ << std_a_*std_a_, 0.0,
    //      0.0, std_yawdd_*std_yawdd_;  
    
    ///* Weights of sigma points
    // set weights
    double weight_0 = (lambda_)/(lambda_ + n_aug_);
    weights_ = VectorXd(n_aug_size_);
    weights_(0) = weight_0;
    for (int i=1; i< n_aug_size_; i++) {  //2n+1 weights
      double weight_ = 0.5/(n_aug_ + lambda_);
      weights_(i) = weight_;
    }
        
    // sigma point matrix for radar measurement
    n_lidar_ = 2;
    n_radar_ = 3;
    // Measurment noise
    

    //if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    //  //Initialize state.  
    //  x_(0) = meas_package.raw_measurements_[0];
    //  x_(1) = meas_package.raw_measurements_[1];
    //}

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
  VectorXd x_aug_ = VectorXd(n_aug_);
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
  
  // Calculate the augmentation matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  
  MatrixXd sigma_pts_ = MatrixXd(n_aug_, n_aug_size_);
  sigma_pts_.col(0) = x_aug_;
  
  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();
  for (int i = 0; i< n_aug_; i++)
  {
    sigma_pts_.col(i+1)       = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    sigma_pts_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  predicted_sigma_pts_.fill(0.0);
  //predict sigma points
  for (int i = 0; i< n_aug_size_; i++)
  {
    //extract values for better readability
    double p_x = sigma_pts_(0,i);
    double p_y = sigma_pts_(1,i);
    double v = sigma_pts_(2,i);
    double yaw = sigma_pts_(3,i);
    double yawd = sigma_pts_(4,i);
    double nu_a = sigma_pts_(5,i);
    double nu_yawdd = sigma_pts_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    predicted_sigma_pts_(0,i) = px_p;
    predicted_sigma_pts_(1,i) = py_p;
    predicted_sigma_pts_(2,i) = v_p;
    predicted_sigma_pts_(3,i) = yaw_p;
    predicted_sigma_pts_(4,i) = yawd_p;
  }
  
  std::cout << "predicted_sigma_pts_" << std::endl;
  std::cout << lambda_ + n_aug_ << std::endl;
  
  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * predicted_sigma_pts_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = predicted_sigma_pts_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  std::cout << "UKF::UpdateLidar Start" << std::endl;
  
  int n_z = 2;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_aug_size_);
  
  //transform sigma points into measurement space
  for (int i = 0; i < n_aug_size_; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = predicted_sigma_pts_(0,i);
    double p_y = predicted_sigma_pts_(1,i);
    
    // measurement model
    Zsig(0,i) = p_x;                        //r
    Zsig(1,i) = p_y;                                 //phi
    
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_aug_size_; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = predicted_sigma_pts_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  std::cout << "UKF::UpdateRadar Start" << std::endl;
  
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_aug_size_);
  
  //transform sigma points into measurement space
  for (int i = 0; i < n_aug_size_; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = predicted_sigma_pts_(0,i);
    double p_y = predicted_sigma_pts_(1,i);
    double v  = predicted_sigma_pts_(2,i);
    double yaw = predicted_sigma_pts_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_aug_size_; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = predicted_sigma_pts_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
