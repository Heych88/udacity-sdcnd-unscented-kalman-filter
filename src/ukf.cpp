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
  
  // State dimension
  n_x_ = 5;
  // Augmented state dimension
  n_aug_ = 7;
  n_aug_size_ = 2 * n_aug_ + 1;
  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.1); // initalise slightly of the origin

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  //create matrix with predicted sigma points as columns
  predicted_sigma_pts_ = MatrixXd(n_x_, n_aug_size_); 
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
    
    // set the sigma points weights
    weights_ = VectorXd(n_aug_size_);
    weights_(0) = (lambda_)/(lambda_ + n_aug_);
    for (int i=1; i< n_aug_size_; i++) {  //2n+1 weights
      weights_(i) = 0.5/(n_aug_ + lambda_);
    }
    
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
  
  // predict the objects position at the current time step 
  UKF::Prediction(dt);
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // Update the predicted position with the measured sensor data
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    UKF::UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UKF::UpdateRadar(meas_package);
  }
}

void UKF::SigmaPointPrediction(MatrixXd &sigma_pts, MatrixXd &pred_sigma_pts,
    double delta_t) {
  
  pred_sigma_pts.fill(0.0);
  //predict sigma points
  for (int i = 0; i< n_aug_size_; i++)
  {
    //extract values for better readability
    double pos_x = sigma_pts(0,i);
    double pos_y = sigma_pts(1,i);
    double vel = sigma_pts(2,i);
    double yaw = sigma_pts(3,i);
    double yaw_dot = sigma_pts(4,i);
    double a_pos = sigma_pts(5,i);
    double a_yaw_dot = sigma_pts(6,i);

    // check for divide by zero
    if (fabs(yaw_dot) > 0.001) {
        pos_x += vel/yaw_dot * (sin(yaw + yaw_dot * delta_t) - sin(yaw));
        pos_y += vel/yaw_dot * (cos(yaw) - cos(yaw + yaw_dot * delta_t));
    } else {
        pos_x += vel * delta_t * cos(yaw);
        pos_y += vel * delta_t * sin(yaw);
    }

    //write predicted sigma point into right column
    pred_sigma_pts(0,i) = pos_x + 0.5 * a_pos * delta_t * delta_t * cos(yaw);
    pred_sigma_pts(1,i) = pos_y + 0.5 * a_pos * delta_t * delta_t * sin(yaw);
    pred_sigma_pts(2,i) = vel + a_pos * delta_t;
    pred_sigma_pts(3,i) = yaw + yaw_dot * delta_t + 0.5 * a_yaw_dot * 
        delta_t * delta_t;
    pred_sigma_pts(4,i) = yaw_dot + a_yaw_dot * delta_t;
  }
}

void UKF::PredictMeanAndCovariance(VectorXd &x, MatrixXd &P, 
    MatrixXd &pred_sigma_pts, const int yaw_pos=-1) {
  //predict state mean
  x.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) {  //iterate over sigma points
    x += weights_(i) * pred_sigma_pts.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = pred_sigma_pts.col(i) - x;
    if(yaw_pos >= 0) {
      //angle normalization
      while (x_diff(yaw_pos)> M_PI) x_diff(yaw_pos) -= 2.*M_PI;
      while (x_diff(yaw_pos)<-M_PI) x_diff(yaw_pos) += 2.*M_PI;
    }

    P += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // update the state augmentation vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;
  
  // Calculate the augmentation matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  
  MatrixXd sigma_pts = MatrixXd(n_aug_, n_aug_size_);
  sigma_pts.col(0) = x_aug;
  
  //create square root matrix
  MatrixXd sqrt_P = P_aug.llt().matrixL();
  
  // Generate sigma points
  for (int i = 0; i< n_aug_; i++)
  {
    // columns 1 -> n_aug_ = x + sqrt((lambda + n_aug_) * P_) 
    sigma_pts.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * sqrt_P.col(i);
    // columns n_aug_+1 -> 2*n_aug_+1 = x - sqrt((lambda + n_aug_) * P_)
    sigma_pts.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * sqrt_P.col(i);
  }
  
  // Predict the sigma point values for this the time step
  UKF::SigmaPointPrediction(sigma_pts, predicted_sigma_pts_, delta_t);
  // Calculate the mean and covariance of the predicted sigma points
  UKF::PredictMeanAndCovariance(x_, P_, predicted_sigma_pts_, 3);
}

void UKF::UpdateState(VectorXd &x, MatrixXd &P, MatrixXd &pred_sigma_pts, 
    MatrixXd &S, MatrixXd &Zsig, VectorXd &z_pred, const int &n_z, 
    MeasurementPackage &meas_package) {
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_aug_size_; i++) { 

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = pred_sigma_pts.col(i) - x;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //update state mean and covariance matrix
  x += K * z_diff;
  P = P - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  const int n_z = 2; // number of sensor output parameters (px, py)
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_aug_size_);
  
  //transform sigma points into measurement space
  Zsig.row(0) = predicted_sigma_pts_.row(0);
  Zsig.row(1) = predicted_sigma_pts_.row(1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  // Calculate the mean and covariance of the predicted sigma points
  UKF::PredictMeanAndCovariance(z_pred, S, Zsig, -1);

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;
  
  UKF::UpdateState(x_, P_, predicted_sigma_pts_, S, Zsig, z_pred, n_z, 
      meas_package);   
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
  int n_z = 3; // number of sensor output parameters. (range, angle, vel)
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_aug_size_);
  
  //transform sigma points into measurement space
  for (int i = 0; i < n_aug_size_; i++) {  

    // extract values for better readibility
    double p_x = predicted_sigma_pts_(0,i);
    double p_y = predicted_sigma_pts_(1,i);
    double v  = predicted_sigma_pts_(2,i);
    double yaw = predicted_sigma_pts_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y); //ramge
    Zsig(1,i) = atan2(p_y,p_x); //angle
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y); //velocity
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  // Calculate the mean and covariance of the predicted sigma points
  UKF::PredictMeanAndCovariance(z_pred, S, Zsig, 1);

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  UKF::UpdateState(x_, P_, predicted_sigma_pts_, S, Zsig, z_pred, n_z, 
      meas_package);
}