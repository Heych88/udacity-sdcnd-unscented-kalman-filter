#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * @param {n_x} - number of state parameters
 * @param {n_aug} - number of state + noise parameters for the augmented state
 * @param {use_laser} - Use lidar data in the calculation
 * @param {use_radar} - Use radar data in the calculation
 * @param {std_a} - process noise due to change in acceleration
 * @param {std_yawdd} - process noise due to change in measured angle
 * @param {std_laspx} - lidar measurement noise in the x plane
 * @param {std_laspy} - lidar measurement noise in the y plane
 * @param {std_radr} - radar noise in the distance measurement
 * @param {std_radphi} - radar noise in the angle measurement
 * @param {std_radrd} - radar noise in the angle acceleration measurement
 */
UKF::UKF(int n_x, int n_aug, bool use_laser, bool use_radar, bool use_nis) {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = use_laser;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = use_radar;
  // if false NIS will not be calculated.
  use_nis_ = use_nis;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.4;
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
  n_x_ = n_x;
  // Augmented state dimension
  n_aug_ = n_aug;
  // Number of sigma points
  n_aug_size_ = 2 * n_aug_ + 1;
  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.1); // initalise slightly of the origin
  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  // create matrix with predicted sigma points as columns
  predicted_sigma_pts_ = MatrixXd(n_x_, n_aug_size_); 
  
  // init count forNormalised Innovation Squared
  update_count_ = 0;
  nis_thresh_count_ = 0;
}

UKF::~UKF() {}

/**
 * Updates the Uncented Kalman Filter with the measured data.
 * Call this function to update the objects position.
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
      x_.head(2) = meas_package.raw_measurements_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
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

/**
 * Predicts the previous sigma points at the current time step.
 * @param {sigma_pts} - matrix containing the previous sigma points
 * @param {pred_sigma_pts} - matrix containing the predicted sigma points 
 * @param {delta_t} - change in time between UKF updates 
 */
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
    double a_pos = sigma_pts(5,i); // acceleration process noise
    double a_yaw_dot = sigma_pts(6,i); // change in angle process noise

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

/**
 * Predicts the mean and covariance from the provided sigma points matrix.
 * @param {x} - the output vector with the predicted mean values
 * @param {P} - the output Matrix with the predicted covariance values
 * @param {pred_sigma_pts} - the matrix containing the sigma points 
 * @param {yaw_pos} - the vector position of the stored yaw value. -1 if yaw 
 *                    not present in the sigma points. 
 */
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
  
  /*****************************************************************************
   *  Generate sigma points
   ****************************************************************************/
  for (int i = 0; i< n_aug_; i++)
  {
    // columns 1 -> n_aug_ = x + sqrt((lambda + n_aug_) * P_) 
    sigma_pts.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * sqrt_P.col(i);
    // columns n_aug_+1 -> 2*n_aug_+1 = x - sqrt((lambda + n_aug_) * P_)
    sigma_pts.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * sqrt_P.col(i);
  }
  
  /*****************************************************************************
   *  Predict the sigma point values for this the time step
   ****************************************************************************/
  UKF::SigmaPointPrediction(sigma_pts, predicted_sigma_pts_, delta_t);
  
  /*****************************************************************************
   *  Calculate the mean and covariance of the predicted sigma points
   ****************************************************************************/
  UKF::PredictMeanAndCovariance(x_, P_, predicted_sigma_pts_, 3);
}

/**
 * Updates the mean and covariance with a proprtion of the predicted and 
 * senor measured position.
 * @param {x} - output vector with the predicted mean values
 * @param {P} - output Matrix with the predicted covariance values
 * @param {S} - measurement covariance matrix 
 * @param {Zsig} - Matrix with the sensor predicted sigma points
 * @param {z_pred} - Vector with the current predicted mean values
 * @param {n_z} - Number of measured sensor parameters, lidar = 2, Radar = 3. 
 * @param {meas_package} - Measurement sensor class
 */
void UKF::UpdateState(VectorXd &x, MatrixXd &P, MatrixXd &pred_sigma_pts, 
    MatrixXd &S, MatrixXd &Zsig, VectorXd &z_pred, const int &n_z, 
    MeasurementPackage &meas_package) {
  // TODO: Check for NaN and inf in the calculations.
  
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
  
  if(use_nis_) {
    // check the error using Normalised Innovation Squared
    UKF::NISState(S, z_diff, n_z);
  }

  //update state mean and covariance matrix
  x += K * z_diff;
  P = P - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  const int n_z = 2; // number of sensor output parameters (pos_x, pos_y)
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_aug_size_);

  /*****************************************************************************
   *  convert predicted sigma points into measurement space
   ****************************************************************************/
  Zsig.row(0) = predicted_sigma_pts_.row(0);
  Zsig.row(1) = predicted_sigma_pts_.row(1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  /*****************************************************************************
   *  Calculate the mean and covariance
   ****************************************************************************/
  UKF::PredictMeanAndCovariance(z_pred, S, Zsig, -1);

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;
  
  /*****************************************************************************
   *  Update the unscented kalman filter
   ****************************************************************************/
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
  
  /*****************************************************************************
   *  transform predicted sigma points into measurement space
   ****************************************************************************/
  for (int i = 0; i < n_aug_size_; i++) {  
    // extract values for better readibility
    double pos_x = predicted_sigma_pts_(0,i);
    double pos_y = predicted_sigma_pts_(1,i);
    double vel  = predicted_sigma_pts_(2,i);
    double yaw = predicted_sigma_pts_(3,i);

    // measurement model
    // range = sqrt(px^2 + py^2)
    Zsig(0,i) = sqrt(fabs(pos_x*pos_x + pos_y*pos_y));
    Zsig(1,i) = atan2(pos_y,pos_x); //angle
    // velocity = (px*cos(ψ)v+py*sin(ψ)v)/sqrt(px^2 + py^2)
    Zsig(2,i) = (pos_x*cos(yaw)*vel + pos_y*sin(yaw)*vel) / Zsig(0,i);
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  /*****************************************************************************
   *  Calculate the mean and covariance
   ****************************************************************************/
  UKF::PredictMeanAndCovariance(z_pred, S, Zsig, 1);

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  /*****************************************************************************
   *  Update the unscented kalman filter
   ****************************************************************************/
  UKF::UpdateState(x_, P_, predicted_sigma_pts_, S, Zsig, z_pred, n_z, 
      meas_package);
}

/**
 * Normalised Innovation Squared for checking noise parameter values,
 * @param {S} - measurement covariance matrix 
 * @param {z_diff} - difference between predicted and measured values
 * @param {n_z} - Number of measured sensor parameters, lidar = 2, Radar = 3.
 */
void UKF::NISState(MatrixXd &S, VectorXd &z_diff, const int &n_z) {
  
  update_count_ += 1;
  
  // compare error with the Chi Squared 0.050 distribution
  if(n_z == 2) { // Lidar measurement check
    if(z_diff.transpose() * S.inverse() * z_diff > 5.991) {
      nis_thresh_count_ += 1;
    }
  } else if(n_z == 3) { // Radar measurement check
    if(z_diff.transpose() * S.inverse() * z_diff > 7.815) {
      nis_thresh_count_ += 1;
    }
  }
  
  std::cout << "Chi Squared 0.05 error = ";
  std::cout << (float(nis_thresh_count_) / float(update_count_)) * 100.0 << "%";
  std::cout << std::endl;  
}