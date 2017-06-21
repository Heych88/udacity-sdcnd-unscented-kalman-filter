#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {  
public:

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;
  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  ///* state covariance matrix
  MatrixXd P_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;
  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;
  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;
  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;
  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;
  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;
  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* State dimension
  int n_x_;
  ///* Augmented state dimension
  int n_aug_;

  /**
 * Unscented Kalman filter Constructor
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
  UKF(int n_x, int n_aug, bool use_laser=true, bool use_radar=true, 
    double std_a=2.1, double std_yawdd=2.1, double std_laspx=0.15, 
    double std_laspy=0.15, double std_radr=0.3, double std_radphi=0.03, 
    double std_radrd=0.3);

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

private:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;
  
  //create matrix with predicted sigma points as columns
  MatrixXd predicted_sigma_pts_;  
  ///* Weights of sigma points
  VectorXd weights_;
  
  ///* time when the state is true, in us
  long long previous_timestamp_;

  ///* sigma points size of the augmented matrix
  int n_aug_size_;
  ///* Sigma point spreading parameter
  double lambda_;
  
  ///* number of updates performed
  long update_count_;
  ///* number of NIS values above the NIS threshold
  long nis_thresh_count_;
  
  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
  
  /**
   * Predicts the previous sigma points at the current time step.
   * @param {sigma_pts} - matrix containing the previous sigma points
   * @param {pred_sigma_pts} - matrix containing the predicted sigma points 
   * @param {delta_t} - change in time between UKF updates 
   */
  void SigmaPointPrediction(MatrixXd &sigma_pts, MatrixXd &pred_sigma_pts, 
        double delta_t);
  
  /**
   * Predicts the mean and covariance from the provided sigma points matrix.
   * @param {x} - the output vector with the predicted mean values
   * @param {P} - the output Matrix with the predicted covariance values
   * @param {pred_sigma_pts} - the matrix containing the sigma points 
   * @param {yaw_pos} - the vector position of the stored yaw value. -1 if yaw 
   *                    not present in the sigma points. 
   */
  void PredictMeanAndCovariance(VectorXd &x, MatrixXd &P, 
        MatrixXd &pred_sigma_pts, const int yaw_pos);
  
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
  void UpdateState(VectorXd &x, MatrixXd &P, MatrixXd &pred_sigma_pts, 
    MatrixXd &S, MatrixXd &Zsig, VectorXd &z_pred, const int &n_z,
    MeasurementPackage &meas_package);
  
  /**
   * Normalised Innovation Squared for checking noise parameter values,
   * @param {S} - measurement covariance matrix 
   * @param {z_diff} - difference between predicted and measured values
   * @param {n_z} - Number of measured sensor parameters, lidar = 2, Radar = 3.
   */
  void NISState(MatrixXd &S, VectorXd &z_diff, const int &n_z);
};

#endif /* UKF_H */
