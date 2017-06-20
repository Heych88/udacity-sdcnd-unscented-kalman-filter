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

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;
  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;
  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  //VectorXd x_aug_;
  ///* state covariance matrix
  MatrixXd P_;
  //MatrixXd P_aug_;
  ///* predicted sigma points matrix
  //MatrixXd sigma_pts_;
  //create matrix with predicted sigma points as columns
  MatrixXd predicted_sigma_pts_;
  //MatrixXd Q_;

  ///* time when the state is true, in us
  long long previous_timestamp_;

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

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;
  ///* Augmented state dimension
  int n_aug_;
  ///* sigma points size of the augmented matrix
  int n_aug_size_;
  ///* Sigma point spreading parameter
  double lambda_;
  
  // measurement parameters for the radar sensor
  int n_lidar_;
  ///* radar measurement noise
  //MatrixXd R_lidar_;

  // measurement parameters for the radar sensor
  int n_radar_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

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
  
  void SigmaPointPrediction(MatrixXd &sigma_pts, MatrixXd &pred_sigma_pts, 
        double delta_t);
  
  void PredictMeanAndCovariance(VectorXd &x, MatrixXd &P, 
        MatrixXd &pred_sigma_pts, const int yaw_pos);
  
  void UpdateState(VectorXd &x, MatrixXd &P, MatrixXd &pred_sigma_pts, 
    MatrixXd &S, MatrixXd &Zsig, VectorXd &z_pred, const int &n_z,
    MeasurementPackage &meas_package);
};

#endif /* UKF_H */
