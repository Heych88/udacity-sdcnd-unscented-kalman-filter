#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}
/*
 * Calculates the Root Mean Square Error of between the the estimated position
 * and the true position.
 * RMSE = sqrt(1/(vector size) * sum(square(estimate-truth)))
 * @param estimations : predicted position
 * @param ground_truth : actual measured position
 * @return : the error between the estimate and ground_truth 
 * below is copied from the Extended Kalman filter project found
 * https://github.com/Heych88/udacity-sdcnd-extended-kalman-filter/blob/master/src/tools.cpp
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse.fill(0.0);
  
  // verify the input data is as expected
  int size = estimations.size();
  if(size != ground_truth.size() || size == 0){
    std::cout << "ERROR: Miss-match in RMSE dimensions. estimations.size() = "
        << size << ",  ground_truth.size() = " << ground_truth.size() 
        << std::endl;
    return rmse;
  }
  
  // RMSE = sqrt(1/(vector size) * sum(square(estimate-truth)))
  for(int i=0; i < size; i++){
    VectorXd difference = estimations[i] - ground_truth[i];
    difference = difference.array() * difference.array();
    rmse = rmse + difference;
  }
  
  rmse = rmse / size;
  rmse = rmse.array().sqrt();
  std::cout << rmse << std::endl;
  return rmse;
}