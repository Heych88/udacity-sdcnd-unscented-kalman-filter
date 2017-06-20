#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

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
  std::cout << "done" << std::endl;
  rmse = rmse / size;
  rmse = rmse.array().sqrt();
  std::cout << rmse << std::endl;
  return rmse;
}