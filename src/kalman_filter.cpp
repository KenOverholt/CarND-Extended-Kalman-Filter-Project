#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	x_ = F_ * x_; //KRO L5s8
	MatrixXd Ft = F_.transpose(); //KRO L5s9
	P_ = F_ * P_ * Ft + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  //L5s7
  VectorXd y = z - H * x;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P * Ht * Si;

  //new state
  x = x + (K * y);
  P = (I - K * H) * P;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //L5s14
  float x = ekf_.x_(0);
  float y = ekf_.x_(1);
  float vx = ekf_.x_(2);
  float vy = ekf_.x_(3);
  
  float rho = sqrt(x*x+y*y);
  float theta = atan2(y, x);
  float rho_dot = (x*vx+y*vy)/rho;
  float z_pred << rho, theta, rho_dot;
  
  VectorXd y = z - z_pred;
  
  //KRO L5s7
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P * Ht * Si;

  //new state
  x = x + (K * y);
  P = (I - K * H) * P;
}
