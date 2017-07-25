#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  Hj_ = tools.CalculateJacobian(ekf_.x_); //KRO set the actual function here OR maybe set this during runtime
  H_laser_ << 1, 0, 0, 0,   // L10, S5
              0, 1, 0, 0;

  //KRO 4x4 state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
	     0, 0, 0, 1;

  //KRO 4x4 matrix:  set real values
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
	    0, 0, 0, 1000;
  
  // set the acceleration noise components
  noise_ax_ = 9; // L5, sect 13, quiz 9
  noise_ay_ = 9; // 
  //KRO finished TODO
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;  //KRO important for RMSE; first two will be overwritten but should play with the last 2

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //measurement_pack.raw_measurements_(0) is rho
      //measurement_pack.raw_measurements_(1) is theta
      ekf_.x_(0) = measurement_pack.raw_measurements_(0) * cos(measurement_pack.raw_measurements_(1)); //KRO 
      ekf_.x_(1) = measurement_pack.raw_measurements_(0) * sin(measurement_pack.raw_measurements_(1)); //KRO
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //measurement_pack.raw_measurements_(0) is x
      //measurement_pack.raw_measurements_(1) is y
      ekf_.x_(0) = measurement_pack.raw_measurements_(0); //KRO
      ekf_.x_(1) = measurement_pack.raw_measurements_(1); //KRO
    }

    ekf_.F_ << 1, 0, 0, 0,  //KRO 1 diagonal matrix
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;
    
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  //Modify the F matrix so that the time is integrated  L5, S8
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // set the process covariance matrix Q as in L5, s9
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << noise_ax_*(dt_4/4), 0, noise_ax_*(dt_3/2), 0,
	     0, noise_ay_*(dt_4/4), 0, noise_ay_*(dt_3/2),
             noise_ax_*(dt_3/2), 0, noise_ax_*dt_2, 0,
             0, noise_ay_*(dt_3/2), 0, noise_ay_*dt_2;
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_); // KRO set to Hj which should be set using the Jacobian funciton in tools.cpp
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);  //KRO
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;  //KRO (is this good or need to copy the values?)
    ekf_.R_ = R_laser_;  //KRO (is this good or need to copy the values?)
    ekf_.Update(measurement_pack.raw_measurements_);  //KRO
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
