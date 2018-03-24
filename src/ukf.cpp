#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static VectorXd processModel(const VectorXd& Xaug, double dt);

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false; // initalize to false

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  Xsig_pred_ = MatrixXd(5, 2*7+1);


  ///* Weights of sigma points
  weights_ = VectorXd(2*7+1);

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {

  }
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {

  }

}

static VectorXd processModel(const VectorXd& Xaug, double dt)
{
    VectorXd xk(5);
    double px = Xaug(0);
    double py = Xaug(1);
    double v = Xaug(2);
    double psi = Xaug(3);
    double psidot = Xaug(4);
    double nu_accel = Xaug(5);
    double nu_psidd = Xaug(6);
    
    xk(0) = Xaug(0) + v/max(psidot,0.0001)*(sin(psi+psidot*dt)-sin(psi)) + 0.5*dt*dt*cos(psi)*nu_accel;
    xk(1) = Xaug(1) + v/max(psidot,0.0001)*(-cos(psi+psidot*dt)+cos(psi)) + 0.5*dt*dt*sin(psi)*nu_accel;
    xk(2) = Xaug(2) + dt*nu_accel;
    xk(3) = Xaug(3) + psidot*dt + 0.5*dt*dt*nu_psidd;
    xk(4) = Xaug(4) + dt*nu_psidd;
    
    return xk;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    //predict sigma points


  // MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // for (int i = 0; i < 2*n_aug+1; i++)
  // {
  //     Xsig_pred.col(i) = processModel(Xsig_aug.col(i), delta_t);
  // }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
