#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "stdlib.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static VectorXd processModel(const VectorXd& Xaug, double dt);
static VectorXd measLidarModel(const VectorXd& x);

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  previous_timestamp_ = 0.f;

  is_initialized_ = false; // initalize to false

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
  
  Xsig_pred_ = MatrixXd(5, 2*7+1);

  ///* Weights of sigma points
  weights_ = VectorXd(2*7+1);

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  //* Radar measurement dimension
  n_z_radar_ = 3;
 
  //* Laser measurement dimension
  n_z_laser_ = 2;
 
  ///* Sigma point spreading parameter
  lambda_ = 3.0;
  
  ///* Measurement Noise Matrices
  R_Radar_ = MatrixXd(3,3);
  R_Radar_ << std_radr_*std_radr_, 0 ,0, 
              0, std_radphi_*std_radphi_, 0, 
              0, 0, std_radrd_*std_radrd_;

  R_Laser_ = MatrixXd(2,2);
  R_Laser_ << std_laspx_*std_laspx_, 0 , 
              0, std_laspy_*std_laspy_;             
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
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state matrix.
      */
      /* Initialize state vector */

      double rho     = meas_package.raw_measurements_[0]; // radial distance to object
      double phi     = meas_package.raw_measurements_[1]; // bearing angle between object and vehicle current heading
      double rho_dot = meas_package.raw_measurements_[2]; // radial velocity to object
      double vx = rho_dot*cos(phi);
      double vy = rho_dot*sin(phi);
      x_ << rho*cos(phi), rho*sin(phi), sqrt(vx*vx+vy*vy), 0, 0; // convert polar to cartesian x = rho*cos(phi), y = rho*sin(phi), assume vx, vy without phidot measurement
      use_radar_ = false;  
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state vector
      */
      double px, py;
      if (meas_package.raw_measurements_[0] < 0.001 and meas_package.raw_measurements_[1] < 0.001)
      {
        px = 0;
        py = 0;
      }
      else
      {
        px = meas_package.raw_measurements_[0];
        py = meas_package.raw_measurements_[1];
      }

      x_ <<  px, py, 0, 0, 0; // x,y distances to object (directly in cartesian coord.)  
      use_radar_ = false;      
    }
  /* Initialize state process covariance matrix */
  for (uint8_t i = 0; i < n_x_; i++)
  {
    P_(i,i) = 1;
  }

  P_(2,2) = 100; // velocity state initially more uncertain

  /* Initialize weights */
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i = 1; i < 2*n_aug_+1; i++)
  {
      weights_(i) = 1/(2*(lambda_+n_aug_));
  }  

  /* Initialize Sigma Points */
  is_initialized_ = true;
  previous_timestamp_ = meas_package.timestamp_;

  }
  else // Run if already initialized
  {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      use_radar_ = true;
      use_laser_ = false;      
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) 
    {
      use_radar_ = false;      
      use_laser_ = true;
    }

    double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
    previous_timestamp_ = meas_package.timestamp_;

    Prediction(dt); // Prediction

    UpdateMeasurement(meas_package); // Update
  }
}

/* Helper function which defines the process model using the augmented state vector */
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
    xk(0) = px; // px
    xk(1) = py; // py  

    if (fabs(psidot) > 0.001) // check divide by zero
    {
      xk(0) += v/psidot*(sin(psi+psidot*dt)-sin(psi));
      xk(1) += v/psidot*(-cos(psi+psidot*dt)+cos(psi));
    }
    else
    {
      xk(0) += v*dt*cos(psi); // simplified model no psi dot term
      xk(1) += v*dt*sin(psi);     
    }

    xk(0) += 0.5*dt*dt*cos(psi)*nu_accel; // add accel noise
    xk(1) += 0.5*dt*dt*sin(psi)*nu_accel; // add accel noise
    xk(2) = v + dt*nu_accel; // add accel noise
    xk(3) = psi + psidot*dt + 0.5*dt*dt*nu_psidd;  // psi    
    xk(4) = psidot + dt*nu_psidd; //psidot
    
    return xk;
}

/* Helper function which creates the augmented state vector, sigma points and covariance matrix  */
void UKF::AugmentedSigmaPoints(VectorXd* pX_aug, MatrixXd* pXsig_aug, MatrixXd* pP_aug)
{
  VectorXd x_aug = *pX_aug; // copy deference for ease of use;
  MatrixXd Xsig_aug = *pXsig_aug; // copy deference for ease of use;
  MatrixXd P_aug = *pP_aug; // copy deference for ease of use;

  //create augmented mean state
  for (int i = 0; i < n_x_; i++)
  { 
     x_aug(i) = x_(i); 
  }
  x_aug(5) = 0;
  x_aug(6) = 0;
 
   //create augmented covariance matrix
  P_aug.block(0,0,n_x_,n_x_) = P_;

  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //calculate square root of P_aug
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug;

  //set remaining sigma points
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)     = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
  }  
  *pX_aug = x_aug;
  *pXsig_aug = Xsig_aug;
  *pP_aug = P_aug;
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
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); 

  Xsig_aug.fill(0.0);

  /* Generate Augmented sigma points from state vector, sigma points, and covariance matrix */
  AugmentedSigmaPoints(&x_aug, &Xsig_aug, &P_aug);

  /* Pass augmented sigma points through process model */
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
      Xsig_pred_.col(i) = processModel(Xsig_aug.col(i), delta_t);
  }

  /* Predicted Mean */
  x_ = Xsig_pred_*weights_;

  P_.fill(0.0);

  /* Predicted Covariance Matrix */
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
     VectorXd state_error = Xsig_pred_.col(i)-x_;
     state_error(3) = tools.unWrapAngle(state_error(3));
     P_ += weights_(i)*state_error*state_error.transpose();
  }

}

VectorXd UKF::measRadarModel(const VectorXd& x)
{
    VectorXd z(3);
    double px = x(0);
    double py = x(1);
    double v = x(2);
    double psi = x(3);
    double psidot = x(4);
    
    z(0) = sqrt(px*px+py*py); // rho
    z(1) = tools.unWrapAngle(atan2(py,px)); // phi
    z(2) = v*(px*cos(psi)+py*sin(psi))/z(0); // rhodot
    
    return z;
}

static VectorXd measLidarModel(const VectorXd& x)
{
    VectorXd z(2);

    z(0) = x(0); // px
    z(1) = x(1); // py

    return z;
}


void UKF::UpdateState(MeasurementPackage meas_package, VectorXd z_pred, MatrixXd Zsig, MatrixXd S)
{
  int n_z_ = 1;
  if (use_radar_)
  {
    n_z_ = n_z_radar_;
  }
  else if(use_laser_)
  {
    n_z_ = n_z_laser_;
  }
 
  VectorXd z_meas(n_z_);  
  if (use_radar_)
  {
    z_meas(0) = meas_package.raw_measurements_[0]; // radial distance to object
    z_meas(1) = meas_package.raw_measurements_[1]; // bearing angle between object and vehicle current heading
    z_meas(2) = meas_package.raw_measurements_[2]; // radial velocity to object
  }
  else if(use_laser_)
  {  
    z_meas(0) = meas_package.raw_measurements_[0]; // position x to object
    z_meas(1) = meas_package.raw_measurements_[1]; // position y to object
  }


  MatrixXd Tc(n_x_, n_z_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
      VectorXd x_diff = Xsig_pred_.col(i)-x_;
      x_diff(3) = tools.unWrapAngle(x_diff(3));
      VectorXd z_diff = Zsig.col(i)-z_meas;
      if (use_radar_)
      {       
        z_diff(1) = tools.unWrapAngle(z_diff(1));
      }
      Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = z_meas-z_pred;
  if (use_radar_)
  {   
    z_diff(1) = tools.unWrapAngle(z_diff(1));
  }

  // state vector and covariance matrix update //
  x_ = x_ + K*z_diff;
  P_ = P_ - K*S*K.transpose();

}

void UKF::UpdateMeasurement(MeasurementPackage meas_package)
{
  int n_z_ = 1;

  if (use_radar_)
  {
    n_z_ = n_z_radar_;
  }
  else if(use_laser_)
  {
    n_z_ = n_z_laser_;
  }

  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S = MatrixXd(n_z_, n_z_);
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);

  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    if (use_radar_)
    {
      Zsig.col(i) = measRadarModel(Xsig_pred_.col(i));
    }
    else if(use_laser_)
    {
      Zsig.col(i) = measLidarModel(Xsig_pred_.col(i));
    } 
  }

  //calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
      z_pred += weights_(i)*Zsig.col(i);
  } 

  //calculate innovation covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
      VectorXd error = Zsig.col(i)-z_pred;
      if (use_radar_)
      {
        error(1) = tools.unWrapAngle(error(1));
      }
      S += weights_(i)*error*error.transpose();        
  }

  if (use_radar_)
  {
    S += R_Radar_;
  }
  else if(use_laser_)
  {
    S += R_Laser_;
  }

  UpdateState(meas_package, z_pred, Zsig, S);
}

