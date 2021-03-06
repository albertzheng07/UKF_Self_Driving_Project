#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

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

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

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

  int n_z_laser_;

  int n_z_radar_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* R Radar Noise Matrix
  MatrixXd R_Radar_;

  ///* R Laser Noise Matrix
  MatrixXd R_Laser_;

  double previous_timestamp_;

  int counter;

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

  // *
  //  * Updates the predicted measurement, measurement sigma points and measurement noise matrix using either the laser or radar measurement
  //  * @param meas_package The measurement at k+1
  void UpdateMeasurement(MeasurementPackage meas_package);

  // *
  //  * Updates the state and the state covariance matrix using either the laser or radar measurement
  //  * @param meas_package The measurement at k+1
  //  * @param zpred The predicted measurement at k+1
  //  * @param Zsig The predicted measurement sigma points at k+1
  //  * @param S The predicted measurement noise matrix at k+1
  void UpdateState(MeasurementPackage meas_package, VectorXd zpred, MatrixXd Zsig, MatrixXd S);

  /* 
  * Generates Augmented Sigma points from state vector, state sigma points and state covariance matrix 
  
  //  * @param x_aug predicted augmented state vector
  //  * @param Xsig_aug predicted augmented state matrix of sigma points
  //  * @param P_aug predicted augmented state covariance matrix
  */
  void AugmentedSigmaPoints(VectorXd* x_aug, MatrixXd* Xsig_aug, MatrixXd* P_aug);

  /* 
  * Returns measurement Radar Model from current state vector
  //  * @param x predicted state vector
  */
  VectorXd measRadarModel(const VectorXd& x);

  /* Tools instance lives here */
  Tools tools;

};

#endif /* UKF_H */
