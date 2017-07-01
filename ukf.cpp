#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // original value
  // std_a_ = 30;
  // new value
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // original value
  // std_yawdd_ = 30;
  // new value
  std_yawdd_ = 0.5;

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

  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  //set state dimension : Number of rows in our state vector
  n_x_ = 5;

  //set augmented dimension : Number of rows in our state vector + 2 rows for the noise processes
  n_aug_ = n_x_ + 2;;

  //Spreading parameter
  lambda_ = 3 - n_x_;

  //Matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //Vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  //NIS for radar
  NIS_radar_ = 0.0;

  //NIS for laser
  NIS_laser_ = 0.0;
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
      if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ ) {
        /**
        CONVERT RADAR FROM POLAR TO CARTESIAN COORDINATES AND INITIALIZE STATE.
        */

        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        double rho_dot = meas_package.raw_measurements_[2];
        x_ << rho*cos(phi), rho*sin(phi), 0, 0, 0;
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        /**
        INITIALIZE STATE.
        */
        x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      }

      //INIT COVARIANCE MATRIX
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      time_us_ = meas_package.timestamp_;

      //DONE INITIALIZING, NO NEED TO PREDICT OR UPDATE
      is_initialized_ = true;
      return;

  }

  //CALCULATE NEW ELASPED TIME AND TIME IS MEASURED IN SECONDS
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;    //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);
  
  //UPDATE STEP BASED ON LIDAR OR RADAR MEASUREMENT
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } 
  else {
    UpdateLidar(meas_package);
  }
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

  // STEP 1 : GENERATE SIGMA POINTS
      
      //SIGMA POINT MATRIX
      VectorXd x_aug = VectorXd(n_aug_);

      //Augmented state covariance
      MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);

      //SIGMA POINT MATRIX
      MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

      //Augmented mean state
      x_aug.head(5) = x_;
      x_aug(5) = 0;
      x_aug(6) = 0;

      //AUGMENTED COVARIANCE MATRIX
      P_aug.fill(0.0);
      P_aug.topLeftCorner(5,5) = P_;
      P_aug(5,5) = std_a_*std_a_;
      P_aug(6,6) = std_yawdd_*std_yawdd_;

      //SQUARE ROOT MATRIX
      MatrixXd L = P_aug.llt().matrixL();

      //IF DECOMPOSITION FAILS, WE HAVE NUMERICAL ISSUES
      if (P_aug.llt().info() == Eigen::NumericalIssue) {
          std::cout << "LLT failed!" << std::endl;
      throw std::range_error("LLT failed");
      }

      //CREATE AUGMENTED SIGMA POINTS
      Xsig_aug.col(0)  = x_aug;
      for (int i = 0; i< n_aug_; i++)
      {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
      }



  // STEP2 : PREDICT SIGMA POINTS
      for (int i = 0; i< 2*n_aug_+1; i++)
      {
        //EXTRACT VALUES FOR BETTER READABILITY
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        //PREDICTED STATE VALUES
        double px_p, py_p;

        //AVOID DIVISION BY ZERO
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
            }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //ADD PROCESS NOISE
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //WRITE PREDICTED SIGMA POINT INTO RIGHT COLUMN
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
      }


  // STEP 3: COMPUTE MEAN AND COVARIANCE MATRIX PREDICTION
        //SET VECTOR FOR WEIGHTS
        double weight_0 = lambda_/(lambda_+n_aug_);
        weights_(0) = weight_0;
        for (int i=1; i<2*n_aug_+1; i++) {  
          double weight = 0.5/(n_aug_+lambda_);
          weights_(i) = weight;
        }
        
        //PREDICTED STATE MEAN
        x_.fill(0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
          x_ = x_ + Xsig_pred_.col(i) * weights_(i);
        }

        //PREDICTED STATE COVARIANCE MATRIX
        P_.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

          //STATE DIFFERENCE
          VectorXd x_diff = Xsig_pred_.col(i) - x_;
          //ANGLE NORMALIZATION
          while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
          while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

          P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
         }

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

  // STEP: LASER UPDATE 
          //EXTRACT MEASUREMENT AS VECTORXD
          VectorXd z = meas_package.raw_measurements_;

          //SET MEASUREMENT DIMENSION, LIDAR CAN MEASURE P_X AND P_Y
          int n_z = 2;

          //CREATE MATRIX FOR SIGMA POINTS IN MEASUREMENT SPACE
          MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

          //TRANSFORM SIGMA POINTS INTO MEASUREMENT SPACE
          for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

            // extract values for better readibility
            double p_x = Xsig_pred_(0, i);
            double p_y = Xsig_pred_(1, i);

            // measurement model
            Zsig(0, i) = p_x;
            Zsig(1, i) = p_y;
          }

          //MEAN PREDICTED MEASUREMENT
          VectorXd z_pred = VectorXd(n_z);
          z_pred.fill(0.0);
          for (int i = 0; i < 2 * n_aug_ + 1; i++) {
            z_pred = z_pred + weights_(i) * Zsig.col(i);
          }

          //MEASUREMENT COVARIANCE MATRIX S
          MatrixXd S = MatrixXd(n_z, n_z);
          S.fill(0.0);
          for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

            //RESIDUAL
            VectorXd z_diff = Zsig.col(i) - z_pred;

            S = S + weights_(i) * z_diff * z_diff.transpose();
          }

          //ADD MEASUREMENT NOISE COVARIANCE MATRIX
          MatrixXd R = MatrixXd(n_z, n_z);
          R << std_laspx_*std_laspx_, 0,
               0, std_laspy_*std_laspy_;
          S = S + R;

          //CREATE MATRIX FOR CROSS CORRELATION Tc
          MatrixXd Tc = MatrixXd(n_x_, n_z);

          //CALCULATE CROSS CORRELATION MATRIX
          Tc.fill(0.0);
          for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

            //residual
            VectorXd z_diff = Zsig.col(i) - z_pred;

            // state difference
            VectorXd x_diff = Xsig_pred_.col(i) - x_;

            Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
          }

          //KALMAN GAIN K;
          MatrixXd K = Tc * S.inverse();

          //RESIDUAL
          VectorXd z_diff = z - z_pred;

          //CALCULATE NIS
          NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
          cout << "NIS LASER: " << NIS_laser_ << endl;

          //UPDATE STATE MEAN AND COVARIANCE MATRIX
          x_ = x_ + K * z_diff;
          P_ = P_ - K*S*K.transpose();
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
  
  // STEPS : RADAR UPDATE
          //EXTRACT MEASUREMENT AS VECTORXD
          VectorXd z = meas_package.raw_measurements_;

          //SET MEASUREMENT DIMENSION, RADAR CAN MEASURE r, phi, and r_dot
          int n_z = 3;

          //CREATE MATRIX FOR SIGMA POINTS IN MEASUREMENT SPACE
          MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

          //TRANSFORM SIGMA POINTS INTO MEASUREMENT SPACE
          for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

            //EXTRACT VALUES FOR BETTER READIBILITY
            double p_x = Xsig_pred_(0, i);
            double p_y = Xsig_pred_(1, i);
            double v   = Xsig_pred_(2, i);
            double yaw = Xsig_pred_(3, i);

            double v1 = cos(yaw)*v;
            double v2 = sin(yaw)*v;

            //MEASUREMENT MODEL
            Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
            
            //AVOID DIVISION BY ZERO
            if(fabs(p_x) < 0.0001 || fabs(p_y) < 0.0001 ){
              cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
            }
            else {
            Zsig(1, i) = atan2(p_y, p_x);                               //phi
            Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
            }
          }

          //MEAN PREDICTED MEASUREMENT
          VectorXd z_pred = VectorXd(n_z);
          z_pred.fill(0.0);
          for (int i = 0; i < 2 * n_aug_ + 1; i++) {
            z_pred = z_pred + weights_(i) * Zsig.col(i);
          }

          //MEASUREMENT COVARIANCE MATRIX S
          MatrixXd S = MatrixXd(n_z, n_z);
          S.fill(0.0);
          for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
            //RESIDUAL
            VectorXd z_diff = Zsig.col(i) - z_pred;

            //ANGLE NORMALIZATION
            while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
            while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

            S = S + weights_(i) * z_diff * z_diff.transpose();
          }

          //ADD MEASUREMENT NOISE COVARIANCE MATRIX
          MatrixXd R = MatrixXd(n_z, n_z);
          R <<  std_radr_*std_radr_, 0, 0,
                0, std_radphi_*std_radphi_, 0,
                0, 0, std_radrd_*std_radrd_;
          S = S + R;

          //CREATE MATRIX FOR CROSS CORRELATION Tc
          MatrixXd Tc = MatrixXd(n_x_, n_z);
          //CALCULATE CROSS CORRELATION MATRIX
          Tc.fill(0.0);
          for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

            //RESIDUAL
            VectorXd z_diff = Zsig.col(i) - z_pred;
            //angle normalization
            while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
            while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

            //STATE DIFFERENCE
            VectorXd x_diff = Xsig_pred_.col(i) - x_;
            //ANGLE NORMALIZATION
            while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
            while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

            Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
          }

          //KALMAN GAIN K;
          MatrixXd K = Tc * S.inverse();

          //RESIDUAL
          VectorXd z_diff = z - z_pred;

          //ANGLE NORMALIZATION
          while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
          while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

          //CALCULATE NIS
          NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
          cout << "NIS RADAR: " << NIS_radar_ << endl;

          //UPDATE STATE MEAN AND COVARIANCE MATRIX
          x_ = x_ + K * z_diff;
          P_ = P_ - K*S*K.transpose();

}