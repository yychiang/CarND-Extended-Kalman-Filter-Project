#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
/////////////////////////////
#include <stdlib.h>
#include <time.h>

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
    R_laser_ = MatrixXd(2, 2);  // laser sensor measurement covariance matrix
    R_radar_ = MatrixXd(3, 3);  // radar sensor measurement covariance matrix
    H_laser_ = MatrixXd(2, 4);  // laser sensor measurement matrix
    Hj_ = MatrixXd(3, 4);       // radar sensor Jacobian matrix

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
            0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;

    /*
    int a1,a5,a9;
    float b1,b5,b9;
    srand(time(NULL));
    a1=rand()%1000+1;
    a5=rand()%1000+1;
    a9=rand()%1000+1;
    b1=a1/1000.0;
    b5=a5/1000.0;
    b9=a9/1000.0;

    cout << b1 <<endl;
    cout << b5 <<endl;
    cout << b9 <<endl;

    R_radar_ << b1,0,0,
            0,b5/10,0,
            0,0,b9/10;


    R_radar_ << 0.361,0,0,
            0,0.0019,0,
            0,0,0.0599;
            */



    /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

    H_laser_ << 1,0,0,0,
            0,1,0,0;

    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd(4,4);
    ekf_.F_ = MatrixXd::Identity(4,4);
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.Q_ = MatrixXd(4,4);

    ekf_.H_ << 1, 0, 0, 0,
            0, 1, 0, 0;

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
        ekf_.x_ << 1, 1, 1, 1;

        double x_position=0;
        double y_position=0;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            x_position = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
            y_position = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);


            //Initialize the state ekf_.x_ with the first measurement
            ekf_.x_ << x_position, y_position, 1, 1;


        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x_position = measurement_pack.raw_measurements_[0];
            y_position = measurement_pack.raw_measurements_[1];
            //Initialize the state ekf_.x_ with the first measurement
            ekf_.x_ << x_position, y_position, 1, 1;
        }


        //Create the covariance matrix.
        ekf_.P_ =  MatrixXd(4,4);

        ekf_.P_ << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1000, 0,
                0, 0, 0, 1000;

        is_initialized_ = true;
        previous_timestamp_ = measurement_pack.timestamp_;
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

    // Update the state transition matrix F according to the new elapsed time.
    // Time is measured in seconds.
    double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;

    // Update the process noise covariance matrix, Q.
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt4 = dt3 * dt;
    double noise_ax=9;
    double noise_ay=9;


    ekf_.Q_ <<  dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
            0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
            dt3/2*noise_ax, 0, dt2*noise_ax, 0,
            0, dt3/2*noise_ay, 0, dt2*noise_ay;


    // Call the Kalman Filter predict() function.
    if ( dt > 0.001 )
    {
        ekf_.Predict();
    }
    //ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

    //measurement update

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // radar updates
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
        //cout << "ekf_.x_[0]=" << ekf_.x_[0] <<endl;
        //cout << "ekf_.x_[1]=" << ekf_.x_[1] <<endl;

        if (ekf_.x_[1]>0.05)
        {
            ekf_.UpdateEKF(measurement_pack.raw_measurements_);
        }
    } else {
        // laser updates
        ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, ekf_.H_, R_laser_, ekf_.Q_);
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    //cout << "x_ = " << ekf_.x_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;

}
