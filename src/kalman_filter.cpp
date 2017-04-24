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
  MatrixXd Ft = F_.transpose();

  x_ = F_ * x_;
  P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  // z is the measured value of the sensors
  H_ = MatrixXd(2, 4);
  H_ << 1, 0, 0, 0,
          0, 1, 0, 0;

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_ ;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);

  long long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  // z is the measured value of the sensors
  //Convert radar from polar to cartesian coordinates and initialize state.

  double rho = sqrt((x_(0)*x_(0) + x_(1)*x_(1)));
  double phi = atan2(x_(1),x_(0));
  double rhorate= (x_(0)*x_(2) + x_(1)*x_[3])/sqrt((x_(0)*x_(0) + x_(1)*x_(1)));

  // z_pred_radar: using above non-linear equations instead of H
    MatrixXd  z_pred_radar(3, 1);
    z_pred_radar << rho, phi, rhorate;

    VectorXd y = z - z_pred_radar;

    /*
    MatrixXd Ht = H_.transpose();
     MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
     */

    //So instead of calculating the transpose of H_ twice we do it just once,
    // and there are only 3 matrix multiplications instead of 4.

    MatrixXd Ht = H_.transpose();
    MatrixXd PHt = P_ * Ht;
    MatrixXd S = H_ * PHt + R_;
    MatrixXd K = PHt * S.inverse();

    //new estimate
    x_ = x_ + (K * y);
    size_t x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
