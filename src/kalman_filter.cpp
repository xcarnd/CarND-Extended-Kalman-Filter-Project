#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in, Eigen::MatrixXd &H_in,
                        Eigen::MatrixXd &R_laser_in, Eigen::MatrixXd &R_radar_in, Eigen::MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
  I_ = MatrixXd(4, 4);
  I_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_laser_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
  double rho_p = sqrt(px * px + py * py);
  double phi_p = atan2(py, px);
  double rho_dp = (px * vx + py * vy) / rho_p;

  VectorXd y(3);
  y << rho_p, phi_p, rho_dp;
  y = z - y;
  while (y[1] < -M_PI) {
    y[1] += 2 * M_PI;
  }
  while (y[1] > M_PI) {
    y[1] -= 2 * M_PI;
  }
  MatrixXd H_j = tools.CalculateJacobian(x_);
  MatrixXd S = H_j * P_ * H_j.transpose() + R_radar_;
  MatrixXd K = P_ * H_j.transpose() * S.inverse();

  x_ = x_ + K * y;
  P_ = (I_ - K * H_j) * P_;
}
