#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in,
	      Eigen::MatrixXd &H_laser_in, Eigen::MatrixXd &H_radar_in,
		  Eigen::MatrixXd &R_laser_in, Eigen::MatrixXd &R_radar_in,
		  Eigen::MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;

  H_laser_ = H_laser_in;

  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

	//Predict the state vector
	x_ = F_ * x_;

	//Predict the state covariance matrix
	MatrixXd Ft_ = F_.transpose();
	P_ = F_ * P_ * Ft_ + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

	//Update
	//Update for linear measurements, ie. lidar

	//Predicted state to measurement space
	VectorXd z_pred = H_laser_ * x_;

	//How differs measured state and predicted state
	VectorXd y = z - z_pred;

	//Calculate the Kalman Gain
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//Estimate based on Kalman Gain, if measurement
	//uncertainty is high, get more weight to prediction
	//if measurement uncertainty is low, get more weight to
	//measurement.
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

	//Update
	//Update for measurements that need non-linear transformation to object's state space

	//First of all compute the Jacobian matrix
	MatrixXd Hj_ = tools.CalculateJacobian(x_);

	//Predicted state to measurement space
	VectorXd z_pred = CartesianToPolar(x_);

	//How differs measured state and predicted state
	VectorXd y = z - z_pred;
	if( y[1] > M_PI )
		y[1] -= 2.f * M_PI;
	if( y[1] < -M_PI )
		y[1] += 2.f * M_PI;

	//Calculate the Kalman Gain
	MatrixXd Hjt = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Hjt + R_radar_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Hjt;
	MatrixXd K = PHt * Si;

	//Estimate based on Kalman Gain, if measurement
	//uncertainty is high, get more weight to prediction
	//if measurement uncertainty is low, get more weight to
	//measurement.
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj_) * P_;
}

Eigen::VectorXd KalmanFilter::CartesianToPolar(VectorXd & cartesian)
{
	//Initialize the matrix for polar coordinates
	Eigen::VectorXd polar = Eigen::VectorXd(3);

	//Get the cartesian position
	float px = cartesian(0);
	float py = cartesian(1);
	float vx = cartesian(2);
	float vy = cartesian(3);

	//Convert from cartesian to polar coordinates
	polar(0) = sqrt(px * px + py * py);
	polar(1) = atan2(py, px);
	polar(2) = (px * vx + py * vy) / sqrt(px * px + py * py);

	return polar;
}

