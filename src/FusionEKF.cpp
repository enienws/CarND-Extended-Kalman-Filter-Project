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

	//laser measurement covariance matrix
	R_laser_ << 0.0225, 0,
		0, 0.0225;

	//radar measurement covariance matrix
	R_radar_ << 0.09, 0, 0,
		0, 0.0009, 0,
		0, 0, 0.09;

	/**
	TODO:
	* Finish initializing the FusionEKF.
	* Set the process and measurement noises
	*/

	//Process noise
	Eigen::MatrixXd Q_ = Eigen::MatrixXd(4,4);
	Q_ << 0, 0, 0, 0,
		      0, 0, 0, 0,
			  0, 0, 0, 0,
			  0, 0, 0, 0;

	//State matrix
	Eigen::VectorXd state = Eigen::VectorXd(4);
	state << 0, 0, 0, 0;

	//State covariance matrix
	Eigen::MatrixXd P_ = Eigen::MatrixXd(4, 4);
	P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

	//Laser measurement matrix
	H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;

	//Transition matrix F_
	Eigen::MatrixXd F_ = Eigen::MatrixXd(4, 4);
	F_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

	//Initialize Jacobian matrix
	Hj_ << 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0;

	//Initialize EKF
	ekf_.Init(state, P_, F_, H_laser_, Hj_, R_laser_, R_radar_, Q_);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

Eigen::VectorXd FusionEKF::PolarToCartesian(Eigen::VectorXd polar)
{
	//Range
	float rho = polar(0);

	//Bearing
	float phi = polar(1);

	//Range velocity
	float rho_dot = polar(2);

	VectorXd cartesian(4);
	cartesian(0) = rho * cos(phi);
	cartesian(1) = rho * sin(phi);
	//Set speed to zero since this method is just used for initialization
	cartesian(2) = rho_dot * cos(phi);
	cartesian(3) = rho_dot * sin(phi);

	return cartesian;
}

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

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
		/**
		Convert radar from polar to cartesian coordinates and initialize state.
		*/
    	ekf_.x_ = PolarToCartesian(measurement_pack.raw_measurements_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
		/**
		Initialize state.
		*/
    	//We can directly use the measurement coming from lidar
    	ekf_.x_(0) = measurement_pack.raw_measurements_(0);
    	ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    	ekf_.x_(2) = ekf_.x_(3) = 0;
    }

    // Saving first timestamp in seconds
    previous_timestamp_ = measurement_pack.timestamp_ ;

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

	static const float noise_ax = 9.0f;
	static const float noise_ay = 9.0f;

	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;

	//Update F
	ekf_.F_(0,2) = dt;
	ekf_.F_(1,3) = dt;

	//Update Q
	float delta_t_2 = dt * dt;
	float delta_t_3 = delta_t_2 * dt;
	float delta_t_4 = delta_t_3 * dt;

	ekf_.Q_(0,0) = delta_t_4 / 4 * noise_ax;
	ekf_.Q_(0,2) = delta_t_3 / 2 * noise_ax;
	ekf_.Q_(1,1) = delta_t_4 / 4 * noise_ay;
	ekf_.Q_(1,3) = delta_t_3 / 2 * noise_ay;
	ekf_.Q_(2,0) = delta_t_3 / 2 * noise_ax;
	ekf_.Q_(2,2) = delta_t_2 * noise_ax;
	ekf_.Q_(3,1) = delta_t_3 / 2 * noise_ay;
	ekf_.Q_(3,3) = delta_t_2 * noise_ay;

	ekf_.Predict();

	/*****************************************************************************
	*  Update
	****************************************************************************/

	/**
	TODO:
	 * Use the sensor type to perform the update step.
	 * Update the state and covariance matrices.
	*/

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
	{
		// Radar updates
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	}
	else
	{
		// Laser updates
		ekf_.Update(measurement_pack.raw_measurements_);

	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
