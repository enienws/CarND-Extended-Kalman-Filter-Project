#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
	/**
	TODO:
	* Calculate the RMSE here.
	*/

	VectorXd rmse(4);
	rmse << 0,0,0,0;

    if(estimations.size() == 0 || ground_truth.size() == 0)
    {
        //error
    }

    if (estimations.size() != ground_truth.size())
    {
        //error
    }

	//accumulate squared residuals
	for(size_t i=0; i < estimations.size(); i++)
	{
        VectorXd difference = estimations[i] - ground_truth[i];
        difference = difference.array() * difference.array();
        rmse += difference;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);


	//pre-compute a set of terms to avoid repeated calculation
	double a1 = px*px+py*py;
	double a2 = sqrt(a1);
	double a3 = (a1*a2);

	//check division by zero
	if(fabs(a1) < 0.0001)
	{
		//Error
		cout << "ERROR - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/a2), (py/a2), 0, 0,
		  -(py/a1), (px/a1), 0, 0,
		  py*(vx*py - vy*px)/a3, px*(px*vy - py*vx)/a3, px/a2, py/a2;

	return Hj;
}
