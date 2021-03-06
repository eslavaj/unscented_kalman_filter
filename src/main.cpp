//============================================================================
// Name        :
// Author      : 
// Version     :
// Copyright   :
// Description :
//============================================================================

#include <iostream>
#include "Eigen/Dense"
#include "ukf.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {

	// Create a UKF instance
	UKF ukf;



	/**
	 * Testing Sigma points generation
	 */
	//MatrixXd Xsig = MatrixXd(5, 11);
	//ukf.GenerateSigmaPoints(&Xsig);


	/**
	 * Testing augmented Sigma points generation
	 */
	//MatrixXd Xsig_aug = MatrixXd(7, 15);
	//ukf.AugmentedSigmaPoints(&Xsig_aug);


	/**
	 * Testing augmented sigma points prediction
	 */
	//MatrixXd Xsig_pred = MatrixXd(15, 5);
	//ukf.SigmaPointPrediction(&Xsig_pred);


	/**
	 * Testing the prediction of mean and variance
	 */
	//VectorXd x_pred = VectorXd(5);
	//MatrixXd P_pred = MatrixXd(5, 5);
	//ukf.PredictMeanAndCovariance(&x_pred, &P_pred);

	/**
	 * Testing the prediction of radar measurements
	 */
	//VectorXd z_out = VectorXd(3);
	//MatrixXd S_out = MatrixXd(3, 3);
	//ukf.PredictRadarMeasurement(&z_out, &S_out);


	/**
	 * Testing update state function
	 */
	VectorXd x_out = VectorXd(5);
	MatrixXd P_out = MatrixXd(5, 5);
	ukf.UpdateState(&x_out, &P_out);


	return 0;
}
