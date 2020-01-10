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
	MatrixXd Xsig_pred = MatrixXd(15, 5);
	ukf.SigmaPointPrediction(&Xsig_pred);


	return 0;
}
