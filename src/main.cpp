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
	MatrixXd Xsig = MatrixXd(5, 11);
	ukf.GenerateSigmaPoints(&Xsig);




	return 0;
}
