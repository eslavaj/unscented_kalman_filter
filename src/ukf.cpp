/*
 * ukf.cpp
 *
 */


#include "ukf.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

UKF::UKF() {
	Init();
}

UKF::~UKF() {

}

void UKF::Init() {

}


/**
 * Generate sigma points
 */
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

	// set state dimension
	int n_x = 5;

	// define spreading parameter
	double lambda = 3 - n_x;

	// set example state
	VectorXd x = VectorXd(n_x);
	x <<   5.7441,
			1.3800,
			2.2049,
			0.5015,
			0.3528;

	// set example covariance matrix
	MatrixXd P = MatrixXd(n_x, n_x);
	P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
			-0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
			0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
			-0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
			-0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

	// create sigma point matrix
	MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

	// calculate square root of P
	MatrixXd A = P.llt().matrixL();

	MatrixXd x_replic = MatrixXd(n_x, n_x);
	x_replic = x.replicate(1, n_x);
	// calculate sigma points ...
	MatrixXd A_scaled = A * sqrt(lambda + n_x);

	/*Adding first sigma point (mean)*/
	Xsig.col(0) = x;
	/*Adding first shifted points*/
	Xsig.block(0, 1, n_x, n_x) = x_replic + A_scaled;
	/*Adding second shifted points*/
	Xsig.block(0, 1 + n_x, n_x, n_x) = x_replic - A_scaled;

	// print result
	// std::cout << "Xsig = " << std::endl << Xsig << std::endl;

	// write result
	*Xsig_out = Xsig;
}

/**
 * expected result:
 * Xsig =
 *  5.7441  5.85768   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441
 *    1.38  1.34566  1.52806     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38
 *  2.2049  2.28414  2.24557  2.29582   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049
 *  0.5015  0.44339 0.631886 0.516923 0.595227   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015
 *  0.3528 0.299973 0.462123 0.376339  0.48417 0.418721 0.405627 0.243477 0.329261  0.22143 0.286879
 */



/**
 * Generate sigma points for the augmented state
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {




	// set state dimension
	//The augmented state includes the state and linear and phase acceleration
	//In this case we add 2 process noise dimensions linear and phase acceleration
	int n_aug = 7;
	int n_x = 5;


	// define spreading parameter
	double lambda = 3 - n_aug;

	// set example state
	VectorXd x = VectorXd(n_aug);
	x <<   5.7441,
			1.3800,
			2.2049,
			0.5015,
			0.3528,
			0,
			0;

	// set example covariance matrix
	MatrixXd P = MatrixXd(n_x, n_x);
	P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
			-0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
			0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
			-0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
			-0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

	//For the augmented state we will use P_aug instead of P
	//P_aug is a matrix that contains P and Q (Q is the process noise covariance matrix)

	// Process noise standard deviation longitudinal acceleration in m/s^2
	double std_a = 0.2;
	// Process noise standard deviation yaw acceleration in rad/s^2
	double std_yawdd = 0.2;
	//Q process noise covariance matrix
	MatrixXd Q = MatrixXd(n_aug - n_x, n_aug - n_x);
	Q << std_a*std_a, 0,
			0, std_yawdd*std_yawdd;

	//P_aug
	MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug);
	P_aug.topLeftCorner(P.rows(), P.cols()) = P;
	P_aug.bottomRightCorner(Q.rows(), Q.cols()) = Q;

	// create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

	// calculate square root of P_aug
	MatrixXd A = P_aug.llt().matrixL();

	MatrixXd x_replic = MatrixXd(n_aug, n_aug);
	x_replic = x.replicate(1, n_aug);
	MatrixXd A_scaled = A * sqrt(lambda + n_aug);

	/*Adding first sigma point (mean)*/
	Xsig_aug.col(0) = x;
	/*Adding first shifted points*/
	Xsig_aug.block(0, 1, n_aug, n_aug) = x_replic + A_scaled;
	/*Adding second shifted points*/
	Xsig_aug.block(0, 1 + n_aug, n_aug, n_aug) = x_replic - A_scaled;

	// print result
	std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

	*Xsig_out = Xsig_aug;
}

/**
 * expected result:
 *  Xsig_aug =
 * 5.7441  5.85768   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441
 *   1.38  1.34566  1.52806     1.38     1.38     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38     1.38     1.38
 * 2.2049  2.28414  2.24557  2.29582   2.2049   2.2049   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049   2.2049   2.2049
 * 0.5015  0.44339 0.631886 0.516923 0.595227   0.5015   0.5015   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015   0.5015   0.5015
 * 0.3528 0.299973 0.462123 0.376339  0.48417 0.418721   0.3528   0.3528 0.405627 0.243477 0.329261  0.22143 0.286879   0.3528   0.3528
 *      0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641        0
 *      0        0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641
 */

