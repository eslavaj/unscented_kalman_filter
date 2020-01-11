/*
 * ukf.cpp
 *
 */


#include "ukf.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

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


void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {

	// set state dimension
	int n_x = 5;

	// set augmented dimension
	int n_aug = 7;

	// create example sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
	Xsig_aug <<
			5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
			1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
			2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
			0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
			0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
			0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
			0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

	// create matrix with predicted sigma points as columns
	MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

	double delta_t = 0.1; // time diff in sec


	// predict sigma points
	for(int j=0; j< Xsig_aug.cols(); j++)
	{
		/*
		 * Reminder: augmented state components:
		 *
		 *   [px py v yaw yaw_vel noise_acc noise_yaw_acc].transpose
		 *
		 *
		 * */

		//double px = Xsig_aug(0,j);
		//double py = Xsig_aug(1,j);
		double v = Xsig_aug(2,j);
		double yaw = Xsig_aug(3,j);
		double yaw_vel = Xsig_aug(4,j);
		double noise_acc = Xsig_aug(5,j);
		double noise_yaw_acc = Xsig_aug(6,j);

		MatrixXd firstTerm = MatrixXd(n_x,1);
		MatrixXd secondTerm = MatrixXd(n_x,1);

		//Check yaw_vel to avoid division by zero
		if(fabs(yaw_vel) <= 0.001)
		{
			firstTerm << v*cos(yaw)*delta_t,
						 v*sin(yaw)*delta_t,
						 0,
						 yaw_vel*delta_t,
						 0;
		}
		else
		{
			firstTerm << (v/yaw_vel)*( sin(yaw + yaw_vel*delta_t) - sin(yaw) ),
						 (v/yaw_vel)*( -cos(yaw + yaw_vel*delta_t) + cos(yaw) ),
						 0,
						 yaw_vel*delta_t,
						 0;
		}


		secondTerm << 0.5*pow(delta_t,2)*cos(yaw)*noise_acc,
					  0.5*pow(delta_t,2)*sin(yaw)*noise_acc,
					  delta_t*noise_acc,
					  0.5*pow(delta_t,2)*noise_yaw_acc,
					  delta_t*noise_yaw_acc;


		Xsig_pred.block(0, j, Xsig_pred.rows(), 1) = Xsig_aug.block(0, j, Xsig_pred.rows(), 1) + firstTerm + secondTerm;


	}


	// print result
	std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

	// write result
	*Xsig_out = Xsig_pred;
}
/*
Xsig_pred =

5.93553 6.06251 5.92217 5.9415 5.92361 5.93516 5.93705 5.93553 5.80832 5.94481 5.92935 5.94553 5.93589 5.93401 5.93553

1.48939 1.44673 1.66484 1.49719 1.508 1.49001 1.49022 1.48939 1.5308 1.31287 1.48182 1.46967 1.48876 1.48855 1.48939

2.2049 2.28414 2.24557 2.29582 2.2049 2.2049 2.23954 2.2049 2.12566 2.16423 2.11398 2.2049 2.2049 2.17026 2.2049

0.53678 0.473387 0.678098 0.554557 0.643644 0.543372 0.53678 0.538512 0.600173 0.395462 0.519003 0.429916 0.530188 0.53678 0.535048

0.3528 0.299973 0.462123 0.376339 0.48417 0.418721 0.3528 0.387441 0.405627 0.243477 0.329261 0.22143 0.286879 0.3528 0.318159
 */



void UKF::calcMeanVariance(MatrixXd & Xsig_in, VectorXd & weights_in, VectorXd & xmean_out, MatrixXd &Mvariance_out, int rowToNormalize2pi)
{

	int nbr_of_sigmapoints = Xsig_in.cols();
	int state_dimension = Xsig_in.rows();

	MatrixXd weigths_mat = MatrixXd(nbr_of_sigmapoints, nbr_of_sigmapoints);
	weigths_mat = weights_in.replicate(1, nbr_of_sigmapoints);

	MatrixXd Xsig_pred_weigthed = Xsig_in*weigths_mat;

	/*Calculating predicted mean*/
	xmean_out = Xsig_pred_weigthed.rowwise().sum();
	xmean_out = xmean_out/nbr_of_sigmapoints;

	MatrixXd x_mat = xmean_out.replicate(1, nbr_of_sigmapoints);
	MatrixXd Xsig_in_centered = Xsig_in - x_mat;

	/*
	 * Normalizing from 0 to 2*pi
	 * Formula to normalize from 0 to 2*pi:  ( offsetValue - ( round( offsetValue / width ) * width ) ) + start
	 * width = 2*pi
	 * start = 0
	 * */
	Xsig_in_centered.block(rowToNormalize2pi, 0, 1, nbr_of_sigmapoints) = Xsig_in_centered.block(rowToNormalize2pi, 0, 1, nbr_of_sigmapoints).array() + \
				( -1*Xsig_in_centered.block(rowToNormalize2pi, 0, 1, nbr_of_sigmapoints)/(2*M_PI) ).array().round()*(2*M_PI);

	/*Calculating predicted covariance*/
	MatrixXd Xsig_in_centered_weighted = (weigths_mat.leftCols(state_dimension).transpose().array() )*Xsig_in_centered.array();
	Mvariance_out = Xsig_in_centered_weighted*( Xsig_in_centered.transpose());

}




/**
 * Predict the mean and covariance of the sigma points:
 */

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

	// set state dimension
	int n_x = 5;
	// set augmented dimension
	int n_aug = 7;
	// define spreading parameter
	double lambda = 3 - n_aug;

	// create example matrix with predicted sigma points
	MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
	Xsig_pred <<
			5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
			1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
			2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
			0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
			0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

	// create vector for weights
	VectorXd weights = VectorXd(2*n_aug+1);

	/*Filling weights vector and convert it weights matrix*/
	for(int i=0; i<weights.rows(); i++)
	{
		if(i==0)
		{
			weights(i) = lambda/(lambda + n_aug);
		}
		else
		{
			weights(i) = 0.5/(lambda + n_aug);
		}
	}

	// create vector for predicted state
	VectorXd x = VectorXd(n_x);

	// create covariance matrix for prediction
	MatrixXd P = MatrixXd(n_x, n_x);

	calcMeanVariance(Xsig_pred, weights, x, P, 3);

	// print result
	std::cout << "Predicted state" << std::endl;
	std::cout << x << std::endl;
	std::cout << "Predicted covariance matrix" << std::endl;
	std::cout << P << std::endl;

	// write result
	*x_out = x;
	*P_out = P;
}


/*

expected result x:
x =

5.93637

1.49035

2.20528

0.536853

0.353577

expected result p:
P =

0.00543425 -0.0024053 0.00341576 -0.00348196 -0.00299378

-0.0024053 0.010845 0.0014923 0.00980182 0.00791091

0.00341576 0.0014923 0.00580129 0.000778632 0.000792973

-0.00348196 0.00980182 0.000778632 0.0119238 0.0112491

-0.00299378 0.00791091 0.000792973 0.0112491 0.0126972


 */


void UKF::TransformToMeasureSpace(MatrixXd & Xsig_in, MatrixXd & Xsig_meaSpace_out)
{
	int nbr_of_sigmapoints = Xsig_in.cols();

	for(int j=0; j< nbr_of_sigmapoints; j++)
	{
		double px = Xsig_in(0, j);
		double py = Xsig_in(1, j);
		double v = Xsig_in(2, j);
		double yaw = Xsig_in(3, j);
		//double yaw_vel = Xsig_in(4, j);

		Xsig_meaSpace_out(0, j) = sqrt(px*px + py*py);
		Xsig_meaSpace_out(1, j) = atan(py/px);
		Xsig_meaSpace_out(2, j) = ( px*cos(yaw)*v + py*sin(yaw)*v )/( Xsig_meaSpace_out(0, j) );

	}
}



/**
 * Predict tadar measurements
 */

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {

  // set state dimension
  int n_x = 5;

  // set augmented dimension
  int n_aug = 7;

  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // define spreading parameter
  double lambda = 3 - n_aug;

  // set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  double weight = 0.5/(lambda+n_aug);
  weights(0) = weight_0;

  for (int i=1; i<2*n_aug+1; ++i) {
    weights(i) = weight;
  }

  // radar measurement noise standard deviation radius in m
  double std_radr = 0.3;

  // radar measurement noise standard deviation angle in rad
  double std_radphi = 0.0175;

  // radar measurement noise standard deviation radius change in m/s
  double std_radrd = 0.1;

  /*Creating measure noise matrix*/
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr*std_radr, 0, 0,
	   0, std_radphi*std_radphi, 0,
	   0, 0, std_radrd*std_radrd;


  // create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  // transform sigma points into measurement space
  TransformToMeasureSpace(Xsig_pred, Zsig);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  // calculate mean predicted and covariance measurement
  calcMeanVariance(Zsig, weights, z_pred, S, 1);

  /*Adding measure noise*/
  S = S + R;

  // print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  // write result
  *z_out = z_pred;
  *S_out = S;
}


/*

expected result z_out:
z_pred =

6.12155

0.245993

2.10313

expected result s_out:
S =

0.0946171 -0.000139448 0.00407016

-0.000139448 0.000617548 -0.000770652

0.00407016 -0.000770652 0.0180917


 * */



