/*
 * ukf.hpp
 *
 *  Created on: Jan 1, 2020
 *      Author: jeslava
 */

#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"

class UKF {
public:
	/**
	 * Constructor
	 */
	UKF();

	/**
	 * Destructor
	 */
	virtual ~UKF();

	/**
	 * Init Initializes Unscented Kalman filter
	 */
	void Init();

	/**
	 * Student assignment functions
	 */
	void GenerateSigmaPoints(Eigen::MatrixXd* Xsig_out);
	void AugmentedSigmaPoints(Eigen::MatrixXd* Xsig_out);
	void SigmaPointPrediction(Eigen::MatrixXd* Xsig_out);
	void PredictMeanAndCovariance(Eigen::VectorXd* x_pred,
			Eigen::MatrixXd* P_pred);
	void PredictRadarMeasurement(Eigen::VectorXd* z_out,
			Eigen::MatrixXd* S_out);
	void UpdateState(Eigen::VectorXd* x_out,
			Eigen::MatrixXd* P_out);


private:
	void calcMeanVariance(Eigen::MatrixXd & Xsig_in, Eigen::VectorXd & weights_in, Eigen::VectorXd & xmean_out, Eigen::MatrixXd &Mvariance_out, int rowToNormalize2pi);
	void TransformToMeasureSpace(Eigen::MatrixXd & Xsig_in, Eigen::MatrixXd & Xsig_meaSpace_out);
	void calcWeightedMeanAndCenter(Eigen::MatrixXd & Xsig_in, Eigen::VectorXd & weights_in, Eigen::VectorXd & xmean_out, Eigen::MatrixXd & Xsig_centered_out);
	void normalizeRow(Eigen::MatrixXd & blockToNormalize, int rowToNormalize, double start, double width);
	void calcWeightedCorrelation(Eigen::MatrixXd & M1_in, Eigen::MatrixXd & M2_in, Eigen::VectorXd & weights_in, Eigen::MatrixXd & Mvariance_out);
};


#endif  // UKF_H
