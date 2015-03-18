#include <stdlib.h>
#include <stdio.h>

#include "ez_ransac.h"

namespace mylib{

	vector<double> rigidRansac(const std::vector<KeyPoint>& skeypoints, const std::vector<KeyPoint>& mkeypoints,
		const vector< DMatch >& allMatches, vector<int>& inliers, double confidence/* = 0.95*/)
	{
		inliers.clear();
		//-- Create input data
		Eigen::MatrixXd dataPoints((int)allMatches.size(), 4);
		for (unsigned int i = 0; i < allMatches.size(); ++i)
		{
			dataPoints(i, 0) = skeypoints[allMatches[i].queryIdx].pt.x;
			dataPoints(i, 1) = skeypoints[allMatches[i].queryIdx].pt.y;
			dataPoints(i, 2) = mkeypoints[allMatches[i].trainIdx].pt.x;
			dataPoints(i, 3) = mkeypoints[allMatches[i].trainIdx].pt.y;
		}
		// RANSAC detect outliers
		auto_ptr< estimators::Solver<Eigen::MatrixXd, Eigen::VectorXd> > ptrSolver(
			new estimators::rigidSolver<Eigen::MatrixXd, Eigen::VectorXd>);
		//for (int i = 0; i < (int)good_matches.size(); i++) inliers.push_back(i);
		vector<Eigen::VectorXd> models;

		ransac::Ransac_Handler ransac_fun_Handler;
		bool result = ransac::Ransac_RobustEstimator
			(
			dataPoints, // the input data
			estimators::rigidSolver<Eigen::MatrixXd, Eigen::VectorXd>::extractor, // How select sampled point from indices
			dataPoints.rows(),  // the number of putatives data
			*(ptrSolver.get()),  // compute the underlying model given a sample set
			estimators::rigidSolver<Eigen::MatrixXd, Eigen::VectorXd>::defaultEvaluator,  // the function to evaluate a given model
			//Ransac Object that contain function:
			// CandidatesSelector, Sampler and TerminationFunction
			ransac_fun_Handler, // the basic ransac object
			1000,  // the maximum rounds for RANSAC routine
			inliers, // inliers to the final solution
			models, // models array that fit input data
			0.9//0.95 // the confidence want to achieve at the end
			);

		if (models.size() < 1)
		{
			return vector<double>();
		}
		vector<double> model(models[0].size());
		for (size_t i = 0; i < models[0].size(); i++)
		{
			model[i] = models[0][i];
		}
		return model;
	}


	vector<double> affineRansac(const std::vector<KeyPoint>& skeypoints, const std::vector<KeyPoint>& mkeypoints,
		const vector< DMatch >& allMatches, vector<int>& inliers, double confidence/* = 0.95*/)
	{
		inliers.clear();
		//-- Create input data
		Eigen::MatrixXd dataPoints((int)allMatches.size(), 4);
		for (unsigned int i = 0; i < allMatches.size(); ++i)
		{
			dataPoints(i, 0) = skeypoints[allMatches[i].queryIdx].pt.x;
			dataPoints(i, 1) = skeypoints[allMatches[i].queryIdx].pt.y;
			dataPoints(i, 2) = mkeypoints[allMatches[i].trainIdx].pt.x;
			dataPoints(i, 3) = mkeypoints[allMatches[i].trainIdx].pt.y;
		}
		// RANSAC detect outliers
		auto_ptr< estimators::Solver<Eigen::MatrixXd, Eigen::VectorXd> > ptrSolver(
			new estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>);
		//for (int i = 0; i < (int)good_matches.size(); i++) inliers.push_back(i);
		vector<Eigen::VectorXd> models;

		ransac::Ransac_Handler ransac_fun_Handler;
		bool result = ransac::Ransac_RobustEstimator
			(
			dataPoints, // the input data
			estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::extractor, // How select sampled point from indices
			dataPoints.rows(),  // the number of putatives data
			*(ptrSolver.get()),  // compute the underlying model given a sample set
			estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::defaultEvaluator,  // the function to evaluate a given model
			//Ransac Object that contain function:
			// CandidatesSelector, Sampler and TerminationFunction
			ransac_fun_Handler, // the basic ransac object
			1000,  // the maximum rounds for RANSAC routine
			inliers, // inliers to the final solution
			models, // models array that fit input data
			confidence//0.95 // the confidence want to achieve at the end
			);

		if (models.size() < 1)
		{
			return vector<double>();
		}
		vector<double> model(models[0].size());
		for (size_t i = 0; i < models[0].size(); i++)
		{
			model[i] = models[0][i];
		}
		return model;
	}
}