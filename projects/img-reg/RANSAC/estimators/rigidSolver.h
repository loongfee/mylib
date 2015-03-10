#ifndef GROUPSAC_ESTIMATORS_RIGID_SOLVER_H
#define GROUPSAC_ESTIMATORS_RIGID_SOLVER_H

#include <cassert>
#include <iostream>
#include <vector>
#include "rigidError.h"
#include "Solver.h"
using namespace std;
#include <Eigen\Eigen>
using namespace Eigen;

namespace groupsac  {
namespace estimators  {

/// Fit a 2D line with to a set of points
/// Specifically, find a and b in the model y = ax + b
///
/// Input data must be typed as follow :
/// X0 Y0
/// X1 Y1
/// X... Y ...
/// Internal extractor function allow to extract sampled data.
template<typename T = Eigen::MatrixXd, typename Model = Eigen::VectorXd>
class rigidSolver : public Solver<T,Model>
{
  /// At least two point are necessary to solve the line equation y = ax+b
  enum { MINIMUM_SAMPLES = 2 };
public :

  int get_MINIMUM_SAMPLES() const {return MINIMUM_SAMPLES;}
  /// See groupsac::estimators::Solver
  bool solve(const T & candidates, vector<Model> & model) const
  {
	  /*
	  | x'|   |  a0  a1 |     | x |    | a2 |
	  |   | = |         |  *  |   |  + |    |
	  | y'|   | -a1  a0 |     | y |    | a3 |	  
	  */
      // Build matrices to solve Ax = b problem:
      Eigen::VectorXd b(candidates.rows() * 2);
      Eigen::MatrixXd A(candidates.rows() * 2, 4);
      //b = candidates.col(6);
	  for(int i = 0;i < (int)candidates.rows();++i)
	  {
		  A(2 * i, 0) = candidates(i, 0);
		  A(2 * i, 1) = candidates(i, 1);
		  A(2 * i, 2) = 1.0;
		  A(2 * i, 3) = 0.0;
		  b(2*i) = candidates(i, 2);

		  A(2 * i + 1, 0) = candidates(i, 1);
		  A(2 * i + 1, 1) = -candidates(i, 0);
		  A(2 * i + 1, 2) = 0.0;
		  A(2 * i + 1, 3) = 1.0;
		  b(2 * i + 1) = candidates(i, 3);
	  }

      // Compute least-squares solution:
      Eigen::VectorXd X;
	  
	  X = (A.transpose()*A).inverse()*A.transpose()*b;
	  model.push_back(X);
      return false;
  }

  /**
  * Return the candidates that are estimated as inliers to the best model
  *
  * \param[in] model (The model(s) that fit the data).
  * \param[in] candidates (The input data).
  * \param[in] threshold (To evaluate if the candidates is an Inlier or Outlier)
  *
  * \return The list of point that are considered as inliers
  */
  static vector<int> defaultEvaluator(vector<Model> & model,
                                      const T & candidates,
                                      double threshold)
  {
	  if(model.size() < 1)
	  {
		  vector<int> inlier;
		  return inlier;
	  }
    assert(model.size() > 0);
    vector< vector<int> > inliers(model.size());
    int bestIndex = 0;
    // For each model compute the number of inliers and the index of the inliers
    // Return the longest inliers vector.
    // must use the pointToLineDist.h file.
    for (size_t i = 0; i < model.size(); ++i)
    {
      const Model & modelToTest = model[i];
      for (size_t j = 0; j < candidates.rows(); ++j)
      {
        double dist = rigidError( modelToTest, candidates.row(j) );

        if ( abs(dist) < threshold)
          inliers[i].push_back(j);
      }
      if ( i > 0 && inliers[bestIndex].size() < inliers[i].size())
      {
        bestIndex = i;
      }
    }
    return inliers[bestIndex];
  }

  /**
    * Extract the sampled indices from the data container.
    *
    * \param[in] data (The input data).
    * \param[in] samples (The indices of data to extract (line or row)).
    *
    * \return The sampled data.
    */
  static T extractor(const T & data, const vector<int> & sampled)
  {
    Eigen::MatrixXd test;
	test.setZero(sampled.size(), data.cols());
    for(size_t i=0; i < sampled.size(); ++i)
      test.row(i) = data.row( sampled[i] );
    return test;
  }
};

}; // namespace estimators
}; // namespace groupsac

#endif //GROUPSAC_ESTIMATORS_LINEFITTINGSOLVER_H
