#ifndef GROUPSAC_ESTIMATORS_RIGID_ERROR_H
#define GROUPSAC_ESTIMATORS_RIGID_ERROR_H

//#include <armadillo>
//using namespace arma;
#include <Eigen\Eigen>
using namespace Eigen;

namespace groupsac  {
namespace estimators  {

/**
* The distance between the given point and the line y=ax+b, i.e. ax-y+b=0
*
* \param[in] ab (ax+b line equation).
* \param[in] pt (The point on which distance will be computed).  
*
* \return The distance from the point to the line.
* http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
*/
	inline double rigidError(const Eigen::VectorXd & model, const Eigen::VectorXd & obs)  {
		/*
		| x'|   |  a0  a1 |     | x |    | a2 |
		|   | = |         |  *  |   |  + |    |
		| y'|   | -a1  a0 |     | y |    | a3 |
		*/
		double d1 = 1.0*model(2) + obs(0)*model(0) + obs(1)*model(1) - obs(2);
		double d2 = 1.0*model(3) - obs(0)*model(1) + obs(1)*model(0) - obs(3);

	//dist += obs(k);
	return std::sqrt(d1*d1 + d2*d2);
  //return -1.0f; ///Todo(pmoulon) wait the matrix framework
}

}; // namespace estimators
}; // namespace groupsac

#endif //GROUPSAC_ESTIMATORS_POINTTOLINEDIST_H
