#pragma once
#ifndef EZ_RANSAC_H
#define EZ_RANSAC_H

#include <ossim/base/ossimString.h>

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <wtypes.h>
#include <atltypes.h>

#include <string>
#include <vector>
#include <iostream>

#include <sstream>
#include <algorithm>
#include <iterator>
#include <RANSAC/ransac/ransac.h>
#include <RANSAC/estimators/Solver.h>
#include "RANSAC/estimators/affineSolver.h"
#include "RANSAC/estimators/affineError.h"
#include <RANSAC/estimators/rigidSolver.h>
#include <RANSAC/estimators/rigidError.h>

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d/features2d.hpp"

using namespace groupsac;
using namespace std;
using namespace cv;
using namespace std;

namespace mylib{

	vector<double> rigidRansac(const std::vector<KeyPoint>& skeypoints, const std::vector<KeyPoint>& mkeypoints,
		const vector< DMatch >& allMatches, vector<int>& inliers, double confidence = 0.95);


	vector<double> affineRansac(const std::vector<KeyPoint>& skeypoints, const std::vector<KeyPoint>& mkeypoints,
		const vector< DMatch >& allMatches, vector<int>& inliers, double confidence = 0.95);
}

#endif