#pragma once
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <gcpUtil.h>
#include <time.h>
#include "gdal_priv.h"
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>

#ifndef _WIN64
#pragma comment(lib, "ossim20.lib")
#pragma comment(lib, "blas_win32_MT.lib")
#pragma comment(lib, "lapack_win32_MT.lib")
#else
#pragma comment(lib, "ossim20x64.lib")
#pragma comment(lib, "blas_win64_MT.lib")
#pragma comment(lib, "lapack_win64_MT.lib")
#endif

#pragma comment(lib, "ossim_plugin.lib")
#pragma comment(lib, "OpenThreads.lib")
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "mlpack.lib")

#include <mprojectdefine.h>
#include <ossim/base/ossimPreferences.h>
#include <ossim_plugin/radi/radiRpcSolver.h>
#include <ossim_plugin/radi/radiRpcModel.h>
#include <func.h>

using namespace std;
using namespace mylib;
using namespace ossimplugins;

namespace ba{
	enum ParamIndex
	{
		ROLL_OFFSET = 0,
		PITCH_OFFSET,
		YAW_OFFSET,
		EPH_SPLIT,	// the above parameters are identical for a scene, and the below parameters are identical for each ephemeris
		ATT_ROLL_OFFSET,
		ATT_PITCH_OFFSET,
		ATT_YAW_OFFSET,
		ATT_VEL_ROLL_OFFSET,
		ATT_VEL_PITCH_OFFSET,
		ATT_VEL_YAW_OFFSET,
		NUM_ADJUSTABLE_PARAMS // not an index
	};   

	NEWMAT::ColumnVector
		getResidue(const vector< HJ1OptStruct >& hj1OptStructList);

	ossimDpt getForwardDeriv(ossimHj1Model* hj1Model, const ossimGpt& gpos, int paramIdx, int ephIdx = 0, double hdelta = 1e-6);

	void buildNormalEquation(const vector< HJ1OptStruct >& hj1OptStructList, 
		const vector<int>& paramList,
		NEWMAT::SymmetricMatrix& A,
		NEWMAT::ColumnVector& residue,
		NEWMAT::ColumnVector& projResidue,
		double pstep_scale = 1e-6);

	void adjustment(const vector< HJ1OptStruct >& hj1OptStructList, const vector<int>& paramList);

	double getParameter(ossimHj1Model* hj1Model, int paramIdx, int ephIdx = 0);
	void setParameter(double v, ossimHj1Model* hj1Model, int paramIdx, int ephIdx = 0);

	int getParameterNum(const vector< HJ1OptStruct >& hj1OptStructList, const vector<int>& paramList);
	int getImageParameterNum(ossimHj1Model* hj1Model, const vector<int>& paramList);

	NEWMAT::ColumnVector getParameters(const vector< HJ1OptStruct >& hj1OptStructList, const vector<int>& paramList);
	void setParameters(NEWMAT::ColumnVector nparm, const vector< HJ1OptStruct >& hj1OptStructList, const vector<int>& paramList);
}