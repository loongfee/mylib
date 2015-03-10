#pragma once
#include "Hj1Optimize.h"
#include <suitesparse/cholmod.h>
#include <splm.h>

namespace Hj1Optimize{
	namespace exteriorOptimize{

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


		struct optimzeDataStruct
		{
			vector< HJ1OptStruct > hj1OptStructList;
			vector<int> paramList;
		};

		NEWMAT::ColumnVector
			getResidue(const vector< HJ1OptStruct >& hj1OptStructList);

		ossimDpt getForwardDeriv(ossimHj1Model* hj1Model, const ossimGpt& gpos, int paramIdx, int ephIdx = 0, double hdelta = 1e-6);
		ossimGpt getInverseDeriv(ossimHj1Model* hj1Model, const ossimDpt& dpos, int paramIdx, int ephIdx = 0, double hdelta = 1e-6);

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

		static void splm_func(double *p, double *hx, int m, int n, void *adata);
		static void splm_jac(double *p, struct splm_crsm *jac, int m, int n, void *adata);
	}
}