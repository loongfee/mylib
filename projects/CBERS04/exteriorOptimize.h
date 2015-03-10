#pragma once
#include "Cbers04Optimize.h"
#include <suitesparse/cholmod.h>
#include <splm.h>

namespace Cbers04Optimize{
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
			vector< Cbers04OptStruct > optStructList;
			vector<int> paramList;
		};

		NEWMAT::ColumnVector
			getResidue(const vector< Cbers04OptStruct >& optStructList);

		ossimDpt getForwardDeriv(radiCbers04Model* cbers04Model, const ossimGpt& gpos, int paramIdx, int ephIdx = 0, double hdelta = 1e-6);
		ossimGpt getInverseDeriv(radiCbers04Model* cbers04Model, const ossimDpt& dpos, int paramIdx, int ephIdx = 0, double hdelta = 1e-6);

		void buildNormalEquation(const vector< Cbers04OptStruct >& optStructList,
			const vector<int>& paramList,
			NEWMAT::SymmetricMatrix& A,
			NEWMAT::ColumnVector& residue,
			NEWMAT::ColumnVector& projResidue,
			double pstep_scale = 1e-6);

		void adjustment(const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList);

		double getParameter(radiCbers04Model* cbers04Model, int paramIdx, int ephIdx = 0);
		void setParameter(double v, radiCbers04Model* cbers04Model, int paramIdx, int ephIdx = 0);

		int getParameterNum(const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList);
		int getImageParameterNum(radiCbers04Model* cbers04Model, const vector<int>& paramList);

		NEWMAT::ColumnVector getParameters(const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList);
		void setParameters(NEWMAT::ColumnVector nparm, const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList);

		static void splm_func(double *p, double *hx, int m, int n, void *adata);
		static void splm_jac(double *p, struct splm_crsm *jac, int m, int n, void *adata);
	}
}