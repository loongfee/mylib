#pragma once

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
}