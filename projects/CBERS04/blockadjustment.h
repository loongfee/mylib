#pragma once
#include "Cbers04Optimize.h"
#include <suitesparse/cholmod.h>
#include <splm.h>

namespace Cbers04Optimize{
	namespace blockadjustment{
		class Cbers04Calibration;

		enum exteriorParamIndex
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
			Cbers04Optimize::blockadjustment::Cbers04Calibration* pCalibration;
			vector< Cbers04OptStruct > optStructList;
			vector<int> innerParamList;
			vector<int> exteriorParamList;
		};

		class Cbers04Calibration
		{
		public:
			Cbers04Calibration(){};
			~Cbers04Calibration(){};
			NEWMAT::ColumnVector
				getResidue(const vector< Cbers04OptStruct >& optStructList);

			ossimDpt getForwardDeriv(radiCbers04Model* cbers04Model, const ossimGpt& gpos, int paramIdx, int ephIdx = 0, double hdelta = 1e-6);
			ossimGpt getInverseDeriv(radiCbers04Model* cbers04Model, const ossimDpt& dpos, int paramIdx, int ephIdx = 0, double hdelta = 1e-6);

			//void adjustment(const vector< Cbers04OptStruct >& optStructList, const vector<int>& innerParamList, const vector<int>& exteriorParamList);
			void adjustment(const vector< Cbers04OptStruct >& optStructList,
				const vector<int>& innerParamList,
				const vector<int>& exteriorParamList);

			void buildErrorEquation(const Cbers04OptStruct& optStruct, const ossimTieGpt& tiePoint,
				const vector<int>& innerParamList,
				const vector<int>& exteriorParamList, NEWMAT::Matrix &A, NEWMAT::Matrix &B,
				NEWMAT::ColumnVector &L, double pstep_scale);

			double getExteriorParameter(radiCbers04Model* cbers04Model, int paramIdx, int ephIdx = 0);
			void setExteriorParameter(double v, radiCbers04Model* cbers04Model, int paramIdx, int ephIdx = 0);

			int getExteriorParameterNum(const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList);
			int getImageExteriorParameterNum(radiCbers04Model* cbers04Model, const vector<int>& paramList);

			NEWMAT::ColumnVector getExteriorParameters(const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList);
			void setExteriorParameters(NEWMAT::ColumnVector nparm, const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList);

			static void splm_func(double *p, double *hx, int m, int n, void *adata);
			static void splm_jac(double *p, struct splm_crsm *jac, int m, int n, void *adata);

		};
	}
}