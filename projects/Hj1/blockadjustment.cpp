#include "blockadjustment.h"

namespace Hj1Optimize{
	namespace blockadjustment{


		ossimDpt Hj1Calibration::getForwardDeriv(ossimHj1Model* hj1Model, const ossimGpt& gpos, int paramIdx,
			int ephIdx/* = 0*/, double hdelta/* = 1e-6*/)
		{
			double den = 0.5 / hdelta;
			ossimDpt res;

			double middle = getExteriorParameter(hj1Model, paramIdx, ephIdx);
			setExteriorParameter(middle + hdelta, hj1Model, paramIdx, ephIdx);
			//res = inverse(gpos);
			res = hj1Model->forward(gpos);
			setExteriorParameter(middle - hdelta, hj1Model, paramIdx, ephIdx);
			//res -= inverse(gpos);
			res -= hj1Model->forward(gpos);
			res = res*den;
			setExteriorParameter(middle, hj1Model, paramIdx, ephIdx);

			return res;
		}

		ossimGpt Hj1Calibration::getInverseDeriv(ossimHj1Model* hj1Model, const ossimDpt& dpos, int paramIdx,
			int ephIdx/* = 0*/, double hdelta/* = 1e-6*/)
		{
			double den = 0.5 / hdelta;
			ossimGpt res, gd;

			double middle = getExteriorParameter(hj1Model, paramIdx, ephIdx);
			setExteriorParameter(middle + hdelta, hj1Model, paramIdx, ephIdx);
			//res = inverse(gpos);
			res = hj1Model->inverse(dpos);
			setExteriorParameter(middle - hdelta, hj1Model, paramIdx, ephIdx);
			//res -= inverse(gpos);
			gd = hj1Model->inverse(dpos);
			setExteriorParameter(middle, hj1Model, paramIdx, ephIdx);

			res.lon = den*(res.lon - gd.lon) * 100000.0; //TBC : approx meters
			res.lat = den*(res.lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
			res.hgt = den*(res.hgt - gd.hgt);

			return res;
		}
		
		// build a error equation according to a gcp
		void Hj1Calibration::buildErrorEquation(const HJ1OptStruct& hjOptStruct, const ossimTieGpt& tiePoint,
			const vector<int>& innerParamList,
			const vector<int>& exteriorParamList, NEWMAT::Matrix &A, NEWMAT::Matrix &B,
			NEWMAT::ColumnVector &L, double pstep_scale)
		{
			int nip = (int)innerParamList.size();
			int nep = exteriorParamList.size();
			int dimObs = 3; //image observation

			A.ReSize(dimObs, nip);
			B.ReSize(dimObs, nep);
			L.ReSize(dimObs);
			//Zeroify matrices that will be accumulated
			A = 0.0;
			B = 0.0;
			L = 0.0;
			
			int c = 0;

			//compute residue
			ossimGpt gd = hjOptStruct.hj1Model->inverse(tiePoint.tie);
			L[0] = (tiePoint.lon - gd.lon) * 100000.0; //approx meters //TBC TBD
			L[1] = (tiePoint.lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
			L[2] = (tiePoint.hgt - gd.hgt); //meters
			//compute all image derivatives regarding parametres for the tie point position
			for (int p = 0; p<nip; p++)
			{
				ossimGpt gdDerpInner;
				gdDerpInner = hjOptStruct.hj1Model->getInverseDeriv(innerParamList[p], tiePoint.tie, pstep_scale);
				A.element(0, p) = gdDerpInner.lon;
				A.element(1, p) = gdDerpInner.lat;
				A.element(2, p) = gdDerpInner.hgt;
			}
			int nEhp = (int)hjOptStruct.hj1Model->theSupportData->theAttitudeSamples.size();
			for (int ii = 0; ii < (int)exteriorParamList.size(); ++ii)
			{
				ossimGpt gdDerpExterior;
				int iParam = exteriorParamList[ii];
				if (iParam < EPH_SPLIT)
				{
					// parameters for a scene
					gdDerpExterior = getInverseDeriv(hjOptStruct.hj1Model, tiePoint.tie, iParam, 0, pstep_scale);
				}
				else
				{
					// parameters for an ephemeris
					for (int iEhp = 0; iEhp < nEhp; ++iEhp)
					{
						gdDerpExterior = getInverseDeriv(hjOptStruct.hj1Model, tiePoint.tie, iParam, 0, pstep_scale);
					}
				}
				B.element(0, ii) = gdDerpExterior.lon;
				B.element(1, ii) = gdDerpExterior.lat;
				B.element(2, ii) = gdDerpExterior.hgt;
			}

		}

		//// build a error equation according to a gcp
		//void buildErrorEquation(const HJ1OptStruct& hjOptStruct, const vector<int>& innerParamList, 
		//	const vector<int>& exteriorParamList, NEWMAT::Matrix &A, NEWMAT::Matrix &B,
		//	NEWMAT::ColumnVector &L, double pstep_scale)
		//{
		//	int nPoints = hjOptStruct.gcpSet->getTiePoints().size();
		//	int nip = (int)innerParamList.size();
		//	int nep = exteriorParamList.size();
		//	int dimObs = 3; //image observation

		//	A.ReSize(dimObs*nPoints, nip);
		//	B.ReSize(dimObs*nPoints, nep);
		//	L.ReSize(dimObs*nPoints);
		//	//Zeroify matrices that will be accumulated
		//	A = 0.0;
		//	B = 0.0;
		//	L = 0.0;

		//	ossimDpt imDerp;
		//	const vector<ossimRefPtr<ossimTieGpt> >& theTPV = hjOptStruct.gcpSet->getTiePoints();
		//	vector<ossimRefPtr<ossimTieGpt> >::const_iterator tit;

		//	int c = 0;
		//	for (tit = theTPV.begin(); tit != theTPV.end(); ++tit)
		//	{
		//		//compute residue
		//		ossimGpt gd = hjOptStruct.hj1Model->inverse((*tit)->tie);
		//		L(c+0) = ((*tit)->lon - gd.lon) * 100000.0; //approx meters //TBC TBD
		//		L(c+1) = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
		//		L(c+2) = (*tit)->hgt - gd.hgt; //meters
		//		//compute all image derivatives regarding parametres for the tie point position
		//		for (int p = 0; p<nip; p++)
		//		{
		//			ossimGpt gdDerpInner;
		//			gdDerpInner = hjOptStruct.hj1Model->getInverseDeriv(innerParamList[p], (*tit)->tie, pstep_scale);
		//			A.element(c + 0, p) = gdDerpInner.lon;
		//			A.element(c + 1, p) = gdDerpInner.lat;
		//			A.element(c + 2, p) = gdDerpInner.hgt;
		//		}
		//		int nEhp = (int)hjOptStruct.hj1Model->theSupportData->theAttitudeSamples.size();
		//		for (int ii = 0; ii < (int)exteriorParamList.size(); ++ii)
		//		{
		//			ossimGpt gdDerpExterior;
		//			int iParam = exteriorParamList[ii];
		//			if (iParam < EPH_SPLIT)
		//			{
		//				// parameters for a scene
		//				gdDerpExterior = getInverseDeriv(hjOptStruct.hj1Model, (*tit)->tie, iParam, 0, pstep_scale);
		//			}
		//			else
		//			{
		//				// parameters for an ephemeris
		//				for (int iEhp = 0; iEhp < nEhp; ++iEhp)
		//				{
		//					gdDerpExterior = getInverseDeriv(hjOptStruct.hj1Model, (*tit)->tie, iParam, 0, pstep_scale);
		//				}
		//			}
		//			B.element(c + 0, ii) = gdDerpExterior.lon;
		//			B.element(c + 1, ii) = gdDerpExterior.lat;
		//			B.element(c + 2, ii) = gdDerpExterior.hgt;
		//		}
		//	}

		//}
		
		void Hj1Calibration::adjustment(const vector< HJ1OptStruct >& hj1OptStructList,
			const vector<int>& innerParamList,
			const vector<int>& exteriorParamList)
		{

			int nImages = (int)hj1OptStructList.size();
			int nip = (int)innerParamList.size();
			int nep = getExteriorParameterNum(hj1OptStructList, exteriorParamList);
			int dimObs = 3;

			int nTotalEquation = 0;
			for (int i = 0; i < nImages; ++i)
			{
				nTotalEquation += dimObs * (int)hj1OptStructList[i].gcpSet->getTiePoints().size();
			}
			if (nImages < 1)
			{
				return;
			}
			int nobs = nTotalEquation / dimObs;

			int Nstate = nip + nep;
			int Nmeasurements = nTotalEquation;
			double *p = new double[Nstate];


			//get current adjustment (between -1 and 1 normally) and convert to ColumnVector
			ossimAdjustmentInfo cadj;
			hj1OptStructList[0].hj1Model->getAdjustment(cadj);
			std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
			NEWMAT::ColumnVector iparm(nip);
			for (int n = 0; n<nip; ++n)
			{
				p[n] = parmlist[innerParamList[n]].getParameter();
			}

			NEWMAT::ColumnVector eparm(nep);
			eparm = getExteriorParameters(hj1OptStructList, exteriorParamList);
			for (size_t i = 0; i < nep; i++)
			{
				p[i+nip] = eparm[i];
			}

			optimzeDataStruct dataStruct;
			dataStruct.pHj1Calibration = this;
			dataStruct.hj1OptStructList = hj1OptStructList;
			dataStruct.innerParamList = innerParamList;
			dataStruct.exteriorParamList = exteriorParamList;

			double opts[SPLM_OPTS_SZ], info[SPLM_INFO_SZ];
			int m, n, ret;
			int nnz;
			m = Nstate;
			n = nTotalEquation;
			nnz = nip * nTotalEquation;
			for (size_t i = 0; i < hj1OptStructList.size(); i++)
			{
				nnz += hj1OptStructList[i].gcpSet->getTiePoints().size() * getImageExteriorParameterNum(hj1OptStructList[i].hj1Model, exteriorParamList);
			}

			opts[0] = SPLM_INIT_MU; opts[1] = SPLM_STOP_THRESH; opts[2] = SPLM_STOP_THRESH;
			opts[3] = SPLM_STOP_THRESH;
			opts[4] = SPLM_DIFF_DELTA; // relevant only if finite difference approximation to Jacobian is used
			opts[5] = SPLM_CHOLMOD; // use CHOLMOD
			//opts[5]=SPLM_PARDISO; // use PARDISO


			ret = sparselm_dercrs(splm_func, splm_jac, p, NULL, m, 0, n, nnz, -1, 1000, opts, info, &dataStruct); // CRS Jacobian


			printf("sparseLM returned %d in %g iter, reason %g\n ", ret, info[5], info[6]);
			double rmse0 = sqrt(info[0] * 2.0 / (double)nTotalEquation);
			double rmse1 = sqrt(info[1] * 2.0 / (double)nTotalEquation);
			printf("Initial RMSE: %g\n Optimized RMSE: %g\n ", rmse0, rmse1);
			//for(i=0; i<m; ++i)
			//	printf("%.7g ", p[i]);
			//printf("\n\nMinimization info:\n");
			//for(i=0; i<SPLM_INFO_SZ; ++i)
			//	printf("%g ", info[i]);
			//printf("\n");

			delete[]p;
		}


		void Hj1Calibration::splm_func(double *p, double *hx, int m, int n, void *adata)
		{
			optimzeDataStruct* pData = (optimzeDataStruct*)adata;
			Hj1Calibration* pThis = pData->pHj1Calibration;
			int nImages = (int)pData->hj1OptStructList.size();
			int nip = (int)pData->innerParamList.size();
			int nep = pThis->getExteriorParameterNum(pData->hj1OptStructList, pData->exteriorParamList);
			int dimObs = 3;

			//double norm2_x = 0.0;

			int nTotalEquation = 0;
			for (int i = 0; i < nImages; ++i)
			{
				nTotalEquation += dimObs * (int)pData->hj1OptStructList[i].gcpSet->getTiePoints().size();
			}


			NEWMAT::ColumnVector iparm(nip);
			for (int image_index = 0; image_index < nImages; ++image_index)
			{
				for (int n = 0; n<nip; ++n)
				{
					pData->hj1OptStructList[image_index].hj1Model->setAdjustableParameter(
						pData->innerParamList[n], iparm(n + 1), false); //do not update now, wait
				}
				pData->hj1OptStructList[image_index].hj1Model->updateModel();
			}

			// update parameters
			NEWMAT::ColumnVector eparm(nep);
			for (size_t i = 0; i < nep; i++)
			{
				eparm[i] = p[i+nip];
			}
			pThis->setExteriorParameters(eparm, pData->hj1OptStructList, pData->exteriorParamList);
			for (int i = 0; i < (int)pData->hj1OptStructList.size(); i++)
			{
				pData->hj1OptStructList[i].hj1Model->updateModel();
			}

			//NEWMAT::ColumnVector newresidue = getResidue(pData->hj1OptStructList);

			double deltap_scale = 1e-4; //step_Scale is 1.0 because we expect parameters to be between -1 and 1
			int nEquationIndex = 0;
			for (std::vector< HJ1OptStruct >::const_iterator iter = pData->hj1OptStructList.begin();
				iter != pData->hj1OptStructList.end(); ++iter)
			{
				for (int j = 0; j < static_cast<int>(iter->gcpSet->getTiePoints().size()); ++j)
				{
					NEWMAT::Matrix A;
					NEWMAT::Matrix B;
					NEWMAT::ColumnVector L;
					// use one gcp to build two error equations
					pThis->buildErrorEquation(*iter, *iter->gcpSet->getTiePoints()[j],
						pData->innerParamList, pData->exteriorParamList, A, B, L, deltap_scale);
					hx[nEquationIndex] = -L.element(0);
					hx[nEquationIndex + 1] = -L.element(1);
					hx[nEquationIndex + 2] = -L.element(2);
					nEquationIndex += 3;
				}
			}
		}

		void Hj1Calibration::splm_jac(double *p, struct splm_crsm *jac, int m, int n, void *adata)
		{
			optimzeDataStruct* pData = (optimzeDataStruct*)adata;
			Hj1Calibration* pThis = pData->pHj1Calibration;
			int nImages = (int)pData->hj1OptStructList.size();
			int nip = (int)pData->innerParamList.size();
			int nep = pThis->getExteriorParameterNum(pData->hj1OptStructList, pData->exteriorParamList);
			int dimObs = 3;

			double norm2_x = 0.0;

			int nTotalEquation = 0;
			for (int i = 0; i < nImages; ++i)
			{
				nTotalEquation += dimObs * (int)pData->hj1OptStructList[i].gcpSet->getTiePoints().size();
			}


			NEWMAT::ColumnVector iparm(nip);
			for (int image_index = 0; image_index < nImages; ++image_index)
			{
				for (int n = 0; n<nip; ++n)
				{
					pData->hj1OptStructList[image_index].hj1Model->setAdjustableParameter(
						pData->innerParamList[n], iparm(n + 1), false); //do not update now, wait
				}
				pData->hj1OptStructList[image_index].hj1Model->updateModel();
			}

			// update parameters
			NEWMAT::ColumnVector eparm(nep);
			for (size_t i = 0; i < nep; i++)
			{
				eparm[i] = p[i + nip];
			}
			pThis->setExteriorParameters(eparm, pData->hj1OptStructList, pData->exteriorParamList);
			for (int i = 0; i < (int)pData->hj1OptStructList.size(); i++)
			{
				pData->hj1OptStructList[i].hj1Model->updateModel();
			}

			int iJacobian = 0;
#define STORE_SPLM_JACOBIAN(col, g)                  \
	do                                      \
			{                                       \
	jac->colidx[ iJacobian ] = col;           \
	jac->val[ iJacobian ] = g;             \
	iJacobian++;                          \
			} while(0)

			double deltap_scale = 1e-4; //step_Scale is 1.0 because we expect parameters to be between -1 and 1
			int nEquationIndex = 0;
			for (int image_index = 0; image_index < nImages; ++image_index)
			{

				for (int j = 0; j < static_cast<int>(pData->hj1OptStructList[image_index].gcpSet->getTiePoints().size()); ++j)
				{
					ossimTieGpt tiePoint = *pData->hj1OptStructList[image_index].gcpSet->getTiePoints()[j];
					NEWMAT::Matrix A;
					NEWMAT::Matrix B;
					NEWMAT::ColumnVector L;
					// use one gcp to build two error equations
					pThis->buildErrorEquation(pData->hj1OptStructList[image_index], tiePoint,
						pData->innerParamList, pData->exteriorParamList, A, B, L, deltap_scale);
					// add the two equations to the block adjustment error equations set

					// In this sample problem, every measurement depends on every element of the
					// state vector, so I loop through all the state vectors here. In practice
					// libdogleg is meant to be applied to sparse problems, where this internal
					// loop would be MUCH shorter than Nstate long
					int nep1 = pThis->getImageExteriorParameterNum(pData->hj1OptStructList[image_index].hj1Model, pData->exteriorParamList);
					//compute all image derivatives regarding parametres for the tie point position
					int pos = 0;
					for (int j = 0; j < image_index; ++j)
					{
						pos += pThis->getImageExteriorParameterNum(pData->hj1OptStructList[image_index].hj1Model, pData->exteriorParamList);
					}
					int nEhp = (int)pData->hj1OptStructList[image_index].hj1Model->theSupportData->theAttitudeSamples.size();
					double val_register;
					jac->rowptr[nEquationIndex] = iJacobian;
					for (int k = 0; k < nip; ++k)
					{
						val_register = A.element(0, k);
						STORE_SPLM_JACOBIAN(k, val_register);
					}
					for (int ii = 0; ii < (int)pData->exteriorParamList.size(); ++ii)
					{
						int iParam = pData->exteriorParamList[ii];
						if (iParam < EPH_SPLIT)
						{
							// parameters for a scene
							val_register = B.element(0, ii);
							STORE_SPLM_JACOBIAN(nip+pos+ii, val_register);
						}
						else
						{
							// parameters for an ephemeris
							for (int iEhp = 0; iEhp < nEhp; ++iEhp)
							{
								val_register = B.element(0, ii);
								STORE_SPLM_JACOBIAN(nip + pos + ii, val_register);
							}
						}
					}
					jac->rowptr[nEquationIndex + 1] = iJacobian;
					for (int k = 0; k < nip; ++k)
					{
						val_register = A.element(1, k);
						STORE_SPLM_JACOBIAN(k, val_register);
					}
					for (int ii = 0; ii < (int)pData->exteriorParamList.size(); ++ii)
					{
						int iParam = pData->exteriorParamList[ii];
						if (iParam < EPH_SPLIT)
						{
							// parameters for a scene
							val_register = B.element(1, ii);
							STORE_SPLM_JACOBIAN(nip + pos + ii, val_register);
						}
						else
						{
							// parameters for an ephemeris
							for (int iEhp = 0; iEhp < nEhp; ++iEhp)
							{
								val_register = B.element(1, ii);
								STORE_SPLM_JACOBIAN(nip + pos + ii, val_register);
							}
						}
					}
					jac->rowptr[nEquationIndex + 2] = iJacobian;
					for (int k = 0; k < nip; ++k)
					{
						val_register = A.element(2, k);
						STORE_SPLM_JACOBIAN(k, val_register);
					}
					for (int ii = 0; ii < (int)pData->exteriorParamList.size(); ++ii)
					{
						int iParam = pData->exteriorParamList[ii];
						if (iParam < EPH_SPLIT)
						{
							// parameters for a scene
							val_register = B.element(2, ii);
							STORE_SPLM_JACOBIAN(nip + pos + ii, val_register);
						}
						else
						{
							// parameters for an ephemeris
							for (int iEhp = 0; iEhp < nEhp; ++iEhp)
							{
								val_register = B.element(2, ii);
								STORE_SPLM_JACOBIAN(nip + pos + ii, val_register);
							}
						}
					}

					nEquationIndex += 3;
				}
			}
			jac->rowptr[nTotalEquation] = iJacobian;
		}

		double Hj1Calibration::getExteriorParameter(ossimHj1Model* hj1Model, int paramIdx, int ephIdx/* = 0*/)
		{
			if (ROLL_OFFSET == paramIdx)
			{
				return hj1Model->theRollOffset;
			}
			else if (PITCH_OFFSET == paramIdx)
			{
				return hj1Model->thePitchOffset;
			}
			else if (YAW_OFFSET == paramIdx)
			{
				return hj1Model->theYawOffset;
			}
			else if (ATT_PITCH_OFFSET == paramIdx)
			{
				return hj1Model->theSupportData->theAttitudeSamples[ephIdx].x;
			}
			else if (ATT_ROLL_OFFSET == paramIdx)
			{
				return hj1Model->theSupportData->theAttitudeSamples[ephIdx].y;
			}
			else if (ATT_YAW_OFFSET == paramIdx)
			{
				return hj1Model->theSupportData->theAttitudeSamples[ephIdx].z;
			}
			else if (ATT_VEL_PITCH_OFFSET == paramIdx)
			{
				return hj1Model->theSupportData->theAttitudeVelSamples[ephIdx].x;
			}
			else if (ATT_VEL_ROLL_OFFSET == paramIdx)
			{
				return hj1Model->theSupportData->theAttitudeVelSamples[ephIdx].y;
			}
			else if (ATT_VEL_YAW_OFFSET == paramIdx)
			{
				return hj1Model->theSupportData->theAttitudeVelSamples[ephIdx].z;
			}

			return 0;
		}

		void Hj1Calibration::setExteriorParameter(double v, ossimHj1Model* hj1Model, int paramIdx, int ephIdx/* = 0*/)
		{
			if (ROLL_OFFSET == paramIdx)
			{
				hj1Model->theRollOffset = v;
			}
			else if (PITCH_OFFSET == paramIdx)
			{
				hj1Model->thePitchOffset = v;
			}
			else if (YAW_OFFSET == paramIdx)
			{
				hj1Model->theYawOffset = v;
			}
			else if (ATT_PITCH_OFFSET == paramIdx)
			{
				hj1Model->theSupportData->theAttitudeSamples[ephIdx].x = v;
			}
			else if (ATT_ROLL_OFFSET == paramIdx)
			{
				hj1Model->theSupportData->theAttitudeSamples[ephIdx].y = v;
			}
			else if (ATT_YAW_OFFSET == paramIdx)
			{
				hj1Model->theSupportData->theAttitudeSamples[ephIdx].z = v;
			}
			else if (ATT_VEL_PITCH_OFFSET == paramIdx)
			{
				hj1Model->theSupportData->theAttitudeVelSamples[ephIdx].x = v;
			}
			else if (ATT_VEL_ROLL_OFFSET == paramIdx)
			{
				hj1Model->theSupportData->theAttitudeVelSamples[ephIdx].y = v;
			}
			else if (ATT_VEL_YAW_OFFSET == paramIdx)
			{
				hj1Model->theSupportData->theAttitudeVelSamples[ephIdx].z = v;
			}
		}

		int Hj1Calibration::getImageExteriorParameterNum(ossimHj1Model* hj1Model, const vector<int>& exteriorParamList)
		{
			int np = 0;
			for (int ii = 0; ii < (int)exteriorParamList.size(); ++ii)
			{
				int iParam = exteriorParamList[ii];
				if (iParam < EPH_SPLIT)
				{
					// parameters for a scene
					np++;
				}
				else
				{
					// parameters for an ephemeris
					np += (int)hj1Model->theSupportData->theAttitudeSamples.size();
				}
			}
			return np;
		}

		int Hj1Calibration::getExteriorParameterNum(const vector< HJ1OptStruct >& hj1OptStructList, const vector<int>& exteriorParamList)
		{
			int nImages = (int)hj1OptStructList.size();
			int np = 0;
			for (int ii = 0; ii < (int)exteriorParamList.size(); ++ii)
			{
				int iParam = exteriorParamList[ii];
				if (iParam < EPH_SPLIT)
				{
					// parameters for a scene
					np += nImages;
				}
				else
				{
					for (int i = 0; i < nImages; ++i)
					{
						// parameters for an ephemeris
						np += (int)hj1OptStructList[i].hj1Model->theSupportData->theAttitudeSamples.size();
					}
				}
			}
			return np;
		}


		NEWMAT::ColumnVector Hj1Calibration::getExteriorParameters(const vector< HJ1OptStruct >& hj1OptStructList, const vector<int>& exteriorParamList)
		{
			int nImages = (int)hj1OptStructList.size();

			int np = getExteriorParameterNum(hj1OptStructList, exteriorParamList);

			NEWMAT::ColumnVector cparm(np);
			int n = 1;
			for (int image_index = 0; image_index < nImages; ++image_index)
			{
				int nEhp = (int)hj1OptStructList[image_index].hj1Model->theSupportData->theAttitudeSamples.size();
				for (int ii = 0; ii < (int)exteriorParamList.size(); ++ii)
				{
					int iParam = exteriorParamList[ii];
					if (iParam < EPH_SPLIT)
					{
						// parameters for a scene
						cparm(n++) = getExteriorParameter(hj1OptStructList[image_index].hj1Model, iParam);
					}
					else
					{
						// parameters for an ephemeris
						for (int iEhp = 0; iEhp < nEhp; ++iEhp)
						{
							cparm(n++) = getExteriorParameter(hj1OptStructList[image_index].hj1Model, iParam, iEhp);
						}
					}
				}
			}
			return cparm;
		}

		void Hj1Calibration::setExteriorParameters(NEWMAT::ColumnVector nparm, const vector< HJ1OptStruct >& hj1OptStructList, const vector<int>& exteriorParamList)
		{
			int nImages = (int)hj1OptStructList.size();

			int np = getExteriorParameterNum(hj1OptStructList, exteriorParamList);

			int n = 1;
			for (int image_index = 0; image_index < nImages; ++image_index)
			{
				int nEhp = (int)hj1OptStructList[image_index].hj1Model->theSupportData->theAttitudeSamples.size();
				for (int ii = 0; ii < (int)exteriorParamList.size(); ++ii)
				{
					int iParam = exteriorParamList[ii];
					if (iParam < EPH_SPLIT)
					{
						// parameters for a scene
						setExteriorParameter(nparm(n++), hj1OptStructList[image_index].hj1Model, iParam);
					}
					else
					{
						// parameters for an ephemeris
						for (int iEhp = 0; iEhp < nEhp; ++iEhp)
						{
							setExteriorParameter(nparm(n++), hj1OptStructList[image_index].hj1Model, iParam, iEhp);
						}
					}
				}
			}
		}
	}
}