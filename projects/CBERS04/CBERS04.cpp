#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <gcpUtil.h>
#include "omp.h"
#include "Cbers04Optimize.h"
#include "innerOptimize.h"
#include "exteriorOptimize.h"
#include "blockadjustment.h"
#include <time.h>
//#include "opencv2/core/core.hpp"
//#include "opencv2/features2d/features2d.hpp"
//#include "opencv2/highgui/highgui.hpp"
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

//#pragma comment(lib, "opencv_calib3d248.lib")
//#pragma comment(lib, "opencv_contrib248.lib")
//#pragma comment(lib, "opencv_core248.lib")
//#pragma comment(lib, "opencv_features2d248.lib")
//#pragma comment(lib, "opencv_flann248.lib")
//#pragma comment(lib, "opencv_gpu248.lib")
//#pragma comment(lib, "opencv_highgui248.lib")
//#pragma comment(lib, "opencv_imgproc248.lib")
//#pragma comment(lib, "opencv_legacy248.lib")
//#pragma comment(lib, "opencv_ml248.lib")
//#pragma comment(lib, "opencv_nonfree248.lib")
//#pragma comment(lib, "opencv_objdetect248.lib")
//#pragma comment(lib, "opencv_photo248.lib")
//#pragma comment(lib, "opencv_stitching248.lib")
//#pragma comment(lib, "opencv_ts248.lib")
//#pragma comment(lib, "opencv_video248.lib")
//#pragma comment(lib, "opencv_videostab248.lib")

#include <mprojectdefine.h>
#include <ossim/base/ossimPreferences.h>
#include <ossim_plugin/radi/radiRpcSolver.h>
#include <ossim_plugin/radi/radiRpcModel.h>
#include <ossim_plugin/radi/radiBlockAdjustment.h>
#include <ossim_plugin/radi/radiBlockTieGptSet.h>
#include <func.h>
using namespace std;
using namespace mylib;
using namespace ossimplugins;


//// build a error equation according to a gcp
//void buildErrorEquation (radiCbers04Model* cbers04Model,
//						 const ossimTieGpt& tiePoint, int nType, NEWMAT::Matrix &A,
//										   NEWMAT::Matrix &B, NEWMAT::ColumnVector &L, double pstep_scale)
//{
//	int np = getNumberOfAdjustableParameters();
//	int no = 2; //image observation
//
//	A.ReSize(no,np);
//	B.ReSize(no,3);
//	L.ReSize(no);
//	//Zeroify matrices that will be accumulated
//	A = 0.0;
//	B = 0.0;
//	L = 0.0;	
//
//	ossimDpt imDerp;
//
//	ossimDpt resIm;
//	resIm = tiePoint.tie - forward(tiePoint);
//	L(1) = resIm.x;
//	L(2) = resIm.y;
//	//compute all image derivatives regarding parametres for the tie point position
//	for(int p=0;p<np;++p)
//	{
//		imDerp = getForwardDeriv( p , tiePoint , pstep_scale);
//		A.element(0, p) = imDerp.x;
//		A.element(1, p) = imDerp.y;
//	}
//
//	if(0 == nType)
//	{// if the unknown corresponding points
//		for(int p=0;p<3;++p)
//		{
//			imDerp = getCoordinateForwardDeriv( p , tiePoint , pstep_scale);
//			B.element(0, p) = imDerp.x;
//			B.element(1, p) = imDerp.y;
//		}
//	}
//}

///** 
// * stable invert stolen from ossimRpcSolver
// */
//NEWMAT::Matrix 
//invert(const NEWMAT::SymmetricMatrix& m)
//{
//	int nCol = m.Ncols();
//	int nRow = m.Nrows();
//	Eigen::Matrix2d em(nRow, nCol);
//	, m.Store(), nCol*nRow*sizeof(double));
//	//for (int i = 0;i < nRow;++i)
//	//{
//	//	for(int j = 0;j < nCol;++j)
//	//	{
//	//		double t = m(i+1,j+1);
//	//		em(i,j) = t;
//	//	}
//	//}
//	Eigen::Matrix2d iem = em.inverse();
//	NEWMAT::Matrix result(nRow, nCol);
//	for (int i = 0;i < nRow;++i)
//	{
//		for(int j = 0;j < nCol;++j)
//		{
//			result.element(i,j) = iem(i,j);
//		}
//	}
//	return result;
//}

ossimDpt getForwardDeriv(ossimplugins::radiCbers04Model* cbers04Model,
						 int iAtt, int iComponent, const ossimGpt& gpos, double hdelta)
{   
	double den = 0.5/hdelta;
	ossimDpt res;
	if (0 == iComponent)
	{
		double middle = cbers04Model->theSupportData->theAttitudeBias[iAtt].x;
		cbers04Model->theSupportData->theAttitudeBias[iAtt].x = middle + hdelta;
		//res = inverse(gpos);
		res = cbers04Model->forward(gpos);
		cbers04Model->theSupportData->theAttitudeBias[iAtt].x = middle - hdelta;
		//res -= inverse(gpos);
		res -= cbers04Model->forward(gpos);
		res = res*den;
		cbers04Model->theSupportData->theAttitudeBias[iAtt].x = middle;
	}
	else if(1 == iComponent)
	{
		double middle = cbers04Model->theSupportData->theAttitudeBias[iAtt].y;
		cbers04Model->theSupportData->theAttitudeBias[iAtt].y = middle + hdelta;
		//res = inverse(gpos);
		res = cbers04Model->forward(gpos);
		cbers04Model->theSupportData->theAttitudeBias[iAtt].y = middle - hdelta;
		//res -= inverse(gpos);
		res -= cbers04Model->forward(gpos);
		res = res*den;
		cbers04Model->theSupportData->theAttitudeBias[iAtt].y = middle;
	}
	else if(2 == iComponent)
	{
		double middle = cbers04Model->theSupportData->theAttitudeBias[iAtt].z;
		cbers04Model->theSupportData->theAttitudeBias[iAtt].z = middle + hdelta;
		//res = inverse(gpos);
		res = cbers04Model->forward(gpos);
		cbers04Model->theSupportData->theAttitudeBias[iAtt].z = middle - hdelta;
		//res -= inverse(gpos);
		res -= cbers04Model->forward(gpos);
		res = res*den;
		cbers04Model->theSupportData->theAttitudeBias[iAtt].z = middle;
	}
	else if(3 == iComponent)
	{
		double middle = cbers04Model->theSupportData->theAttitudeVelBias[iAtt].x;
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].x = middle + hdelta;
		//res = inverse(gpos);
		res = cbers04Model->forward(gpos);
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].x = middle - hdelta;
		//res -= inverse(gpos);
		res -= cbers04Model->forward(gpos);
		res = res*den;
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].x = middle;
	}
	else if(4 == iComponent)
	{
		double middle = cbers04Model->theSupportData->theAttitudeVelBias[iAtt].y;
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].y = middle + hdelta;
		//res = inverse(gpos);
		res = cbers04Model->forward(gpos);
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].y = middle - hdelta;
		//res -= inverse(gpos);
		res -= cbers04Model->forward(gpos);
		res = res*den;
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].y = middle;
	}
	else if(5 == iComponent)
	{
		double middle = cbers04Model->theSupportData->theAttitudeVelBias[iAtt].z;
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].z = middle + hdelta;
		//res = inverse(gpos);
		res = cbers04Model->forward(gpos);
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].z = middle - hdelta;
		//res -= inverse(gpos);
		res -= cbers04Model->forward(gpos);
		res = res*den;
		cbers04Model->theSupportData->theAttitudeVelBias[iAtt].z = middle;
	}
	return res;
}

void buildNormalEquation(radiCbers04Model* cbers04Model,
	const ossimTieGptSet& tieSet,
	NEWMAT::SymmetricMatrix& A,
	NEWMAT::ColumnVector& residue,
	NEWMAT::ColumnVector& projResidue,
	double pstep_scale)
{
	//goal:       build Least Squares system
	//constraint: never store full Jacobian matrix in memory (can be huge)
	//            so we build the matrices incrementally
	// the system can be built using forward() or inverse() depending on the projection capabilities : useForward()
	//
	//TBD : add covariance matrix for each tie point	
	//init
	int np = (int)cbers04Model->theSupportData->theAttitudeBias.size()*6;
	int dimObs = 2;
	int no = dimObs * tieSet.size(); //number of observations

	A.ReSize(np);
	residue.ReSize(no);
	projResidue.ReSize(np);
	//Zeroify matrices that will be accumulated
	A           = 0.0;
	projResidue = 0.0;

	const vector<ossimRefPtr<ossimTieGpt> >& theTPV = tieSet.getTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::const_iterator tit;
	unsigned long c=1;

	//image observations 
	std::vector<ossimDpt> imDerp(np);
	ossimDpt resIm;
	// loop on tie points
	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		//compute residue
		resIm = (*tit)->tie - cbers04Model->forward(*(*tit));
		residue(c++) = resIm.x;
		residue(c++) = resIm.y;

		//compute all image derivatives regarding parametres for the tie point position
		for(int p=0;p<np;p+=6)
		{
			int iAtt = p / 6;
			imDerp[p] = getForwardDeriv( cbers04Model, iAtt , 0, *(*tit) , pstep_scale);
			imDerp[p+1] = getForwardDeriv( cbers04Model, iAtt , 1, *(*tit) , pstep_scale);
			imDerp[p+2] = getForwardDeriv( cbers04Model, iAtt , 2, *(*tit) , pstep_scale);
			imDerp[p+3] = getForwardDeriv( cbers04Model, iAtt , 3, *(*tit) , pstep_scale);
			imDerp[p+4] = getForwardDeriv( cbers04Model, iAtt , 4, *(*tit) , pstep_scale);
			imDerp[p+5] = getForwardDeriv( cbers04Model, iAtt , 5, *(*tit) , pstep_scale);
		}

		//compute influence of tie point on all sytem elements
		for(int p1=0;p1<np;++p1)
		{        
			//proj residue: J * residue
			projResidue.element(p1) += imDerp[p1].x * resIm.x + imDerp[p1].y * resIm.y;

			//normal matrix A = transpose(J)*J
			for(int p2=p1;p2<np;++p2)
			{
				A.element(p1,p2) += imDerp[p1].x * imDerp[p2].x + imDerp[p1].y * imDerp[p2].y;
			}
		}
	}
}

NEWMAT::ColumnVector
	getResidue(radiCbers04Model* cbers04Model, const ossimTieGptSet& tieSet, bool useImageObs = true)
{
	NEWMAT::ColumnVector residue;
	int dimObs;
	//bool useImageObs = useForward(); //caching
	if (useImageObs)
	{
		dimObs = 2; //image observation
	} else {
		dimObs = 3; //ground observations
	}
	int no = dimObs * tieSet.size(); //number of observations
	residue.ReSize(no);
	const vector<ossimRefPtr<ossimTieGpt> >& theTPV = tieSet.getTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::const_iterator tit;
	unsigned long c=1;
	if(!cbers04Model->m_proj) return residue;
	if (useImageObs)
	{ 
		ossimDpt resIm;
		for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
		{
			resIm = (*tit)->tie - cbers04Model->forward(**tit);
			residue(c++) = resIm.x;
			residue(c++) = resIm.y;
		}
	} else {
		ossimGpt gd;
		ossimDpt tmplast,tmpnew;	// loong
		for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
		{	
			//loong
			//gd = inverse((*tit)->tie);
			gd = cbers04Model->inverse((*tit)->tie);
			tmplast=cbers04Model->m_proj->forward((*tit)->getGroundPoint());
			tmpnew=cbers04Model->m_proj->forward(gd );
			residue(c++) =  tmplast.lon-tmpnew.lon;
			residue(c++) = tmplast.lat-tmpnew.lat;
			residue(c++) =  (*tit)->hgt - gd.hgt; //TBD : normalize to meters?
			//
			//residue(c++) = ((*tit)->lon - gd.lon) * 100000.0; //approx meters //TBC TBD
			//residue(c++) = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
			//residue(c++) = (*tit)->hgt - gd.hgt; //meters
		}
	} //end of if (useImageObs)
	return residue;
}

double custom_opt(radiCbers04Model* cbers04Model, const ossimTieGptSet& tieSet)
{
	//use a simple Levenberg-Marquardt non-linear optimization
	//note : please limit the number of tie points
	//
	//INPUTS: requires Jacobian matrix (partial derivatives with regards to parameters)
	//OUPUTS: will also compute parameter covariance matrix
	//
	//TBD: use targetVariance!
	
	int np = (int)cbers04Model->theSupportData->theAttitudeBias.size()*6;
	int nobs = tieSet.size();
	
	//setup initail values
	int iter=0;
	int iter_max = 200;
	double minResidue = 1e-10; //TBC
	double minDelta = 1e-10; //TBC
	
	//build Least Squares initial normal equation
	// don't waste memory, add samples one at a time
	NEWMAT::SymmetricMatrix A;
	NEWMAT::ColumnVector residue;
	NEWMAT::ColumnVector projResidue;
	double deltap_scale = 1e-4; //step_Scale is 1.0 because we expect parameters to be between -1 and 1
	buildNormalEquation(cbers04Model, tieSet, A, residue, projResidue, deltap_scale);
	double ki2=residue.SumSquare();
	
	//get current adjustment (between -1 and 1 normally) and convert to ColumnVector
	NEWMAT::ColumnVector cparm(np), nparm(np);
	for(int n=0;n<np;n+=6)
	{
		int iAtt = n / 6;
		cparm(n+1) = cbers04Model->theSupportData->theAttitudeBias[iAtt].x;
		cparm(n+2) = cbers04Model->theSupportData->theAttitudeBias[iAtt].y;
		cparm(n+3) = cbers04Model->theSupportData->theAttitudeBias[iAtt].z;
		cparm(n+4) = cbers04Model->theSupportData->theAttitudeVelBias[iAtt].x;
		cparm(n+5) = cbers04Model->theSupportData->theAttitudeVelBias[iAtt].y;
		cparm(n+6) = cbers04Model->theSupportData->theAttitudeVelBias[iAtt].z;
	}
	
	double damping_speed = 2.0;
	//find max diag element for A
	double maxdiag=0.0;
	for(int d=1;d<=np;++d) {
		if (maxdiag < A(d,d)) maxdiag=A(d,d);
	}
	double damping = 1e-3 * maxdiag;
	double olddamping = 0.0;
	bool found = false;
	
	//DEBUG TBR
	// cout<<"rms="<<sqrt(ki2/nobs)<<" ";
	// cout.flush();
	
	while ( (!found) && (iter < iter_max) ) //non linear optimization loop
	{
		bool decrease = false;
	
		do
		{
			//add damping update to normal matrix
			for(int d=1;d<=np;++d) A(d,d) += damping - olddamping;
			olddamping = damping;
	
			//NEWMAT::ColumnVector deltap = invert(A)*projResidue;
			NEWMAT::ColumnVector deltap = cbers04Model->solveLeastSquares(A, projResidue);
	
			if (deltap.NormFrobenius() <= minDelta) 
			{
				found = true;
			} else {
				//update adjustment
				nparm = cparm + deltap;
				for(int n=0;n<np;n+=6)
				{
					int iAtt = n / 6;
					cbers04Model->theSupportData->theAttitudeBias[iAtt].x = nparm(n+1);
					cbers04Model->theSupportData->theAttitudeBias[iAtt].y = nparm(n+2);
					cbers04Model->theSupportData->theAttitudeBias[iAtt].z = nparm(n+3);
					cbers04Model->theSupportData->theAttitudeVelBias[iAtt].x = nparm(n+4);
					cbers04Model->theSupportData->theAttitudeVelBias[iAtt].y = nparm(n+5);
					cbers04Model->theSupportData->theAttitudeVelBias[iAtt].z = nparm(n+6);
				}
				cbers04Model->updateModel();
	
				//check residue is reduced
				NEWMAT::ColumnVector newresidue = getResidue(cbers04Model, tieSet);
				double newki2=newresidue.SumSquare();
				double res_reduction = (ki2 - newki2) / (deltap.t()*(deltap*damping + projResidue)).AsScalar();
				//DEBUG TBR
				cout<<sqrt(newki2/nobs)<<" ";
				cout.flush();
	
				if (res_reduction > 0)
				{
					//accept new parms
					cparm = nparm;
					ki2=newki2;
	
					deltap_scale = max(1e-15, deltap.NormInfinity()*1e-4);
	
					buildNormalEquation(cbers04Model, tieSet, A, residue, projResidue, deltap_scale);
					olddamping = 0.0;
	
					found = ( projResidue.NormInfinity() <= minResidue );
					//update damping factor
					damping *= max( 1.0/3.0, 1.0-std::pow((2.0*res_reduction-1.0),3));
					damping_speed = 2.0;
					decrease = true;
				} else {
					//cancel parameter update
					for(int n=0;n<np;n+=6)
					{
						int iAtt = n / 6;
						cbers04Model->theSupportData->theAttitudeBias[iAtt].x = nparm(n+1);
						cbers04Model->theSupportData->theAttitudeBias[iAtt].y = nparm(n+2);
						cbers04Model->theSupportData->theAttitudeBias[iAtt].z = nparm(n+3);
						cbers04Model->theSupportData->theAttitudeVelBias[iAtt].x = nparm(n+4);
						cbers04Model->theSupportData->theAttitudeVelBias[iAtt].y = nparm(n+5);
						cbers04Model->theSupportData->theAttitudeVelBias[iAtt].z = nparm(n+6);
					}
					cbers04Model->updateModel();
	
					damping *= damping_speed;
					damping_speed *= 2.0;
				}
			}
		} while (!decrease && !found);
		++iter;
	}
	
	//DEBUG TBR
	cout<<endl;

	for(int n=0;n<np;n+=3)
	{
		int iAtt = n / 6;
		cout<<cbers04Model->theSupportData->theAttitudeBias[iAtt].x
			<<" "<<cbers04Model->theSupportData->theAttitudeBias[iAtt].y
			<<" "<<cbers04Model->theSupportData->theAttitudeBias[iAtt].z
			<<" "<<cbers04Model->theSupportData->theAttitudeVelBias[iAtt].x
			<<" "<<cbers04Model->theSupportData->theAttitudeVelBias[iAtt].y
			<<" "<<cbers04Model->theSupportData->theAttitudeVelBias[iAtt].z
			<<endl;
	}
	cout<<endl;
	
	//compute parameter correlation
	// use normal matrix inverse
	//TBD
	
	return ki2/nobs;
}

//void Cbers04CalibrationTest()
//{
//	ossimInit::instance()->initialize();
//	ossimFilename workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320";
//	workfold = "E:\\HJ1\\HJ-1A_CCD-1_MYC_201403260030_201403260036\\Scene03";
//	MyProject prj;
//	ossimFilename sourceFile = workfold + "\\IMAGE.tif";
//
//	ossimFilename gcpFile = workfold + "\\gcp.txt";
//	ossimFilename demPath = "D:\\workspace\\dem\\srtm90";
//
//	ossimFilename outFile = workfold + "\\rect.tif";
//	ossimFilename reportFile = workfold + "\\report.txt";
//	ossimFilename residualFile = workfold + "\\residual.txt";
//	ModelType modelType = OrbitType;
//
//	prj.m_ImgFileNameUnc = sourceFile;
//	prj.m_ImgFileNameout = outFile;
//	prj.m_ModelType = modelType;
//	prj.ReadGcpAndProjection(ossimFilename(gcpFile));
//
//	prj.theMgr = ossimElevManager::instance();
//	prj.theMgr->loadElevationPath(ossimFilename(demPath));
//	prj.m_DemPath=ossimFilename(demPath);
//
//	prj.GetElevations(prj.m_CtrlGptSet);
//	prj.GetElevations(prj.m_ChkGptSet);
//
//	prj.SavePointToFile(workfold + "\\gcp-hgt.txt", prj.m_CtrlGptSet, NULL);
//
//	prj.InitiateSensorModel();
//	prj.m_sensorModel->m_proj = prj.m_MapPar;
//
//
//	ossimDpt imagepoint,cimagepoint;
//	ossimGpt goundpoint,tGpt;
//
//	int num = static_cast<int>(prj.m_CtrlGptSet->size());
//
//	vector<ossimRefPtr<ossimTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();
//
//	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
//
//
//	//setup initail values
//	int iter=0;
//	int iter_max = 200;
//	double minDelta = 1e-5;
//	double optimizer_delta = 1e10;
//	vector<int> exteriorList;
//	vector<int> innerList;
//	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
//	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
//	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::ROLL_A1_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::PITCH_A1_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::YAW_A1_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::ROLL_A2_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::PITCH_A2_OFFSET);
//	////exteriorList.push_back(radiCbers04Model::AdjustParamIndex::YAW_A2_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
//	////////exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCDLEN_OFFSET);
//	////exteriorList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A0_OFFSET);
//	////exteriorList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A1_OFFSET);
//	////exteriorList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A2_OFFSET);
//	////exteriorList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A3_OFFSET);
//	////exteriorList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A4_OFFSET);
//	//////////exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AX0_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY6_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::X_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::Y_OFFSET);
//	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::Z_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCDLEN_OFFSET);
//	innerList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A0_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A1_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A2_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A3_OFFSET);
//	////innerList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A4_OFFSET);
//	////innerList.push_back(radiCbers04Model::AdjustParamIndex::AX0_OFFSET);
//	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
//	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
//	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
//	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
//	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
//	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY6_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::LINE_A0_OFFSET);
//	//innerList.push_back(radiCbers04Model::AdjustParamIndex::LINE_OFFSET);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::ROLL_A0_OFFSET, 34);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::PITCH_A0_OFFSET, -150);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::YAW_A0_OFFSET, -442);
//	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCDLEN_OFFSET, 34.518);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_A0_OFFSET, 1835);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, -2.206);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, 357.707);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 591.811);
//	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 12.3884);
//	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -1.28087);
//
//
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::ROLL_A0_OFFSET, 33.7765);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::PITCH_A0_OFFSET, -134.602);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::YAW_A0_OFFSET, -441.351);
//	////////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCDLEN_OFFSET, 34.518);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_A0_OFFSET, 1835);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, -2.33876);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, 361.788);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 650.208);
//
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 113.617);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::ROLL_A0_OFFSET, 1.23165);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::PITCH_A0_OFFSET, -947.914);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::YAW_A0_OFFSET, -677.247);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 3.03778);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, 75.352);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 12.4454);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 477.775);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, 2.9194);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 371.127);
//	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY6_OFFSET, 387.375);
//
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 67.6818);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::ROLL_A0_OFFSET, 1.23165);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::PITCH_A0_OFFSET, -947.914);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::YAW_A0_OFFSET, -677.247);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 3.03778);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, 42.7842);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 12.7577);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 427.358);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -12.8224);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 456.429);
//
//	//// 1B-CCD1
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, -94.5678);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::ROLL_A0_OFFSET, 0.791816);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::PITCH_A0_OFFSET, 793.844);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::YAW_A0_OFFSET, 635.167);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 3.04783);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -71.931);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 10.2021);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 300.976);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -8.41909);
//	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 259.129);
//
//	// 1A-CCD1
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 132.172);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, -0.178295);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 702.947);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, 663.609);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 2.57239);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, 88.9769);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 4.56039);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 488.326);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, 93.6162);
//	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 322.148);
//	prj.m_sensorModel->updateModel();
//	prj.m_sensorModel->saveState(prj.geom);
//	ossimAdjustmentInfo cadj;
//	prj.m_sensorModel->getAdjustment(cadj);
//	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
//	int np = (int)parmlist.size();
//	NEWMAT::ColumnVector old_parm(np), new_parm(np);
//	for(int n=0;n<np;++n)
//	{
//		old_parm(n+1) = parmlist[n].getParameter();
//	}
//	while(optimizer_delta > minDelta && iter < iter_max)
//	{
//
//		for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
//		{
//			imagepoint=(*tit)->getImagePoint();
//			goundpoint=(*tit)->getGroundPoint();
//			tGpt = prj.m_MapPar->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
//			tGpt.hgt = (*tit)->hgt;
//			(*tit)->setGroundPoint(tGpt);
//		}
//		
//		//// do inner orientation
//		//prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, innerList);
//		//prj.m_sensorModel->updateModel();
//		//prj.m_sensorModel->saveState(prj.geom);
//
//		//// do exterior orientation
//		//prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, exteriorList);
//		//prj.m_sensorModel->updateModel();
//		//prj.m_sensorModel->saveState(prj.geom);
//
//		//custom_opt((radiCbers04Model*)(prj.m_sensorModel), *prj.m_CtrlGptSet);
//		////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, 0);
//		//prj.m_sensorModel->updateModel();
//		//prj.m_sensorModel->saveState(prj.geom);
//
//		ossimAdjustmentInfo cadj;
//		prj.m_sensorModel->getAdjustment(cadj);
//		std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
//		for(int n=0;n<np;++n)
//		{
//			new_parm(n+1) = parmlist[n].getParameter();
//			cout<<new_parm(n+1)<<endl;
//		}
//		// then calculate the change value of the optimizers
//		optimizer_delta = (new_parm - old_parm).NormInfinity();
//		old_parm = new_parm;
//		cout<<"iteration "<<1+iter++<<" :"<<optimizer_delta<<endl;
//
//		for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
//		{
//			ossimDpt dpt = prj.m_sensorModel->m_proj->forward(*(*tit));
//			ossimGpt gpt(dpt.x,dpt.y);
//			(*tit)->setGroundPoint(ossimGpt(dpt.x,dpt.y,(*tit)->hgt));
//		}
//
//		fstream fs;
//		fs.open(reportFile.c_str(), ios_base::out);
//		fs.setf(ios::fixed, ios::floatfield);
//		fs.precision(2);
//		prj.OutputReport(fs, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
//		fs.close();
//
//		fstream ofs;
//		ofs.open(residualFile.c_str(), ios_base::out);
//		for (int i = 0;i < prj.m_CtrlGptSet->getTiePoints().size();++i)
//		{
//			ossimGpt ll;
//			prj.m_sensorModel->lineSampleToWorld(prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint(), ll);
//			ossimDpt gpt = prj.m_MapPar->forward(ll);
//			double delta_lat = gpt.x - prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint().lat;
//			double delta_lon = gpt.y - prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint().lon;
//			if (0 != i)
//			{
//				ofs<<endl;
//			}
//			ofs<<i+1<<"\t"
//				<<prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().x<<"\t"
//				<<prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().y<<"\t"
//				//<<gptSet->getTiePoints()[i]->getGroundPoint().lat<<"\t"
//				//<<gptSet->getTiePoints()[i]->getGroundPoint().lon<<"\t"
//				//<<gptSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
//				<<delta_lat<<"\t"
//				<<delta_lon;
//		}
//		ofs.close();
//	}// end while
//	//while(optimizer_delta > minDelta && iter < iter_max)
//	//{
//
//	//	// do inner orientation
//	//	eigen_levmar_optimization(prj.m_sensorModel, prj.m_CtrlGptSet, exteriorList);
//	//	prj.m_sensorModel->updateModel();
//	//	prj.m_sensorModel->saveState(prj.geom);
//
//	//	// do exterior orientation
//	//	alglib_optimization(prj.m_sensorModel, prj.m_CtrlGptSet, exteriorList);
//	//	prj.m_sensorModel->updateModel();
//	//	prj.m_sensorModel->saveState(prj.geom);
//
//
//
//	//	ossimAdjustmentInfo cadj;
//	//	prj.m_sensorModel->getAdjustment(cadj);
//	//	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
//	//	for(int n=0;n<np;++n)
//	//	{
//	//		new_parm(n+1) = parmlist[n].getParameter();
//	//		cout<<new_parm(n+1)<<endl;
//	//	}
//	//	// then calculate the change value of the optimizers
//	//	optimizer_delta = (new_parm - old_parm).NormInfinity();
//	//	old_parm = new_parm;
//	//	cout<<"iteration "<<1+iter++<<" :"<<optimizer_delta<<endl;
//	//}// end while
//
//
//
//	prj.m_OutBandList.clear();
//	prj.m_OutBandList.push_back(0);
//	system("pause");
//
//	//prj.starline	=	00000;
//	//prj.starpixel	=	0000;
//	//prj.endline		=	11999;
//	//prj.endpixel	=	11999;
//	prj.Orthograph(outFile);
//}

//void create3DGridPoints(const ossimDrect& imageBounds,
//						const ossimSensorModel& proj,
//						double height,
//						ossimTieGptSet*& pTieGptSet,
//						ossim_uint32 xSamples,
//						ossim_uint32 ySamples,
//						bool latlon,
//						bool shiftTo0Flag)
//{
//	if (NULL == pTieGptSet)
//	{
//		pTieGptSet = new ossimTieGptSet;
//	}
//	ossim_uint32 x,y;
//	ossim_float64 w = imageBounds.width();
//	ossim_float64 h = imageBounds.height();
//	ossimGpt gpt;
//	ossimGpt defaultGround;
//	if(ySamples < 1) ySamples = 12;
//	if(xSamples < 1) xSamples = 12;
//	srand(time(0));
//	double xnorm;
//	double ynorm;
//	ossimDpt ul = imageBounds.ul();
//	ossimDpt shiftTo0(-ul.x,
//		-ul.y);
//	for(y = 0; y < ySamples; ++y)
//	{
//		for(x = 0; x < xSamples; ++x)
//		{
//			ossimDpt imagePoint;
//			if(ySamples > 1)
//			{
//				ynorm = (double)y/(double)(ySamples - 1);
//			}
//			else
//			{
//				ynorm = 0.0;
//			}
//			if(xSamples > 1)
//			{
//				xnorm = (double)x/(double)(xSamples - 1);
//			}
//			else
//			{
//				xnorm = 0.0;
//			}
//
//			ossimDpt dpt((w-1)*xnorm + ul.x,
//				(h-1)*ynorm + ul.y);
//
//			proj.lineSampleHeightToWorld(dpt, height, gpt);
//			//proj.worldToLineSample(gpt, dpt);
//
//			ossimDpt dpt1;
//			ossimDpt dpt2;
//			proj.worldToLineSample(proj.m_proj->inverse(ossimDpt(844063.7962,5300035.258)), dpt1);
//			proj.worldToLineSample(proj.m_proj->inverse(ossimDpt(846431.7615,5299320.665)), dpt2);
//			ossimGpt gpt1;
//			ossimGpt gpt2;
//			proj.lineSampleHeightToWorld(dpt1, height, gpt1);
//			proj.lineSampleHeightToWorld(dpt2, height, gpt2);
//
//
//			proj.lineSampleHeightToWorld(ossimDpt(3999.666667,43984.72889), height, gpt1);
//			proj.lineSampleHeightToWorld(ossimDpt(3999.666667,43985.72963), height, gpt2);
//
//			//proj.lineSampleHeightToWorld(dpt, height, gpt);
//			//ossimGpt gptTemp;
//			//proj.lineSampleToWorld(dpt,
//			//	gptTemp);
//			//ossimDpt dptTemp;
//			//proj.worldToLineSample(gptTemp, dptTemp);
//			//ossimDpt res = dptTemp - dpt;
//			//gpt.changeDatum(defaultGround.datum());
//			if(shiftTo0Flag)
//			{
//				imagePoint = dpt + shiftTo0;
//			}
//			else
//			{
//				imagePoint = dpt;
//			}
//			//gpt.hgt = height;
//
//			if (latlon)
//			{
//				dpt = ossimDpt(gpt.lon, gpt.lat);
//			}
//			else
//			{
//				dpt = proj.m_proj->forward(gpt);
//			}
//			ossimGpt groundPoint(dpt.x,dpt.y, gpt.hgt);
//
//			ossimString strId;
//			char tmpStr[256];
//			sprintf(tmpStr, "%d", y * xSamples + x + 1);
//			strId = tmpStr;
//			ossimTieGpt *aTiePt = new ossimTieGpt(groundPoint, imagePoint, 1.0, strId);
//			pTieGptSet->addTiePoint(aTiePt);
//		}
//	}
//}

void HJ1_create3DGridPoints(ossimFilename workfold = "")
{
	//ossimFilename workfold = "E:\\beijing";
	//ossimFilename workfold = "E:\\HJ1\\test";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401050432_201401050443";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_SYC_201401050300_201401050309";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401070348_201401070354";
	if (workfold.empty())
	{
		workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320";
	}
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\IMAGE.tif";
	//ossimFilename sourceFile = workfold + "\\IMAGE_1.tif";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";

	ossimFilename gcpFile = workfold + "\\gcp.txt";
	ossimFilename demPath = "D:\\workspace\\dem\\srtm90";

	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report.txt";

	ossimFilename virtual_gcpfile = workfold + "\\virtual_gcp.txt";
	ossimFilename virtual_chkfile = workfold + "\\virtual_chk.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpFile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	//prj.SavePointToFile(gcp_hgt_File, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);

	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 67.6818);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 1.23165);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, -947.914);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, -677.247);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 3.03778);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, 42.7842);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 12.7577);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 427.358);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -12.8224);
	prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 456.429);
	prj.m_sensorModel->updateModel();
	prj.m_sensorModel->saveState(prj.geom);

	prj.OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);


	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourceFile);
	if(!handler) return;   //应该弹出警告对话框
	ossimTieGptSet* gptSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	double max_height = 4000.0;
	double min_height = 0.0;

	//ossimGpt gpt;
	//ossimDpt dpt;
	//prj.m_sensorModel->lineSampleToWorld(ossimDpt(0.0000000000, 34000.0000000000), gpt);
	//prj.m_sensorModel->worldToLineSample(gpt, dpt);
	//int nLevels = 5;
	//for (int i=0;i < nLevels;++i)
	//{
	//	double hgt = min_height + i*(max_height-min_height)/(nLevels-1);
	//	create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 15, 15, true, false);
	//}

	//ossimTieGptSet* chkSet = new ossimTieGptSet;
	//nLevels = 10;
	//for (int i=0;i < nLevels;++i)
	//{
	//	double hgt = min_height + i*(max_height-min_height)/(nLevels-1);
	//	create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, chkSet, 30, 30, true, false);
	//}

	int nLevels = 1;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height;// + i*(max_height-min_height)/(nLevels-1);
		//create3DGridPoints(ossimIrect(0, 43000, 11999, 44351), *prj.m_sensorModel, hgt, gptSet, 12000, 3, false, false);
		create3DGridPoints(ossimIrect(0, 40000, 11999, 44351), *prj.m_sensorModel, hgt, gptSet, 1, 4351, false, false);
	}

	nLevels = 6;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height + i*(max_height-min_height)/(nLevels-1);
		create3DGridPoints(ossimIrect(0, 30000, 11999, 33000), *prj.m_sensorModel, hgt, chkSet, 21, 21, false, false);
	}
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 0.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 1000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 2000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 3000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 4000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 5000.0, gptSet, 10, 10, false);
	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 20, 20, true);
	//gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 11, 11, false);
	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 21, 21, false);
	handler->close();
	delete handler;

	prj.SavePointToFile(virtual_gcpfile, gptSet, NULL);
	//prj.OutputReport(reportFile, prj.m_sensorModel, gptSet, prj.m_ChkGptSet, true);
	prj.SavePointToFile(virtual_chkfile, chkSet, NULL);
}

ossimplugins::radiRpcModel* Cbers04CreateRpcs(ossimSensorModel* sensorModel, ossimIrect theImageClipRect, ossimFilename workfold)
{
	ossimFilename projectionFile = workfold + "\\gcp.geom";
	ossimFilename reportFile = workfold + "\\report0.txt";
	ossimFilename rpcFile = workfold + "\\rpcStruct.txt";

	// statistic the max and min Height
	ossim_uint32 xSamples = 20;
	ossim_uint32 ySamples = 20;

	double max_Height = -1.0e10;
	double min_Height = 1.0e10;
	ossim_uint32 x,y;
	ossim_float64 w = theImageClipRect.width();
	ossim_float64 h = theImageClipRect.height();
	ossimGpt gpt;
	ossimGpt defaultGround;
	if(ySamples < 1) ySamples = 12;
	if(xSamples < 1) xSamples = 12;

	double xnorm;
	double ynorm;
	ossimDpt ul = theImageClipRect.ul();
	ossimDpt shiftTo0(-ul.x, -ul.y);
	for(y = 0; y < ySamples; ++y)
	{
		for(x = 0; x < xSamples; ++x)
		{
			ossimDpt imagePoint;
			if(ySamples > 1)
			{
				ynorm = (double)y/(double)(ySamples - 1);
			}
			else
			{
				ynorm = 0.0;
			}
			if(xSamples > 1)
			{
				xnorm = (double)x/(double)(xSamples - 1);
			}
			else
			{
				xnorm = 0.0;
			}

			ossimDpt dpt((w-1)*xnorm + ul.x,
				(h-1)*ynorm + ul.y);

			sensorModel->lineSampleToWorld(dpt, gpt);
			if (!gpt.hasNans())
			{
				if (gpt.height() > max_Height)
				{
					max_Height = gpt.height();
				}
				if (gpt.height() < min_Height)
				{
					min_Height = gpt.height();
				}
			}
		}
	}

	min_Height -= 500.0;
	max_Height += 500.0;
	ossimTieGptSet *gptSet = NULL;
	int nLevels = 5;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_Height + i*(max_Height-min_Height)/(nLevels-1);
		create3DGridPoints(theImageClipRect, *sensorModel, hgt, gptSet, 15, 15, true, false);
	}
	ossimTieGptSet *chkSet = NULL;
	nLevels = 7;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_Height + i*(max_Height-min_Height)/(nLevels-1);
		create3DGridPoints(theImageClipRect, *sensorModel, hgt, chkSet, 21, 21, true, false);
	}

	ossimplugins::radiRpcSolver *solver = new ossimplugins::radiRpcSolver(true, false);
	int num = (int)gptSet->getTiePoints().size();
	vector < ossimDpt > imagePoints(num);
	vector < ossimGpt > groundControlPoints(num);
	for(int i = 0;i < num;i++)
	{
		groundControlPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
		imagePoints[i] = gptSet->getTiePoints()[i]->getImagePoint();
	}
	//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LS, 1e-8);
	//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::RIDGE, 1e-6);
	//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LASSO, 1e-5);
	clock_t  clockBegin1, clockEnd1;
	clockBegin1 = clock();
	solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LASSO, 1e-4);
	clockEnd1 = clock();
	printf("l1-norm regularized least squares Time consuming: %lf s\n", (clockEnd1 - clockBegin1)*1e-3);
	clock_t  clockBegin2, clockEnd2;
	clockBegin2 = clock();
	solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LS, 1e-8);
	clockEnd2 = clock();
	printf("Least squares Time consuming: %lf s\n", (clockEnd2 - clockBegin2)*1e-3);

	ossimImageGeometry *imageGeom = solver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimplugins::radiRpcModel* rpcModel = new ossimplugins::radiRpcModel();
	rpcModel->loadState(geom, "projection.");

	ossimKeywordlist kwl;
	kwl.addFile(projectionFile);
	rpcModel->m_proj = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(kwl));

	fstream rpcStructFile;
	rpcStructFile.open(rpcFile.c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStructFile);
	rpcStructFile.close();
	//OutputReport1(reportFile0, prj.m_sensorModel, gptSet, chkSet);
	//mylib::OutputReport(reportFile, rpcModel, gptSet, chkSet, false, true);
	mylib::OutputReport(reportFile, rpcModel, gptSet, chkSet, true, true);
	return rpcModel;
}

void Cbers04CreateRpcs()
{
	ossimInit::instance()->initialize();
	ossimFilename workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320";
	//workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01";
	workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02";
	//workfold = "E:\\HJ1\\HJ-1A_CCD-1_MYC_201403260030_201403260036";
	ossimFilename sourceFile = workfold + "\\IMAGE.tif";
	ossimFilename demPath = "D:\\workspace\\dem\\srtm90";
	ossimElevManager::instance()->loadElevationPath(ossimFilename(demPath));
	ossimRefPtr<ossimImageHandler> handler   = ossimImageHandlerRegistry::instance()->open(sourceFile);

	string strPath(sourceFile.c_str());
	strPath = SBeforeLast(strPath, '\\') + "\\";
	ossimString descfile = (strPath + "desc.xml");

	radiCbers04Model* sensorModel = new radiCbers04Model();
	sensorModel->setupOptimizer(descfile);
	//ossimKeywordlist kwl;
	//kwl.addFile(projectionFile);
	//sensorModel->m_proj = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(kwl));
	ossimIrect theImageClipRect = handler->getBoundingRect();
	Cbers04CreateRpcs(sensorModel, theImageClipRect, workfold);
	////workfold = "I:\\wuhan\\spot\\2005\\2.5";
	////sourceFile = workfold + "\\SCENE01\\imagery.tif";
	//ossimFilename projectionFile = workfold + "\\projection.txt";
	//ossimFilename reportFile0 = workfold + "\\report0.txt";
	//ossimFilename reportFile1 = workfold + "\\report1.txt";
	//ossimFilename rpcFile = workfold + "\\rpcStruct.txt";


	//// statistic the max and min Height
	//ossim_uint32 xSamples = 20;
	//ossim_uint32 ySamples = 20;

	//double max_Height = -1.0e10;
	//double min_Height = 1.0e10;
	//ossim_uint32 x,y;
	//ossim_float64 w = theImageClipRect.width();
	//ossim_float64 h = theImageClipRect.height();
	//ossimGpt gpt;
	//ossimGpt defaultGround;
	//if(ySamples < 1) ySamples = 12;
	//if(xSamples < 1) xSamples = 12;

	//double xnorm;
	//double ynorm;
	//ossimDpt ul = theImageClipRect.ul();
	//ossimDpt shiftTo0(-ul.x, -ul.y);
	//for(y = 0; y < ySamples; ++y)
	//{
	//	for(x = 0; x < xSamples; ++x)
	//	{
	//		ossimDpt imagePoint;
	//		if(ySamples > 1)
	//		{
	//			ynorm = (double)y/(double)(ySamples - 1);
	//		}
	//		else
	//		{
	//			ynorm = 0.0;
	//		}
	//		if(xSamples > 1)
	//		{
	//			xnorm = (double)x/(double)(xSamples - 1);
	//		}
	//		else
	//		{
	//			xnorm = 0.0;
	//		}

	//		ossimDpt dpt((w-1)*xnorm + ul.x,
	//			(h-1)*ynorm + ul.y);

	//		sensorModel->lineSampleToWorld(dpt, gpt);
	//		if (!gpt.hasNans())
	//		{
	//			if (gpt.height() > max_Height)
	//			{
	//				max_Height = gpt.height();
	//			}
	//			if (gpt.height() < min_Height)
	//			{
	//				min_Height = gpt.height();
	//			}
	//		}
	//	}
	//}

	//min_Height -= 500.0;
	//max_Height += 500.0;
	//ossimTieGptSet *gptSet = NULL;
	//int nLevels = 5;
	//for (int i=0;i < nLevels;++i)
	//{
	//	double hgt = min_Height + i*(max_Height-min_Height)/(nLevels-1);
	//	create3DGridPoints(theImageClipRect, *sensorModel, hgt, gptSet, 15, 15, true, false);
	//}
	//ossimTieGptSet *chkSet = NULL;
	//nLevels = 7;
	//for (int i=0;i < nLevels;++i)
	//{
	//	double hgt = min_Height + i*(max_Height-min_Height)/(nLevels-1);
	//	create3DGridPoints(theImageClipRect, *sensorModel, hgt, chkSet, 21, 21, true, false);
	//}

	//ossimplugins::radiRpcSolver *solver = new ossimplugins::radiRpcSolver(true, false);
	//int num = (int)gptSet->getTiePoints().size();
	//vector < ossimDpt > imagePoints(num);
	//vector < ossimGpt > groundControlPoints(num);
	//for(int i = 0;i < num;i++)
	//{
	//	groundControlPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
	//	imagePoints[i] = gptSet->getTiePoints()[i]->getImagePoint();
	//}
	////solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LS, 1e-8);
	////solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::RIDGE, 1e-5);
	//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LASSO, 1e-5);

	//ossimImageGeometry *imageGeom = solver->createRpcModel();
	//ossimKeywordlist geom;
	//imageGeom->saveState(geom);
	//ossimRpcModel* rpcModel = new ossimRpcModel();
	//rpcModel->loadState(geom, "projection.");


	//fstream rpcStructFile;
	//rpcStructFile.open(rpcFile.c_str(), ios_base::out);
	//rpcModel->saveRpcModelStruct(rpcStructFile);
	//rpcStructFile.close();
	////OutputReport1(reportFile0, prj.m_sensorModel, gptSet, chkSet);
	//mylib::OutputReport(reportFile1, rpcModel, gptSet, chkSet, true, true);
}

void LineOffsetTest()
{
	ossimInit::instance()->initialize();
	ossimFilename workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320";
	workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01";
	workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02";
	workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01\\L2\\L1";
	//workfold = "E:\\HJ1\\HJ-1A_CCD-1_MYC_201403260030_201403260036";
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\IMAGE.tif";

	ossimFilename gcpFile = workfold + "\\gcp.txt";
	ossimFilename demPath = "D:\\workspace\\dem\\srtm90";

	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report.txt";
	ossimFilename residualFile = workfold + "\\residual.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpFile));

	prj.theMgr = ossimElevManager::instance();
	//prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	//prj.GetElevations(prj.m_CtrlGptSet);
	//prj.GetElevations(prj.m_ChkGptSet);
	
	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;


	ossimDpt imagepoint,cimagepoint;
	ossimGpt goundpoint,tGpt;

	int num = static_cast<int>(prj.m_CtrlGptSet->size());

	//setup initail values
	int iter=0;
	int iter_max = 200;
	double minDelta = 1e-1;
	double optimizer_delta = 1e10;

	vector<int> exteriorList;
	vector<int> innerList;
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);

	innerList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AX0_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AX1_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AX2_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AX3_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);

	vector<Cbers04Optimize::Cbers04OptStruct> optStructList;
	Cbers04Optimize::Cbers04OptStruct cbers04Struct((radiCbers04Model*)prj.m_sensorModel, prj.m_CtrlGptSet);
	optStructList.push_back(cbers04Struct);

	vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[0].gcpSet->refTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
	// loop on tie points
	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	int step = 100;
	int minLineOffset = 2000;
	int maxLineOffset = 3000;
	for (int iLineOffset = minLineOffset;iLineOffset <= maxLineOffset;iLineOffset += step)
	{
		string strPath(sourceFile.c_str());
		strPath = SBeforeLast(strPath, '\\') + "\\";
		ossimString descfile = (strPath + "desc.xml");

		prj.m_sensorModel = new ossimplugins::radiCbers04Model();
		prj.m_sensorModel->setupOptimizer(descfile);
		prj.m_sensorModel->m_proj = prj.m_MapPar;
		prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, iLineOffset);
		optStructList[0].cbers04Model = (radiCbers04Model*)prj.m_sensorModel;

		cout<<"iLineOffset = "<<iLineOffset<<":\n";
		//Cbers04Optimize::innerOptimize::adjustment(optStructList, exteriorList);
		Cbers04Optimize::innerOptimize::adjustment(optStructList, innerList);

		prj.m_sensorModel->saveState(prj.geom);

		ossimRefPtr<ossimImageHandler> handler   = ossimImageHandlerRegistry::instance()->open(sourceFile);
		ossimIrect theImageClipRect = handler->getBoundingRect();
		char buf[2048];
		sprintf_s(buf, "%s\\test_results\\%d", workfold.c_str(), iLineOffset);
		string newFold(buf);
		if (!ossimFilename(buf).exists())
		{
			_mkdir(buf);
		}
		Cbers04CreateRpcs(prj.m_sensorModel, theImageClipRect, newFold);
	}
}

//using namespace cv;
//int writeMatToFile(const cv::Mat &I, string path) {
//
//	//load the matrix size
//	int matWidth = I.size().width, matHeight = I.size().height;
//
//	//read type from Mat
//	int type = I.type();
//
//	//declare values to be written
//	float fvalue;
//	double dvalue;
//	cv::Vec3f vfvalue;
//	cv::Vec3d vdvalue;
//
//	//create the file stream
//	fstream file;
//	file.open(path.c_str(), ios::out | ios::binary );
//	if (!file)
//		return -1;
//
//	//write type and size of the matrix first
//	file.write((const char*) &type, sizeof(type));
//	file.write((const char*) &matWidth, sizeof(matWidth));
//	file.write((const char*) &matHeight, sizeof(matHeight));
//
//	//write data depending on the image's type
//	switch (type)
//	{
//	default:
//		cout << "Error: wrong Mat type: must be CV_32F, CV_64F, CV_32FC3 or CV_64FC3" << endl;
//		break;
//		// FLOAT ONE CHANNEL
//	case CV_32F:
//		cout << "Writing CV_32F image" << endl;
//		for (int i=0; i < matWidth*matHeight; ++i) {
//			fvalue = I.at<float>(i);
//			file.write((const char*) &fvalue, sizeof(fvalue));
//		}
//		break;
//		// DOUBLE ONE CHANNEL
//	case CV_64F:
//		cout << "Writing CV_64F image" << endl;
//		for (int i=0; i < matWidth*matHeight; ++i) {
//			dvalue = I.at<double>(i);
//			file.write((const char*) &dvalue, sizeof(dvalue));
//		}
//		break;
//
//		// FLOAT THREE CHANNELS
//	case CV_32FC3:
//		cout << "Writing CV_32FC3 image" << endl;
//		for (int i=0; i < matWidth*matHeight; ++i) {
//			vfvalue = I.at<cv::Vec3f>(i);
//			file.write((const char*) &vfvalue, sizeof(vfvalue));
//		}
//		break;
//
//		// DOUBLE THREE CHANNELS
//	case CV_64FC3:
//		cout << "Writing CV_64FC3 image" << endl;
//		for (int i=0; i < matWidth*matHeight; ++i) {
//			vdvalue = I.at<cv::Vec3d>(i);
//			file.write((const char*) &vdvalue, sizeof(vdvalue));
//		}
//		break;
//
//	}
//
//	//close file
//	file.close();
//
//	return 0;
//}

void Cbers04CalibrationTest(ossimFilename workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320")
{
	ossimInit::instance()->initialize();
	//ossimFilename workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320";
	//workfold = "E:\\HJ1\\HJ-1B_CCD-1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01";
	//workfold = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201405250250_201405250300";
	//workfold = "E:\\HJ1\\HJ-1B_CCD-1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02";
	//workfold = "E:\\HJ1\H\J-1A_CCD-1_MYC_201403260030_201403260036";
	//workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01\\L2\\L1";
	//workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201405190109_201405190117";
	//workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250";
	//workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140";
	//workfold = "E:\\HJ1\\HJ-1A_CCD-1_MYC_201403240114_201403240125";
	//workfold = "E:\\HJ1\\HJ-1A_CCD-1_MYC_201403260030_201403260036";
	//workfold = "E:\\HJ1\\HJ-1A_CCD-1_MYC_201404110217_201404110229";
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\IMAGE_B3.tif";

	ossimFilename gcpFile = workfold + "\\gcp.txt";
	ossimFilename demPath = "D:\\workspace\\dem\\srtm90";

	ossimFilename rpcFile0 = sourceFile;
	rpcFile0.setExtension("rpb0");
	ossimFilename rpcFile = sourceFile;
	rpcFile.setExtension("rpb");
	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report.txt";
	ossimFilename residualFile = workfold + "\\residual.txt";
	ossimFilename rpc_reportFile = workfold + "\\rpc_report.txt";
	ossimFilename l1_reportFile = workfold + "\\l1_report.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpFile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.SavePointToFile(workfold + "\\gcp-hgt.txt", prj.m_CtrlGptSet, NULL);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	//int step = 50;
	//int nRows = ceil(12000./step);
	//int nCols = ceil(12000./step);

	//GDALDatasetH hDstDS = new GDALDatasetH;
	//GDALDriverH hDriver = GDALGetDriverByName( "GTiff" );
	//GDALDataType eDT = GDALDataType::GDT_Float32;
	//int dataTypeSize = GDALGetDataTypeSize(eDT) / 8;

	//int bandMap[] = {1, 2, 3};
	//hDstDS = GDALCreate( hDriver, "imagingRays.tif", nCols, nRows, 3, eDT, NULL );

	////cv::Mat imagimgRayMat(cv::Size(nCols, nRows), CV_64FC3);
	//float *pData = new float [nCols * nRows * 3];

	////cv::Mat imagimgRayMat(nRows, nCols, CV_64FC3);
	////double *input = (double*)(imagimgRayMat.data);
	////fstream attfs;
	////attfs.open("att.txt", ios_base::out);
	//for (int i = 0;i < nRows;i += 1)
	//{
	//	for (int j = 0;j < nCols;j += 1)
	//	{
	//		//ossimColumnVector3d image_ray;
	//		//((radiCbers04Model*)prj.m_sensorModel)->imagingRay(ossimDpt(j*step, i*step), image_ray);
	//		////imagimgRayMat.at<cv::Vec3b>(i,j)[0] = image_ray[0];
	//		////imagimgRayMat.at<cv::Vec3b>(i,j)[1] = image_ray[1];
	//		////imagimgRayMat.at<cv::Vec3b>(i,j)[2] = image_ray[2];
	//		////attfs<<i+1<<" "<<j+1<<" "<<image_ray[0]<<" "<<image_ray[1]<<" "<<image_ray[2]<<endl;
	//		////pData[3 * nCols * i + j] = image_ray[0];
	//		////pData[3 * nCols * i + j + 1] = image_ray[1];
	//		////pData[3 * nCols * i + j + 2] = image_ray[2];
	//		//pData[0 * nCols * nRows + i * nCols + j] = image_ray[0];
	//		//pData[1 * nCols * nRows + i * nCols + j] = image_ray[1];
	//		//pData[2 * nCols * nRows + i * nCols + j] = image_ray[2];

	//		ossimGpt gpt;
	//		((radiCbers04Model*)prj.m_sensorModel)->lineSampleToWorld(ossimDpt(j*step, i*step), gpt);
	//		pData[0 * nCols * nRows + i * nCols + j] = gpt.lon;
	//		pData[1 * nCols * nRows + i * nCols + j] = gpt.lat;
	//		pData[2 * nCols * nRows + i * nCols + j] = gpt.hgt;
	//	}
	//}
	//((GDALDataset*)hDstDS)->RasterIO(GF_Write, 0, 0, nCols, nRows, pData, nCols, nRows, eDT, 3,
	//	bandMap, 0, 0, 0);
	//GDALClose( hDstDS );
	//CPLFree( pData );
	////attfs.close();	

	ossimDpt imagepoint,cimagepoint;
	ossimGpt goundpoint,tGpt;

	int num = static_cast<int>(prj.m_CtrlGptSet->size());
	
	//setup initail values
	int iter=0;
	int iter_max = 200;
	double minDelta = 1e-1;
	double optimizer_delta = 1e10;
	vector<int> exteriorList;
	vector<int> innerList;
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::LINE_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::X_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::Y_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::Z_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);

	//innerList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
	////innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_OFFSET);
	//////innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROT_CENTER_OFFSET);
	//////innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROT_ANGLE);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_K1);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_K2);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_P1);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_P2);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX1_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX2_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX3_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AZ0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AZ1_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AZ2_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AZ3_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);


	vector<int> customList;
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::YAW_OFFSET);
	customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_ROLL_OFFSET);
	customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_PITCH_OFFSET);
	customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_YAW_OFFSET);
	customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_ROLL_OFFSET);
	customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_PITCH_OFFSET);
	customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_YAW_OFFSET);

	// 1B-CCD1
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 0.286288);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 7.92975);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, 6.32533);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 2.91934);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 1.00579);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -2.98177);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 2.20805);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 29.1069);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -0.000378363);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 1.54369);

	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, -11.2791);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 799.159);
	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, 5.43226);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 30.7768);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 15.4054);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -31.2659);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 27.4321);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 296.461);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -4.00537);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 15.2944);

	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 13.5843);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, -0.992082);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, -682.535);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::X_OFFSET, 0.00206263);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::Y_OFFSET, 0.07996);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::Z_OFFSET, -0.0013782);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 0.0);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, -23.0567);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -19.4154);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 4.01698);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 280.183);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -5.02263);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 24.814);

	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, 1835);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 3.61222);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 788.033);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, -662.256);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 27.9111);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, -37.822);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -32.6087);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 21.3682);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 288.869);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -6.27305);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 19.7977);

	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, 1835);
	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, 900);
	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 1.97166);
	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 742.698);
	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, -642.169);
	////prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 19.9917);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 10.2219);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -23.8694);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 22.274);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 289.748);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -0.426784);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 18.8512);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AZ0_OFFSET, 24.6023);

	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 11.7077);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, -4.65194);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, -2.37805);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, -57.374);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_OFFSET, -18.2822);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, -76.5686);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -26.2697);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, -2.08596);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 281.399);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -3.85236);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 26.6545);


	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 84.8855);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_OFFSET, 43.9206);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AX0_OFFSET, 14.3323);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AX1_OFFSET, 5.11701);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AX2_OFFSET, -5.09833);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AX3_OFFSET, -9.27952);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 123.857);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, 43.0508);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, -0.555577);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 306.981);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, 15.5828);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 9.49014);

	//// 1A-CCD1
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 0.228391);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 7.01682);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, 6.5892);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 2.21632);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 0.807303);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -3.12348);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 1.56282);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 28.0006);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, 0.427014);
	//prj.m_sensorModel->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 2.8663);
	prj.m_sensorModel->updateModel();
	prj.m_sensorModel->saveState(prj.geom);
	ossimAdjustmentInfo cadj;
	prj.m_sensorModel->getAdjustment(cadj);
	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	int np = (int)parmlist.size();
	NEWMAT::ColumnVector old_parm(np), new_parm(np);
	for(int n=0;n<np;++n)
	{
		old_parm(n+1) = parmlist[n].getParameter();
	}
	vector<Cbers04Optimize::Cbers04OptStruct> optStructList;
	Cbers04Optimize::Cbers04OptStruct cbers04Struct((radiCbers04Model*)prj.m_sensorModel, prj.m_CtrlGptSet);
	optStructList.push_back(cbers04Struct);
	int nImages = (int)optStructList.size();
	while(optimizer_delta > minDelta && iter < iter_max)
	{
		for(int image_index = 0;image_index < nImages;++image_index)
		{
			vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
			vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
			// loop on tie points
			for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
			{
				imagepoint=(*tit)->getImagePoint();
				goundpoint=(*tit)->getGroundPoint();
				tGpt = prj.m_MapPar->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
				tGpt.hgt = (*tit)->hgt;
				(*tit)->setGroundPoint(tGpt);
			}
		}

		//ossimplugins::radiRpcSolver *solver = new ossimplugins::radiRpcSolver(true, false);
		//vector < ossimDpt > imagePoints;
		//vector < ossimGpt > groundControlPoints;
		//prj.m_CtrlGptSet->getSlaveMasterPoints(imagePoints, groundControlPoints);
		//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LASSO, 1e-5);

		//ossimImageGeometry *imageGeom = solver->createRpcModel();
		//ossimKeywordlist geom;
		//imageGeom->saveState(geom);
		//ossimRpcModel* rpcModel = new ossimRpcModel();
		//rpcModel->loadState(geom, "projection.");

		//rpcModel->m_proj = prj.m_sensorModel->m_proj;
		//prj.m_sensorModel = rpcModel;

		//// do exterior orientation
		////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, exteriorList);
		////prj.m_sensorModel->updateModel();
		//Cbers04Optimize::innerOptimize::adjustment(optStructList, exteriorList);

		//// do inner orientation
		////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, innerList);
		////prj.m_sensorModel->updateModel();
		//Cbers04Optimize::innerOptimize::adjustment(optStructList, innerList);

		//Cbers04Optimize::innerOptimize::adjustment(optStructList, vector<int>(radiCbers04Model::AdjustParamIndex::LINE_OFFSET));

		//custom_opt((radiCbers04Model*)(prj.m_sensorModel), *prj.m_CtrlGptSet);
		//////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, 0);
		////prj.m_sensorModel->updateModel();
		//Cbers04Optimize::exteriorOptimize::adjustment(optStructList, customList);


		////((radiCbers04Model*)prj.m_sensorModel)->updateApproximateModel();
		prj.m_sensorModel->saveState(prj.geom);
		ossimAdjustmentInfo cadj;
		prj.m_sensorModel->getAdjustment(cadj);
		std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
		cout<<endl;
		for(int n=0;n<np;++n)
		{
			new_parm(n+1) = parmlist[n].getParameter();
			cout<<parmlist[n].getDescription()<<":"<<new_parm(n+1)<<endl;
		}
		cout<<endl;
		// then calculate the change value of the optimizers
		optimizer_delta = (new_parm - old_parm).NormInfinity();
		old_parm = new_parm;
		cout<<"iteration "<<1+iter++<<" :"<<optimizer_delta<<endl;

		for(int image_index = 0;image_index < nImages;++image_index)
		{
			vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
			vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
			// loop on tie points
			for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
			{
				ossimDpt dpt = prj.m_sensorModel->m_proj->forward(*(*tit));
				ossimGpt gpt(dpt.x,dpt.y);
				(*tit)->setGroundPoint(ossimGpt(dpt.x,dpt.y,(*tit)->hgt));
			}
		}

		//fstream fs;
		//fs.open(reportFile.c_str(), ios_base::out);
		//fs.setf(ios::fixed, ios::floatfield);
		//fs.precision(2);
		//prj.OutputReport(fs, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
		//fs.close();
		//OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, true);
		//OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, true);
		//prj.OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);
		OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);

		fstream ofs;
		ofs.open(residualFile.c_str(), ios_base::out);
		for (int i = 0;i < prj.m_CtrlGptSet->getTiePoints().size();++i)
		{
			ossimGpt ll;
			prj.m_sensorModel->lineSampleToWorld(prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint(), ll);
			ossimDpt dpt;
			optStructList[0].cbers04Model->worldToLineSample(ll, dpt);
			ossimDpt gpt = prj.m_MapPar->forward(ll);
			double delta_lat = gpt.x - prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint().lat;
			double delta_lon = gpt.y - prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint().lon;
			//double delta_lat = dpt.x - prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().x;
			//double delta_lon = dpt.y - prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().y;
			if (0 != i)
			{
				ofs<<endl;
			}
			ofs<<i+1<<"\t"
				<<prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().x<<"\t"
				<<prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().y<<"\t"
				<<delta_lat<<"\t"
				<<delta_lon;
		}
		ofs.close();
	}// end while


	ossimRefPtr<ossimImageHandler> handler   = ossimImageHandlerRegistry::instance()->open(sourceFile);
	ossimIrect theImageClipRect = handler->getBoundingRect();
	//theImageClipRect.set_ul(ossimIpt( 0, 27000));
	//theImageClipRect.set_lr(ossimIpt( 12000, 42000));
	//theImageClipRect.set_ll(ossimIpt(0, 6000));
	//theImageClipRect.set_ul(ossimIpt( 0, 0));
	//theImageClipRect.set_lr(ossimIpt( 12000, 3000));
	ossimplugins::radiRpcModel* rpcModel = Cbers04CreateRpcs(prj.m_sensorModel, theImageClipRect, workfold);
	rpcModel->writeRpcFile(rpcFile0);


	OutputReport(rpc_reportFile, (ossimSensorModel*&)rpcModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, (ossimSensorModel*&)rpcModel, prj.geom);
	rpcModel->writeRpcFile(rpcFile);
	OutputReport(l1_reportFile, (ossimSensorModel*&)rpcModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);


	//prj.m_OutBandList.clear();
	//prj.m_OutBandList.push_back(0);
	//system("pause");
	//((radiCbers04Model*)prj.m_sensorModel)->updateApproximateModel();

	//prj.starline	=	5000;
	//prj.starpixel	=	5000;
	//prj.endline		=	6000;
	//prj.endpixel	=	6000;
	////prj.m_sensorModel = rpcModel;
	//prj.Orthograph(outFile);
}

bool addCbers04Scene(ossimFilename workfold, vector<Cbers04Optimize::Cbers04OptStruct>& optStructList)
{
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\IMAGE_B3.tif";

	ossimFilename gcpFile = workfold + "\\gcp.txt";
	ossimFilename demPath = "D:\\workspace\\dem\\srtm90";

	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report.txt";
	ossimFilename residualFile = workfold + "\\residual.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpFile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.SavePointToFile(workfold + "\\gcp-hgt.txt", prj.m_CtrlGptSet, NULL);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	Cbers04Optimize::Cbers04OptStruct cbers04Struct((radiCbers04Model*)prj.m_sensorModel, prj.m_CtrlGptSet);
	optStructList.push_back(cbers04Struct);

	return true;
}


bool addCbers04Scene(ossimFilename workfold, vector<Cbers04Optimize::Cbers04OptStruct>& optStructList, ossimFilename gcpfilename)
{
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\IMAGE.tif";
	ossimString descfile = workfold + "\\desc.xml";

	ossimFilename gcpFile;
	gcpFile = workfold + "\\" + gcpfilename;
	ossimFilename demPath = "D:\\workspace\\dem\\srtm90";

	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report.txt";
	ossimFilename residualFile = workfold + "\\residual.txt";


	ossimTieGptSet* gptSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;

	ossimplugins::radiCbers04Model* sensorModel = new ossimplugins::radiCbers04Model();
	sensorModel->setupOptimizer(descfile);
	mylib::readGcpFile(gcpFile, gptSet, chkSet, NULL);

	Cbers04Optimize::Cbers04OptStruct cbers04Struct((radiCbers04Model*)sensorModel, gptSet);
	optStructList.push_back(cbers04Struct);

	return true;
}

void MultiCbers04CalibrationTest()
{
	ossimInit::instance()->initialize();
	vector<Cbers04Optimize::Cbers04OptStruct> optStructList;
	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250", optStructList);
	////addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-1_MYC_201405190109_201405190117", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140", optStructList);
	////addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01", optStructList);
	////addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-1_MYC_201403240114_201403240125", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-1_MYC_201403260030_201403260036", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-1_MYC_201404110217_201404110229", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-2\\HJ-1B_CCD-2_MYC_201405270251_201405270301\\Scene07", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-2\\HJ-1B_CCD-2_MYC_201406040119_201406040127\\Scene04", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-2\\HJ-1B_CCD-2_MYC_201406030231_201406030241\\Scene12", optStructList);

	//addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201404110217_201404110229\\Scene04", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201404110217_201404110229\\Scene06", optStructList);
	//addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201405250250_201405250300\\Scene07", optStructList);

	addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201405250250_201405250300\\Scene07", optStructList);
	addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329\\Scene08", optStructList);
	addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201407280332_201407280335\\Scene03", optStructList);
	addCbers04Scene("E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201407280153_201407280201\\Scene07", optStructList);

	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329", optStructList);
	//setup initail values
	int iter=0;
	int iter_max = 200;
	double minDelta = 1e-5;
	double optimizer_delta = 1e10;
	vector<int> exteriorList;
	vector<int> innerList;
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	exteriorList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::X_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::Y_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::Z_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::LINE_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
	//exteriorList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);

	//innerList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX1_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX2_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX3_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AZ0_OFFSET);
	////innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROT_CENTER_OFFSET);
	////innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROT_ANGLE);
	////innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET);
	////innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_K1);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_K2);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_P1);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_P2);


	vector<int> customList;
	customList.push_back(Cbers04Optimize::exteriorOptimize::ROLL_OFFSET);
	customList.push_back(Cbers04Optimize::exteriorOptimize::PITCH_OFFSET);
	customList.push_back(Cbers04Optimize::exteriorOptimize::YAW_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_YAW_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_YAW_OFFSET);

	int nImages = (int)optStructList.size();
	for(int image_index = 0;image_index < nImages;++image_index)
	{
		////optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, 1835);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 30.7768);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 15.4054);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -31.2659);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 27.4321);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 296.461);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -4.00537);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 15.2944);

		//// 1B-CCD1
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 0.286288);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 7.92975);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, 6.32533);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 2.91934);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 1.00579);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -2.98177);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 2.20805);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 29.1069);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -0.000378363);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 1.54369);

		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, 1835);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 3.61222);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 788.033);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, -662.256);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 27.9111);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, -37.822);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -32.6087);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 21.3682);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 288.869);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -6.27305);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 19.7977);

		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 67.2059);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, 33.9899);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 26.5904);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 289.017);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -7.69588);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 21.9504);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AZ0_OFFSET, 121.722);

		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 10.2219);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -23.8694);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 22.274);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 289.748);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -0.426784);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 18.8512);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AZ0_OFFSET, 24.6023);

		//// 1A-CCD1
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET, 0.228391);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET, 7.01682);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET, 6.5892);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 2.21632);
		//optStructList[image_index].cbers04Model>setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 0.807303);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -3.12348);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 1.56282);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 28.0006);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, 0.427014);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 2.8663);

		optStructList[image_index].cbers04Model->updateModel();
	}
	ossimAdjustmentInfo cadj;
	optStructList[0].cbers04Model->getAdjustment(cadj);
	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	int np = (int)parmlist.size();
	NEWMAT::ColumnVector old_parm(np), new_parm(np);
	for(int n=0;n<np;++n)
	{
		old_parm(n+1) = parmlist[n].getParameter();
	}

	ossimDpt imagepoint,cimagepoint;
	ossimGpt goundpoint,tGpt;
	while(optimizer_delta > minDelta && iter < iter_max)
	{
		for(int image_index = 0;image_index < nImages;++image_index)
		{
			vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
			vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
			// loop on tie points
			for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
			{
				imagepoint=(*tit)->getImagePoint();
				goundpoint=(*tit)->getGroundPoint();
				tGpt = optStructList[image_index].cbers04Model->m_proj->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
				tGpt.hgt = (*tit)->hgt;
				(*tit)->setGroundPoint(tGpt);
			}
		}

		//// do exterior orientation
		////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, exteriorList);
		////prj.m_sensorModel->updateModel();
		Cbers04Optimize::innerOptimize::adjustment(optStructList, exteriorList);

		//// do inner orientation
		////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, innerList);
		////prj.m_sensorModel->updateModel();
		Cbers04Optimize::innerOptimize::adjustment(optStructList, innerList);

		//////custom_opt((radiCbers04Model*)(prj.m_sensorModel), *prj.m_CtrlGptSet);
		////////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, 0);
		//////prj.m_sensorModel->updateModel();
		//Cbers04Optimize::exteriorOptimize::adjustment(optStructList, customList);





		ossimAdjustmentInfo cadj;
		optStructList[0].cbers04Model->getAdjustment(cadj);
		std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
		cout<<endl;
		for(int n=0;n<np;++n)
		{
			new_parm(n+1) = parmlist[n].getParameter();
			cout<<parmlist[n].getDescription()<<":"<<new_parm(n+1)<<endl;
		}
		cout<<endl;
		// then calculate the change value of the optimizers
		optimizer_delta = (new_parm - old_parm).NormInfinity();
		old_parm = new_parm;
		cout<<"iteration "<<1+iter++<<" :"<<optimizer_delta<<endl;

		for(int image_index = 0;image_index < nImages;++image_index)
		{
			vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
			vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
			// loop on tie points
			for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
			{
				ossimDpt dpt = optStructList[image_index].cbers04Model->m_proj->forward(*(*tit));
				ossimGpt gpt(dpt.x,dpt.y);
				(*tit)->setGroundPoint(ossimGpt(dpt.x,dpt.y,(*tit)->hgt));
			}
		}

		//fstream fs;
		//fs.open(reportFile.c_str(), ios_base::out);
		//fs.setf(ios::fixed, ios::floatfield);
		//fs.precision(2);
		//prj.OutputReport(fs, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
		//fs.close();

		//fstream ofs;
		//ofs.open(residualFile.c_str(), ios_base::out);
		//for (int i = 0;i < prj.m_CtrlGptSet->getTiePoints().size();++i)
		//{
		//	ossimGpt ll;
		//	prj.m_sensorModel->lineSampleToWorld(prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint(), ll);
		//	ossimDpt gpt = prj.m_MapPar->forward(ll);
		//	double delta_lat = gpt.x - prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint().lat;
		//	double delta_lon = gpt.y - prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint().lon;
		//	if (0 != i)
		//	{
		//		ofs<<endl;
		//	}
		//	ofs<<i+1<<"\t"
		//		<<prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().x<<"\t"
		//		<<prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().y<<"\t"
		//		<<delta_lat<<"\t"
		//		<<delta_lon;
		//}
		//ofs.close();
	}// end while
	system("pause");
}


void calibration()
{
	ossimInit::instance()->initialize();
	vector<Cbers04Optimize::Cbers04OptStruct> optStructList;
	ossimFilename foldName = "E:\\HJ1\\CBERS04\\CBERS04_MUX_MYC_201501230418_201501230428\\Scene31";
	//foldName = "E:\\HJ1\\CBERS04\\CBERS04_P5M_MYC_201501250132_201501250137\\Scene14";
	//foldName = "E:\\HJ1\\CBERS04\\CBERS04_MUX_MYC_201501250132_201501250137\\Scene07";
	//foldName = "E:\\HJ1\\CBERS04\\CBERS04_MUX_MYC_201501270340_201501270352\\Scene36";
	//foldName = "E:\\HJ1\\CBERS04\\CBERS04_MUX_MYC_201501280446_201501280454\\Scene18";
	//foldName = "E:\\HJ1\\CBERS04\\CBERS04_WFI_MYC_201501270340_201501270352\\Scene04";
	//foldName = "E:\\HJ1\\CBERS04\\test\\Scene31";
	ossimFilename reportFile = foldName + "\\cbers04Report.txt";
	ossimFilename attFile = foldName + "\\attList.txt";
	addCbers04Scene(foldName, optStructList, "gcp.txt");
	addCbers04Scene("E:\\HJ1\\CBERS04\\CBERS04_MUX_MYC_201501250132_201501250137\\Scene07", optStructList, "gcp.txt");
	addCbers04Scene("E:\\HJ1\\CBERS04\\CBERS04_MUX_MYC_201501270340_201501270352\\Scene36", optStructList, "gcp.txt");
	//addCbers04Scene("E:\\HJ1\\CBERS04\\CBERS04_P5M_MYC_201501250132_201501250137\\Scene14", optStructList, "gcp.txt");
	int nImages = (int)optStructList.size();

	for (int image_index = 0; image_index < nImages; ++image_index)
	{
		mylib::get_elevation(optStructList[image_index].gcpSet);
	}

	int iter = 0;
	int iter_max = 200;
	double minDelta = 1e-5;
	double optimizer_delta = 1e10;
	vector<int> innerList;
	innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::X_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::Y_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::Z_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET);

	vector<int> exteriorList;
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::ROLL_OFFSET);
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::PITCH_OFFSET);
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::YAW_OFFSET);


	vector<int> customList;
	customList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	customList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	customList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);
	//customList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET);
	//customList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_YAW_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_YAW_OFFSET);

	for (int image_index = 0; image_index < nImages; ++image_index)
	{
		//	//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, 1835);
		//	optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 30.7768);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 15.4054);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -31.2659);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 27.4321);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 296.461);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -4.00537);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 15.2944);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET, -4);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET, 39);

		optStructList[image_index].cbers04Model->updateModel();
	}
	ossimAdjustmentInfo cadj;
	optStructList[0].cbers04Model->getAdjustment(cadj);
	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	int np = (int)parmlist.size();
	NEWMAT::ColumnVector old_parm(np), new_parm(np);
	for (int n = 0; n<np; ++n)
	{
		old_parm(n + 1) = parmlist[n].getParameter();
	}

	ossimDpt imagepoint, cimagepoint;
	ossimGpt goundpoint, tGpt;
	
	while (optimizer_delta > minDelta && iter < iter_max)
	{

		//Cbers04Optimize::innerOptimize::adjustment(optStructList, customList);
		////// do inner orientation
		//////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, innerList);
		//////prj.m_sensorModel->updateModel();
		Cbers04Optimize::innerOptimize::adjustment(optStructList, innerList);

		//////Cbers04Optimize::innerOptimize::adjustment(optStructList, customList);
		//////// do exterior orientation
		////////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, exteriorList);
		////////prj.m_sensorModel->updateModel();
		//Cbers04Optimize::exteriorOptimize::adjustment(optStructList, exteriorList);


		ossimAdjustmentInfo cadj;
		optStructList[0].cbers04Model->getAdjustment(cadj);
		std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
		cout << endl;
		for (int n = 0; n<np; ++n)
		{
			new_parm(n + 1) = parmlist[n].getParameter();
			cout << parmlist[n].getDescription() << ":" << new_parm(n + 1) << endl;
		}
		cout << endl;
		// then calculate the change value of the optimizers
		optimizer_delta = (new_parm - old_parm).NormInfinity();
		old_parm = new_parm;
		cout << "iteration " << 1 + iter++ << " :" << optimizer_delta << endl;


		FILE* pf = fopen(reportFile.c_str(), "w+");
		int icount = 0;
		for (int image_index = 0; image_index < nImages; ++image_index)
		{
			vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
			vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
			// loop on tie points
			for (tit = theTPV.begin(); tit != theTPV.end(); ++tit)
			{
				////compute residue
				//ossimGpt gd = optStructList[image_index].cbers04Model->inverse((*tit)->tie);
				//ossimGpt res;
				//res.lon = ((*tit)->lon - gd.lon) * 100000.0; //approx meters //TBC TBD
				//res.lat = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
				//res.hgt = (*tit)->hgt - gd.hgt; //meters

				//fprintf(pf, "%6d%10.2lf%10.2lf%15.6lf%15.6lf%15.6lf\n", icount + 1, (*tit)->tie.x,
				//	(*tit)->tie.y, res.lon, res.lat, res.hgt);

				ossimDpt dpt = optStructList[image_index].cbers04Model->forward((*tit)->getGroundPoint());
				ossimDpt res;
				res.x = dpt.x - (*tit)->getImagePoint().x;
				res.y = dpt.y - (*tit)->getImagePoint().y;
				fprintf(pf, "%6d%10.2lf%10.2lf%15.6lf%15.6lf%15.6lf\n", icount + 1, (*tit)->tie.x,
					(*tit)->tie.y, res.x, res.y);

				icount++;
			}
		}

		fclose(pf);

		//pf = fopen(attFile.c_str(), "w+");
		//fprintf(pf, "%10s%15s%15s%15s%15s%15s%15s\n", "id", "att0.x", "att.x", "att0.y", "att.y", "att0.z", "att.z");
		//double step = 12000 / (double)nl;
		//for (size_t i = 0; i < nl; i++)
		//{
		//	double line = i * step;
		//	double t_line;
		//	optStructList[i].cbers04Model->theSupportData->getLineTime(line, t_line);
		//	ossimDpt3d att(0.0, 0.0, 0.0);
		//	ossimDpt3d att0(0.0, 0.0, 0.0);
		//	optStructList[i].cbers04Model->theSupportData->getAttitude(t_line, att0);
		//	att.x = att0.x + optStructList[i].cbers04Model->thePitchOffset*1e-3;
		//	att.y = att0.y + optStructList[i].cbers04Model->theRollOffset*1e-3;
		//	att.z = att0.z + optStructList[i].cbers04Model->theYawOffset*1e-3;
		//	fprintf(pf, "%10.1lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf\n", line, att0.x, att.x, att0.y, att.y, att0.z, att.z);
		//}
		//fclose(pf);
	}// end while

}

void calcAttitude()
{
	ossimInit::instance()->initialize();
	vector<Cbers04Optimize::Cbers04OptStruct> optStructList;

	ossimFilename reportFile = "E:\\HJ1\\Scene08\\attReport.txt";
	ossimFilename attFile = "E:\\HJ1\\Scene08\\attList.txt";

	int nl = 30;
	for (size_t i = 0; i < nl; i++)
	{
		char buf[256];
		sprintf_s(buf, "gcp\\%d.txt", i);
		addCbers04Scene("E:\\HJ1\\Scene08", optStructList, buf);
	}

	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329", optStructList);
	//setup initail values
	int iter = 0;
	int iter_max = 200;
	double minDelta = 1e-5;
	double optimizer_delta = 1e10;
	vector<int> innerList;
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::X_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::Y_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::Z_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET);

	vector<int> exteriorList;
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::ROLL_OFFSET);
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::PITCH_OFFSET);
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::YAW_OFFSET);


	//vector<int> customList;
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::YAW_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_YAW_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_YAW_OFFSET);

	int nImages = (int)optStructList.size();
	for (int image_index = 0; image_index < nImages; ++image_index)
	{
	//	//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, 1835);
	//	optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 30.7768);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 15.4054);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -31.2659);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 27.4321);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 296.461);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -4.00537);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 15.2944);
		optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET, -4);
		optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET, 39);

		optStructList[image_index].cbers04Model->updateModel();
	}
	ossimAdjustmentInfo cadj;
	optStructList[0].cbers04Model->getAdjustment(cadj);
	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	int np = (int)parmlist.size();
	NEWMAT::ColumnVector old_parm(np), new_parm(np);
	for (int n = 0; n<np; ++n)
	{
		old_parm(n + 1) = parmlist[n].getParameter();
	}

	ossimDpt imagepoint, cimagepoint;
	ossimGpt goundpoint, tGpt;

	//Cbers04Optimize::blockadjustment::sparselm_adjustment(optStructList, innerList, exteriorList);

	////ossimAdjustmentInfo cadj;
	////optStructList[0].cbers04Model->getAdjustment(cadj);
	////std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	////cout << endl;
	////for (int n = 0; n<np; ++n)
	////{
	////	new_parm(n + 1) = parmlist[n].getParameter();
	////	cout << parmlist[n].getDescription() << ":" << new_parm(n + 1) << endl;
	////}
	////cout << endl;
	////// then calculate the change value of the optimizers
	////optimizer_delta = (new_parm - old_parm).NormInfinity();
	////old_parm = new_parm;
	////cout << "iteration " << 1 + iter++ << " :" << optimizer_delta << endl;


	//fstream ofs;
	//ofs.open(reportFile.c_str(), ios_base::out);
	//int icount = 0;
	//for (int image_index = 0; image_index < nImages; ++image_index)
	//{
	//	vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
	//	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
	//	// loop on tie points
	//	for (tit = theTPV.begin(); tit != theTPV.end(); ++tit)
	//	{
	//		//ossimDpt dpt = optStructList[image_index].cbers04Model->forward(*(*tit));
	//		//imagepoint = (*tit)->getImagePoint();
	//		//ossimDpt res = dpt - imagepoint;

	//		//ofs << icount + 1 << "\t"
	//		//	<< (*tit)->tie.x << "\t"
	//		//	<< (*tit)->tie.y << "\t"
	//		//	<< res.x << "\t"
	//		//	<< res.y
	//		//	<< endl;
	//		//icount++;


	//		//compute residue
	//		ossimGpt gd = optStructList[image_index].cbers04Model->inverse((*tit)->tie);
	//		ossimGpt res;
	//		res.lon = ((*tit)->lon - gd.lon) * 100000.0; //approx meters //TBC TBD
	//		res.lat = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
	//		res.hgt = (*tit)->hgt - gd.hgt; //meters

	//		ofs << icount + 1 << "\t"
	//			<< (*tit)->tie.x << "\t"
	//			<< (*tit)->tie.y << "\t"
	//			<< res.lon << "\t"
	//			<< res.lat << "\t"
	//			<< res.hgt << "\t"
	//			<< endl;
	//		icount++;
	//	}
	//}

	//ofs.close();

	//FILE* pf = fopen(attFile.c_str(), "w+");
	//double step = 12000 / (double)nl;
	//for (size_t i = 0; i < nl; i++)
	//{
	//	double line = i * step;
	//	double t_line;
	//	optStructList[i].cbers04Model->theSupportData->getLineTime(line, t_line);
	//	ossimDpt3d att(0.0, 0.0, 0.0);
	//	optStructList[i].cbers04Model->theSupportData->getAttitude(t_line, att);
	//	att.x += optStructList[i].cbers04Model->thePitchOffset;
	//	att.y += optStructList[i].cbers04Model->theRollOffset;
	//	att.z += optStructList[i].cbers04Model->theYawOffset;
	//	fprintf(pf, "%7.1lf%15.6lf%15.6lf%15.6lf\n", line, att.x, att.y, att.z);
	//}
	//fclose(pf);

	while (optimizer_delta > minDelta && iter < iter_max)
	{
		//for (int image_index = 0; image_index < nImages; ++image_index)
		//{
		//	vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
		//	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
		//	// loop on tie points
		//	for (tit = theTPV.begin(); tit != theTPV.end(); ++tit)
		//	{
		//		imagepoint = (*tit)->getImagePoint();
		//		goundpoint = (*tit)->getGroundPoint();
		//		tGpt = optStructList[image_index].cbers04Model->m_proj->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
		//		tGpt.hgt = (*tit)->hgt;
		//		(*tit)->setGroundPoint(tGpt);
		//	}
		//}

		//// do exterior orientation
		////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, exteriorList);
		////prj.m_sensorModel->updateModel();
		Cbers04Optimize::exteriorOptimize::adjustment(optStructList, exteriorList);

		////// do inner orientation
		//////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, innerList);
		//////prj.m_sensorModel->updateModel();
		Cbers04Optimize::innerOptimize::adjustment(optStructList, innerList);

		ossimAdjustmentInfo cadj;
		optStructList[0].cbers04Model->getAdjustment(cadj);
		std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
		cout << endl;
		for (int n = 0; n<np; ++n)
		{
			new_parm(n + 1) = parmlist[n].getParameter();
			cout << parmlist[n].getDescription() << ":" << new_parm(n + 1) << endl;
		}
		cout << endl;
		// then calculate the change value of the optimizers
		optimizer_delta = (new_parm - old_parm).NormInfinity();
		old_parm = new_parm;
		cout << "iteration " << 1 + iter++ << " :" << optimizer_delta << endl;


		FILE* pf = fopen(reportFile.c_str(), "w+");
		int icount = 0;
		for (int image_index = 0; image_index < nImages; ++image_index)
		{
			vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
			vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
			// loop on tie points
			for (tit = theTPV.begin(); tit != theTPV.end(); ++tit)
			{
				//ossimDpt dpt = optStructList[image_index].cbers04Model->forward(*(*tit));
				//imagepoint = (*tit)->getImagePoint();
				//ossimDpt res = dpt - imagepoint;

				//ofs << icount + 1 << "\t"
				//	<< (*tit)->tie.x << "\t"
				//	<< (*tit)->tie.y << "\t"
				//	<< res.x << "\t"
				//	<< res.y
				//	<< endl;
				//icount++;


				//compute residue
				ossimGpt gd = optStructList[image_index].cbers04Model->inverse((*tit)->tie);
				ossimGpt res;
				res.lon = ((*tit)->lon - gd.lon) * 100000.0; //approx meters //TBC TBD
				res.lat = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
				res.hgt = (*tit)->hgt - gd.hgt; //meters

				fprintf(pf, "%6d%10.2lf%10.2lf%15.6lf%15.6lf%15.6lf\n", icount + 1, (*tit)->tie.x,
					(*tit)->tie.y, res.lon, res.lat, res.hgt);
				icount++;
			}
		}

		fclose(pf);

		pf = fopen(attFile.c_str(), "w+");
		fprintf(pf, "%10s%15s%15s%15s%15s%15s%15s\n", "id", "att0.x", "att.x", "att0.y", "att.y", "att0.z", "att.z");
		double step = 12000 / (double)nl;
		for (size_t i = 0; i < nl; i++)
		{
			double line = i * step;
			double t_line;
			optStructList[i].cbers04Model->theSupportData->getLineTime(line, t_line);
			ossimDpt3d att(0.0, 0.0, 0.0);
			ossimDpt3d att0(0.0, 0.0, 0.0);
			optStructList[i].cbers04Model->theSupportData->getAttitude(t_line, att0);
			att.x = att0.x + optStructList[i].cbers04Model->thePitchOffset*1e-3;
			att.y = att0.y + optStructList[i].cbers04Model->theRollOffset*1e-3;
			att.z = att0.z + optStructList[i].cbers04Model->theYawOffset*1e-3;
			fprintf(pf, "%10.1lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf%15.6lf\n", line, att0.x, att.x, att0.y, att.y, att0.z, att.z);
		}
		fclose(pf);
	}// end while

}

void calcAttitudeSparse()
{
	ossimInit::instance()->initialize();
	vector<Cbers04Optimize::Cbers04OptStruct> optStructList;

	ossimFilename reportFile = "E:\\HJ1\\Scene08\\attReport.txt";
	ossimFilename attFile = "E:\\HJ1\\Scene08\\attList.txt";

	int nl = 10;
	for (size_t i = 0; i < nl; i++)
	{
		char buf[256];
		sprintf_s(buf, "gcp\\%d.txt", i);
		addCbers04Scene("E:\\HJ1\\Scene08", optStructList, buf);
	}

	//addCbers04Scene("E:\\HJ1\\HJ-1B_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329", optStructList);
	//setup initail values
	int iter = 0;
	int iter_max = 200;
	double minDelta = 1e-5;
	double optimizer_delta = 1e10;
	vector<int> innerList;
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_YAW_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AX0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY0_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY1_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY2_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY3_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY4_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::AY5_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::X_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::Y_OFFSET);
	//innerList.push_back(radiCbers04Model::AdjustParamIndex::Z_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET);
	innerList.push_back(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET);

	vector<int> exteriorList;
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::ROLL_OFFSET);
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::PITCH_OFFSET);
	exteriorList.push_back(Cbers04Optimize::exteriorOptimize::YAW_OFFSET);


	//vector<int> customList;
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::YAW_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_YAW_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_ROLL_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_PITCH_OFFSET);
	//customList.push_back(Cbers04Optimize::exteriorOptimize::ATT_VEL_YAW_OFFSET);

	int nImages = (int)optStructList.size();
	for (int image_index = 0; image_index < nImages; ++image_index)
	{
		//	//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::LINE_OFFSET, 1835);
		//	optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::FOCAL_OFFSET, 30.7768);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY0_OFFSET, 15.4054);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY1_OFFSET, -31.2659);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY2_OFFSET, 27.4321);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY3_OFFSET, 296.461);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY4_OFFSET, -4.00537);
		//optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::AY5_OFFSET, 15.2944);
		optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_CENTER_X_OFFSET, -4);
		optStructList[image_index].cbers04Model->setAdjustableParameter(radiCbers04Model::AdjustParamIndex::CCD_CENTER_Y_OFFSET, 39);

		optStructList[image_index].cbers04Model->updateModel();
	}
	ossimAdjustmentInfo cadj;
	optStructList[0].cbers04Model->getAdjustment(cadj);
	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	int np = (int)parmlist.size();
	NEWMAT::ColumnVector old_parm(np), new_parm(np);
	for (int n = 0; n<np; ++n)
	{
		old_parm(n + 1) = parmlist[n].getParameter();
	}

	ossimDpt imagepoint, cimagepoint;
	ossimGpt goundpoint, tGpt;

	Cbers04Optimize::blockadjustment::Cbers04Calibration *pCalibration = new Cbers04Optimize::blockadjustment::Cbers04Calibration();
	pCalibration->adjustment(optStructList, innerList, exteriorList);

	delete pCalibration;

	//ossimAdjustmentInfo cadj;
	//optStructList[0].cbers04Model->getAdjustment(cadj);
	//std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	//cout << endl;
	//for (int n = 0; n<np; ++n)
	//{
	//	new_parm(n + 1) = parmlist[n].getParameter();
	//	cout << parmlist[n].getDescription() << ":" << new_parm(n + 1) << endl;
	//}
	//cout << endl;
	//// then calculate the change value of the optimizers
	//optimizer_delta = (new_parm - old_parm).NormInfinity();
	//old_parm = new_parm;
	//cout << "iteration " << 1 + iter++ << " :" << optimizer_delta << endl;


	FILE* pf = fopen(reportFile.c_str(), "w+");
	int icount = 0;
	for (int image_index = 0; image_index < nImages; ++image_index)
	{
		vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->refTiePoints();
		vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
		// loop on tie points
		for (tit = theTPV.begin(); tit != theTPV.end(); ++tit)
		{
			//ossimDpt dpt = optStructList[image_index].cbers04Model->forward(*(*tit));
			//imagepoint = (*tit)->getImagePoint();
			//ossimDpt res = dpt - imagepoint;

			//ofs << icount + 1 << "\t"
			//	<< (*tit)->tie.x << "\t"
			//	<< (*tit)->tie.y << "\t"
			//	<< res.x << "\t"
			//	<< res.y
			//	<< endl;
			//icount++;


			//compute residue
			ossimGpt gd = optStructList[image_index].cbers04Model->inverse((*tit)->tie);
			ossimGpt res;
			res.lon = ((*tit)->lon - gd.lon) * 100000.0; //approx meters //TBC TBD
			res.lat = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
			res.hgt = (*tit)->hgt - gd.hgt; //meters

			fprintf(pf, "%6d%10.2lf%10.2lf%15.6lf%15.6lf%15.6lf\n", icount + 1, (*tit)->tie.x,
				(*tit)->tie.y, res.lon, res.lat, res.hgt);
			icount++;
		}
	}

	fclose(pf);

	pf = fopen(attFile.c_str(), "w+");
	double step = 12000 / (double)nl;
	for (size_t i = 0; i < nl; i++)
	{
		double line = i * step;
		double t_line;
		optStructList[i].cbers04Model->theSupportData->getLineTime(line, t_line);
		ossimDpt3d att(0.0, 0.0, 0.0);
		optStructList[i].cbers04Model->theSupportData->getAttitude(t_line, att);
		att.x += optStructList[i].cbers04Model->thePitchOffset*1e-4;
		att.y += optStructList[i].cbers04Model->theRollOffset*1e-4;
		att.z += optStructList[i].cbers04Model->theYawOffset*1e-4;
		fprintf(pf, "%10.1lf%15.6lf%15.6lf%15.6lf\n", line, att.x, att.y, att.z);
	}
	fclose(pf);
}

void Cbers04L2Rpc()
{
	fstream fs;
	//fs.open("D:\\workspace\\HJ\\test\\HJ1B-CCD2-449-56-20100917-L20000393774-SatAngle.txt", ios_base::in);
	//fs.open("D:\\workspace\\HJ\\test\\HJ1B-CCD1-456-84-20140416-L20001144407-SatAngle.txt", ios_base::in);
	fs.open("D:\\workspace\\HJ\\test\\HJ1B-CCD1-12-44-20140410-L20001141942-SatAngle.txt", ios_base::in);
	char buf[1024];
	int step = 600; // 20 points
	int threshold = 17;
	ossimTieGptSet *gcpSet = new ossimTieGptSet;
	int iCount = 0;
	while (fs.getline(buf, 1024))
	{
		string strLine(buf);
		vector<string> strList;
		mylib::SplitString(strLine, " ", strList, false);
		if (strList.size() != 8)
		{
			continue;
		}
		if (atoi(strList[0].c_str())%step > threshold && atoi(strList[0].c_str())%step < (600 -threshold))
		{
			continue;
		}
		if (atoi(strList[1].c_str())%step > threshold && atoi(strList[1].c_str())%step < (600 -threshold))
		{
			continue;
		}

		ossimDpt dpt(atoi(strList[1].c_str()), atoi(strList[0].c_str()));
		ossimGpt gpt(atof(strList[5].c_str()), atof(strList[4].c_str()));
		gpt.hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(gpt);
		gpt = ossimGpt(gpt.lon, gpt.lat, gpt.hgt);

		char strId[10];
		_itoa(++iCount, strId, 10);
		gcpSet->addTiePoint(new ossimTieGpt(gpt, dpt, 0.0, ossimString(strId)));
	}

	mylib::saveGcpFile("D:\\workspace\\HJ\\test\\gcp.txt", gcpSet);
	fs.close();
}

void rpcFromPointsTest()
{
	ossimFilename workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01";
	workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02";
	ossimFilename gcpFile = workfold + "\\gcp.txt";
	ossimFilename rpcReportFile = workfold + "\\rpc_report.txt";
	ossimRpcModel* rpcModel = mylib::createRpcModelFromPoints(gcpFile);

	mylib::projection2ll(workfold + "\\gcp-hgt.txt", workfold + "\\ll.txt");

	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	mylib::readGcpFile(gcpFile, gcpSet, chkSet, &prjKwl);
	rpcModel->m_proj = PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(prjKwl));
	
	fstream rpcStructFile;
	rpcStructFile.open((workfold + "\\rpc.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStructFile);
	rpcStructFile.close();

	mylib::get_elevation(gcpSet, prjKwl, 0.0);
	mylib::get_elevation(chkSet, prjKwl, 0.0);
	//mylib::OutputReport(rpcReportFile, rpcModel, gcpSet, chkSet, false, false);
	mylib::OutputReport(rpcReportFile, rpcModel, gcpSet, chkSet, true, false);
}

void ll2projection()
{
	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	ossimFilename folder = "E:\\HJ1\\HJ-1A_CCD-1_MYC_201404110217_201404110229";
	mylib::readGcpFile(folder+"\\ll.txt", gcpSet, chkSet, &prjKwl);
	prjKwl.clear();
	prjKwl.addFile(ossimFilename(folder+"\\gcp.geom"));
	ossimRefPtr<ossimMapProjection> theMapProjection = PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(prjKwl));

	for (vector<ossimRefPtr<ossimTieGpt> >::iterator iter = gcpSet->refTiePoints().begin() ; iter != gcpSet->refTiePoints().end() && iter->valid();)
	{
		ossimGpt oldGpt = (*iter)->getGroundPoint();
		ossimDpt dpt = theMapProjection->forward(oldGpt);
		(*iter)->setGroundPoint(ossimGpt(dpt.x, dpt.y, oldGpt.hgt));
		iter++;
	}

	for (vector<ossimRefPtr<ossimTieGpt> >::iterator iter = chkSet->refTiePoints().begin() ; iter != chkSet->refTiePoints().end() && iter->valid();)
	{
		ossimGpt oldGpt = (*iter)->getGroundPoint();
		ossimDpt dpt = theMapProjection->forward(oldGpt);
		(*iter)->setGroundPoint(ossimGpt(dpt.x, dpt.y, oldGpt.hgt));
		iter++;
	}

	mylib::saveGcpFile(folder+"\\gcp.txt", gcpSet, chkSet, &prjKwl);
}

void ll2projection(ossimFilename folder)
{
	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	mylib::readGcpFile(folder+"\\ll.txt", gcpSet, chkSet, &prjKwl);
	prjKwl.clear();
	prjKwl.addFile(ossimFilename(folder+"\\gcp.geom"));
	ossimRefPtr<ossimMapProjection> theMapProjection = PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(prjKwl));

	for (vector<ossimRefPtr<ossimTieGpt> >::iterator iter = gcpSet->refTiePoints().begin() ; iter != gcpSet->refTiePoints().end() && iter->valid();)
	{
		ossimGpt oldGpt = (*iter)->getGroundPoint();
		ossimDpt dpt = theMapProjection->forward(oldGpt);
		(*iter)->setGroundPoint(ossimGpt(dpt.x, dpt.y, oldGpt.hgt));
		iter++;
	}

	for (vector<ossimRefPtr<ossimTieGpt> >::iterator iter = chkSet->refTiePoints().begin() ; iter != chkSet->refTiePoints().end() && iter->valid();)
	{
		ossimGpt oldGpt = (*iter)->getGroundPoint();
		ossimDpt dpt = theMapProjection->forward(oldGpt);
		(*iter)->setGroundPoint(ossimGpt(dpt.x, dpt.y, oldGpt.hgt));
		iter++;
	}

	mylib::saveGcpFile(folder+"\\gcp.txt", gcpSet, chkSet, &prjKwl);
}

int main( int argc, char* argv[] )
{
	string dropbox = getenv("DROPBOX");
	//string preferences_file = dropbox + "\\Programs\\ossimTest_2.0\\preference_office.txt";
	//string preferences_file = dropbox + "\\Programs\\batch_orth\\batch_orth\\preference.txt";
	string preferences_file = string(getenv("OSSIM2_DIR")) + "\\preference.txt";

	GDALAllRegister();
	//OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持


	ossimPreferences::instance()->loadPreferences(ossimFilename(preferences_file));
	ossimInit::instance()->initialize();
	//ossimElevManager::instance()->loadElevationPath("D:\\workspace\\dem\\srtm90");

	clock_t  clockBegin, clockEnd;
	clockBegin = clock();
	//ll2projection("E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201403240251_201403240259");
	//ll2projection("E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201405250250_201405250300");
	//ll2projection("E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329");
	//Cbers04L2Rpc();
	//LineOffsetTest();
	//rpcFromPointsTest();
	//Cbers04CalibrationTest();

	ossimFilename folder = "E:\\HJ1\\HJ-1B_CCD-1\\HJ-1B_CCD-1_MYC_201405190109_201405190117";
	////Cbers04CalibrationTest(folder+"");
	//Cbers04CalibrationTest(folder+"\\Scene01");
	//Cbers04CalibrationTest(folder+"\\Scene02");
	//Cbers04CalibrationTest(folder+"\\Scene03");
	//Cbers04CalibrationTest(folder+"\\Scene04");
	//Cbers04CalibrationTest(folder+"\\Scene05");
	//Cbers04CalibrationTest(folder+"\\Scene06");
	//Cbers04CalibrationTest(folder+"\\Scene07");
	//Cbers04CalibrationTest(folder+"\\Scene08");
	//Cbers04CalibrationTest(folder+"\\Scene09");

	//folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329";
	////Cbers04CalibrationTest(folder+"");
	//Cbers04CalibrationTest(folder+"\\Scene01");
	//Cbers04CalibrationTest(folder+"\\Scene02");
	//Cbers04CalibrationTest(folder+"\\Scene03");
	//Cbers04CalibrationTest(folder+"\\Scene04");
	//Cbers04CalibrationTest(folder+"\\Scene05");
	//Cbers04CalibrationTest(folder+"\\Scene06");
	//Cbers04CalibrationTest(folder+"\\Scene07");
	//Cbers04CalibrationTest(folder+"\\Scene08");
	//Cbers04CalibrationTest(folder+"\\Scene09");
	//Cbers04CreateRpcs();
	//MultiCbers04CalibrationTest();
	//HJ1_create3DGridPoints();

	//projection2ll("E:\\HJ1\\HJ-1B_CCD-1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\gcp.txt", "E:\\HJ1\\HJ-1B_CCD-1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\gcp_ll.txt");
	//folder = "E:\\HJ1\\HJ-1B_CCD-2\\HJ-1B_CCD-2_MYC_201406030231_201406030241\\Scene12";
	//folder = "E:\\HJ1\\HJ-1B_CCD-2\\HJ-1B_CCD-2_MYC_201405270251_201405270301\\Scene07";
	//folder = "E:\\HJ1\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201404110217_201404110229\\Scene04";
	//folder = "E:\\HJ1\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201404110217_201404110229\\Scene06";
	//folder = "E:\\HJ1\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201405250250_201405250300\\Scene07";

	folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201405250250_201405250300\\Scene07";
	//folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329\\Scene08";
	//folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201406040208_201406040220\\Scene13";
	folder = "E:\\HJ1\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201405250250_201405250300\\Scene07";
	folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201407280332_201407280335\\Scene03";
	folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201407280153_201407280201\\Scene07";
	folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329\\Scene08";
	folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201405250250_201405250300\\Scene07";
	folder = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201405250250_201405250300\\Scene07";
	//Cbers04CalibrationTest(folder);
	//calcAttitude();
	//calcAttitudeSparse();
	calibration();

	clockEnd = clock();
	printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
}