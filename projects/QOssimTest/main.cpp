#include <iostream>
#include <iterator>
#include <fstream>
using namespace std;

#include <stdlib.h>
#include <stdio.h>

#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#include <cv.h>
#include <cxcore.h>
#include <highgui.h>

#include "omp.h"
#include <ossim/parallel/ossimMpi.h>
#include <ossim/base/ossimPreferences.h>

//#include "opencv-ossim/ossimOpenCVCannyFilter.h"


//#include "feature_detection.h"
#include "test\RpcTest.hpp"
#include "test\featureTest.hpp"
#include <gcpUtil.h>

#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")
#pragma comment(lib, "opencv_ts300.lib")
#pragma comment(lib, "opencv_world300.lib")
using namespace cv;


void get_elevation(ossimFilename infile, ossimFilename outfile, ossimFilename elevationPath="D:\\workspace\\ossimdem")
{
	FILE *pfIn = fopen(infile.c_str(), "r+" );
	FILE *pfOut = fopen(outfile.c_str(), "w+");

	ossimElevManager* theElevManager = ossimElevManager::instance();
	theElevManager->loadElevationPath(elevationPath);
	int id;
	double L2Line, L2Sample;
	double Lon, Lat;
	double Hgt;
	int index = 1;
	while (EOF != fscanf(pfIn, "%d%lf%lf%lf%lf%lf", &id, &L2Line, &L2Sample, &Lon, &Lat, &Hgt))
	{
		//ossimGpt gpt(Lon, Lat, Hgt);
		ossimGpt gpt(Lat, Lon, Hgt);
		Hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(gpt);
		fprintf(pfOut, "%6d%20.9lf%20.9lf%20.9lf%20.9lf%20.9lf\n", index++, L2Line, L2Sample, Lon, Lat, Hgt);
	}
	fclose(pfIn);
	fclose(pfOut);
}


void landsat_optimization(ossimFilename gcpFile, ossimFilename headerFile, ossimFilename elevationpath)
{

	ossimFilename workfold = headerFile.path();
	ossimFilename sourcefile = workfold + "\\header.dat";
	ossimFilename gcpfile = gcpFile;
	ossimFilename reportfile = workfold + "\\report.txt";

	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));

	prj.theMgr = ossimElevManager::instance();
	//prj.theMgr->loadElevationPath(ossimFilename(elevationpath));//
	prj.m_DemPath=ossimFilename(elevationpath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.InitiateSensorModel(sourcefile);
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.m_sensorModel->loadState(prj.geom);
	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
}

//void batch_gcp_optimization()
//{
//	ossimFilename inGcpDir = "F:\\testdata\\landsat\\gcp_whole\\L5\\2008\\0";
//	ossimFilename inSrcDir = "F:\\testdata\\landsat\\gcp_whole\\L5\\2008";
//	ossimFilename outDir = "F:\\testdata\\landsat\\gcp_whole\\L5\\2008\\0_height";
//	ossimFilename elevationpath="F:\\testdata\\ossimdem";
//
//	if (!outDir.isDir())
//	{
//		_mkdir(outDir.c_str());
//	}
//
//	//ossimFilename elevationPath="E:\\l5\\srtm_57_06";
//	//get_elevation(infile, outfile, elevationPath);
//	
//	std::vector<ossimFilename> allGcpFiles;
//	//ossimString strReg = "/\.txt$/i ";
//	//ossimString strReg = "^(L5-TM-)+(.txt)$/i";
//	//ossimString strReg = "(L5-TM)[\\s\\S]*(.txt)$";
//	//ossimString strReg = "L5-TM[\s\S]*.txt";
//	ossimString strReg = "^L5-TM.*.txt$";
//	ossimDirectory(inGcpDir).findAllFilesThatMatch(allGcpFiles, strReg);
//
//	QStringList allHeaderFiles;
//	QFindFile(inSrcDir.c_str(), allHeaderFiles, QStringList("header.dat"));
//	//std::vector<ossimFilename> allHeaderFiles;
//	//for (int i = 0; i < (int)tempList.size(); i++)
//	//{
//	//	allHeaderFiles.push_back(ossimFilename(tempList[i].c_str()));
//	//}
//	//std::vector<ossimFilename> allHeaderFiles;
//	//strReg = "^(header.dat)$";
//	//strReg = "LD\\d{10}";
//	//ossimDirectory(inSrcDir).findAllFilesThatMatch(allHeaderFiles, strReg);
//
//	typedef std::pair<ossimFilename, ossimFilename> gcp_header_pair;
//	std::vector<gcp_header_pair> pairList;
//	for (int i = 0; i < (int)allGcpFiles.size(); i++)
//	{
//		// L5-TM-114-030-20080803-LD2010004377
//		// 0    1    2     3      4               5
//		vector<ossimString> strList = ossimString(allGcpFiles[i].fileNoExtension()).split('-');
//		ossimString productId = strList[5];
//		for (int j = 0; j < (int)allHeaderFiles.size(); j++)
//		{
//			QString foldName = QAfterLast(QBeforeLast(allHeaderFiles[j], '\\'), '\\');
//			if (foldName.contains(productId.c_str()))
//			{
//				gcp_header_pair pair;
//				pair.first = allGcpFiles[i];
//				pair.second = ossimFilename(allHeaderFiles[j].toLatin1());
//				pairList.push_back(pair);
//				allHeaderFiles.erase(allHeaderFiles.begin()+j);
//				printf("%s\n", pair.first);
//				landsat_optimization(pair.first, pair.second, elevationpath);
//				break;
//			}
//		}
//	}
//
//	//for (int i = 0;i < (int)allHeaderFiles.size();++i)
//	//{
//	//	ossimFilename infile = allHeaderFiles[i];
//	//	ossimFilename outfile = outDir + "\\" + infile.file();
//	//	MyProject prj;
//	//	prj.theMgr = ossimElevManager::instance();
//	//	if(!prj.theMgr->loadElevationPath(elevationPath))
//	//	{
//	//		cout<<"warning: 加载DEM失败！"<<endl;
//	//		return;
//	//	}
//	//	prj.ReadGcpAndProjection(infile);
//	//	prj.theMgr = ossimElevManager::instance();
//	//	prj.theMgr->loadElevationPath(elevationPath);
//	//	prj.m_DemPath = elevationPath;
//	//	prj.GetElevations(prj.m_CtrlGptSet);
//	//	prj.GetElevations(prj.m_ChkGptSet);
//	//	prj.SavePointToFile(outfile, prj.m_CtrlGptSet, prj.m_ChkGptSet);
//	//	
//	//}
//}

bool computer_error(ossimTieGptSet* gptset,ossimTieGptSet* chkset, std::vector<ossimGpt>& ctrErrorList, std::vector<ossimGpt>& chkErrorList)
{
	ctrErrorList.clear();
	chkErrorList.clear();
	ossimLeastSquaresBilin theLatFit;
	ossimLeastSquaresBilin theLonFit;
	ossimLeastSquaresBilin theXFit;
	ossimLeastSquaresBilin theYFit;
	theLatFit.clear();theYFit.clear();
	theLonFit.clear();theXFit.clear();
	const ossim_uint32 SIZE = gptset->size();
	if (SIZE < 4)
	{
		for (ossim_uint32 i = 0; i < SIZE; ++i) 
		{
			ossimGpt p;
			p.lat=0;p.lon=0;p.hgt=0;
			ctrErrorList.push_back(p);
		}

		for (size_t i = 0 ; i < chkset->getTiePoints().size(); ++i)
		{
			ossimGpt p;
			p.lat=0;p.lon=0;p.hgt=0;
			chkErrorList.push_back(p);
		}

		return true;
	}
	vector<ossimRefPtr<ossimTieGpt> >& theTPV = gptset->refTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
	ossimDpt imagepoint;
	ossimGpt goundpoint;
	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		theLatFit.addSample(imagepoint.x,
			imagepoint.y,
			goundpoint.lat);
		theLonFit.addSample(imagepoint.x,
			imagepoint.y,
			goundpoint.lon);
		theXFit.addSample(goundpoint.lon,
			goundpoint.lat,
			imagepoint.y);
		theYFit.addSample(goundpoint.lon,
			goundpoint.lat,
			imagepoint.x);
	}
	theLatFit.solveLS();
	theLonFit.solveLS();
	theXFit.solveLS();
	theYFit.solveLS();
	double pa,  pb_x,  pc_y, pd_xy;
	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		ossimGpt error1;
		theLatFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lat=(pa+pb_x*imagepoint.x+pc_y*imagepoint.y+pd_xy*imagepoint.x*imagepoint.y)-goundpoint.lat;
		theLonFit.getLSParms(pa,pb_x, pc_y, pd_xy);
		error1.lon=(pa+pb_x*imagepoint.x+pc_y*imagepoint.y+pd_xy*imagepoint.x*imagepoint.y)-goundpoint.lon;
		error1.hgt=sqrt(error1.lat*error1.lat+error1.lon*error1.lon);
		//cout<<error1<<endl;
		theXFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lon=(pa+pb_x*goundpoint.lon+pc_y*goundpoint.lat+pd_xy*goundpoint.lon*goundpoint.lat)-imagepoint.y;
		theYFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lat=(pa+pb_x*goundpoint.lon+pc_y*goundpoint.lat+pd_xy*goundpoint.lon*goundpoint.lat)-imagepoint.x;
		error1.hgt=sqrt(error1.lat*error1.lat+error1.lon*error1.lon);
		//cout<<error1<<endl;
		ctrErrorList.push_back(error1);
	}
	for (size_t i = 0 ; i < chkset->getTiePoints().size(); ++i)
	{
		imagepoint = chkset->refTiePoints()[i]->getImagePoint();
		goundpoint = chkset->refTiePoints()[i]->getGroundPoint();
		ossimGpt error1;
		theLatFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lat=(pa+pb_x*imagepoint.x+pc_y*imagepoint.y+pd_xy*imagepoint.x*imagepoint.y)-goundpoint.lat;
		theLonFit.getLSParms(pa,pb_x, pc_y, pd_xy);
		error1.lon=(pa+pb_x*imagepoint.x+pc_y*imagepoint.y+pd_xy*imagepoint.x*imagepoint.y)-goundpoint.lon;
		error1.hgt=sqrt(error1.lat*error1.lat+error1.lon*error1.lon);
		//cout<<error1<<endl;
		theXFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lon=(pa+pb_x*goundpoint.lon+pc_y*goundpoint.lat+pd_xy*goundpoint.lon*goundpoint.lat)-imagepoint.y;
		theYFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lat=(pa+pb_x*goundpoint.lon+pc_y*goundpoint.lat+pd_xy*goundpoint.lon*goundpoint.lat)-imagepoint.x;
		error1.hgt=sqrt(error1.lat*error1.lat+error1.lon*error1.lon);
		//cout<<error1<<endl;
		chkErrorList.push_back(error1);
	}
}

bool pm_computer_error(ossimTieGptSet* gptset,ossimTieGptSet* chkset, std::vector<ossimGpt>& ctrErrorList, std::vector<ossimGpt>& chkErrorList)
{
	ctrErrorList.clear();
	chkErrorList.clear();
	ossimLeastSquaresBilin theLatFit;
	ossimLeastSquaresBilin theLonFit;
	ossimLeastSquaresBilin theXFit;
	ossimLeastSquaresBilin theYFit;
	theLatFit.clear();theYFit.clear();
	theLonFit.clear();theXFit.clear();
	const ossim_uint32 SIZE = gptset->size();
	if (SIZE < 4)
	{
		for (ossim_uint32 i = 0; i < SIZE; ++i) 
		{
			ossimGpt p;
			p.lat=0;p.lon=0;p.hgt=0;
			ctrErrorList.push_back(p);
		}

		for (size_t i = 0 ; i < chkset->getTiePoints().size(); ++i)
		{
			ossimGpt p;
			p.lat=0;p.lon=0;p.hgt=0;
			chkErrorList.push_back(p);
		}

		return true;
	}
	vector<ossimRefPtr<ossimTieGpt> >& theTPV = gptset->refTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;
	ossimDpt imagepoint;
	ossimGpt goundpoint;
	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		theLatFit.addSample(imagepoint.x,
			imagepoint.y,
			goundpoint.lat);
		theLonFit.addSample(imagepoint.x,
			imagepoint.y,
			goundpoint.lon);
		theXFit.addSample(goundpoint.lon,
			goundpoint.lat,
			imagepoint.y);
		theYFit.addSample(goundpoint.lon,
			goundpoint.lat,
			imagepoint.x);
	}
	theLatFit.solveLS();
	theLonFit.solveLS();
	theXFit.solveLS();
	theYFit.solveLS();
	double pa,  pb_x,  pc_y, pd_xy;
	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		ossimGpt error1;
		theLatFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lat=(pa+pb_x*imagepoint.x+pc_y*imagepoint.y+pd_xy*imagepoint.x*imagepoint.y)-goundpoint.lat;
		theLonFit.getLSParms(pa,pb_x, pc_y, pd_xy);
		error1.lon=(pa+pb_x*imagepoint.x+pc_y*imagepoint.y+pd_xy*imagepoint.x*imagepoint.y)-goundpoint.lon;
		error1.hgt=sqrt(error1.lat*error1.lat+error1.lon*error1.lon);
		//cout<<error1<<endl;
		theXFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lon=(pa+pb_x*goundpoint.lon+pc_y*goundpoint.lat+pd_xy*goundpoint.lon*goundpoint.lat)-imagepoint.y;
		theYFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lat=(pa+pb_x*goundpoint.lon+pc_y*goundpoint.lat+pd_xy*goundpoint.lon*goundpoint.lat)-imagepoint.x;
		error1.hgt=sqrt(error1.lat*error1.lat+error1.lon*error1.lon);
		//cout<<error1<<endl;
		ctrErrorList.push_back(error1);
	}
	for (size_t i = 0 ; i < chkset->getTiePoints().size(); ++i)
	{
		imagepoint = chkset->refTiePoints()[i]->getImagePoint();
		goundpoint = chkset->refTiePoints()[i]->getGroundPoint();
		ossimGpt error1;
		theLatFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lat=(pa+pb_x*imagepoint.x+pc_y*imagepoint.y+pd_xy*imagepoint.x*imagepoint.y)-goundpoint.lat;
		theLonFit.getLSParms(pa,pb_x, pc_y, pd_xy);
		error1.lon=(pa+pb_x*imagepoint.x+pc_y*imagepoint.y+pd_xy*imagepoint.x*imagepoint.y)-goundpoint.lon;
		error1.hgt=sqrt(error1.lat*error1.lat+error1.lon*error1.lon);
		//cout<<error1<<endl;
		theXFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lon=(pa+pb_x*goundpoint.lon+pc_y*goundpoint.lat+pd_xy*goundpoint.lon*goundpoint.lat)-imagepoint.y;
		theYFit.getLSParms(pa,pb_x,pc_y,pd_xy);
		error1.lat=(pa+pb_x*goundpoint.lon+pc_y*goundpoint.lat+pd_xy*goundpoint.lon*goundpoint.lat)-imagepoint.x;
		error1.hgt=sqrt(error1.lat*error1.lat+error1.lon*error1.lon);
		//cout<<error1<<endl;
		chkErrorList.push_back(error1);
	}
}

bool RfmComputerError(ossimTieGptSet* gptset,ossimTieGptSet* chkset, std::vector<ossimGpt>& ctrErrorList, std::vector<ossimGpt>& chkErrorList)
{
	int nCtrl = (int)gptset->getTiePoints().size();
	ossimRpcModel* rpcModel = new ossimRpcModel();
	ossimRpcSolver *solver = new ossimRpcSolver(true, false);
	vector < ossimDpt > imagePoints(nCtrl);
	vector < ossimGpt > groundControlPoints(nCtrl);
	ossimGpt gpt;
	for(int i = 0;i < nCtrl;i++)
	{
		groundControlPoints[i] = gptset->getTiePoints()[i]->getGroundPoint();
		imagePoints[i] = gptset->getTiePoints()[i]->getImagePoint();
	}
	solver->solveCoefficients(imagePoints, groundControlPoints);
	ossimImageGeometry *imageGeom = solver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	rpcModel->loadState(geom, "projection.");

	ctrErrorList.clear();
	chkErrorList.clear();
	for (unsigned int i = 0;i < gptset->getTiePoints().size();++i)
	{
		ossimGpt gpt = gptset->getTiePoints()[i]->getGroundPoint();
		ossimDpt residue = gptset->getTiePoints()[i]->getImagePoint() - rpcModel->forward(gpt);
		ctrErrorList.push_back(ossimGpt(residue.x, residue.y));
	}

	for (unsigned int i = 0;i < chkset->getTiePoints().size();++i)
	{
		ossimGpt gpt = chkset->getTiePoints()[i]->getGroundPoint();
		ossimDpt residue = chkset->getTiePoints()[i]->getImagePoint() - rpcModel->forward(gpt);
		chkErrorList.push_back(ossimGpt(residue.x, residue.y));
	}

	return true;
}

void Gcp_Auto_Filter(ossimTieGptSet* ctrSet,ossimTieGptSet* chkSet, double threshold=2.5, int ntimes = 10)
{
	int nCtrl = ctrSet->getTiePoints().size();
	int nChk = chkSet->getTiePoints().size();
	int i, t;
	std::vector<ossimGpt> ctrErrorList;
	std::vector<ossimGpt> chkErrorList;
	computer_error(ctrSet,chkSet, ctrErrorList, chkErrorList);

	for(t = ntimes;t > 0; t--)
	{
		double threshold_new = threshold * t;
		int num = 0;

		for(i = 0;i < nCtrl;i++)
		{
			double error = sqrt(ctrErrorList[i].lat * ctrErrorList[i].lat + ctrErrorList[i].lon * ctrErrorList[i].lon);
			if(threshold_new < error)
			{
				chkSet->refTiePoints().push_back(ctrSet->getTiePoints()[i]);
				ctrSet->refTiePoints().erase(ctrSet->getTiePoints().begin() + i);
				chkErrorList.push_back(ctrErrorList[i]);
				ctrErrorList.erase(ctrErrorList.begin() + i);
				num++;
				nCtrl--;
				nChk++;
				i--;
			}
		}

		computer_error(ctrSet,chkSet, ctrErrorList, chkErrorList);
		for(;;)
		{
			for(;;)
			{
				num = 0;
				for(i = 0;i < nChk;i++)
				{
					double error = sqrt(chkErrorList[i].lat * chkErrorList[i].lat + chkErrorList[i].lon * chkErrorList[i].lon);
					if(threshold >= error)
					{
						ctrSet->refTiePoints().push_back(chkSet->getTiePoints()[i]);
						chkSet->refTiePoints().erase(chkSet->getTiePoints().begin() + i);
						ctrErrorList.push_back(chkErrorList[i]);
						chkErrorList.erase(chkErrorList.begin() + i);
						num++;
						nCtrl++;
						nChk--;
						i--;
					}
				}
				if(0 == num)
					break;
				computer_error(ctrSet,chkSet, ctrErrorList, chkErrorList);
			}

			num = 0;
			for(i = 0;i < nCtrl;i++)
			{
				double error = sqrt(ctrErrorList[i].lat * ctrErrorList[i].lat + ctrErrorList[i].lon * ctrErrorList[i].lon);
				if(threshold_new < error)
				{
					chkSet->refTiePoints().push_back(ctrSet->getTiePoints()[i]);
					ctrSet->refTiePoints().erase(ctrSet->getTiePoints().begin() + i);
					chkErrorList.push_back(ctrErrorList[i]);
					ctrErrorList.erase(ctrErrorList.begin() + i);
					num++;
					nCtrl--;
					nChk++;
					i--;
				}
			}
			if(0 == num)
				break;
			computer_error(ctrSet,chkSet, ctrErrorList, chkErrorList);
		}
	}
}


void distribution_optimization(const ossimTieGptSet*& allSet, int num, ossimTieGptSet* selectedSet, ossimTieGptSet* unselectedSet, int nSelect = 30)
{
	if (!selectedSet)
	{
		return;
	}
	int nTotal = (int)allSet->getTiePoints().size();
	if (nSelect > nTotal)
	{
		nSelect = nTotal;
	}
	if (0 == nSelect)
	{
		selectedSet->clearTiePoints();
		*unselectedSet = *allSet;
	}

	list<int> candidateList;
	for (int i = 0; i < nTotal; i++)
	{
		candidateList.push_back(i);
	}

	srand(unsigned(time(0)));
	int start_index = rand() % nTotal;

	selectedSet->addTiePoint(allSet->getTiePoints()[start_index]);
	candidateList.remove(start_index);
	
	for (int i = 0; i < num-1; i++)
	{
		double max_minDist = 0.0;
		double max_minIndex = 0;
		for (list<int>::iterator iter = candidateList.begin();iter != candidateList.end(); iter++)
		{
			double minDist = DBL_MAX;
			for (int j = 0; j < (int)selectedSet->getTiePoints().size(); j++)
			{
				ossimDpt allDpt = allSet->getTiePoints()[*iter]->getImagePoint();
				ossimDpt selectedDpt = selectedSet->getTiePoints()[j]->getImagePoint();
				double dist = (allDpt.x - selectedDpt.x)*(allDpt.x - selectedDpt.x);
				dist += (allDpt.y - selectedDpt.y)*(allDpt.y - selectedDpt.y);
				if (minDist > dist)
				{
					minDist = dist;
				}
			}
			if (max_minDist < minDist)
			{
				max_minDist = minDist;
				max_minIndex = *iter;
			}
		}
		selectedSet->addTiePoint(allSet->getTiePoints()[max_minIndex]);
	}

	if(unselectedSet)
	{
		unselectedSet->clearTiePoints();
		for (list<int>::iterator iter = candidateList.begin();iter != candidateList.end(); iter++)
		{
			unselectedSet->addTiePoint(allSet->getTiePoints()[*iter]);
		}
	}
}

void LandsatRpc(ossimFilename gcpFile, ossimFilename headerFile, ossimFilename elevationpath)
{
	ossimFilename workfold = headerFile.path();
	ossimFilename sourcefile = workfold + "\\header.dat";
	ossimFilename gcpfile = gcpFile;
	ossimFilename reportfile = workfold + "\\report.txt";

	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));
	prj.theMgr = ossimElevManager::instance();
	//prj.theMgr->loadElevationPath(ossimFilename(elevationpath));//
	prj.m_DemPath=ossimFilename(elevationpath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(6);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);
	prj.InitiateSensorModel(sourcefile);

	int num = static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());
	ossimRpcSolver *solver = new ossimRpcSolver(true, false);
	vector < ossimDpt > imagePoints(num);
	vector < ossimGpt > groundControlPoints(num);
	ossimGpt gpt;
	for(int i = 0;i < num;i++)
	{
		groundControlPoints[i] = prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint();
		ossimGpt gpt = prj.m_MapPar->inverse(ossimDpt(groundControlPoints[i].lat, groundControlPoints[i].lon));
		gpt.hgt = groundControlPoints[i].hgt;
		groundControlPoints[i] = gpt;
		imagePoints[i] = prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint();
	}
	solver->solveCoefficients(imagePoints, groundControlPoints);

	ossimImageGeometry *imageGeom = solver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimRpcModel* rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom, "projection.");

	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	//// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;fstream rpcStruct;
	fstream rpcStruct;
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.m_sensorModel->saveState(prj.geom);
	rpcModel->loadState(prj.geom);
	rpcStruct.open((workfold+"\\rpcStruct0.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	//prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);
	//return;

	//double tmp;
	//for(int i = 0;i < num;i++)
	//{
	//	groundControlPoints.push_back(prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint());
	//	imagePoints.push_back(prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint());
	//}
	//prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, &tmp);		//用优化控制点进行模型优化
	//prj.m_sensorModel->updateModel();
	//prj.m_sensorModel->saveState(prj.geom);
	//prj.OutputReport(reportfileall, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	//for(int i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	//{
	//	ossimDpt dpt = prj.m_MapPar->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
	//	ossimGpt gpt(dpt.x,dpt.y);
	//	prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(ossimGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	//}
	
	//// 全参数优化
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.m_sensorModel->saveState(prj.geom);
	//rpcModel->loadState(prj.geom);
	// 像方仿射变换
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);
	return;
	//prj.Orthograph(outfile);
}

void Gcp_Auto_Filter(ossimFilename inGcpFile, ossimFilename outGcpFile, double threshold=2.5, int ntimes = 10)
{
	MyProject prj;
	prj.ReadGcpAndProjection(inGcpFile);

	Gcp_Auto_Filter(prj.m_CtrlGptSet, prj.m_ChkGptSet, threshold, ntimes);

	//prj.SavePointToFile(outGcpFile, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	prj.SavePointToFile(outGcpFile, prj.m_CtrlGptSet, NULL);
}

void batch_gcp_filter()
{
	ossimFilename inDir = "I:\\testdata\\landsat\\gcp\\gcp_all";
	ossimFilename outDir = "I:\\testdata\\landsat\\gcp\\gcp_filtered";
	std::vector<ossimFilename> allGcpFiles;
	ossimString strReg = "^.*-GCP-ALL.txt$";
	ossimDirectory(inDir).findAllFilesThatMatch(allGcpFiles, strReg);
	if (!outDir.exists())
	{
		_mkdir(outDir.c_str());
	}

	int nTotal = (int)allGcpFiles.size();
	for (int i = 0;i < nTotal;++i)
	{
		ossimFilename outFile = outDir + "\\" + allGcpFiles[i].file();
		Gcp_Auto_Filter(allGcpFiles[i], outFile, 3.5, 10);
	}
}

//ossimProjection* newUtmView(const ossimGpt& centerGround,
//						   const ossimDpt& metersPerPixel)
//{
//	ossimUtmProjection* utm = new ossimUtmProjection;
//
//	// we will make it a square pixel in meters
//	double averageGsd = (metersPerPixel.x + metersPerPixel.y)*.5;
//	utm->setZone(centerGround);
//	utm->setMetersPerPixel(ossimDpt(metersPerPixel));
//
//	return utm;
//}

void reprojection()
{
	MyProject prj;
	ossimTieGptSet *gptSet = new ossimTieGptSet;
	ossimTieGptSet *chkSet = new ossimTieGptSet;
	prj.ReadGcpAndProjection("F:\\testdata\\alos_batch\\multiple\\gcp_hgt.txt", gptSet, chkSet);
	ossimProjection *utm = newUtmView(ossimGpt(0.0,123.0,0), ossimDpt(10.0, 10.0));

	vector<ossimRefPtr<ossimTieGpt> >& theTPV = gptSet->refTiePoints();

	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		ossimDpt imagepoint=(*tit)->getImagePoint();
		ossimGpt goundpoint=(*tit)->getGroundPoint();
		ossimGpt tGpt = prj.m_MapPar->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
		tGpt.hgt = (*tit)->hgt;

		ossimDpt utmDpt = utm->forward(tGpt);
		ossimGpt utmGpt(utmDpt.lon, utmDpt.lat, tGpt.hgt);
		(*tit)->setGroundPoint(utmGpt);
	}

	prj.SavePointToFile("F:\\testdata\\alos_batch\\multiple\\gcp_utm.txt", gptSet, NULL);

}

void initialProjection()
{
	ProjectionParameters pp;
	pp.DatumName = "西安80";
	pp.TrueOriginLongitude = 123;
	pp.TrueOriginLatitude = 0.0;
	pp.EastingFalse = 21500000.0;
	pp.NorthingFalse = 0.0;
	pp.PixelSizeX = 1.0;
	pp.PixelSizeY = 1.0;

	ossimKeywordlist kwl;
	ossimMapProjection* pMapProjection = CreateProjection(pp,kwl);
	fstream fs;
	fs.open("D:\\workspace\\Radarsat2\\projection.txt", std::ios_base::out);
	fs<<kwl;
	fs.close();
}

void Rfm_Outlier(ossimTieGptSet*& gcpSet, ossimTieGptSet*& chkSet)
{
	std::vector<ossimGpt> ctrErrorList;
	std::vector<ossimGpt> chkErrorList;
	int nCtrl = gcpSet->getTiePoints().size();
	int nChk = chkSet->getTiePoints().size();
	int i, t;
	int ntimes = 5;
	double threshold = 2.0; // pixel

	RfmComputerError(gcpSet, chkSet, ctrErrorList, chkErrorList);

	for(t = ntimes;t > 0; t--)
	{
		double threshold_new = threshold * t;
		int num = 0;

		for(i = 0;i < nCtrl;i++)
		{
			double error = sqrt(ctrErrorList[i].lat * ctrErrorList[i].lat + ctrErrorList[i].lon * ctrErrorList[i].lon);
			if(threshold_new < error)
			{
				chkSet->refTiePoints().push_back(gcpSet->getTiePoints()[i]);
				gcpSet->refTiePoints().erase(gcpSet->getTiePoints().begin() + i);
				chkErrorList.push_back(ctrErrorList[i]);
				ctrErrorList.erase(ctrErrorList.begin() + i);
				num++;
				nCtrl--;
				nChk++;
				i--;
			}
		}
		RfmComputerError(gcpSet, chkSet, ctrErrorList, chkErrorList);
		for(;;)
		{
			for(;;)
			{
				num = 0;
				for(i = 0;i < nChk;i++)
				{
					double error = sqrt(chkErrorList[i].lat * chkErrorList[i].lat + chkErrorList[i].lon * chkErrorList[i].lon);
					if(threshold >= error)
					{
						gcpSet->refTiePoints().push_back(chkSet->getTiePoints()[i]);
						chkSet->refTiePoints().erase(chkSet->getTiePoints().begin() + i);
						ctrErrorList.push_back(chkErrorList[i]);
						chkErrorList.erase(chkErrorList.begin() + i);
						num++;
						nCtrl++;
						nChk--;
						i--;
					}
				}
				if(0 == num)
					break;
				RfmComputerError(gcpSet, chkSet, ctrErrorList, chkErrorList);
			}

			num = 0;
			for(i = 0;i < nCtrl;i++)
			{
				double error = sqrt(ctrErrorList[i].lat * ctrErrorList[i].lat + ctrErrorList[i].lon * ctrErrorList[i].lon);
				if(threshold_new < error)
				{
					chkSet->refTiePoints().push_back(gcpSet->getTiePoints()[i]);
					gcpSet->refTiePoints().erase(gcpSet->getTiePoints().begin() + i);
					chkErrorList.push_back(ctrErrorList[i]);
					ctrErrorList.erase(ctrErrorList.begin() + i);
					num++;
					nCtrl--;
					nChk++;
					i--;
				}
			}
			if(0 == num)
				break;
			RfmComputerError(gcpSet, chkSet, ctrErrorList, chkErrorList);
		}
	}
}

void GcpsReformat(ossimFilename infile, ossimFilename outfile)
{
	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	readGcpFile(infile, gcpSet, chkSet, &prjKwl);

	ossimElevManager::instance()->loadElevationPath("D:\\workspace\\dem\\srtm90");
	vector<ossimRefPtr<ossimTieGpt> >& theGcp = gcpSet->refTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::iterator iter,tit;
	ossimGpt tGpt;
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(ossimFilename("E:\\HJ1\\test1\\masterOut1.tif"));
	//ossimProjection *prj = handler->getImageGeometry()->getProjection();
	for (iter = theGcp.begin() ; iter != theGcp.end() && iter->valid() ;)
	{

		tGpt = (*iter)->getGroundPoint();
		//ossimDpt dpt = prj->forward(tGpt);	// ll to projection

		double hgt =ossimElevManager::instance()->getHeightAboveEllipsoid(tGpt);
		//(*iter)->setGroundPoint(ossimGpt(dpt.x, dpt.y, hgt));
		(*iter)->hgt = hgt;
		if(ossim::isnan((*iter)->hgt))// || (*iter)->hgt == 0.0)
			//if(ossim::isnan((*iter)->hgt) || (*iter)->hgt == 0.0)
		{
			// 如果为空，则去掉该点
			theGcp.erase(iter);
			// 如果为空，则赋为平均值
			//(*iter)->hgt = default_hgt;
			++iter;
		}
		else{
			++iter;
		}
	}

	Rfm_Outlier(gcpSet, chkSet);

	saveGcpFile(outfile, gcpSet, chkSet, NULL);
}


//bool correctEnvisatASAR()
//{	
//	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossim_plugin.dll");
//	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimgdal_plugin.dll");
//	ossimInit::instance()->initialize();
//
//	//ossimFilename workfold = "D:\\workspace\\Envisat\\beijing";
//	//ossimFilename sourceFile = workfold + "\\ASA_IMP_1PNBEI20040114_022837_000000172023_00218_09791_1828.N1";
//	ossimFilename workfold = "D:\\workspace\\Radarsat2";
//	ossimFilename sourceFile = workfold + "\\imagery_HH.tif";
//	//ossimFilename demPath = "E:\\dem_bj";
//	//ossimFilename demPath = "D:\\workspace\\dem";
//	//ossimFilename demPath = "D:\\workspace\\beijingDem";
//	ossimFilename demPath = "D:\\workspace\\dem-hgt";
//	ossimFilename gcpfile = workfold + "\\gcps.txt";
//	ossimFilename projection_file = workfold + "\\projection.txt";
//	//ossimFilename gcpfile = workfold + "\\gcp_80_12.txt";
//	//ossimFilename chkFile = workfold + "\\gcp.txt";
//	ossimFilename outFile = workfold + "\\rect.tif";
//	ossimFilename reportFile0 = workfold + "\\report.txt";
//
//	
//	ossimElevManager::instance()->loadElevationPath(demPath);
//
//	//ossimSensorModel* pSensorModel =  new ossimEnvisatAsarModel();
//	//ossimEnvisatAsarModel* pEnvisatAsarModel = new ossimEnvisatAsarModel();
//	ossimRadarSat2RPCModel* pRadarSat2RPCModel = new ossimRadarSat2RPCModel();
//	pRadarSat2RPCModel->open(sourceFile);
//	ossimRpcModel* pRpcModel = new ossimRpcModel;
//	ossimRpcModel::rpcModelStruct rpcStruct;
//	pRadarSat2RPCModel->getRpcParameters(rpcStruct);
//	pRpcModel->setAttributes(rpcStruct);
//	//ossimGeometricSarSensorModel* pSensorModel = pEnvisatAsarModel;
//
//	ossimKeywordlist model_kwl;
//	std::list<ossimGpt> groundCoordinates;
//    std::list<ossimDpt> imageCoordinates;
//	//pRadarSat2RPCModel->optimizeModel(groundCoordinates, imageCoordinates);
//	//pRadarSat2RPCModel->optimizeFit();
//	//pRadarSat2RPCModel->updateModel();
//	pRpcModel->saveState(model_kwl);
//	pRpcModel->loadState(model_kwl);
//	//cout<<model_kwl<<endl;
//
//	
//	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(ossimFilename(sourceFile));
//	if(!handler)
//	{
//		//	cout << "Unable to open input image: "<< endl;
//		return false;
//	}
//	ossimKeywordlist in_geom_kwl,out_geom_kwl,tt_geom_kwl;
//	handler->getImageGeometry()->saveState(in_geom_kwl);
//	//cout<<in_geom_kwl<<endl;
//	
//	ossimImageGeometry imageGeom;
//	imageGeom.loadState(model_kwl);
//	handler->setImageGeometry(&imageGeom);
//	//handler->loadState(model_kwl);
//
//	// 指定输出投影
//	out_geom_kwl.addFile(projection_file);
//	ossimRefPtr<ossimProjection> proj;
//	proj = ossimProjectionFactoryRegistry::instance()->createProjection(out_geom_kwl);
//	ossimMapProjection* pNewProjection = PTR_CAST(ossimMapProjection, proj.get());
//
//	
//	// 选择输出波段	
//	ossimBandSelector* theBandSelector;
//	theBandSelector = new ossimBandSelector;
//	theBandSelector->connectMyInputTo(0, handler);	
//	int nbands = handler->getNumberOfInputBands();
//	vector<ossim_uint32> outBandList(nbands);
//	for (int i = 0; i < nbands; i++)
//	{
//		outBandList[i] = i;
//	}
//	theBandSelector->setOutputBandList(outBandList);
//	
//	// 输出指定范围图像
//	vector<ossimDpt> polygon;
//	ossimIrect bound,boundw;
//	ossimPolyCutter* theCutter = new ossimPolyCutter;
//	ossimDpt imagesize(handler->getImageRectangle().width(), handler->getImageRectangle().height());
//	int starline,starpixel,endpixel,endline;
//	starline	=	0;
//	starpixel	=	0;
//	endpixel	=	0;
//	endline		=	0;
//	if (0 == endpixel)
//	{
//		endpixel = imagesize.x - 1;
//	}
//	if (0 == endline)
//	{
//		endline = imagesize.y - 1;
//	}/////////////////////////////////////////////////////////////////////
//	ossimDpt ps(starpixel,starline),p2(endpixel,starline),p3(endpixel,endline),p4(starpixel,endline),p5(starpixel,starline);
//	polygon.push_back(ps);
//	polygon.push_back(p2);
//	polygon.push_back(p3);
//	polygon.push_back(p4);
//	polygon.push_back(p5);
//	theCutter->connectMyInputTo(theBandSelector);
//	theCutter->setPolygon(polygon);
//	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
//	theCutter->setNumberOfPolygons(1);
//	
//
//	ossimImageRenderer* renderer = new ossimImageRenderer;	
//	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
//	renderer->connectMyInputTo(theCutter);
//	renderer->setView(pNewProjection);
//	
//	
//	ossimImageFileWriter* writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
//	if(!writer)
//	{
//		//	cout << "Unable to writer input image: " << endl;
//		return 0;
//	}
//	writer->saveState(tt_geom_kwl);
//	tt_geom_kwl.add("pixel_type","area",true);
//	writer->loadState(tt_geom_kwl);
//	// cout<<"writer  "<<endl;
//	//cout<<tt_geom;
//	writer->setFilename(outFile);
//	//ossimStdOutProgress progress(0,true);
//	//writer->addListener(&progress);
//	// 开始监听	
//	myOutProgress* progress = new myOutProgress(0,true);
//	writer->addListener(progress);
//
//	// 写图像
//	writer->connectMyInputTo(0,renderer);
//	writer->execute();
//
//	// 结束监听
//	writer->disableListener();
//	writer->removeListener(progress);
//
//	// 清理
//	delete progress;
//	delete writer;
//	handler->close();
//	delete theBandSelector;
//	in_geom_kwl.clear();
//	out_geom_kwl.clear();
//	tt_geom_kwl.clear();
//	return true;
//}

void getResidues()
{
	string folder = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401070348_201401070354\\check";
	string strImageFile = folder + "\\1792bias.tif";
	string gcpFile = folder + "\\check1792bias.txt";
	string residueFile = folder + "\\residue1792bias.txt";

	ossimTieGptSet *gptSet = ReadGptFromFile(gcpFile);
	ossimRefPtr<ossimImageHandler> handler = ossimImageHandlerRegistry::instance()->open(strImageFile);
	ossimRefPtr<ossimProjection> prj = handler->getImageGeometry()->getProjection();

	fstream ofs;
	ofs.open(residueFile.c_str(), ios_base::out);
	for (int i = 0;i < gptSet->getTiePoints().size();++i)
	{
		ossimGpt ll;
		prj->lineSampleToWorld(gptSet->getTiePoints()[i]->getImagePoint(), ll);
		ossimDpt gpt = prj->forward(ll);
		double delta_lat = gpt.x - gptSet->getTiePoints()[i]->getGroundPoint().lat;
		double delta_lon = gpt.y - gptSet->getTiePoints()[i]->getGroundPoint().lon;
		if (0 != i)
		{
			ofs<<endl;
		}
		ofs<<i+1<<"\t"
			<<gptSet->getTiePoints()[i]->getImagePoint().x<<"\t"
			<<gptSet->getTiePoints()[i]->getImagePoint().y<<"\t"
			//<<gptSet->getTiePoints()[i]->getGroundPoint().lat<<"\t"
			//<<gptSet->getTiePoints()[i]->getGroundPoint().lon<<"\t"
			//<<gptSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
			<<delta_lat<<"\t"
			<<delta_lon;
	}
	ofs.close();
}

int main(int argc, char* argv[])
{
	string dropbox = getenv("DROPBOX");
	//string preferences_file = dropbox + "\\Programs\\ossimTest_2.0\\preference_office.txt";
	string preferences_file = "D:\\opensource\\ossim\\preference.txt";

	ossimPreferences::instance()->loadPreferences(ossimFilename(preferences_file));
	ossimInit::instance()->initialize();
	double t = (double)getTickCount();

	//landsat_optimization("D:\\workspace\\LD2010003816\\200.txt", "D:\\workspace\\LD2010003816\\header.dat", "D:\\workspace\\dem\\aster30");
	feature_landsat();
	//Spot_create3DGridPoints();
	//getResidues();
	//Rcp_Generic();
	//spotTest();
	//line_match();
	//ossimFilename strSource = "E:\\testdata\\source.TIF";
	//ossimFilename strRefer = "E:\\testdata\\ref.TIF";
	//ossimFilename strSource = "D:\\workspace\\source.TIF";
	//ossimFilename strRefer = "D:\\workspace\\ref.TIF";
	//ossimFilename strSource = "E:\\Dropbox\\Programs\\data\\ding\\1.tif";
	//ossimFilename strRefer = "E:\\Dropbox\\Programs\\data\\ding\\3.tif";
	//openSurMatch(strRefer, strSource);
	//opencvPointMatch(strSource);
	//linesegment_detection(argc, argv);
	//get_elevation();
	//batch_get_elevation();
	//landsat_optimization();
	//batch_gcp_optimization();
	//landsat_optimization("I:\\testdata\\landsat\\l5\\gcp1\\1\\gcp_300.txt", "I:\\testdata\\landsat\\l5\\gcp1\\1\\L2\\header.dat", "I:\\testdata\\ossimdem");
	//LandsatRpc("I:\\testdata\\landsat\\l5\\gcp1\\1\\gcp_80+200.txt", "I:\\testdata\\landsat\\l5\\gcp1\\1\\L2\\header.dat", "I:\\testdata\\ossimdem");
	//AlosCreateVirtualGpcs("D:\\workspace\\alos\\d1001883-009_alpsmw254402725_o1b2r_uw_rpc", "D:\\workspace\\alos\\d1001883-009_alpsmw254402725_o1b2r_uw_rpc\\virtualGcps.txt");
	//AlosCreateVirtualGpcs("I:\\testdata\\alos\\beijing\\d1002016-107_ALPSMW253092805_O1B2R_UW_RPC", "D:\\workspace\\alos\\d1001883-009_alpsmw254402725_o1b2r_uw_rpc\\virtualGcps.txt");
	//batch_AlosCreateVirtualGpcs("I:\\testdata\\alos\\jl_source");
	//batch_AlosCreateVirtualGpcs("F:\\testdata\\alos\\ln_jl_mul");
	//batch_AlosCreateVirtualGpcs("F:\\testdata\\alos_batch\\multiple");
	//Gcp_Auto_Filter("F:\\testdata\\landsat\\l5\\gcp1\\1\\L5-TM-136-032-20090528-LD2010002092-L4-GCP-ALL.txt", "F:\\testdata\\landsat\\l5\\gcp1\\1\\filtered_gcp.txt", 2.5, 10);
	//batch_gcp_filter();
	//Rcp_Generic();
	//straightline_landsat();
	//AlosRpc();
	//reprojection();
	//AlosBatchOrth("I:\\testdata\\alos_batch\\source", "I:\\testdata\\alos_batch\\output");
	//spotTest();
	//initialProjection();
	//feature_radarsat2();
	//correctEnvisatASAR();
	//feature_landsat();
	//feature_gf1();
	//AlosRpc_createTieGptSet();
	//if(ossimMpi::instance()->getRank() == 0)
	//{
		//Hj1Test();
		//Spot_create3DGridPoints();
	//}
	//Hj1Rpc();

	//reprojectionPoints("E:\\projection\\gcp_hgt.txt", "E:\\projection\\outProjection.txt", "E:\\projection\\wgs84_hgt.txt");
	//projection2ll("E:\\projection\\gcp_hgt.txt", "E:\\projection\\ll_hgt.txt");
	//Spot_create3DGridPoints();
	//TheosTest();
	//Theos_create3DGridPoints();

	ossimElevManager::instance()->loadElevationPath("D:\\workspace\\dem\\srtm90");
	//matchPoint("D:\\workspace\\Landsat\\L5146030_03020070817.TM-GLS2005\\L5146030_03020070817_B10.TIF", \
		"D:\\workspace\\Landsat\\146030_UTM.tif", "D:\\workspace\\Landsat\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\Landsat\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat", \
		"D:\\workspace\\Landsat\\146030_UTM.tif", "D:\\workspace\\Landsat\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\Landsat\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat", \
	"I:\\testdata\\2000_UTM_std\\si.shp", "D:\\workspace\\Landsat\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\HJ\\db\\HJ1B-CCD2-449-56-20100917-L20000393774.tif", \
		"I:\\testdata\\std\\index.shp", "D:\\workspace\\HJ\\db\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\HJ\\db1\\hj_tm_119027.tif", \
		"D:\\workspace\\HJ\\db1\\LC81190272013280BJC00.TIF", "D:\\workspace\\HJ\\db1\\ossim_matching.txt");
	
	//matchPoint("D:\\workspace\\Landsat\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat", \
		"D:\\workspace\\Landsat\\si.shp", "D:\\workspace\\Landsat\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\Landsat\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat", \
		"D:\\workspace\\Landsat\\L5146030_03020070817.TM-GLS2005\\L5146030_03020070817_B10.TIF", "D:\\workspace\\Landsat\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\Spot5\\beijing\\281268-20040602-2.5\\scene01\\imagery.tif", \
		"D:\\workspace\\Spot5\\beijing\\281268_20040602.tif", "D:\\workspace\\Spot5\\beijing\\ossim-matching.txt");
	//SpotRpcTest1();

	//matchPoint("D:\\workspace\\HJ\\db1\\hj_tm_119027.tif", \
		"D:\\workspace\\HJ\\db1\\tm_119027.tif", "D:\\workspace\\HJ\\db1\\ossim_matching.txt");
	//matchPoint("D:\\workspace\\HJ\\db\\HJ1B-CCD2-449-56-20100917-L20000393774.TIF", \
		"D:\\workspace\\HJ\\db1\\refs-117\\index.shp", "D:\\workspace\\HJ\\db1\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\HJ\\db\\HJ1B-CCD2-449-56-20100917-L20000393774.tif", \
		"I:\\testdata\\tmchina\\china.vrt", "D:\\workspace\\HJ\\db\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\HJ\\db1\\hj_tm_119027.tif", \
		"D:\\workspace\\HJ\\db1\\index.shp", "D:\\workspace\\HJ\\db1\\ossim_matching.txt");

	//matchPoint("E:\\HJ1\\test1\\ccd1.tif", \
		"E:\\HJ1\\test1\\121037-l8.tif", "E:\\HJ1\\test1\\ossim_matching.txt");

	//matchPoint("D:\\workspace\\LD2010003816\\header.dat", \
		"D:\\workspace\\LD2010003816\\TM_123036.tif", "D:\\workspace\\LD2010003816\\ossim_matching.txt");

	//GcpsReformat("E:\\HJ1\\test1\\gcp1.xml", "E:\\HJ1\\test1\\gcp.txt");

	//matchPoint("I:\\testdata\\landsat\\test\\ezhou.tif", \
	"I:\\testdata\\tmchina\\china.vrt", "I:\\testdata\\landsat\\test\\ossim_matching.txt");

	//Landsat_create3DGridPoints();
	//AlosRpc_createTieGptSet();

//#pragma omp parallel num_threads(8)
//	{
//		cout<<"Hello world!\n";
//	}

	t = ((double)getTickCount() - t)/getTickFrequency();
	std::cout << "detection time [s]: " << t/1.0 << std::endl;
	cout<<"OK!"<<endl;
	//system("pause");
	return 0;
}