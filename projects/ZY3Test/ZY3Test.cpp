#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <gcpUtil.h>
#include <func.h>
#include "omp.h"
#include <ossim/base/ossimPreferences.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>
#include <ossim/projection/ossimSensorModel.h>
#include <ossim_plugin/radi/ossimHj1Model.h>
#include <time.h>

using namespace std;
using namespace mylib;

#ifndef _WIN64
#pragma comment(lib, "ossim20.lib")
#pragma comment(lib, "blas_win32_MT.lib")
#pragma comment(lib, "lapack_win32_MT.lib")
#else
#pragma comment(lib, "ossim20x64.lib")
#pragma comment(lib, "blas_win64_MT.lib")
#pragma comment(lib, "lapack_win64_MT.lib")
#endif
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "ossim_plugin.lib")
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")

void ZY3Test()
{
	ossimFilename workfold = "E:\\HJ1\\ZY-3\\MUX\\Scene02";
	//workfold = "E:\\HJ1\\ZY-3\\MUX\\Scene03";
	ossimFilename inputFile = workfold + "\\IMAGE.TIF";
	//inputFile = workfold + "\\IMAGE.jp2";
	ossimFilename virtual_gcpfile = workfold + "\\virtual_gcp.txt";
	ossimFilename virtual_chkfile = workfold + "\\virtual_chk.txt";


	ossimRefPtr<ossimImageHandler> handler = ossimImageHandlerRegistry::instance()->open(inputFile);
	if (!handler)
	{
		printf("Failed to open image file: %s\n", inputFile);
		return;
	}

	ossimRefPtr<ossimSensorModel> sensorModel;

	sensorModel = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());
	
	ossimTieGptSet* gptSet = new ossimTieGptSet;
	double max_height = 1100.0;
	double min_height = 0.0;

	ossimIrect imageRect(handler->getImageRectangle().ul(),
		handler->getImageRectangle().lr());
	int nLevels = 5;
	nLevels = 1;
	for (int i = 0; i < nLevels; ++i)
	{
		double hgt = min_height + i*(max_height - min_height) / max(1, (nLevels - 1));
		mylib::create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 15, 15, true, false);
		//create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 8821, 1, true, false);
		//create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 1, 9215, true, false);
	}

	ossimTieGptSet* chkSet = new ossimTieGptSet;
	nLevels = 10;
	for (int i = 0; i < nLevels; ++i)
	{
		double hgt = min_height + i*(max_height - min_height) / max(1, (nLevels - 1));
		mylib::create3DGridPoints(imageRect, *sensorModel, hgt, chkSet, 30, 30, true, false);
	}
	handler->close();

	mylib::saveGcpFile(virtual_gcpfile, gptSet, NULL, false);
	mylib::saveGcpFile(virtual_chkfile, chkSet, NULL, false);
}


void Cbers04Test()
{
	//ossimFilename workfold = "E:\\HJ1\\CBERS04\\P5\\Scene09";
	ossimFilename workfold = "E:\\HJ1\\CBERS04\\test\\Scene03";
	workfold = "E:\\HJ1\\CBERS04\\CBERS04_MUX_MYC_201501230418_201501230428\\Scene31";
	//workfold = "E:\\HJ1\\CBERS04\\CBERS04_MUX_MYC_201501250132_201501250137\\Scene07";
	//workfold = "E:\\HJ1\\CBERS04\\CBERS04_P5M_MYC_201501250132_201501250137\\Scene14";
	//workfold = "E:\\HJ1\\CBERS04\\P5\\Scene07";
	//workfold = "E:\\HJ1\\ZY-3\\MUX\\Scene03";
	ossimFilename inputFile = workfold + "\\IMAGE.TIF";
	//inputFile = workfold + "\\IMAGE.jp2";
	ossimFilename virtual_gcpfile = workfold + "\\virtual_gcp.txt";
	ossimFilename virtual_chkfile = workfold + "\\virtual_chk.txt";


	ossimRefPtr<ossimImageHandler> handler = ossimImageHandlerRegistry::instance()->open(inputFile);
	if (!handler)
	{
		printf("Failed to open image file: %s\n", inputFile);
		return;
	}

	ossimRefPtr<ossimSensorModel> sensorModel;

	sensorModel = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());

	ossimTieGptSet* gptSet = new ossimTieGptSet;
	double max_height = 1100.0;
	double min_height = 0.0;

	ossimIrect imageRect(handler->getImageRectangle().ul(),
		handler->getImageRectangle().lr());
	int nLevels = 5;
	//nLevels = 1;
	for (int i = 0; i < nLevels; ++i)
	{
		double hgt = min_height + i*(max_height - min_height) / max(1, (nLevels - 1));
		mylib::create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 15, 15, true, false);
		//create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 8821, 1, true, false);
		//create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 1, 9215, true, false);
	}

	ossimTieGptSet* chkSet = new ossimTieGptSet;
	nLevels = 10;
	for (int i = 0; i < nLevels; ++i)
	{
		double hgt = min_height + i*(max_height - min_height) / max(1, (nLevels - 1));
		mylib::create3DGridPoints(imageRect, *sensorModel, hgt, chkSet, 30, 30, true, false);
	}
	handler->close();

	mylib::saveGcpFile(virtual_gcpfile, gptSet, NULL, false);
	mylib::saveGcpFile(virtual_chkfile, chkSet, NULL, false);
}


void HJ1Test()
{
	//ossimFilename workfold = "E:\\HJ1\\CBERS04\\P5\\Scene09";
	ossimFilename workfold = "E:\\HJ1\\Scene08";
	//workfold = "E:\\HJ1\\ZY-3\\MUX\\Scene03";
	ossimFilename inputFile = workfold + "\\IMAGE.TIF";
	//inputFile = workfold + "\\IMAGE.jp2";
	ossimFilename virtual_gcpfile = workfold + "\\virtual_gcp.txt";
	ossimFilename virtual_chkfile = workfold + "\\virtual_chk.txt";


	ossimRefPtr<ossimImageHandler> handler = ossimImageHandlerRegistry::instance()->open(inputFile);
	if (!handler)
	{
		printf("Failed to open image file: %s\n", inputFile);
		return;
	}
	ossimRefPtr<ossimSensorModel> sensorModel;

	//ossimRefPtr<ossimplugins::ossimHj1Model> sensorModel = new ossimplugins::ossimHj1Model();

	//ossimFilename strPath = inputFile.path();
	//ossimString descfile = strPath + "\\desc.xml";
	//sensorModel->setupOptimizer(descfile);


	//// HJ model
	//ossimRefPtr<ossimImageHandler> handler = ossimImageHandlerRegistry::instance()->open(inputFile);
	//if (!handler)
	//{
	//	printf("Failed to open image file: %s\n", inputFile);
	//	return;
	//}

	//ossimRefPtr<ossimSensorModel> sensorModel;

	//sensorModel = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());


	// rpc model
	ossimplugins::radiRpcModel *rpcModel = new ossimplugins::radiRpcModel;
	if (!rpcModel->parseRpcFile(inputFile))
	{
		//
		cout << "Failed to find rpc file for rpc model." << endl;
		return;
	}
	rpcModel->setUseL1(false);
	sensorModel = PTR_CAST(ossimSensorModel, rpcModel);

	//sensorModel = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());

	double max_height = 1100.0;
	double min_height = 0.0;

	ossimIrect imageRect(handler->getImageRectangle().ul(),
		handler->getImageRectangle().lr());
	//ossimIrect imageRect(handler->getImageRectangle().ul()+ossimDpt(100,100),
	//	handler->getImageRectangle().lr() + ossimDpt(-100, -100));

	//int np = 100;
	//double step = imageRect.height() / (double)np;
	//for (size_t i = 0; i < np; i++)
	//{
	//	double line = i * step;
	//	double t_line;
	//	sensorModel->theSupportData->getLineTime(line, t_line);
	//	ossimDpt3d att(0.0, 0.0, 0.0);
	//	sensorModel->theSupportData->getAttitude(t_line, att);
	//	printf("%7.1lf%15.6lf%15.6lf%15.6lf\n", line, att.x, att.y, att.z);
	//}

	int nl = 30;
	int ns = 30;
	double xnorm, ynorm;

	ossim_float64 w = imageRect.width();
	ossim_float64 h = imageRect.height();
	ossimDpt ul = imageRect.ul();
	for (size_t y = 0; y < nl; y++)
	{
		ossimTieGptSet* gptSet = new ossimTieGptSet;
		for (size_t x = 0; x < ns; ++x)
		{
			ossimDpt imagePoint;
			if (nl > 1)
			{
				ynorm = (double)y / (double)(nl - 1);
			}
			else
			{
				ynorm = 0.0;
			}
			if (ns > 1)
			{
				xnorm = (double)x / (double)(ns - 1);
			}
			else
			{
				xnorm = 0.0;
			}

			ossimDpt dpt((int)((w - 1)*xnorm + ul.x),
				int((h - 1)*ynorm + ul.y));
			ossimGpt gpt;
			sensorModel->lineSampleToWorld(dpt, gpt);
			imagePoint = dpt;
			dpt = ossimDpt(gpt.lat, gpt.lon);
			ossimGpt groundPoint(dpt.x, dpt.y, gpt.hgt);

			ossimString strId;
			char tmpStr[256];
			sprintf(tmpStr, "%d", x+1);
			strId = tmpStr;
			ossimTieGpt *aTiePt = new ossimTieGpt(groundPoint, imagePoint, 1.0, strId);
			gptSet->addTiePoint(aTiePt);
		}

		char buf[256];
		sprintf_s(buf, "%s\\gcp\\%d.txt", workfold.c_str(), y);
		mylib::saveGcpFile(ossimFilename(buf), gptSet, NULL, NULL, true);
	}

	int nLevels = 5;
	nLevels = 1;
	//for (int i = 0; i < nLevels; ++i)
	//{
	//	double hgt = min_height + i*(max_height - min_height) / max(1, (nLevels - 1));
		//create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 1, 100, true, false);
	//	//create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 8821, 1, true, false);
	//	//create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 1, 9215, true, false);
	//}

	ossimTieGptSet* chkSet = new ossimTieGptSet;
	//nLevels = 10;
	//for (int i = 0; i < nLevels; ++i)
	//{
	//	double hgt = min_height + i*(max_height - min_height) / max(1, (nLevels - 1));
	//	create3DGridPoints(imageRect, *sensorModel, hgt, chkSet, 30, 30, true, false);
	//}
	handler->close();

	//mylib::saveGcpFile(virtual_gcpfile, gptSet, NULL, false);
	//mylib::saveGcpFile(virtual_chkfile, chkSet, NULL, false);
}

void main()
{

	const char* pszPreferenceFile = "";
	ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
	if (preferences_file.exists())
	{
		ossimPreferences::instance()->loadPreferences(preferences_file);
	}
	else
	{
		preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
		if (preferences_file.exists())
		{
			ossimPreferences::instance()->loadPreferences(preferences_file);
		}
	}
	ossimInit::instance()->initialize();

	//ZY3Test();
	Cbers04Test();
	//HJ1Test();
}

