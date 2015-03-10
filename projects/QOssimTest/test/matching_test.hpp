#pragma once
#include "base/fileUtil.h"
#include "base/strUtil.h"

#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimObjectFactoryRegistry.h>
#include <ossim/base/ossimObjectFactoryRegistry.h>
#include <ossim/imaging/ossimImageChain.h>
#include <ossim/base/ossimPreferences.h>
#include <ossim/init/ossimInit.h>
#include <ossim/base/ossimStdOutProgress.h>
#include <ossim/elevation/ossimElevManager.h>
#include <ossim/imaging/ossimImageHandler.h>
#include <ossim/imaging/ossimImageHandlerRegistry.h>
#include <ossim/base/ossimCommon.h>
#include <ossim/base/ossimRefPtr.h>
#include <ossim/base/ossimOutputSource.h>
#include <ossim/imaging/ossimImageSourceSequencer.h>
#include <ossim/base/ossimProcessInterface.h>
#include <ossim/imaging/ossimImageSourceFilter.h>
#include <ossim/imaging/ossimBandSelector.h>
#include <ossim/imaging/ossimImageFileWriter.h>
#include <ossim/imaging/ossimImageWriterFactoryRegistry.h>

#include "cv.h"
#include "highgui.h"
#include "opencv2/nonfree/nonfree.hpp"

#include "base\quick_selection.h"
#include "line_match/lineConstant.h"
#include "line_match/line_func.h"
#include "line_match/ecmlrAffine.h"

void openSurMatch(ossimFilename masterFile, ossimFilename slaverFile)
{
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimsurf_plugin.dll");
	ossimInit::instance()->initialize();
	ossimString ascii_file = "D:\\workspace\\test.txt";
	ossimRefPtr<ossimObject> icObject = ossimObjectFactoryRegistry::instance()->createObject(ossimString("ossimSurfMatch"));//自定义自动选取同名点的类，为dll加载形式。
	ossimOutputSource* icSource = PTR_CAST(ossimOutputSource, icObject.get());

	ossimImageSourceSequencer* theImageSourceSequencer = new ossimImageSourceSequencer;
	theImageSourceSequencer->setTileSize(ossimIpt(1024,1024));

	ossimStdOutProgress progress(0,true);
	icSource->addListener(&progress);
	ossimProcessInterface* icProcessInterface = PTR_CAST(ossimProcessInterface, icObject.get());
	ossimPropertyInterface* icPropertyInterface = PTR_CAST(ossimPropertyInterface, icObject.get());

	icPropertyInterface->setProperty("master_band", "0");
	icPropertyInterface->setProperty("slave_band", "0");
	icPropertyInterface->setProperty("hessian_threshold", "500");
	icPropertyInterface->setProperty("master_filename", masterFile);//参考影像
	icPropertyInterface->setProperty("slaver_filename", slaverFile);//待校正影像
	icPropertyInterface->setProperty("output_filename", ascii_file); 
	// if (!file_tmp.exists()) 

	bool result = icProcessInterface->execute();

	icObject=NULL;
}

void opencvPointMatch(ossimFilename strFile)
{
	ossimFilename strOutFile = "D:\\workspace\\test.tif";
	ossimString ascii_file = "D:\\workspace\\test.txt";

	ossimKeywordlist kwl; 

	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimopencv_plugin.dll");
	ossimInit::instance()->initialize();

	ossimImageHandler *handler = ossimImageHandlerRegistry::instance()->open(strFile);

	if(!handler)
	{
		cout << "Unable to open input image: "<< strFile << endl;
		return;
	}

	ossimRefPtr<ossimObject> opencvObj = ossimObjectFactoryRegistry::instance()->createObject(ossimString("ossimOpenCVSURFFeatures"));
	ossimImageSourceFilter* filter = PTR_CAST(ossimImageSourceFilter, opencvObj.get());
	ossimPropertyInterface* PropertyInterface = PTR_CAST(ossimPropertyInterface, opencvObj.get());

	PropertyInterface->setProperty("hessian_threshold", "50");

	//ossimOpenCVCannyFilter * filter = new ossimOpenCVCannyFilter();
	//ossimOpenCVSURFFeatures * filter = new ossimOpenCVSURFFeatures();
	//ossimOpenCVPyrSegmentation * filter = new ossimOpenCVPyrSegmentation();

	//if(kwl.addFile(strOutFile)){
	//	filter->loadState(kwl);
	//}

	vector<ossim_uint32> outBandList;
	outBandList.push_back(1);
	ossimBandSelector* theBandSelector = new ossimBandSelector;
	theBandSelector->connectMyInputTo(0, handler);
	theBandSelector->setOutputBandList(outBandList);

	filter->connectMyInputTo(0,theBandSelector);

	ossimImageFileWriter* writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
	writer->setFilename(strOutFile);
	writer->connectMyInputTo(filter);

	ossimImageSourceSequencer* theImageSourceSequencer = new ossimImageSourceSequencer;
	theImageSourceSequencer->setTileSize(ossimIpt(1024,1024));
	writer->changeSequencer(theImageSourceSequencer);

	ossimStdOutProgress progress(0, true);
	writer->addListener(&progress);

	//writer->setTileSize(ossimIpt(1024,1024));

	cv::initModule_nonfree();
	writer->execute();
	writer->removeListener(&progress);
}

void opencvTest(ossimFilename strFile, vector<cvline_polar> &vec_lines)
{
	//ossimFilename strFile = "E:\\testdata\\QB\\Src\\2\\005554608050_01_P001_MUL\\05JUL04180116-M2AS-005554608050_01_P001.TIL";
	//ossimFilename strFile = "D:\\workspace\\LD2010003816\\header.dat";
	//ossimFilename strFile = "D:\\workspace\\HP\\HP0014.jpg";
	ossimFilename strOutFile = "D:\\workspace\\test.tif";
	ossimString ascii_file = "D:\\workspace\\test.txt";

	ossimKeywordlist kwl; 

	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimopencv_plugin.dll");
	ossimInit::instance()->initialize();

	ossimImageHandler *handler = ossimImageHandlerRegistry::instance()->open(strFile);

	if(!handler)
	{
		cout << "Unable to open input image: "<< strFile << endl;
		return;
	}

	ossimRefPtr<ossimObject> opencvObj = ossimObjectFactoryRegistry::instance()->createObject(ossimString("ossimOpenCVLsdFilter"));
	ossimImageSourceFilter* filter = PTR_CAST(ossimImageSourceFilter, opencvObj.get());
	ossimPropertyInterface* PropertyInterface = PTR_CAST(ossimPropertyInterface, opencvObj.get());

	PropertyInterface->setProperty("scale", "0.2");
	PropertyInterface->setProperty("saveImage", "1");
	PropertyInterface->setProperty("outASCII", ascii_file);

	//ossimOpenCVCannyFilter * filter = new ossimOpenCVCannyFilter();
	//ossimOpenCVSURFFeatures * filter = new ossimOpenCVSURFFeatures();
	//ossimOpenCVPyrSegmentation * filter = new ossimOpenCVPyrSegmentation();

	//if(kwl.addFile(strOutFile)){
	//	filter->loadState(kwl);
	//}

	vector<ossim_uint32> outBandList;
	outBandList.push_back(1);
	ossimBandSelector* theBandSelector = new ossimBandSelector;
	theBandSelector->connectMyInputTo(0, handler);
	theBandSelector->setOutputBandList(outBandList);

	filter->connectMyInputTo(0,theBandSelector);

	ossimImageFileWriter* writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
	writer->setFilename(strOutFile);
	writer->connectMyInputTo(filter);

	ossimImageSourceSequencer* theImageSourceSequencer = new ossimImageSourceSequencer;
	theImageSourceSequencer->setTileSize(ossimIpt(1024,1024));
	writer->changeSequencer(theImageSourceSequencer);

	ossimStdOutProgress progress(0, true);
	writer->addListener(&progress);

	//writer->setTileSize(ossimIpt(1024,1024));

	cv::initModule_nonfree();
	writer->execute();
	writer->removeListener(&progress);


	struct line_data 
	{
		vector<double> lsd_data;
		cvline_polar polar_data;
	};
	std::vector< line_data > lines_all;
	FILE *pf = fopen(ascii_file.c_str(), "r+");
	vector<double> tmp(7);
	while(EOF != fscanf(pf, "%lf%lf%lf%lf%lf%lf%lf\n", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6]))
	{
		cvline_polar line_polar = line2polar(fPoint(tmp[0], tmp[1]), fPoint(tmp[2], tmp[3]));
		line_data ld;
		ld.lsd_data = tmp;
		ld.polar_data = line_polar;
		lines_all.push_back(ld);
	}
	fclose(pf);


	std::vector< line_data > lines_merged;
	// merge the colinear straight lines
	double merge_threshold = 3.0;
	for (int i = 0;i < (int)lines_all.size();++i)
	{
		bool is_exist = false;
		for (int j = 0;j < (int)lines_merged.size();++j)
		{
			fPoint dis = segment2segment(lines_all[i].polar_data, lines_merged[j].polar_data);
			if((dis.x*dis.x + dis.y*dis.y) < merge_threshold)
			{
				// merge
				cvline_polar new_line = mergeline(lines_merged[j].polar_data, lines_all[i].polar_data);
				vector<double> tmp(7);
				tmp[0] = new_line.pt1.x;
				tmp[1] = new_line.pt1.y;
				tmp[2] = new_line.pt2.x;
				tmp[3] = new_line.pt2.y;
				// width
				tmp[4] = std::sqrt((tmp[0] - tmp[2])*(tmp[0] - tmp[2])+(tmp[1] - tmp[3])*(tmp[1] - tmp[3]));
				// p
				tmp[5] = std::min(lines_merged[j].lsd_data[5], lines_all[i].lsd_data[5]);
				// -log10(NFA)
				tmp[6] = std::max(lines_merged[j].lsd_data[6], lines_all[i].lsd_data[6]);

				line_data ld;
				ld.lsd_data = tmp;
				ld.polar_data = new_line;
				lines_merged[j] = ld;
				is_exist = true;
				break;
			}
		}
		if(!is_exist)
		{
			lines_merged.push_back(lines_all[i]);
		}
	}

	double minLen, minNFA;
	int seletion_num = 500;
	vector<double> nfas;
	minLen = 5.0;
	int n = (int)lines_merged.size();
	for(int i=0;i<n;i++)
	{
		cvline_polar outLine;
		fPoint line[2];
		line[0].x = lines_merged[i].lsd_data[0];
		line[0].y = lines_merged[i].lsd_data[1];
		line[1].x = lines_merged[i].lsd_data[2];
		line[1].y = lines_merged[i].lsd_data[3];

		outLine = line2polar(line);
		double len = segment_length(outLine);

		if(len >= minLen) nfas.push_back(-lines_merged[i].lsd_data[6]);
	}

	minNFA = -quick_select(&nfas[0], 0, (int)nfas.size()-1, seletion_num);
	for(int i=0;i<(int)nfas.size();i++)
	{
		{
			if(lines_merged[i].lsd_data[6] >= minNFA)
			{
				vec_lines.push_back(lines_merged[i].polar_data);
			}
		}
		//printf("\n");
	}

	int x = handler->getImageRectangle().width();
	int y = handler->getImageRectangle().height();
	IplImage* black= cvCreateImage(cvSize(x,y), IPL_DEPTH_8U, 1);
	cvZero(black);
	for(int i=0;i<(int)vec_lines.size();++i)
	{
		//if(vec[i].x1==vec[i].x2||vec[i].y1==vec[i].y2)
		cvLine(black,cvPoint(vec_lines[i].pt1.x,vec_lines[i].pt1.y),cvPoint(vec_lines[i].pt2.x,vec_lines[i].pt2.y),CV_RGB(255,255,255),1, CV_AA);
	}
	cvSaveImage(strOutFile.c_str(), black);
}

void line_match()
{
	//ossimFilename strSource = "D:\\workspace\\LD2010003816\\source.TIF";
	ossimFilename strSource = "D:\\workspace\\source.TIF";
	//ossimFilename strSource = "E:\\Dropbox\\Programs\\data\\sitongqiao\\2.jpg";
	//ossimFilename strSource = "D:\\workspace\\HP\\HP0014.jpg";
	vector<cvline_polar> lines_Source;
	opencvTest(strSource, lines_Source);

	//ossimFilename strRefer = "D:\\workspace\\LD2010003816\\TM_123036.TIF";
	//ossimFilename strRefer = "D:\\workspace\\ref.TIF";
	ossimFilename strRefer = "E:\\Dropbox\\Programs\\data\\sitongqiao\\1.jpg";
	vector<cvline_polar> lines_Refer;
	opencvTest(strRefer, lines_Refer);

	ecmlrAffine* ecmlr = new ecmlrAffine(lines_Source, lines_Refer);
	ecmlr->m_modelImageFile = strSource.c_str();
	ecmlr->m_observeImageFile = strRefer.c_str();
	//ecmlr->setDrawResult(true);
	vector<double> modelParameters;
	ecmlr->solve(modelParameters);

	//affine_warp_Image< vector<double> >(strSource.c_str(), strRefer.c_str(), modelParameters);


	for (int i = 0;i < (int)modelParameters.size();++i)
	{
		cout<<modelParameters[i]<<endl;
	}
	cout<<"finished!"<<endl;
}