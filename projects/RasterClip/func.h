#pragma once

/////////////// ossim
#include <ossim/base/ossimConstants.h>
#include <ossim/base/ossimKeywordNames.h>
#include <ossim/base/ossimTrace.h>
#include <ossim/base/ossimStdOutProgress.h>
#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimKeywordlist.h>
#include <ossim/base/ossimDrect.h>
#include <ossim/base/ossimImageTypeLut.h>
#include <ossim/imaging/ossimJpegWriter.h>
#include <ossim/imaging/ossimImageHandler.h>
#include <ossim/imaging/ossimRLevelFilter.h>
#include <ossim/imaging/ossimImageSource.h>
#include <ossim/imaging/ossimImageHandlerRegistry.h>
#include <ossim/imaging/ossimImageWriterFactoryRegistry.h>
#include <ossim/imaging/ossimImageWriterFactory.h>
#include <ossim/imaging/ossimImageFileWriter.h>
#include <ossim/imaging/ossimCacheTileSource.h>
#include <ossim/imaging/ossimBandSelector.h>
#include <ossim/init/ossimInit.h>
#include <ossim/imaging/ossimImageSourceFilter.h>
#include <ossim/base/ossimPolygon.h>
#include <ossim/imaging/ossimImageDataHelper.h>
#include <ossim/imaging/ossimPolyCutter.h>
#include <ossim/projection/ossimProjection.h>
#include <ossim/projection/ossimMapProjection.h>
#include <ossim/projection/ossimProjectionFactoryRegistry.h>
#include <ossim/projection/ossimTransMercatorProjection.h>
#include <ossim/base/ossimDirectory.h>
#include <ossim/plugin/ossimSharedPluginRegistry.h>
#include <ossim/projection/ossimMapProjectionFactory.h>
#include <ossim/projection/ossimProjectionFactoryRegistry.h>
#include <ossim/imaging/ossimImageRenderer.h>
#include <ossim/imaging/ossimBandSelector.h>
#include <ossim/imaging/ossimMaskFilter.h>
#include <ossim/imaging/ossimSingleImageChain.h>
#include <ossim/base/ossimStringProperty.h>
#include <ossim/base/ossimTimer.h>
//#include <ossim/imaging/ossimImageGeometry.h>
#include <ossim/projection/ossimImageViewProjectionTransform.h>
#include <ossim/projection/ossimImageViewTransform.h>
#include <ossim/ossimConfig.h>
#include <ossim/base/ossimCommon.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>
#include <ossim/base/ossimFilename.h>
#include <ossim/init/ossimInit.h>
#include <ossim/base/ossimNotify.h>
#include <ossim/base/ossimKeywordNames.h>
#include <ossim/base/ossimStdOutProgress.h>
#include <ossim/imaging/ossimBitMaskWriter.h>
#include <ossim/imaging/ossimImageHandlerRegistry.h>
#include <ossim/imaging/ossimImageSourceSequencer.h>
#include <ossim/imaging/ossimPixelFlipper.h> // for its keywords
#include <ossim/imaging/ossimBitMaskTileSource.h>
#include <ossim/imaging/ossimMaskFilter.h>

/////////////// gdal
#include <ogrsf_frmts.h>
#include <gdal.h>
#include <ogr_api.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>

///////////// stl
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <exception>
#include <math.h>
#include <vector>

using namespace std;

struct Point2f 
{
	double x;
	double y;
	Point2f(const double& ox, const double& oy)
	{
		x = ox;
		y = oy;
	}
	Point2f()
	{
		x = 0.0;
		y = 0.0;
	}
};


ossimGpt linesample2Wolrd(const Point2f linesample, double *adfThisGeoTransform, OGRCoordinateTransformation *poCT);

vector<ossimGpt> getRasterBoundary(const char* szFilename);

vector<ossim_uint32> getOutputBandList(ossimImageHandler* handler, string strList="");

ossimDrect calcRegion(vector<ossimDpt> clipRect, ossimMapProjection* inmapinfo, ossimImageViewTransform*  IVT);

/////////////////////////////////////////////////
// 
//////////////////////////////////////////////////

//// 根据任意区域裁剪
bool cut_by_region(const char* inputImageFile, vector<ossimGpt> boundaryList, const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList = vector<ossim_uint32>(), int offset = 0);

//// 根据矩形区域进行裁剪
bool cut_by_rect(const char* inputImageFile, vector<ossimGpt> boundaryList, const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList = vector<ossim_uint32>(), int offset = 0);

/// 根据shapefile矢量边界裁剪
bool cut_by_shp(const char* inputImageFile, const char* shpfileName, const char* fieldName, const char* searchValue,
	const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList = vector<ossim_uint32>(), int offset = 0);

/// 根据shapefile矢量边界裁剪，输出矩形
bool cut_by_shp_rect(const char* inputImageFile, const char* shpfileName, const char* fieldName, const char* searchValue,
	const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList = vector<ossim_uint32>(), int offset = 0);

/// 根据shapefile矢量边界蒙板，不改变大小
bool mask_by_shp(const char* inputImageFile, const char* shpfileName, const char* fieldName, const char* searchValue,
	const char* strPrjFile, const char* strOutFile);
//
//// 强制转为方形区域
//void toSquare(double &x1, double &y1, double &x2, double &y2);
//
//// 自定义输出矩形比例,ratio是规定矩形的长宽比
//void toCustomRect(double &x1, double &y1, double &x2, double &y2, double ratio = 1.0);
//
//void custum_cutbyshpfile(ossimString strProvince, ossimString strPrj, ossimString strOut);
//
//void square_cutbyshpfile(ossimString strProvince, ossimString strPrj, ossimString strOut);
//
//void cutbyshpfile_band(ossimString strProvince, ossimString strPrj, ossimString strOut, int band);
//
//
//void t0()
//{
//	ossimInit::instance()->initialize();
//	OGRRegisterAll();
//	GDALAllRegister();
//
//	std::vector<ossimDrect> mapDrect;
//	ossimFilename shpfilename, imgfilePath;
//	ossimFilename outfilename;
//
//	//shpfilename="C:\\data\\2005_administrative_divisions_result\\sichuan\\sichuan_sp2.shp";
//	//outfilename="C:\\data\\2005_administrative_divisions_result\\sichuan\\";
//	//shpfilename="C:\\xiaolu\\china_100wan.shp";
//	//outfilename="C:\\data\\2005_administrative_divisions_result\\test\\";
//
//	//ossimString strProvince = "香港特别行政区";
//	//ossimString strPrj = "guangdonglambertcc.txt";
//	//ossimString strOut = "xianggang";
//	//ossimString strProvince = "北京市";
//	ossimString strTask = "Z:\\loong\\province\\taskFile.txt";
//	vector<ossimString> strProvinceList;
//	vector<ossimString> strPrjList;
//	vector<ossimString> strOutList;
//	fstream fs;
//	fs.open(strTask, std::ios_base::in);
//
//	ossimString strProvince;
//	ossimString strPrj;
//	ossimString strOut;
//	while (fs >> strProvince >> strPrj >> strOut)
//	{
//		strProvinceList.push_back(strProvince);
//		strPrjList.push_back(strPrj);
//		strOutList.push_back(strOut);
//	}
//	fs.close();
//
//	int num = strProvinceList.size();
//	for (int i = 0; i < num; i++)
//	{
//		cout << "正在完成" << num << "个中的第" << i + 1 << "个：" << strProvinceList[i] << "...\n";
//
//		ossimFilename outFilename = outPath + strOutList[i] + ".TIF";	//输出文件名
//		if (!outFilename.exists())
//		{
//			//square_cutbyshpfile(strProvinceList[i], strPrjList[i], strOutList[i]);
//			custum_cutbyshpfile(strProvinceList[i], strPrjList[i], strOutList[i]);
//		}
//		/*for(int j = 0;j < 7;j++)
//		{
//		cutbyshpfile_band(strProvinceList[i], strPrjList[i], strOutList[i], j + 1);
//		}*/
//		cout << endl;
//	}
//}
//
//void t1()
//{
//	ossimInit::instance()->initialize();
//	GDALAllRegister();
//
//	ossimString strTask = "S:\\tflong\\province\\8.tsk";
//
//	vector<ossimString> strProvinceList;
//	vector<ossimString> strPrjList;
//	vector<ossimDrect> rectList;
//
//	fstream fs;
//	fs.open(strTask, std::ios_base::in);
//
//	ossimString strId;
//	ossimString strProvince;
//	ossimString strPrj;
//	double minLong, minLat, maxLong, maxLat;
//	while (fs >> strId >> strProvince >> strPrj >> minLong >> minLat >> maxLong >> maxLat)
//	{
//		strProvinceList.push_back(strProvince);
//		strPrjList.push_back(strPrj);
//		rectList.push_back(ossimDrect(minLong, minLat, maxLong, maxLat));
//	}
//	fs.close();
//
//	//ossimDrect rect(115.5,29.0,121.0,22.5);
//	//cutbyRegion(rect,"fujianlambertcc.txt","福建");
//
//	ossimString strAdd = "";
//	int num = strProvinceList.size();
//	for (int i = 0; i < num; i++)
//	{
//		cout << "正在完成" << num << "个中的第" << i + 1 << "个：" << strProvinceList[i] << "...\n";
//
//		ossimFilename outFilename = outPath + strProvinceList[i] + strAdd + ".tif";	//输出文件名
//		if (!outFilename.exists())
//		{
//			cutbyRegion(rectList[i], strPrjList[i], outFilename);
//		}
//		cout << endl;
//	}
//}