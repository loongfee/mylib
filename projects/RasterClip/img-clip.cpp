#include <stdlib.h>
#include <math.h>


/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"
#include <ossim/base/ossimPreferences.h>
#include "gdalwarp.h"

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
#include <ossim/parallel/ossimMultiThreadSequencer.h>
#include <ossim/parallel/ossimMtDebug.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>

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

#include <strUtil.h>

#include "arguments.h"
#include <gflags/gflags.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#pragma comment(lib, "gflags.lib")
#pragma comment(lib, "gflags_nothreads.lib")
#pragma comment(lib, "ws2_32.lib")
namespace fs = boost::filesystem;
namespace po = boost::program_options;

using namespace std;

#ifndef _WIN64
#pragma comment(lib, "ossim20.lib")
#else
#pragma comment(lib, "ossim20x64.lib")
#endif
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "ossim_plugin.lib")
#pragma comment(lib, "OpenThreads.lib")


#include "global.h"

const char *pszPreferenceFile = "preference.txt";

static string strSrcfile = "china.vrt";
static string strDstfile = "cliped.tif";
static const char* shpfileName = "bnd1.shp";
static const char* strRefImageFile = "";
static const char* fieldName = "";
static const char* searchValue = "";
static const char* strPrjFile = "";
static char *pszScaleType = "DEFAULT";
//static char *pszSampleType = "BICUBIC";
static char *pszOutFileType = "TIFF";
static int offset = 10;
static bool bBoundary = false;
static bool bShapefile = false;
static bool bRefImage = false;
static bool bRect = false;
static vector<ossimGpt> boundaryList;
static vector<ossim_uint32> bandList;
static vector<ossim_float64> nullValueList;
static int nThreads = 0;
static string resampleType = "cubic";

const char* pszApplicationName = "";

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


ossimGpt linesample2Wolrd(const Point2f linesample, double *adfThisGeoTransform, OGRCoordinateTransformation *poCT)
{
	Point2f latlng;
	double P = linesample.x;
	double L = linesample.y;
	latlng.x = adfThisGeoTransform[0] + P*adfThisGeoTransform[1] + L*adfThisGeoTransform[2];
	latlng.y = adfThisGeoTransform[3] + P*adfThisGeoTransform[4] + L*adfThisGeoTransform[5];

	if( poCT == NULL || !poCT->Transform( 1, &latlng.x, &latlng.y ) )
	{
		printf( "Transformation failed.\n" );
	}
	return ossimGpt(latlng.y, latlng.x);
}

void setResampleType(ossimImageRenderer* renderer)
{
	if (_stricmp(resampleType.c_str(), "nearest") == 0)
	{
		renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_NEAREST_NEIGHBOR);
	}
	else if (_stricmp(resampleType.c_str(), "bilinear") == 0)
	{
		renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	}
	else if (_stricmp(resampleType.c_str(), "cubic") == 0)
	{
		renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
	}
	else
	{
		renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
	}
}

vector<ossimGpt> getRasterBoundary(const char* szFilename)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDatasetH hSrcDs = GDALOpen(szFilename, GA_ReadOnly );
	if (NULL == hSrcDs)
	{
		cout<<"Failed to open file: \""<<szFilename<<"\"."<<endl;
		exit(0);
	}
	GDALDataset  *hSrcDataset = (GDALDataset  *)hSrcDs;
	int nWidth = hSrcDataset->GetRasterXSize();
	int nHeight = hSrcDataset->GetRasterYSize();

	Point2f ul(0, 0);
	Point2f ur(nWidth, 0);
	Point2f lr(nWidth, nHeight);
	Point2f ll(0, nHeight);


	double adfThisGeoTransform[6];
	GDALGetGeoTransform( hSrcDs, adfThisGeoTransform );
	char* pszProjection = (char *)GDALGetProjectionRef( hSrcDs );


	OGRSpatialReference geoSRS;
	OGRSpatialReference prjSRS;
	prjSRS.importFromWkt(&pszProjection);
	geoSRS.SetWellKnownGeogCS( "WGS84" );

	OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation( &prjSRS,
		&geoSRS );

	vector<ossimGpt> boundary(4);

	boundary[0] = linesample2Wolrd(ul, adfThisGeoTransform, poCT);
	boundary[1] = linesample2Wolrd(ur, adfThisGeoTransform, poCT);
	boundary[2] = linesample2Wolrd(lr, adfThisGeoTransform, poCT);
	boundary[3] = linesample2Wolrd(ll, adfThisGeoTransform, poCT);

	GDALClose( hSrcDataset );
	return boundary;
}

vector<ossim_uint32> getOutputBandList(ossimImageHandler* handler, string strList/* = ""*/)
{
	vector<ossim_uint32> outList;
	if (!handler)
	{
		return outList;
	}

	size_t nbands = handler->getNumberOfInputBands();
	if (strList.empty())
	{
		// 默认输出所有波段
		for (size_t i = 0; i < nbands; i++)
		{
			outList.push_back(i);
		}
	}
	else
	{
		size_t nList = strList.size();
		for (size_t i = 0; i < nList; i++)
		{
			size_t iBand = strList[i] - '0';
			if (iBand >= 1 && iBand <= nbands)
			{
				// 如果波段号合法则加入
				outList.push_back(iBand - 1);
			}
		}
	}

	return outList;
}


void toSquare(double &x1, double &y1, double &x2, double &y2)
{
	double xSize = fabs(x2 - x1);
	double ySize = fabs(y2 - y1);
	double offset = fabs(xSize - ySize) * 0.5;
	if (xSize >= ySize)
	{
		if (y1 < y2)
		{
			y1 -= offset;
			y2 += offset;
		}
		else
		{
			y1 += offset;
			y2 -= offset;
		}
	}
	else
	{
		if (x1 < x2)
		{
			x1 -= offset;
			x2 += offset;
		}
		else
		{
			x1 += offset;
			x2 -= offset;
		}
	}
}

//自定义输出矩形比例,ratio是规定矩形的长宽比
void toCustomRect(double &x1, double &y1, double &x2, double &y2, double ratio = 1.0)
{
	double xSize = fabs(x2 - x1);
	double ySize = fabs(y2 - y1);
	double offset = 0.0;
	if (ySize < xSize * ratio)
	{
		offset = (xSize * ratio - ySize) * 0.5;
		if (y1 < y2)
		{
			y1 -= offset;
			y2 += offset;
		}
		else
		{
			y1 += offset;
			y2 -= offset;
		}
	}
	else if (ySize > xSize * ratio)
	{
		offset = (ySize / ratio - xSize) * 0.5;
		if (x1 < x2)
		{
			x1 -= offset;
			x2 += offset;
		}
		else
		{
			x1 += offset;
			x2 -= offset;
		}
	}
}

ossimDrect calcRegion(vector<ossimDpt> clipRect, ossimMapProjection* inmapinfo, ossimImageViewTransform*  IVT)
{
	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, IVT);


	ossimDpt dpt = IVT->imageToView(clipRect[0]);
	double xmin = dpt.x;
	double xmax = dpt.x;
	double ymin = dpt.y;
	double ymax = dpt.y;
	for (int i = 0; i < (int)clipRect.size(); i++)
	{
		ossimDpt dpt = IVT->imageToView(clipRect[i]);
		xmin = (dpt.x < xmin) ? dpt.x : xmin;
		ymin = (dpt.y < ymin) ? dpt.y : ymin;
		xmax = (dpt.x > xmax) ? dpt.x : xmax;
		ymax = (dpt.y > ymax) ? dpt.y : ymax;
	}

	ossimDrect viewRegion(xmin, ymin, xmax, ymax);
	return viewRegion;
}

bool mask_by_shp(const char* inputImageFile, const char* shpfileName, const char* fieldName, const char* searchValue,
				 const char* strPrjFile, const char* strOutFile)
{	
	int argc = 0;
	char* argv[20];
	argv[argc++] = "gdalwarp";

	argv[argc++] = "-cutline";
	argv[argc++] = const_cast<char*>(shpfileName);
	argv[argc++] = "-cwhere";
	argv[argc] = new char[1024];
	sprintf(argv[argc++], "%s=\'%s\'", fieldName, searchValue);
	argv[argc++] = "-multi";
	argv[argc++] = const_cast<char*>(inputImageFile);
	argv[argc++] = const_cast<char*>(strOutFile);

	gdalwarp::gdalwarp( argc, argv);
	return true;
}


///////////////////////////////////////////
// 根据矩形区域进行裁剪
/////////////////////////////////////////////

bool cut_by_rect(const char* inputImageFile, vector<ossimGpt> boundaryList, const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList/* = vector<ossim_uint32>()*/, int offset/* = 0*/)
{
	ossimImageHandler* handler;
	handler = ossimImageHandlerRegistry::instance()->open(ossimFilename(inputImageFile));
	if (!handler)
	{
		cout << "图像文件" << inputImageFile << "打开失败！" << endl;
		system("Pause");
		return false;
	}

	ossimRefPtr<ossimProjection> proj;
	// 输入投影
	ossimKeywordlist in_geom_kwl;
	handler->getImageGeometry()->saveState(in_geom_kwl);
	proj = handler->getImageGeometry()->getProjection();
	ossimMapProjection* inmapinfo = PTR_CAST(ossimMapProjection, proj.get());

	// 输出投影
	ossimKeywordlist out_geom_kwl;
	if (ossimFilename(strPrjFile).exists())
	{
		out_geom_kwl.addFile(strPrjFile);
	}
	else
	{
		handler->getImageGeometry()->getProjection()->saveState(out_geom_kwl);
	}
	proj = ossimProjectionFactoryRegistry::instance()->createProjection(out_geom_kwl);
	ossimMapProjection* outmapinfo = PTR_CAST(ossimMapProjection, proj.get());
	ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
	imageGeom->setProjection(outmapinfo);

	vector<ossimDpt>  polygon;
	polygon.clear();

	for (size_t i = 0; i < boundaryList.size(); i++)
	{
		ossimDpt pt = inmapinfo->worldToLineSample(boundaryList[i]);
		polygon.push_back(pt);
	}

	//统计外接矩形
	double xmin = polygon[0].x;
	double xmax = polygon[0].x;
	double ymin = polygon[0].y;
	double ymax = polygon[0].y;

	for (int i = 0; i < (int)polygon.size(); i++)
	{
		xmin = (polygon[i].x < xmin) ? polygon[i].x : xmin;
		ymin = (polygon[i].y < ymin) ? polygon[i].y : ymin;
		xmax = (polygon[i].x > xmax) ? polygon[i].x : xmax;
		ymax = (polygon[i].y > ymax) ? polygon[i].y : ymax;
	}


	//ossimImageViewProjectionTransform IVP(imageGeom.get(), handler->getImageGeometry().get());
	ossimImageViewProjectionTransform IVP(handler->getImageGeometry().get(), imageGeom.get());
	//ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, renderer->getImageViewTransform());

	//统计外接矩形
	ossimDpt dpt;
	IVP.imageToView(polygon[0], dpt);
	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
	for (int i = 1; i < (int)polygon.size(); i++)
	{
		IVP.imageToView(polygon[i], dpt);
		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
	}

	IVP.viewToImage(UL1, dpt);
	ossimDpt UL(dpt.x - offset, dpt.y - offset);
	IVP.viewToImage(LR1, dpt);
	ossimDpt LR(dpt.x + offset, dpt.y + offset);
	IVP.imageToView(UL, UL1);
	IVP.imageToView(LR, LR1);
	//ossimDrect viewRegion = calcRegion(clipRect, inmapinfo, renderer->getImageViewTransform());
	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
	ossimDpt UR;
	//ossimDpt UR = renderer->getImageViewTransform()->viewToImage(ossimDpt(LR1.lon, UL1.lat));
	IVP.viewToImage(ossimDpt(LR1.lon, UL1.lat), UR);
	ossimDpt LL;
	//ossimDpt LL = renderer->getImageViewTransform()->viewToImage(ossimDpt(UL1.lon, LR1.lat));
	IVP.viewToImage(ossimDpt(UL1.lon, LR1.lat), LL);

	vector<ossimDpt> clipRect;
	clipRect.clear();
	clipRect.push_back(UL);
	clipRect.push_back(UR);
	clipRect.push_back(LR);
	clipRect.push_back(LL);
	clipRect.push_back(UL);

	//handler->setImageGeometry(imageGeom.get());//1128
	//handler->loadState(out_geom_kwl);

	//选择输出波段 
	ossimBandSelector *theBandSelector = new ossimBandSelector;
	if (bandList.size() < 1)
	{
		size_t nbands = handler->getNumberOfInputBands();
		for (size_t i = 0; i < nbands; i++)
		{
			bandList.push_back(i);
		}
	}
	// 设置空值
	if (nullValueList.size() > 0 && nullValueList.size() == bandList.size())
	{
		for (size_t i = 0; i < bandList.size(); i++)
		{
			handler->setNullPixelValue(bandList[i], nullValueList[i]);
		}
	}

	theBandSelector->connectMyInputTo(0, handler);
	theBandSelector->setOutputBandList(bandList);

	ossimPolyCutter* theCutter = new ossimPolyCutter;

	theCutter->setNumberOfPolygons(0);
	theCutter->addPolygon(clipRect);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);
	theCutter->connectMyInputTo(theBandSelector);

	ossimImageRenderer* renderer = new ossimImageRenderer;
	setResampleType(renderer);
	//renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	//renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
	renderer->connectMyInputTo(theCutter);
	renderer->setView(outmapinfo);

	//ossimImageFileWriter* writer =
	//	ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_gtiff"));
	ossimImageFileWriter* writer;
	if (0 == strcmp(ossimFilename(strOutFile).ext().upcase(), "TIF"))
	{
		//writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_gtiff"));
		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));

	}
	else
	{
		writer = ossimImageWriterFactoryRegistry::instance()->
			createWriterFromExtension(ossimFilename(strOutFile).ext());
	}

	writer->setFilename(ossimFilename(strOutFile));

	ossimStdOutProgress progress(0, true);
	writer->addListener(&progress);

	if (nThreads < 1)
	{
		nThreads = OpenThreads::GetNumberOfProcessors() * 2;
	}
	cout << "using " << nThreads << " threads:" << endl;
	//ossimMtDebug::instance()->seqDebugEnabled = true;
	//ossimMtDebug::instance()->seqTimedBlocksDt = 50000;
	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
	writer->changeSequencer(sequencer.get());
	writer->connectMyInputTo(0, renderer);

	//writer->initialize();
	writer->setAreaOfInterest(viewRegion);
	writer->execute();

	//delete writer;
	//delete theCutter;
	//delete renderer;
	//polygon.clear();
	handler->close();

	//########## DEBUG CODE FOR TIMING MULTI THREAD LOCKS ##############
	if (sequencer.valid() && (nThreads != 0))
	{
		ossimMultiThreadSequencer* mts = dynamic_cast<ossimMultiThreadSequencer*>(sequencer.get());
		if (mts != NULL)
		{

			double jgtt = mts->d_jobGetTileT;
			ossim_uint32 num_threads = mts->getNumberOfThreads();
			double jgttpj = jgtt / num_threads;
			cout << setprecision(3) << endl;
			cout << "Multi-threading metrics ---" << endl;
			cout << "   Number of threads:      " << num_threads << endl;
			cout << "   Max cache used:         " << mts->d_maxCacheUsed << endl;
			cout << "   Cache emptied count:    " << ossimString::toString(mts->d_cacheEmptyCount) << endl;
			cout << "   Time waiting on jobs:   " << mts->d_idleTime2 << " s" << endl;
			cout << "   Time waiting on cache:  " << mts->d_idleTime5 << " s" << endl;
			cout << "   Handler getTile T:      " << mts->handlerGetTileT() << " s" << endl;
			cout << "   Job getTile T:          " << jgtt << " s" << endl;
			cout << "   Average getTile T/job:  " << jgttpj << " s\n" << endl;
		}
	}

	return true;
}

//bool cut_by_rect(const char* inputImageFile, vector<ossimGpt> boundaryList, const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList/* = vector<ossim_uint32>()*/, int offset/* = 0*/)
//{
//	ossimImageHandler* handler;
//	handler = ossimImageHandlerRegistry::instance()->open(ossimFilename(inputImageFile));
//	if (!handler)
//	{
//		cout << "图像文件" << inputImageFile << "打开失败！" << endl;
//		system("Pause");
//		return false;
//	}
//
//	ossimRefPtr<ossimProjection> proj;
//	// 输入投影
//	ossimKeywordlist in_geom_kwl;
//	handler->getImageGeometry()->saveState(in_geom_kwl);
//	proj = handler->getImageGeometry()->getProjection();
//	ossimMapProjection* inmapinfo = PTR_CAST(ossimMapProjection, proj.get());
//
//	// 输出投影
//	ossimKeywordlist out_geom_kwl;
//	if (ossimFilename(strPrjFile).exists())
//	{
//		out_geom_kwl.addFile(strPrjFile);
//	}
//	else
//	{
//		handler->getImageGeometry()->getProjection()->saveState(out_geom_kwl);
//	}
//	proj = ossimProjectionFactoryRegistry::instance()->createProjection(out_geom_kwl);
//	ossimMapProjection* outmapinfo = PTR_CAST(ossimMapProjection, proj.get());
//
//
//	vector<ossimDpt>  polygon;
//	polygon.clear();
//
//	for (size_t i = 0; i < boundaryList.size(); i++)
//	{
//		ossimDpt pt = inmapinfo->worldToLineSample(boundaryList[i]);
//		polygon.push_back(pt);
//	}
//
//	//统计外接矩形
//	double xmin = polygon[0].x;
//	double xmax = polygon[0].x;
//	double ymin = polygon[0].y;
//	double ymax = polygon[0].y;
//
//	for (int i = 0; i < (int) polygon.size(); i++)
//	{
//		xmin = (polygon[i].x < xmin) ? polygon[i].x : xmin;
//		ymin = (polygon[i].y < ymin) ? polygon[i].y : ymin;
//		xmax = (polygon[i].x > xmax) ? polygon[i].x : xmax;
//		ymax = (polygon[i].y > ymax) ? polygon[i].y : ymax;
//	}
//
//	ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
//	imageGeom->setProjection(outmapinfo);
//	//handler->setImageGeometry(imageGeom.get());//1128
//	//handler->loadState(out_geom_kwl);
//
//	//选择输出波段 
//	ossimBandSelector *theBandSelector = new ossimBandSelector;
//	if (bandList.size() < 1)
//	{
//		size_t nbands = handler->getNumberOfInputBands();
//		for (size_t i = 0; i < nbands; i++)
//		{
//			bandList.push_back(i);
//		}
//	}
//	theBandSelector->connectMyInputTo(0, handler);
//	theBandSelector->setOutputBandList(bandList);
//
//	ossimPolyCutter* theCutter = new ossimPolyCutter;
//
//	vector<ossimDpt> clipRect;
//	//clipRect.push_back(ossimDpt(xmin, ymin));
//	//clipRect.push_back(ossimDpt(xmax, ymin));
//	//clipRect.push_back(ossimDpt(xmax, ymax));
//	//clipRect.push_back(ossimDpt(xmin, ymin));
//	//clipRect.push_back(ossimDpt(xmin, ymin));
//
//
//	//theCutter->connectMyInputTo(theBandSelector);
//
//	//theCutter->setNumberOfPolygons(0);
//	//theCutter->addPolygon(clipRect);
//	//theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
//	//theCutter->setNumberOfPolygons(1);
//
//	//ossimIrect boundw1 = theCutter->getBoundingRect();
//	//ossimIrect boundw = handler->getBoundingRect();
//
//
//	ossimImageRenderer* renderer = new ossimImageRenderer;
//
//	//renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
//	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
//
//	renderer->setView(outmapinfo);
//	renderer->connectMyInputTo(theBandSelector);
//	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, renderer->getImageViewTransform());
//
//	//统计外接矩形
//	ossimDpt dpt = renderer->getImageViewTransform()->imageToView(polygon[0]);
//	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
//	for (int i = 1; i < (int) polygon.size(); i++)
//	{
//		dpt = renderer->getImageViewTransform()->imageToView(polygon[i]);
//		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
//		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
//		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
//		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
//	}
//
//	dpt = renderer->getImageViewTransform()->viewToImage(UL1);
//	ossimDpt UL(dpt.x - offset, dpt.y - offset);
//	dpt = renderer->getImageViewTransform()->viewToImage(LR1);
//	ossimDpt LR(dpt.x + offset, dpt.y + offset);
//	UL1 = renderer->getImageViewTransform()->imageToView(UL);
//	LR1 = renderer->getImageViewTransform()->imageToView(LR);
//	//ossimDrect viewRegion = calcRegion(clipRect, inmapinfo, renderer->getImageViewTransform());
//	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
//	ossimDpt UR = renderer->getImageViewTransform()->viewToImage(ossimDpt(LR1.lon, UL1.lat));
//	ossimDpt LL = renderer->getImageViewTransform()->viewToImage(ossimDpt(UL1.lon, LR1.lat));
//
//	clipRect.clear();
//	clipRect.push_back(UL);
//	clipRect.push_back(UR);
//	clipRect.push_back(LR);
//	clipRect.push_back(LL);
//	clipRect.push_back(UL);
//
//	theCutter->setNumberOfPolygons(0);
//	theCutter->addPolygon(clipRect);
//	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
//	theCutter->setNumberOfPolygons(1);
//
//	theCutter->connectMyInputTo(theBandSelector);
//
//	renderer->setView(outmapinfo);
//	renderer->connectMyInputTo(theCutter);
//	
//	//ossimImageFileWriter* writer =
//	//	ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_gtiff"));
//	ossimImageFileWriter* writer;
//	if (0 == strcmp(ossimFilename(strOutFile).ext().upcase(), "TIF"))
//	{
//		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_gtiff"));
//
//	}
//	else
//	{
//		writer = ossimImageWriterFactoryRegistry::instance()->
//			createWriterFromExtension( ossimFilename(strOutFile).ext() );
//	}
//
//	writer->setFilename(ossimFilename(strOutFile));
//
//	ossimStdOutProgress progress(0, true);
//	writer->addListener(&progress);
//
//	if (nThreads < 1)
//	{
//		nThreads = OpenThreads::GetNumberOfProcessors() * 2;
//	}
//	cout << "using " << nThreads << " threads:" << endl;
//
//	//ossimMtDebug::instance()->seqDebugEnabled = true;
//	//ossimMtDebug::instance()->seqTimedBlocksDt = 50000;
//	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
//	writer->changeSequencer(sequencer.get());
//	writer->connectMyInputTo(0, renderer);
//
//	//writer->initialize();
//	writer->setAreaOfInterest(viewRegion);
//	writer->execute();
//
//	//delete writer;
//	//delete theCutter;
//	//delete renderer;
//	//polygon.clear();
//	handler->close();
//
//	//########## DEBUG CODE FOR TIMING MULTI THREAD LOCKS ##############
//	if (sequencer.valid() && (nThreads != 0))
//	{
//		ossimMultiThreadSequencer* mts = dynamic_cast<ossimMultiThreadSequencer*>(sequencer.get());
//		if (mts != NULL)
//		{
//
//			double jgtt = mts->d_jobGetTileT;
//			ossim_uint32 num_threads = mts->getNumberOfThreads();
//			double jgttpj = jgtt / num_threads;
//			cout << setprecision(3) << endl;
//			cout << "Multi-threading metrics ---" << endl;
//			cout << "   Number of threads:      " << num_threads << endl;
//			cout << "   Max cache used:         " << mts->d_maxCacheUsed << endl;
//			cout << "   Cache emptied count:    " << ossimString::toString(mts->d_cacheEmptyCount) << endl;
//			cout << "   Time waiting on jobs:   " << mts->d_idleTime2 << " s" << endl;
//			cout << "   Time waiting on cache:  " << mts->d_idleTime5 << " s" << endl;
//			cout << "   Handler getTile T:      " << mts->handlerGetTileT() << " s" << endl;
//			cout << "   Job getTile T:          " << jgtt << " s" << endl;
//			cout << "   Average getTile T/job:  " << jgttpj << " s\n" << endl;
//		}
//	}
//
//	return true;
//}


///////////////////////////////////////////
// 根据任意边界进行裁剪
///////////////////////////////////////////
bool cut_by_region(const char* inputImageFile, vector<ossimGpt> boundaryList, const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList/* = vector<ossim_uint32>()*/, int offset/* = 0*/)
{
	ossimImageHandler* handler;
	handler = ossimImageHandlerRegistry::instance()->open(ossimFilename(inputImageFile));
	if (!handler)
	{
		cout << "图像文件" << inputImageFile << "打开失败！" << endl;
		return false;
	}

	ossimRefPtr<ossimProjection> proj;
	// 输入投影
	ossimKeywordlist in_geom_kwl;
	handler->getImageGeometry()->saveState(in_geom_kwl);
	proj = handler->getImageGeometry()->getProjection();
	ossimMapProjection* inmapinfo = PTR_CAST(ossimMapProjection, proj.get());

	// 输出投影
	ossimKeywordlist out_geom_kwl;
	if (ossimFilename(strPrjFile).exists())
	{
		out_geom_kwl.addFile(strPrjFile);
	}
	else
	{
		handler->getImageGeometry()->getProjection()->saveState(out_geom_kwl);
	}
	proj = ossimProjectionFactoryRegistry::instance()->createProjection(out_geom_kwl);
	ossimMapProjection* outmapinfo = PTR_CAST(ossimMapProjection, proj.get());

	ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
	imageGeom->setProjection(outmapinfo);

	vector<ossimDpt>  polygon;
	polygon.clear();

	ossimDpt last_pt;
	for (size_t i = 0; i < boundaryList.size(); i++)
	{
		ossimDpt pt = inmapinfo->worldToLineSample(boundaryList[i]);
		ossimDpt delta_pt = last_pt - pt;
		if (i > 0 && (last_pt - pt).length() < 60.0)
		{
			// 忽略距离小于1的冗余点
			continue;
		}

		polygon.push_back(pt);
		last_pt = pt;
	}


	ossimImageViewProjectionTransform IVP(handler->getImageGeometry().get(), imageGeom.get());
	//统计外接矩形
	ossimDpt dpt;
	IVP.imageToView(polygon[0], dpt);
	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
	for (int i = 1; i < (int)polygon.size(); i++)
	{
		IVP.imageToView(polygon[i], dpt);
		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
	}
	IVP.viewToImage(UL1, dpt);
	ossimDpt UL(dpt.x - offset, dpt.y - offset);
	IVP.viewToImage(LR1, dpt);
	ossimDpt LR(dpt.x + offset, dpt.y + offset);
	IVP.imageToView(UL, UL1);
	IVP.imageToView(LR, LR1);
	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
	ossimDpt UR;
	IVP.viewToImage(ossimDpt(LR1.lon, UL1.lat), UR);
	ossimDpt LL;
	IVP.viewToImage(ossimDpt(UL1.lon, LR1.lat), LL);

	//选择输出波段 
	ossimBandSelector *theBandSelector = new ossimBandSelector;
	if (bandList.size() < 1)
	{
		size_t nbands = handler->getNumberOfInputBands();
		for (size_t i = 0; i < nbands; i++)
		{
			bandList.push_back(i);
		}
	}
	// 设置空值
	if (nullValueList.size() > 0 && nullValueList.size() == bandList.size())
	{
		for (size_t i = 0; i < bandList.size(); i++)
		{
			handler->setNullPixelValue(bandList[i], nullValueList[i]);
		}
	}
	theBandSelector->connectMyInputTo(handler);
	theBandSelector->setOutputBandList(bandList);

	ossimPolyCutter* theCutter = new ossimPolyCutter;
	theCutter->connectMyInputTo(theBandSelector);

	theCutter->setNumberOfPolygons(0);
	theCutter->addPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	ossimImageRenderer* renderer = new ossimImageRenderer;
	setResampleType(renderer);
	//renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	//renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
	renderer->connectMyInputTo(theCutter);
	renderer->setView(outmapinfo);

	//ossimImageFileWriter* writer =
	//	ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
	ossimImageFileWriter* writer;
	if (0 == strcmp(ossimFilename(strOutFile).ext().upcase(), "TIF"))
	{
		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_gtiff"));

	}
	else
	{
		writer = ossimImageWriterFactoryRegistry::instance()->
			createWriterFromExtension(ossimFilename(strOutFile).ext());
	}


	writer->setFilename(ossimFilename(strOutFile));


	ossimStdOutProgress progress(0, true);
	writer->addListener(&progress);
	if (nThreads < 1)
	{
		nThreads = OpenThreads::GetNumberOfProcessors() * 2;
	}
	cout << "using " << nThreads << " threads:" << endl;
	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
	writer->changeSequencer(sequencer.get());
	writer->connectMyInputTo(0, renderer);

	writer->setAreaOfInterest(viewRegion);
	writer->execute();

	polygon.clear();
	return true;
}

//bool cut_by_region(const char* inputImageFile, vector<ossimGpt> boundaryList, const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList/* = vector<ossim_uint32>()*/, int offset/* = 0*/)
//{
//	ossimImageHandler* handler;
//	handler = ossimImageHandlerRegistry::instance()->open(ossimFilename(inputImageFile));
//	if (!handler)
//	{
//		cout << "图像文件" << inputImageFile << "打开失败！" << endl;
//		system("Pause");
//		return false;
//	}
//
//	ossimRefPtr<ossimProjection> proj;
//	// 输入投影
//	ossimKeywordlist in_geom_kwl;
//	handler->getImageGeometry()->saveState(in_geom_kwl);
//	proj = handler->getImageGeometry()->getProjection();
//	ossimMapProjection* inmapinfo = PTR_CAST(ossimMapProjection, proj.get());
//
//	// 输出投影
//	ossimKeywordlist out_geom_kwl;
//	if (ossimFilename(strPrjFile).exists())
//	{
//		out_geom_kwl.addFile(strPrjFile);
//	}
//	else
//	{
//		handler->getImageGeometry()->getProjection()->saveState(out_geom_kwl);
//	}
//	proj = ossimProjectionFactoryRegistry::instance()->createProjection(out_geom_kwl);
//	ossimMapProjection* outmapinfo = PTR_CAST(ossimMapProjection, proj.get());
//
//	vector<ossimDpt>  polygon;
//	polygon.clear();
//
//	ossimDpt last_pt;
//	for (size_t i = 0; i < boundaryList.size(); i++)
//	{
//		ossimDpt pt = inmapinfo->worldToLineSample(boundaryList[i]);
//		ossimDpt delta_pt = last_pt - pt;
//		if (i > 0 && (last_pt - pt).length() < 60.0)
//		{
//			// 忽略距离小于1的冗余点
//			continue;
//		}
//
//		polygon.push_back(pt);
//		last_pt = pt;
//	}
//
//	ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
//	imageGeom->setProjection(outmapinfo);
//	handler->setImageGeometry(imageGeom.get());//1128
//
//	//选择输出波段 
//	ossimBandSelector *theBandSelector = new ossimBandSelector;
//	if (bandList.size() < 1)
//	{
//		size_t nbands = handler->getNumberOfInputBands();
//		for (size_t i = 0; i < nbands; i++)
//		{
//			bandList.push_back(i);
//		}
//	}
//	theBandSelector->connectMyInputTo(handler);
//
//	ossimPolyCutter* theCutter = new ossimPolyCutter;
//	theCutter->connectMyInputTo(theBandSelector);
//	theBandSelector->setOutputBandList(bandList);
//
//	theCutter->setNumberOfPolygons(0);
//	theCutter->addPolygon(polygon);
//	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
//	theCutter->setNumberOfPolygons(1);
//
//	//ossimIrect boundw1 = theCutter->getBoundingRect();
//	//ossimIrect boundw = handler->getBoundingRect();
//
//
//	//ossimImageFileWriter* writer =
//	//	ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
//	ossimImageFileWriter* writer;
//	if (0 == strcmp(ossimFilename(strOutFile).ext().upcase(), "TIF"))
//	{
//		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_gtiff"));
//
//	}
//	else
//	{
//		writer = ossimImageWriterFactoryRegistry::instance()->
//			createWriterFromExtension(ossimFilename(strOutFile).ext());
//	}
//
//
//	writer->setFilename(ossimFilename(strOutFile));
//
//	ossimImageRenderer* renderer = new ossimImageRenderer;
//
//	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
//	//renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
//
//	renderer->setView(outmapinfo);
//	renderer->connectMyInputTo(theCutter);
//	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, renderer->getImageViewTransform());
//
//	//统计外接矩形
//	ossimDpt dpt = renderer->getImageViewTransform()->imageToView(polygon[0]);
//	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
//	for (int i = 1; i < (int)polygon.size(); i++)
//	{
//		dpt = renderer->getImageViewTransform()->imageToView(polygon[i]);
//		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
//		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
//		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
//		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
//	}
//
//	dpt = renderer->getImageViewTransform()->viewToImage(UL1);
//	ossimDpt UL(dpt.x - offset, dpt.y - offset);
//	dpt = renderer->getImageViewTransform()->viewToImage(LR1);
//	ossimDpt LR(dpt.x + offset, dpt.y + offset);
//	UL1 = renderer->getImageViewTransform()->imageToView(UL);
//	LR1 = renderer->getImageViewTransform()->imageToView(LR);
//	//ossimDrect viewRegion = calcRegion(clipRect, inmapinfo, renderer->getImageViewTransform());
//	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
//	ossimDpt UR = renderer->getImageViewTransform()->viewToImage(ossimDpt(LR1.lon, UL1.lat));
//	ossimDpt LL = renderer->getImageViewTransform()->viewToImage(ossimDpt(UL1.lon, LR1.lat));
//
//	ossimStdOutProgress progress(0, true);
//	writer->addListener(&progress);
//	if (nThreads < 1)
//	{
//		nThreads = OpenThreads::GetNumberOfProcessors() * 2;
//	}
//	cout<<"using "<<nThreads<<" threads:"<<endl;
//	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
//	writer->changeSequencer(sequencer.get());
//	writer->connectMyInputTo(0, renderer);
//
//	writer->setAreaOfInterest(viewRegion);
//	writer->execute();
//
//	polygon.clear();
//	return true;
//}

/////////////////////////////////////////////////
// 根据shapefile矢量边界裁剪
//////////////////////////////////////////////////
bool cut_by_shp(const char* inputImageFile, const char* shpfileName, const char* fieldName, const char* searchValue,
				const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList = vector<ossim_uint32>(), int offset = 0)
{
	OGRDataSource       *poDS;
	poDS = OGRSFDriverRegistrar::Open(shpfileName, FALSE);
	if (poDS == NULL)
	{
		cout << "矢量文件" << shpfileName << "打开失败！" << endl;
		system("Pause");
		return false;
	}

	OGRLayer  *poLayer;
	poLayer = poDS->GetLayer(0);//GetLayerByName( "point" );
	OGRFeature *poFeature;
	poLayer->ResetReading();

	OGRSpatialReference *inOSRS = poLayer->GetSpatialRef();
	OGRSpatialReference poLatLong;
	poLatLong.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *poTransform = OGRCreateCoordinateTransformation(inOSRS, &poLatLong);
	if (poTransform == NULL)
	{
		cout << "Projection Error!" << endl;
		system("Pause");
	}

	vector<ossimGpt> boundaryList;
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		int polygonnum = 0;
		OGRGeometry *poGeometry;
		poGeometry = poFeature->GetGeometryRef();
		ossimString strFieldName = poFeature->GetFieldAsString(fieldName);

		if (poGeometry != NULL &&
			wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon &&
			0 == string(searchValue).compare(string(strFieldName)))
		{
			OGRPolygon* pol = (OGRPolygon *)poGeometry;
			OGRLinearRing* ring = pol->getExteriorRing();

			polygonnum = ring->getNumPoints();
			for (int i = 0; i<polygonnum; i++) {
				OGRPoint poPoint;
				ring->getPoint(i, &poPoint);
				double x = poPoint.getX();
				double y = poPoint.getY();
				if (!poTransform->Transform(1, &x, &y))
				{
					continue;
				}
				ossimGpt gpt(y, x);
				boundaryList.push_back(gpt);
			}
		}
		OGRFeature::DestroyFeature(poFeature);
	}////////end while

	OGRDataSource::DestroyDataSource(poDS);

	ossimFilename tmp_file = ossimFilename(strOutFile).path() + "\\tmp.tif";
	//cutbyRegion(inputImageFile, boundaryList, strPrjFile, tmp_file);
	cut_by_rect(inputImageFile, boundaryList, strPrjFile, tmp_file, bandList, offset);

	if (ossimFilename(strOutFile).exists())
	{
		// 如果结果文件已经存在，先将其删除
		remove(strOutFile);
	}
	mask_by_shp(tmp_file, shpfileName, fieldName, searchValue,	strPrjFile, strOutFile);

	remove(tmp_file);

	return true;
}

/////////////////////////////////////////////////
// 根据shapefile矢量边界裁剪，输出矩形
//////////////////////////////////////////////////
bool cut_by_shp_rect(const char* inputImageFile, const char* shpfileName, const char* fieldName, const char* searchValue,
					 const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList = vector<ossim_uint32>(), int offset = 0)
{
	OGRDataSource       *poDS;
	poDS = OGRSFDriverRegistrar::Open(shpfileName, FALSE);
	if (poDS == NULL)
	{
		cout << "矢量文件" << shpfileName << "打开失败！" << endl;
		system("Pause");
		return false;
	}

	OGRLayer  *poLayer;
	poLayer = poDS->GetLayer(0);//GetLayerByName( "point" );
	OGRFeature *poFeature;
	poLayer->ResetReading();

	OGRSpatialReference *inOSRS = poLayer->GetSpatialRef();
	OGRSpatialReference poLatLong;
	poLatLong.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *poTransform = OGRCreateCoordinateTransformation(inOSRS, &poLatLong);
	if (poTransform == NULL)
	{
		cout << "Projection Error!" << endl;
		system("Pause");
	}

	vector<ossimGpt> boundaryList;
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		int polygonnum = 0;
		OGRGeometry *poGeometry;
		poGeometry = poFeature->GetGeometryRef();
		ossimString strFieldName = poFeature->GetFieldAsString(fieldName);

		if (poGeometry != NULL &&
			wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon &&
			0 == string(searchValue).compare(string(strFieldName)))
		{
			OGRPolygon* pol = (OGRPolygon *)poGeometry;
			OGRLinearRing* ring = pol->getExteriorRing();

			polygonnum = ring->getNumPoints();
			for (int i = 0; i<polygonnum; i++) {
				OGRPoint poPoint;
				ring->getPoint(i, &poPoint);
				double x = poPoint.getX();
				double y = poPoint.getY();
				if (!poTransform->Transform(1, &x, &y))
				{
					continue;
				}
				ossimGpt gpt(y, x);
				boundaryList.push_back(gpt);
			}
		}
		OGRFeature::DestroyFeature(poFeature);
	}////////end while

	OGRDataSource::DestroyDataSource(poDS);

	return cut_by_rect(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
}



/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: %s "
		"[-shp shapefile -field fieldname -value value]\n"
		"[-boundary npt [lon lat lon lat ...] [-nv \"value [value...]\"]\n"
		"\t[-prj projectionfile] [-b nb b1 b2 ...] [-format outfiletype]\n"
		"\t[-scale scaletype] [-sample sampletype] [-nt nThreads]\n"
		"\t[-rect bRect] [-ref refImage] [-pref preferencefile]\n"
		"\t <inputfile> <outputfile>\n"

		"\n"
		"-boundary npt [lon lat lon lat ...]: use boundary points, and input n boundary "
			"points (no less than 2). When 2 points are used, the raster will be cutted by a rect region whose"
			" whose ul and lr points are indicated by these 2 points.\n"
		"  -i inputfile\t: input raster file\n"
		"  -o outputfile\t: output file name\n"
		"  -shp shapefile: cut by shapefile\n"
		"  -field fieldname: query the filedname in shapefile\n"
		"  -value value: query value of filedname within the shapefile\n"
		"  -rect: to clip by a rect region\n"
		"  -b nb b1 b2 ...: output band list, start from band 1\n"
		"  -nv \"value [value...]\": null value of each band\n"
		"  -nt nThreads: number of used threads\n"
		"  -prj projectionfile: output projection\n"
		"  -pref preferencefile: preference file\n"
		"  -scale scaletype: scaletype, e.g. UCHAR, UINT16, UINT8, FLOAT32, FLOAT64...\n"
		"  -sample sampletype: sampletype, e.g. NEAREST_NEIGHBOR, BILINEAR, BICUBIC\n"
		"  -format outfiletype: scaletype, e.g. TIFF, PIX, IMG\n",
		pszApplicationName);

	printf("Available resampling methods:\n"
		"\tnearest, bilinear, cubic (default)");

	if( pszErrorMsg != NULL )
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

//#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
//	do { if (i + nExtraArg >= argc) \
//	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)
void checkHasEnoughAdditionalArgs(int argc, char** argv, int currectIdx, int requiredNum)
{
	if (currectIdx + requiredNum >= argc)
	{
		char buf[2048];
		sprintf(buf, "%s option requires %d argument(s)", argv[currectIdx], requiredNum);
		Usage(buf);
	}
}

int parse1(int argc, char** argv)
{
	std::vector<string> unlabeledArgs;
	if (argc > 1)
	{
		/* -------------------------------------------------------------------- */
		/*      Parse arguments.                                                */
		/* -------------------------------------------------------------------- */
		for( int i = 1; i < argc; i++ )
		{
			if( 0 == _stricmp(argv[i],"-shp") )
			{
				bShapefile = true;
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				shpfileName = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-field") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				fieldName = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-value") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				searchValue = argv[++i] ;
			}
			else if(0 == _stricmp(argv[i],"-nt") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				nThreads = atoi( argv[++i] );	
			}
			else if (0 == _stricmp(argv[i], "-nv"))
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				string strNullValues = argv[++i];
				nullValueList.clear();
				vector<string> strList;
				mylib::splitString(strNullValues, { " " }, strList);
				for (size_t ib = 0; ib < strList.size(); ib++)
				{
					nullValueList.push_back(stod(strList[ib]));
				}
			}
			else if( 0 == _stricmp(argv[i],"-boundary") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				int nBoundaryPoints = atoi(argv[++i]);
				boundaryList.clear();
				for (int icount = 0;icount < nBoundaryPoints;++icount)
				{
					checkHasEnoughAdditionalArgs(argc, argv, i, 2);
					double lon = CPLAtofM(argv[++i]);
					double lat = CPLAtofM(argv[++i]);
					boundaryList.push_back(ossimGpt(lat, lon));
				}
				bBoundary = true;
			}
			else if( 0 == _stricmp(argv[i],"-b") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 2);
				int nBands = atoi(argv[++i]);
				if (nBands < 1)
				{
					printf("output band number cannot be less than 1\n");
					exit(0);
				}
				checkHasEnoughAdditionalArgs(argc, argv, i, nBands);
				for (int ib = 0;ib < nBands;++ib)
				{
					bandList.push_back(atoi(argv[++i])-1);
				}
			}
			else if( 0 == _stricmp(argv[i],"-prj") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				strPrjFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-scale") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				pszScaleType = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-rect") )
			{
				bRect = true ;
			}
			else if( 0 == _stricmp(argv[i],"-sample") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				resampleType = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-format") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				pszOutFileType = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-offset") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				offset = atoi(argv[++i]) ;
			}
			else if( 0 == _stricmp(argv[i],"-pref") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				pszPreferenceFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-ref") )
			{
				checkHasEnoughAdditionalArgs(argc, argv, i, 1);
				bRefImage = true;
				strRefImageFile = argv[++i] ;
			}
			else if (argv[i][0] == '-')
			{
				printf("unknown option %s\n", argv[i]);
				Usage();
			}
			else
			{
				unlabeledArgs.push_back(argv[i]);
			}
		}

		if (unlabeledArgs.size() < 2)
		{
			printf("\nFAILURE: No target filename specified.\n");
			Usage();
		}
		else
		{
			strSrcfile = unlabeledArgs[0];
			strDstfile = unlabeledArgs.back();
		}

		// 检查兼容性
		if ( bBoundary && bShapefile)
		{
			//
			Usage();
			printf("\"-boundary\" and \"-shp\" are not compatible, you need to choose only one of them.\n");
		}
		else if ( bBoundary && bRefImage)
		{
			//
			Usage();
			printf("\"-boundary\" and \"-ref\" are not compatible, you need to choose only one of them.\n");
		}
		else if ( bShapefile && bRefImage)
		{
			//
			Usage();
			printf("\"-shp\" and \"-ref\" are not compatible, you need to choose only one of them.\n");
		}
		else if (!bBoundary && !bShapefile && !bRefImage)
		{
			Usage();
			printf("you should set one of \"-boundary\", \"-shp\" and \"-ref\"to perform raster cutting.\n");
		}
		else if (bBoundary)
		{
			if(boundaryList.size() < 2)
			{
				Usage();
				printf("at least 2 boundary points are needed.\n");
			}
			else if (boundaryList.size() == 2)
			{
				bRect = true;
			}
	
	
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
			std::string	tempString;
			ossimArgumentParser::ossimParameter	stringParam(tempString);
			ossimArgumentParser argumentParser(&argc, argv);
			ossimInit::instance()->addOptions(argumentParser);
			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
			ossimInit::instance()->initialize(argumentParser);
	
			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
			if (bRect)
			{
				cut_by_rect(strSrcfile.c_str(), boundaryList, strPrjFile, strDstfile.c_str(), bandList, offset);
			}
			else
			{
				cut_by_region(strSrcfile.c_str(), boundaryList, strPrjFile, strDstfile.c_str(), bandList, offset);
			}
			clockEnd = clock();
			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
		else if (bShapefile)
		{
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
			std::string	tempString;
			ossimArgumentParser::ossimParameter	stringParam(tempString);
			ossimArgumentParser argumentParser(&argc, argv);
			ossimInit::instance()->addOptions(argumentParser);
			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
			ossimInit::instance()->initialize(argumentParser);
			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
	
			if (bRect)
			{
				cut_by_shp_rect(strSrcfile.c_str(), shpfileName, fieldName, searchValue, strPrjFile, strDstfile.c_str(), bandList, offset);
			}
			else{
				cut_by_shp(strSrcfile.c_str(), shpfileName, fieldName, searchValue, strPrjFile, strDstfile.c_str(), bandList, offset);
			}
			clockEnd = clock();
			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
	
		else if (bRefImage)
		{
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
			std::string	tempString;
			ossimArgumentParser::ossimParameter	stringParam(tempString);
			ossimArgumentParser argumentParser(&argc, argv);
			ossimInit::instance()->addOptions(argumentParser);
			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
			ossimInit::instance()->initialize(argumentParser);
			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
	
			boundaryList = getRasterBoundary(strRefImageFile);
	
			cut_by_rect(strSrcfile.c_str(), boundaryList, strPrjFile, strDstfile.c_str(), bandList, offset);
			clockEnd = clock();
			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
	}
}


//DEFINE_string(shape_file, "", "shape file");
//DEFINE_string(field_name, "", "field name");
//DEFINE_string(field_value, "", "field value");
//DEFINE_string(boundary, "", "boundary points (lon1,lat1,lon2,lat2...)");
//DEFINE_string(band, "", "bands (band1,band2...)");
//DEFINE_string(prj, "", "projection file");
//DEFINE_string(pref, "", "preference file");
//DEFINE_string(ref, "", "reference file");
//DEFINE_string(resample, "", "resample type");
//DEFINE_int32(offset, 0, "offset from the boundary");
//DEFINE_bool(rect, false, "output rect region");
//DEFINE_int32(nt, 0, "number of used threads");
//DEFINE_string(input_file, "", "input file");
//DEFINE_string(output_file, "", "output file");
#include <cmdline.h>
#include <optionparser.h>

#include <regex>
#include <tclap/CmdLine.h>
using namespace TCLAP;

int parse2(int argc, char** argv)
{

	// Wrap everything in a try block.  Do this every time, 
	// because exceptions will be thrown for problems. 
	try {

		// Define the command line object.
		CmdLine cmd("Command description message", ' ', "0.9");

		// Define a value argument and add it to the command line.
		UnlabeledValueArg<string> infile("infile", "input file", true, "", "infile");
		cmd.add(infile);
		UnlabeledValueArg<string> outfile("outfile", "output file", false, "", "outfile");
		cmd.add(outfile);
		ValueArg<string> pref("", "pref", "preference file", false, "", "pref");
		cmd.add(pref);
		ValueArg<int> nthread("", "nt", "number of used threads", false, 0, "nThreads");
		cmd.add(nthread);
		ValueArg<string> band("b", "band", "output bands, start from band 1", false, "", "\"band1 [band2...]\"");
		cmd.add(band);
		SwitchArg rect("", "rect", "output rect region", false);
		cmd.add(rect);
		ValueArg<int> offset("", "offset", "offset from the boundary", false, 0, "offset");
		cmd.add(offset);

		vector<string> resample_allowed;
		resample_allowed.push_back("nearest");
		resample_allowed.push_back("bilinear");
		resample_allowed.push_back("cubic");
		ValuesConstraint<string> resample_allowedVals(resample_allowed);
		ValueArg<string> resample("", "resample", "resample method", false, "cubic", &resample_allowedVals);
		cmd.add(resample);
		ValueArg<string> prj("", "prj", "output projection file", false, "", "projection_file");
		cmd.add(prj);
		ValueArg<string> ref("", "ref", "reference file", false, "", "ref");
		cmd.add(ref);
		ValueArg<string> boundary("", "boundary", "boundary points", false, "", "\"lon1 lat1 lon2 lat2 [lon3 lat3...]\"");
		cmd.add(boundary);
		ValueArg<string> field_value("", "field_value", "field value", false, "", "field_value");
		cmd.add(field_value);
		ValueArg<string> field_name("", "field_name", "query the filedname in shapefile", false, "", "field_name");
		cmd.add(field_name);
		ValueArg<string> shape_file("", "shape_file", "shape file", false, "", "shp_file");
		cmd.add(shape_file);

		// Parse the args.
		cmd.parse(argc, argv);


		string _infile = infile.getValue();
		string _outfile = outfile.getValue();
		string _pref = pref.getValue();
		nThreads = nthread.getValue();
		int _offset = offset.getValue();
		string _band = band.getValue();
		resampleType = resample.getValue();
		string _prj = prj.getValue();
		string _ref = ref.getValue();
		string _boundary = boundary.getValue();
		string _field_value = field_value.getValue();
		string _field_name = field_name.getValue();
		string _shape_file = shape_file.getValue();

		// 检查兼容性
		if (!_boundary.empty() && !_shape_file.empty())
		{
			//
			Usage();
			printf("\"-boundary\" and \"-shp\" are not compatible, you need to choose only one of them.\n");
		}
		else if (!_boundary.empty() && !_ref.empty())
		{
			//
			Usage();
			printf("\"-boundary\" and \"-ref\" are not compatible, you need to choose only one of them.\n");
		}
		else if (!_shape_file.empty() && !_ref.empty())
		{
			//
			Usage();
			printf("\"-shp\" and \"-ref\" are not compatible, you need to choose only one of them.\n");
		}
		else if (!_boundary.empty() && !_shape_file.empty() && !_ref.empty())
		{
			Usage();
			printf("you should set one of \"-boundary\", \"-shp\" and \"-ref\"to perform raster cutting.\n");
		}
		else if (!_boundary.empty())
		{
			std::vector<string> ptList;
			mylib::splitString(_boundary.c_str(), { " " }, ptList);
			if (0 != ptList.size() % 2)
			{
				printf("each point should contain two coordinates.\n");
				return 0;
			}
			boundaryList.clear();
			for (size_t ip = 0; ip < ptList.size(); ip += 2)
			{
				boundaryList.push_back(ossimGpt(std::stod(ptList[ip + 1]), std::stod(ptList[ip])));
			}
			if (boundaryList.size() < 2)
			{
				Usage();
				printf("at least 2 boundary points are needed.\n");
			}
			else if (boundaryList.size() == 2)
			{
				bRect = true;
			}

			std::vector<string> bList;// = mylib::split(_band, " ");
			mylib::splitString(_band.c_str(), { " " }, bList);
			bandList.clear();
			for (size_t ip = 0; ip < bList.size(); ip += 1)
			{
				bandList.push_back(std::stoi(bList[ip]) - 1);
			}

			ossimFilename preferences_file = ossimFilename(_pref);
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

			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
			if (bRect)
			{
				cut_by_rect(_infile.c_str(), boundaryList, _prj.c_str(), _outfile.c_str(), bandList, _offset);
			}
			else
			{
				cut_by_region(_infile.c_str(), boundaryList, _prj.c_str(), _outfile.c_str(), bandList, _offset);
			}
			clockEnd = clock();
			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
		else if (!_shape_file.empty())
		{
			ossimFilename preferences_file = ossimFilename(_pref);
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
			clock_t  clockBegin, clockEnd;
			clockBegin = clock();

			if (bRect)
			{
				cut_by_shp_rect(_infile.c_str(), _shape_file.c_str(), _field_name.c_str(),
					_field_value.c_str(), _prj.c_str(), _outfile.c_str(), bandList, _offset);
			}
			else{
				cut_by_shp(_infile.c_str(), _shape_file.c_str(), _field_name.c_str(),
					_field_value.c_str(), _prj.c_str(), _outfile.c_str(), bandList, _offset);
			}
			clockEnd = clock();
			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
		}

		else if (bRefImage)
		{
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
			std::string	tempString;
			ossimArgumentParser::ossimParameter	stringParam(tempString);
			ossimArgumentParser argumentParser(&argc, argv);
			ossimInit::instance()->addOptions(argumentParser);
			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
			ossimInit::instance()->initialize(argumentParser);
			clock_t  clockBegin, clockEnd;
			clockBegin = clock();

			boundaryList = getRasterBoundary(strRefImageFile);

			cut_by_rect(_infile.c_str(), boundaryList, strPrjFile, _outfile.c_str(), bandList, _offset);
			clockEnd = clock();
			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
	}
	catch (ArgException &e)  // catch any exceptions
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		return 0;
	}
}

int main( int argc, char** argv )
{
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持
	// 程序名
	pszApplicationName = argv[0];
	parse1(argc, argv);

	//cmdline::parser a;
	//a.add<string>("shape_file", "", "shape file", false, "");
	//a.add<string>("field_name", "", "field name", false, "");
	//a.add<string>("field_value", "", "field value", false, "");
	//a.add<string>("prj", "", "projection file", false, "");
	//a.add<string>("pref", "", "preference file", false, "");
	//a.add<string>("ref", "", "reference file", false, "");
	//a.add<string>("boundary", "", "\"lon lat [lon lat...]\" boundary points", false, "");
	//a.add<string>("resample", "", "resample method", false, "cubic", cmdline::oneof<string>("nearest", "bilinear", "cubic"));
	//a.add<int>("offset", "", "offset from the boundary", false, 0);
	//a.add("rect", "", "output rect region");
	//a.add<int>("nthread", "nt", "number of used threads", false, 0);
	//a.add("help", "?", "print this message");
	//a.add<string>("band", "b", "output band", false, "");
	////a.add< string >("band", 'b', "output band", false,  "" );
	//a.footer("<infile> [<outfile>]");
	//a.set_program_name("img-clip");

	//bool ok = a.parse(argc, argv);
	//a.parse_check(argc, argv);

	//if (argc == 1 || a.exist("help")){
	//	cerr << a.usage();
	//	return 0;
	//}

	//if (a.exist("nthread")){
	//	cout << "nThread = " << a.get<int>("nthread")<<endl;
	//	return 0;
	//}

	//if (!ok){
	//	cerr << a.error() << endl << a.usage();
	//	return 0;
	//}


	////cout << a.get< string >("band") << endl;
	////vector<string> bandList = a.get< vector<string> >("band");
	////cout << bandList.size() << endl;
	////for (size_t i = 0; i < bandList.size(); i++)
	////{
	////	cout << bandList[i] << " ";
	////}
	////cout << endl;

	////cout << a.get<string>("host") << ":" << a.get<int>("port") << endl;

	//for (size_t i = 0; i<a.rest().size(); i++)
	//	cout << "- " << a.rest()[i] << endl;

	//gflags::ParseCommandLineFlags(&argc, &argv, true);

	//// Declare the supported options.
	//po::options_description desc("Allowed options");
	//desc.add_options()
	//	("help,h", "produce help message")
	//	("shp", po::value<string>(), "shape file")
	//	("field", po::value<string>(), "field name")
	//	("value", po::value<string>(), "field value")
	//	("boundary", po::value< vector<string> >(), "boundary point lon1,lat1")
	//	("band,b", po::value< vector<string> >(), "band")
	//	("prj", po::value<string>(), "projection file")
	//	("resample", po::value<string>(), "resample type")
	//	("offset", po::value<int>(), "offset from the boundary")
	//	("pref", po::value<string>(), "preference file")
	//	("ref", po::value<string>(), "reference file")
	//	("rect", "output rect region")
	//	("nt", po::value<int>(), "number of used threads")
	//	("input-file", po::value<string>(), "input file")
	//	("output-file", po::value<string>(), "output file");

	//po::positional_options_description p;
	//p.add("input-file", 1);
	//p.add("output-file", 1);

	//po::variables_map vm;
	//try{
	//	po::store(po::command_line_parser(argc, argv).
	//		options(desc).positional(p).run(), vm);
	//	//po::store(po::parse_command_line(argc, argv, desc), vm);
	//}
	//catch (std::exception const&  ex)
	//{
	//	printf("%s\n", ex.what());
	//}
	//po::notify(vm);

	//if (vm.count("help")) {
	//	cout << desc << "\n";
	//	return 1;
	//}

	//if (vm.count("boundary")) {
	//	vector<string> pointList = vm["boundary"].as< vector<string> >();
	//	cout << "boundary points are: ";
	//	//cout << vm["boundary"].as< vector<string> >() << ".\n";
	//	for (size_t ip = 0; ip < pointList.size(); ip++)
	//	{
	//		cout << pointList[ip] << " ";
	//	}
	//	cout << ".\n";
	//	//cout << vm["boundary"].as< vector<string> >() << ".\n";
	//}
	//else {
	//	cout << "Compression level was not set.\n";
	//}

	//if (vm.count("input-file"))
	//{
	//	cout << "input-file is: " << vm["input-file"].as< string >() << "\n";
	//}

	//if (vm.count("output-file"))
	//{
	//	cout << "output-file is: " << vm["output-file"].as< string >() << "\n";
	//}

	////system("pause");


	//struct arg_file *shape_file = arg_file0(NULL, "shp,shapefile", NULL, "shape file");
	//struct arg_str *field_name = arg_str0(NULL, "field", NULL, "field name");
	//struct arg_str *field_value = arg_str0(NULL, "value", NULL, "field value");
	//struct arg_dbl *boundary = arg_dbln(NULL, "boundary", NULL, 0, argc + 2, "boundary points (lon1 lat1 lon2 lat2 ...)");
	//struct arg_file *preference_file = arg_file0(NULL, "pref,preference", NULL, "projection file");
	//struct arg_int *out_bands = arg_intn("b", NULL, NULL, 0, argc + 3, "output bands (begin with 1)");
	//struct arg_int *use_threads = arg_int0(NULL, "nt,nthreads", NULL, "the number of threads (if 0, automatically choose upon the processors)");
	//struct arg_int *offset = arg_int0(NULL, "offset", NULL, "offset from the given boundary");
	//struct arg_file *reference_file = arg_file0(NULL, "ref,referencefile", NULL, "reference file");
	//struct arg_lit *sample = arg_lit0(NULL, "rpc", "mandatorily use rpc ");
	//struct arg_lit *rect = arg_lit0(NULL, "l1", "output rect region");
	//struct arg_lit *overwrite = arg_lit0("ow", "overwrite", "overwrite the existing file");
	//struct arg_lit  *help = arg_lit0("h", "help", "print this help and exit");
	//struct arg_lit  *version = arg_lit0("v", "version", "print version information and exit");
	//struct arg_file *infile = arg_file1(NULL, NULL, "<input file>", "input file");
	//struct arg_file *outfile = arg_file1(NULL, NULL, "<output file>", "output file");
//	struct arg_end  *end = arg_end(20);
//
//
//	void* argtable[] = { shape_file, field_name, field_value, boundary,
//		preference_file, out_bands, use_threads, offset, reference_file, sample, rect, overwrite, help, version, infile, outfile, end };
//	const char* progname = "img-clip";
//	int nerrors;
//	int exitcode = 0;
//
//	/* verify the argtable[] entries were allocated sucessfully */
//	if (arg_nullcheck(argtable) != 0)
//	{
//		/* NULL entries were detected, some allocations must have failed */
//		printf("%s: insufficient memory\n", progname);
//		exitcode = 1;
//		goto exit;
//	}
//
//	// default value
//	use_threads->ival[0] = 0;
//
//	/* Parse the command line as defined by argtable[] */
//	nerrors = arg_parse(argc, argv, argtable);
//
//	/* special case: '--help' takes precedence over error reporting */
//	if (help->count > 0)
//	{
//		printf("Usage: %s", progname);
//		arg_print_syntax(stdout, argtable, "\n");
//		printf("This program ortho-rectifies remotely sensed image.\n");
//		//printf("for parsing command line arguments. Argtable accepts integers\n");
//		//printf("in decimal (123), hexadecimal (0xff), octal (0o123) and binary\n");
//		//printf("(0b101101) formats. Suffixes KB, MB and GB are also accepted.\n");
//		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
//		exitcode = 0;
//		goto exit;
//	}
//
//	/* special case: '--version' takes precedence error reporting */
//	if (version->count > 0)
//	{
//		printf("'%s' Ortho-rectification.\n", progname);
//		printf("December 2014, Long Tengfei.\n");
//		exitcode = 0;
//		goto exit;
//	}
//
//	/* If the parser returned any errors then display them and exit */
//	if (nerrors > 0)
//	{
//		/* Display the error details contained in the arg_end struct.*/
//		arg_print_errors(stdout, end, progname);
//		printf("Try '%s -h or --help' for more information.\n", progname);
//		exitcode = 1;
//		goto exit;
//	}
//
//	/* special case: uname with no command line options induces brief help */
//	if (argc == 1)
//	{
//		printf("Try '%s -h or --help' for more information.\n", progname);
//		exitcode = 0;
//		goto exit;
//	}
//
//	if (boundary->count > 0 && boundary)
//	{
//	}
//
//	// 检查兼容性
//	if (bBoundary && bShapefile)
//	{
//		//
//		Usage();
//		printf("\"-boundary\" and \"-shp\" are not compatible, you need to choose only one of them.\n");
//	}
//	else if (bBoundary && bRefImage)
//	{
//		//
//		Usage();
//		printf("\"-boundary\" and \"-ref\" are not compatible, you need to choose only one of them.\n");
//	}
//	else if (bShapefile && bRefImage)
//	{
//		//
//		Usage();
//		printf("\"-shp\" and \"-ref\" are not compatible, you need to choose only one of them.\n");
//	}
//	else if (!bBoundary && !bShapefile && !bRefImage)
//	{
//		Usage();
//		printf("you should set one of \"-boundary\", \"-shp\" and \"-ref\"to perform raster cutting.\n");
//	}
//	else if (bBoundary)
//	{
//		if (boundaryList.size() < 2)
//		{
//			Usage();
//			printf("at least 2 boundary points are needed.\n");
//		}
//		else if (boundaryList.size() == 2)
//		{
//			bRect = true;
//		}
//
//
//		ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//		if (preferences_file.exists())
//		{
//			ossimPreferences::instance()->loadPreferences(preferences_file);
//		}
//		else
//		{
//			preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//			if (preferences_file.exists())
//			{
//				ossimPreferences::instance()->loadPreferences(preferences_file);
//			}
//		}
//		std::string	tempString;
//		ossimArgumentParser::ossimParameter	stringParam(tempString);
//		ossimArgumentParser argumentParser(&argc, argv);
//		ossimInit::instance()->addOptions(argumentParser);
//		argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//		ossimInit::instance()->initialize(argumentParser);
//
//		clock_t  clockBegin, clockEnd;
//		clockBegin = clock();
//		if (bRect)
//		{
//			cut_by_rect(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
//		}
//		else
//		{
//			cut_by_region(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
//		}
//		clockEnd = clock();
//		printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
//	}
//	else if (bShapefile)
//	{
//		ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//		if (preferences_file.exists())
//		{
//			ossimPreferences::instance()->loadPreferences(preferences_file);
//		}
//		else
//		{
//			preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//			if (preferences_file.exists())
//			{
//				ossimPreferences::instance()->loadPreferences(preferences_file);
//			}
//		}
//		std::string	tempString;
//		ossimArgumentParser::ossimParameter	stringParam(tempString);
//		ossimArgumentParser argumentParser(&argc, argv);
//		ossimInit::instance()->addOptions(argumentParser);
//		argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//		ossimInit::instance()->initialize(argumentParser);
//		clock_t  clockBegin, clockEnd;
//		clockBegin = clock();
//
//		if (bRect)
//		{
//			cut_by_shp_rect(inputImageFile, shpfileName, fieldName, searchValue, strPrjFile, strOutFile, bandList, offset);
//		}
//		else{
//			cut_by_shp(inputImageFile, shpfileName, fieldName, searchValue, strPrjFile, strOutFile, bandList, offset);
//		}
//		clockEnd = clock();
//		printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
//	}
//
//	else if (bRefImage)
//	{
//		ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//		if (preferences_file.exists())
//		{
//			ossimPreferences::instance()->loadPreferences(preferences_file);
//		}
//		else
//		{
//			preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//			if (preferences_file.exists())
//			{
//				ossimPreferences::instance()->loadPreferences(preferences_file);
//			}
//		}
//		std::string	tempString;
//		ossimArgumentParser::ossimParameter	stringParam(tempString);
//		ossimArgumentParser argumentParser(&argc, argv);
//		ossimInit::instance()->addOptions(argumentParser);
//		argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//		ossimInit::instance()->initialize(argumentParser);
//		clock_t  clockBegin, clockEnd;
//		clockBegin = clock();
//
//		boundaryList = getRasterBoundary(strRefImageFile);
//
//		cut_by_rect(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
//		clockEnd = clock();
//		printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
//	}
//
//	/* normal case: take the command line options at face value */
//	exitcode = mymain(gcp_file, report_file, projection_file, preference_file,
//		model_file, out_bands, use_threads, overview, rpc, l1, overwrite, help, version, infile, outfile);
//
//exit:
//	/* deallocate each non-null entry in argtable[] */
//	arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
//
//	return exitcode;
//
//
//	if (argc > 1)
//	{
//		/* -------------------------------------------------------------------- */
//		/*      Parse arguments.                                                */
//		/* -------------------------------------------------------------------- */
//		for( int i = 1; i < argc; i++ )
//		{
//			if( 0 == _stricmp(argv[i],"-i") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				inputImageFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-o") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				strOutFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-shp") )
//			{
//				bShapefile = true;
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				shpfileName = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-field") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				fieldName = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-value") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				searchValue = argv[++i] ;
//			}
//			else if(0 == _stricmp(argv[i],"-nt") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				nThreads = atoi( argv[++i] );	
//			}
//			else if( 0 == _stricmp(argv[i],"-boundary") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				int nBoundaryPoints = atoi(argv[++i]);
//				boundaryList.clear();
//				for (int icount = 0;icount < nBoundaryPoints;++icount)
//				{
//					CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(2);
//					double lon = CPLAtofM(argv[++i]);
//					double lat = CPLAtofM(argv[++i]);
//					boundaryList.push_back(ossimGpt(lat, lon));
//				}
//				bBoundary = true;
//			}
//			else if( 0 == _stricmp(argv[i],"-b") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				int nBands = atoi(argv[++i]);
//				if (nBands < 1)
//				{
//					printf("output band number cannot less than 1\n");
//					exit(0);
//				}
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nBands);
//				for (int ib = 0;ib < nBands;++ib)
//				{
//					bandList.push_back(atoi(argv[++i])-1);
//				}
//			}
//			else if( 0 == _stricmp(argv[i],"-prj") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				strPrjFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-scale") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszScaleType = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-rect") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
//				bRect = true ;
//			}
//			else if( 0 == _stricmp(argv[i],"-sample") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszSampleType = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-format") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszOutFileType = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-offset") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				offset = atoi(argv[++i]) ;
//			}
//			else if( 0 == _stricmp(argv[i],"-pref") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszPreferenceFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-ref") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				bRefImage = true;
//				strRefImageFile = argv[++i] ;
//			}
//			else
//			{
//				Usage();
//			}
//		}
//
//		// 检查兼容性
//		if ( bBoundary && bShapefile)
//		{
//			//
//			Usage();
//			printf("\"-boundary\" and \"-shp\" are not compatible, you need to choose only one of them.\n");
//		}
//		else if ( bBoundary && bRefImage)
//		{
//			//
//			Usage();
//			printf("\"-boundary\" and \"-ref\" are not compatible, you need to choose only one of them.\n");
//		}
//		else if ( bShapefile && bRefImage)
//		{
//			//
//			Usage();
//			printf("\"-shp\" and \"-ref\" are not compatible, you need to choose only one of them.\n");
//		}
//		else if (!bBoundary && !bShapefile && !bRefImage)
//		{
//			Usage();
//			printf("you should set one of \"-boundary\", \"-shp\" and \"-ref\"to perform raster cutting.\n");
//		}
//		else if (bBoundary)
//		{
//			if(boundaryList.size() < 2)
//			{
//				Usage();
//				printf("at least 2 boundary points are needed.\n");
//			}
//			else if (boundaryList.size() == 2)
//			{
//				bRect = true;
//			}
//
//
//			ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//			if (preferences_file.exists())
//			{
//				ossimPreferences::instance()->loadPreferences(preferences_file);
//			}
//			else
//			{
//				preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//				if (preferences_file.exists())
//				{
//					ossimPreferences::instance()->loadPreferences(preferences_file);
//				}
//			}
//			std::string	tempString;
//			ossimArgumentParser::ossimParameter	stringParam(tempString);
//			ossimArgumentParser argumentParser(&argc, argv);
//			ossimInit::instance()->addOptions(argumentParser);
//			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//			ossimInit::instance()->initialize(argumentParser);
//
//			clock_t  clockBegin, clockEnd;
//			clockBegin = clock();
//			if (bRect)
//			{
//				cut_by_rect(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
//			}
//			else
//			{
//				cut_by_region(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
//			}
//			clockEnd = clock();
//			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
//		}
//		else if (bShapefile)
//		{
//			ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//			if (preferences_file.exists())
//			{
//				ossimPreferences::instance()->loadPreferences(preferences_file);
//			}
//			else
//			{
//				preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//				if (preferences_file.exists())
//				{
//					ossimPreferences::instance()->loadPreferences(preferences_file);
//				}
//			}
//			std::string	tempString;
//			ossimArgumentParser::ossimParameter	stringParam(tempString);
//			ossimArgumentParser argumentParser(&argc, argv);
//			ossimInit::instance()->addOptions(argumentParser);
//			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//			ossimInit::instance()->initialize(argumentParser);
//			clock_t  clockBegin, clockEnd;
//			clockBegin = clock();
//
//			if (bRect)
//			{
//				cut_by_shp_rect(inputImageFile, shpfileName, fieldName, searchValue, strPrjFile, strOutFile, bandList, offset);
//			}
//			else{
//				cut_by_shp(inputImageFile, shpfileName, fieldName, searchValue, strPrjFile, strOutFile, bandList, offset);
//			}
//			clockEnd = clock();
//			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
//		}
//
//		else if (bRefImage)
//		{
//			ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//			if (preferences_file.exists())
//			{
//				ossimPreferences::instance()->loadPreferences(preferences_file);
//			}
//			else
//			{
//				preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//				if (preferences_file.exists())
//				{
//					ossimPreferences::instance()->loadPreferences(preferences_file);
//				}
//			}
//			std::string	tempString;
//			ossimArgumentParser::ossimParameter	stringParam(tempString);
//			ossimArgumentParser argumentParser(&argc, argv);
//			ossimInit::instance()->addOptions(argumentParser);
//			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//			ossimInit::instance()->initialize(argumentParser);
//			clock_t  clockBegin, clockEnd;
//			clockBegin = clock();
//
//			boundaryList = getRasterBoundary(strRefImageFile);
//
//			cut_by_rect(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
//			clockEnd = clock();
//			printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
//		}
//	}
//	else
//	{
//		boundaryList.push_back(ossimGpt(41, 116, 0));
//		boundaryList.push_back(ossimGpt(40, 117, 0));
//		pszPreferenceFile = "D:\\opensource\\ossim\\preference.txt";
//		inputImageFile = "G:\\testdata\\tmchina\\china.vrt";
//		//strPrjFile = "D:\\workspace\\dem\\prj.geom";
//		strOutFile = "D:\\workspace\\dem\\dem.tif";
//		bandList.clear();
//		bandList.push_back(0);
//		//nThreads = 1;
//
//		if (boundaryList.size() < 2)
//		{
//			Usage();
//			printf("at least 2 boundary points are needed.\n");
//		}
//		else if (boundaryList.size() == 2)
//		{
//			bRect = true;
//		}
//
//		ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//		if (preferences_file.exists())
//		{
//			ossimPreferences::instance()->loadPreferences(preferences_file);
//		}
//		else
//		{
//			preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//			if (preferences_file.exists())
//			{
//				ossimPreferences::instance()->loadPreferences(preferences_file);
//			}
//		}
//		std::string	tempString;
//		ossimArgumentParser::ossimParameter	stringParam(tempString);
//		ossimArgumentParser argumentParser(&argc, argv);
//		ossimInit::instance()->addOptions(argumentParser);
//		argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//		ossimInit::instance()->initialize(argumentParser);
//
//		clock_t  clockBegin, clockEnd;
//		clockBegin = clock();
//		if (bRect)
//		{
//			cut_by_rect(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
//		}
//		else
//		{
//			cut_by_region(inputImageFile, boundaryList, strPrjFile, strOutFile, bandList, offset);
//		}
//		clockEnd = clock();
//		printf("Time elapsed : %lf s\n", (clockEnd - clockBegin)*1e-3);
//
//		Usage(0);
//	}
	return 0;
}