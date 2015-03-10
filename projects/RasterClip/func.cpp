#include "func.h"
#include "gdalwarp.h"
#include "global.h"

/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"
#include "boost/shared_ptr.hpp"

#include <ossim/parallel/ossimMultiThreadSequencer.h>


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

/////////////////////////////////////////////////
// 根据shapefile矢量边界裁剪
//////////////////////////////////////////////////
bool cut_by_shp(const char* inputImageFile, const char* shpfileName, const char* fieldName, const char* searchValue,
	const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList/* = vector<ossim_uint32>()*/, int offset/* = 0*/)
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
	const char* strPrjFile, const char* strOutFile, vector<ossim_uint32> bandList/* = vector<ossim_uint32>()*/, int offset/* = 0*/)
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

///////////////////////////////////////////
// 根据矩形区域进行裁剪
///////////////////////////////////////////
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
		proj = ossimProjectionFactoryRegistry::instance()->createProjection(out_geom_kwl);
	}
	else
	{
		out_geom_kwl = in_geom_kwl;
	}
	ossimMapProjection* outmapinfo = PTR_CAST(ossimMapProjection, proj.get());

		
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

	for (int i = 0; i < (int) polygon.size(); i++)
	{
		xmin = (polygon[i].x < xmin) ? polygon[i].x : xmin;
		ymin = (polygon[i].y < ymin) ? polygon[i].y : ymin;
		xmax = (polygon[i].x > xmax) ? polygon[i].x : xmax;
		ymax = (polygon[i].y > ymax) ? polygon[i].y : ymax;
	}
	
	handler->loadState(out_geom_kwl);

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
	theBandSelector->connectMyInputTo(0, handler);
	theBandSelector->setOutputBandList(bandList);

	ossimPolyCutter* theCutter = new ossimPolyCutter;

	vector<ossimDpt> clipRect;
	//clipRect.push_back(ossimDpt(xmin, ymin));
	//clipRect.push_back(ossimDpt(xmax, ymin));
	//clipRect.push_back(ossimDpt(xmax, ymax));
	//clipRect.push_back(ossimDpt(xmin, ymin));
	//clipRect.push_back(ossimDpt(xmin, ymin));


	//theCutter->connectMyInputTo(theBandSelector);

	//theCutter->setNumberOfPolygons(0);
	//theCutter->addPolygon(clipRect);
	//theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	//theCutter->setNumberOfPolygons(1);

	//ossimIrect boundw1 = theCutter->getBoundingRect();
	//ossimIrect boundw = handler->getBoundingRect();


	ossimImageRenderer* renderer = new ossimImageRenderer;

	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);

	renderer->connectMyInputTo(theBandSelector);
	renderer->setView(outmapinfo);
	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, renderer->getImageViewTransform());

	//统计外接矩形
	ossimDpt dpt = renderer->getImageViewTransform()->imageToView(polygon[0]);
	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
	for (int i = 1; i < (int) polygon.size(); i++)
	{
		dpt = renderer->getImageViewTransform()->imageToView(polygon[i]);
		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
	}

	dpt = renderer->getImageViewTransform()->viewToImage(UL1);
	ossimDpt UL(dpt.x - offset, dpt.y - offset);
	dpt = renderer->getImageViewTransform()->viewToImage(LR1);
	ossimDpt LR(dpt.x + offset, dpt.y + offset);
	UL1 = renderer->getImageViewTransform()->imageToView(UL);
	LR1 = renderer->getImageViewTransform()->imageToView(LR);
	//ossimDrect viewRegion = calcRegion(clipRect, inmapinfo, renderer->getImageViewTransform());
	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
	ossimDpt UR = renderer->getImageViewTransform()->viewToImage(ossimDpt(LR1.lon, UL1.lat));
	ossimDpt LL = renderer->getImageViewTransform()->viewToImage(ossimDpt(UL1.lon, LR1.lat));

	clipRect.clear();
	clipRect.push_back(UL);
	clipRect.push_back(UR);
	clipRect.push_back(LR);
	clipRect.push_back(LL);
	clipRect.push_back(UL);

	theCutter->setNumberOfPolygons(0);
	theCutter->addPolygon(clipRect);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	theCutter->connectMyInputTo(theBandSelector);

	renderer->connectMyInputTo(theBandSelector);

	ossimImageFileWriter* writer =
		ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));

	writer->setFilename(ossimFilename(strOutFile));

	ossimStdOutProgress progress(0, true);
	writer->addListener(&progress);

	if (nThreads < 1)
	{
		nThreads = OpenThreads::GetNumberOfProcessors() * 2;
	}
	cout<<"using "<<nThreads<<" threads:"<<endl;
	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
	writer->changeSequencer(sequencer.get());
	writer->connectMyInputTo(0, renderer);

	writer->setAreaOfInterest(viewRegion);
	writer->execute();

	//delete writer;
	//delete theCutter;
	//	delete renderer;
	polygon.clear();
	handler->close();

	return true;
}


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
		proj = ossimProjectionFactoryRegistry::instance()->createProjection(out_geom_kwl);
	}
	else
	{
		out_geom_kwl = in_geom_kwl;
	}
	ossimMapProjection* outmapinfo = PTR_CAST(ossimMapProjection, proj.get());


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
	
	ossimPolyCutter* theCutter = new ossimPolyCutter;
	handler->loadState(out_geom_kwl);

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
	theCutter->connectMyInputTo(theBandSelector);
	theBandSelector->setOutputBandList(bandList);

	theCutter->setNumberOfPolygons(0);
	theCutter->addPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	//ossimIrect boundw1 = theCutter->getBoundingRect();
	//ossimIrect boundw = handler->getBoundingRect();


	ossimImageFileWriter* writer =
		ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));


	writer->setFilename(ossimFilename(strOutFile));

	ossimImageRenderer* renderer = new ossimImageRenderer;

	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);

	renderer->connectMyInputTo(theCutter);
	renderer->setView(outmapinfo);
	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, renderer->getImageViewTransform());

	//统计外接矩形
	ossimDpt dpt = renderer->getImageViewTransform()->imageToView(polygon[0]);
	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
	for (int i = 1; i < (int)polygon.size(); i++)
	{
		dpt = renderer->getImageViewTransform()->imageToView(polygon[i]);
		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
	}

	dpt = renderer->getImageViewTransform()->viewToImage(UL1);
	ossimDpt UL(dpt.x - offset, dpt.y - offset);
	dpt = renderer->getImageViewTransform()->viewToImage(LR1);
	ossimDpt LR(dpt.x + offset, dpt.y + offset);
	UL1 = renderer->getImageViewTransform()->imageToView(UL);
	LR1 = renderer->getImageViewTransform()->imageToView(LR);
	//ossimDrect viewRegion = calcRegion(clipRect, inmapinfo, renderer->getImageViewTransform());
	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
	ossimDpt UR = renderer->getImageViewTransform()->viewToImage(ossimDpt(LR1.lon, UL1.lat));
	ossimDpt LL = renderer->getImageViewTransform()->viewToImage(ossimDpt(UL1.lon, LR1.lat));

	ossimStdOutProgress progress(0, true);
	writer->addListener(&progress);
	if (nThreads < 1)
	{
		nThreads = OpenThreads::GetNumberOfProcessors() * 2;
	}
	cout<<"using "<<nThreads<<" threads:"<<endl;
	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
	writer->changeSequencer(sequencer.get());
	writer->connectMyInputTo(0, renderer);

	writer->setAreaOfInterest(viewRegion);
	writer->execute();

	polygon.clear();
	return true;
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
	argv[argc++] = const_cast<char*>(inputImageFile);
	argv[argc++] = const_cast<char*>(strOutFile);

	gdalwarp::gdalwarp( argc, argv);
	return true;
}
