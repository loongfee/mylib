
/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"
#include "boost/shared_ptr.hpp"

#include <stdlib.h>
#include <math.h>

#include "strUtil.h"
#include "fileUtil.h"
#include <QTextStream>
#include <QObject>
#include <QtCore/QCoreApplication>

#include <fstream>

using namespace std;

#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "OpenThreads.lib")


#include <OpenThreads/Thread>
#include <OpenThreads/Mutex>
#include <OpenThreads/Barrier>

namespace rasterBoundary{

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

COORD GetConsoleCursorPosition(HANDLE hHandle)
{  
	CONSOLE_SCREEN_BUFFER_INFO info={0};  
	GetConsoleScreenBufferInfo( hHandle , &info );  
	return info.dwCursorPosition;  
}

enum EdgeType
{
	top = 0,
	right,
	left,
	bottom,
};

void findTop(const char* szFilename, Point2f &pt, int tile_size = 1024, int null_value = 0)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDatasetH hSrcDs = GDALOpen(szFilename, GA_ReadOnly );

	CPLErr eErr = CE_None;
	int nRasterXSizeRead,nRasterYSizeRead;

	GDALDataset  *hSrcDataset = (GDALDataset  *)hSrcDs;

	nRasterXSizeRead=hSrcDataset->GetRasterXSize();
	nRasterYSizeRead=hSrcDataset->GetRasterYSize();
	int nBands = hSrcDataset->GetRasterCount();
	GDALDataType eDT = hSrcDataset->GetRasterBand(1)->GetRasterDataType();

	const int nLines	= tile_size;

	int nTiles = (int)ceil(nRasterYSizeRead / (double)nLines);

	GByte *pData;
	for (int i=0;i<nTiles;++i)
	{
		int BufHeight = nLines;
		//列末尾小块处理
		if (i == nTiles-1)
		{
			BufHeight = nRasterYSizeRead - (nTiles-1) * nLines;
			BufHeight = min(BufHeight, nLines);
		}

		pData=(GByte *) CPLMalloc (nRasterXSizeRead * BufHeight * nBands * GDALGetDataTypeSize( eDT ) / 8);
		eErr = hSrcDataset->RasterIO( GF_Read, 0, i*nLines, nRasterXSizeRead , BufHeight, pData, nRasterXSizeRead, BufHeight, eDT, nBands, 0, 0 ,0, 0);
		for (int ii = 0; ii < BufHeight;++ii)
		{
			for (int jj = 0;jj < nRasterXSizeRead;++jj)
			{
				for (int iBand = 0;iBand < nBands;++iBand)
				{
					//int pos = ii * BufWidth * nBands + jj * nBands + iBand;
					int pos = nRasterXSizeRead*BufHeight*iBand + ii * nRasterXSizeRead + jj;					
					if(SRCVAL(pData, eDT, pos) != null_value)
					{
						pt.x = jj;
						pt.y = i * nLines + ii;
						CPLFree( pData );
						pData = NULL;
						GDALClose( hSrcDataset );
						return;
					}
				}
			}
		}
		CPLFree( pData );
		pData = NULL;
	}
	GDALClose( hSrcDataset );
}

void findRight(const char* szFilename, Point2f &pt, int tile_size = 1024, int null_value = 0)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDatasetH hSrcDs = GDALOpen(szFilename, GA_ReadOnly );

	CPLErr eErr = CE_None;
	int nRasterXSizeRead,nRasterYSizeRead;

	GDALDataset  *hSrcDataset = (GDALDataset  *)hSrcDs;

	nRasterXSizeRead=hSrcDataset->GetRasterXSize();
	nRasterYSizeRead=hSrcDataset->GetRasterYSize();
	int nBands = hSrcDataset->GetRasterCount();
	GDALDataType eDT = hSrcDataset->GetRasterBand(1)->GetRasterDataType();

	const int nSamples	= tile_size;

	int nTiles = (int)ceil(nRasterXSizeRead / (double)nSamples);

	GByte *pData;
	for (int i=nTiles-1;i>=0;--i)
	{
		int BufWidth = nSamples;
		//列末尾小块处理
		if (i == nTiles-1)
		{
			BufWidth = nRasterXSizeRead - (nTiles-1) * nSamples;
			BufWidth = min(BufWidth, nSamples);
		}

		pData=(GByte *) CPLMalloc (BufWidth * nRasterYSizeRead * nBands * GDALGetDataTypeSize( eDT ) / 8);
		eErr = hSrcDataset->RasterIO( GF_Read, i*nSamples, 0, BufWidth , nRasterYSizeRead, pData, BufWidth, nRasterYSizeRead, eDT, nBands, 0, 0 ,0, 0);

		for (int jj = BufWidth-1;jj >= 0;--jj)
		{
			for (int ii = 0; ii < nRasterYSizeRead;++ii)
			{
				for (int iBand = 0;iBand < nBands;++iBand)
				{
					//int pos = ii * BufWidth * nBands + jj * nBands + iBand;
					int pos = BufWidth*nRasterYSizeRead*iBand + ii * BufWidth + jj;
					if(SRCVAL(pData, eDT, pos) != null_value)
					{
						pt.x = i * nSamples + jj;
						pt.y = ii;
						CPLFree( pData );
						pData = NULL;
						GDALClose( hSrcDataset );
						return;
					}
				}
			}
		}
		CPLFree( pData );
		pData = NULL;
	}
	GDALClose( hSrcDataset );
}

void findLeft(const char* szFilename, Point2f &pt, int tile_size = 1024, int null_value = 0)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDatasetH hSrcDs = GDALOpen(szFilename, GA_ReadOnly );

	CPLErr eErr = CE_None;
	int nRasterXSizeRead,nRasterYSizeRead;

	GDALDataset  *hSrcDataset = (GDALDataset  *)hSrcDs;

	nRasterXSizeRead=hSrcDataset->GetRasterXSize();
	nRasterYSizeRead=hSrcDataset->GetRasterYSize();
	int nBands = hSrcDataset->GetRasterCount();
	GDALDataType eDT = hSrcDataset->GetRasterBand(1)->GetRasterDataType();

	const int nSamples	= tile_size;

	int nTiles = (int)ceil(nRasterXSizeRead / (double)nSamples);

	GByte *pData;
	for (int i=0;i<nTiles;++i)
	{
		int BufWidth = nSamples;
		//列末尾小块处理
		if (i == nTiles-1)
		{
			BufWidth = nRasterXSizeRead - (nTiles-1) * nSamples;
			BufWidth = min(BufWidth, nSamples);
		}

		pData=(GByte *) CPLMalloc (BufWidth * nRasterYSizeRead * nBands * GDALGetDataTypeSize( eDT ) / 8);
		eErr = hSrcDataset->RasterIO( GF_Read, i*nSamples, 0, BufWidth , nRasterYSizeRead, pData, BufWidth, nRasterYSizeRead, eDT, nBands, 0, 0 ,0, 0);

		for (int jj = 0;jj < BufWidth;++jj)
		{
			for (int ii = nRasterYSizeRead-1; ii >= 0;--ii)
			{
				for (int iBand = 0;iBand < nBands;++iBand)
				{
					int pos = BufWidth*nRasterYSizeRead*iBand + ii * BufWidth + jj;
					if(SRCVAL(pData, eDT, pos) != null_value)
					{
						pt.x = i * nSamples + jj;
						pt.y = ii;
						CPLFree( pData );
						pData = NULL;
						GDALClose( hSrcDataset );
						return;
					}
				}
			}
		}
		CPLFree( pData );
		pData = NULL;
	}
	GDALClose( hSrcDataset );
}


void findBottom(const char* szFilename, Point2f &pt, int tile_size = 1024, int null_value = 0)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDatasetH hSrcDs = GDALOpen(szFilename, GA_ReadOnly );

	CPLErr eErr = CE_None;
	int nRasterXSizeRead,nRasterYSizeRead;

	GDALDataset  *hSrcDataset = (GDALDataset  *)hSrcDs;

	nRasterXSizeRead=hSrcDataset->GetRasterXSize();
	nRasterYSizeRead=hSrcDataset->GetRasterYSize();
	int nBands = hSrcDataset->GetRasterCount();
	GDALDataType eDT = hSrcDataset->GetRasterBand(1)->GetRasterDataType();

	const int nLines	= tile_size;

	int nTiles = (int)ceil(nRasterYSizeRead / (double)nLines);

	GByte *pData;
	for (int i=nTiles-1;i>=0;--i)
	{
		int BufHeight = nLines;
		//列末尾小块处理
		if (i == nTiles-1)
		{
			BufHeight = nRasterYSizeRead - (nTiles-1) * nLines;
			BufHeight = min(BufHeight, nLines);
		}

		pData=(GByte *) CPLMalloc (nRasterXSizeRead * BufHeight * nBands * GDALGetDataTypeSize( eDT ) / 8);
		eErr = hSrcDataset->RasterIO( GF_Read, 0, i*nLines, nRasterXSizeRead , BufHeight, pData, nRasterXSizeRead, BufHeight, eDT, nBands, 0, 0 ,0, 0);

		for (int ii = BufHeight-1; ii >= 0;--ii)
		{
			for (int jj = 0;jj < nRasterXSizeRead;++jj)
			{
				for (int iBand = 0;iBand < nBands;++iBand)
				{
					//int pos = ii * BufWidth * nBands + jj * nBands + iBand;
					int pos = nRasterXSizeRead*BufHeight*iBand + ii * nRasterXSizeRead + jj;
					if(SRCVAL(pData, eDT, pos) != null_value)
					{
						pt.x = jj;
						pt.y = i * nLines + ii;
						CPLFree( pData );
						pData = NULL;
						GDALClose( hSrcDataset );
						return;
					}
				}
			}
		}
		CPLFree( pData );
		pData = NULL;
	}
	GDALClose( hSrcDataset );
}

static OpenThreads::Barrier bar;
static int GLOBAL_NUM_THREADS;

class processThread: public OpenThreads::Thread
{
public:
	processThread()
		: OpenThreads::Thread(){};
	void setData(EdgeType edgeType, const char* szFilename, int tile_size = 1024, int null_value = 0)
	{
		m_EdgeType = edgeType;
		m_szFilename = szFilename;
		m_TileSsize = tile_size;
		m_NullValue = null_value;
	};
	virtual  ~processThread()
	{
		//m_pRegistration = NULL;
	};
	virtual void run()
	{
		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.
		if (EdgeType::top == m_EdgeType)
		{
			findTop(m_szFilename, m_EdgePoint, m_TileSsize, m_NullValue);
		}
		else if (EdgeType::right == m_EdgeType)
		{
			findRight(m_szFilename, m_EdgePoint, m_TileSsize, m_NullValue);
		}
		else if(EdgeType::bottom == m_EdgeType)
		{
			findBottom(m_szFilename, m_EdgePoint, m_TileSsize, m_NullValue);
		}
		else if(EdgeType::left == m_EdgeType)
		{
			findLeft(m_szFilename, m_EdgePoint, m_TileSsize, m_NullValue);
		}

		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.

		//---------------------------------------------------------------------
		// Now that we've done our work, wait for a sign that we should quit.
		//
		while (true) {

			_quitmutex.lock();
			if(_quitflag == true) break;
			_quitmutex.unlock();

			OpenThreads::Thread::YieldCurrentThread();
		}
	};
	void quit() {
		_quitmutex.lock();
		_quitflag = true;
		_quitmutex.unlock();
	};

	void setId(const int &id)
	{
		myId = id;
	}
	Point2f getEdgePoint(){return m_EdgePoint;};
private:
	int myId;
	EdgeType m_EdgeType;
	Point2f m_EdgePoint;
	const char* m_szFilename;
	int m_TileSsize;
	int m_NullValue;
	volatile bool _quitflag;
	OpenThreads::Mutex _quitmutex;
};

Point2f linesample2Wolrd(const Point2f linesample, double *adfThisGeoTransform, OGRCoordinateTransformation *poCT)
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
	return latlng;
}

vector<Point2f> getClipedBoundary(const char* szFilename, int null_value = 0, int tile_size = 1024)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDatasetH hSrcDs = GDALOpen(szFilename, GA_ReadOnly );
	GDALDataset  *hSrcDataset = (GDALDataset  *)hSrcDs;
	
	Point2f ul;
	Point2f ur;
	Point2f lr;
	Point2f ll;/*
			   #pragma omp parallel sections
			   {
			   #pragma omp section
			   {
			   findTop(hSrcDs, ul, tile_size, null_value);
			   }
			   #pragma omp section
			   {
			   findRight(hSrcDs, ur, tile_size, null_value);
			   }
			   #pragma omp section
			   {
			   findBottom(hSrcDs, lr, tile_size, null_value);
			   }
			   #pragma omp section
			   {
			   findLeft(hSrcDs, ll, tile_size, null_value);
			   }
			   }*/

	GLOBAL_NUM_THREADS = 4 + 1;
	OpenThreads::Thread::SetConcurrency(4);
	OpenThreads::Thread::Init();

	processThread * topThread = new processThread;
	topThread->setData(EdgeType::top, szFilename, tile_size, null_value);
	topThread->start();

	processThread * rightThread = new processThread;
	rightThread->setData(EdgeType::right, szFilename, tile_size, null_value);
	rightThread->start();

	processThread * bottomThread = new processThread;
	bottomThread->setData(EdgeType::bottom, szFilename, tile_size, null_value);
	bottomThread->start();

	processThread * leftThread = new processThread;
	leftThread->setData(EdgeType::left, szFilename, tile_size, null_value);
	leftThread->start();

	bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished
	ul = topThread->getEdgePoint();
	ur = rightThread->getEdgePoint();
	lr = bottomThread->getEdgePoint();
	ll = leftThread->getEdgePoint();
	//findTop(szFilename, ul, tile_size, null_value);
	//findRight(szFilename, ur, tile_size, null_value);
	//findBottom(szFilename, lr, tile_size, null_value);
	//findLeft(szFilename, ll, tile_size, null_value);


	double adfThisGeoTransform[6];
	GDALGetGeoTransform( hSrcDs, adfThisGeoTransform );
	char* pszProjection = (char *)GDALGetProjectionRef( hSrcDs );


	OGRSpatialReference geoSRS;
	OGRSpatialReference prjSRS;
	prjSRS.importFromWkt(&pszProjection);
	geoSRS.SetWellKnownGeogCS( "WGS84" );

	OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation( &prjSRS,
		&geoSRS );

	vector<Point2f> boundary(4);

	boundary[0] = linesample2Wolrd(ul, adfThisGeoTransform, poCT);
	boundary[1] = linesample2Wolrd(ur, adfThisGeoTransform, poCT);
	boundary[2] = linesample2Wolrd(lr, adfThisGeoTransform, poCT);
	boundary[3] = linesample2Wolrd(ll, adfThisGeoTransform, poCT);
	
	GDALClose( hSrcDataset );
	return boundary;
}

bool CreateClipedSpatialIndexOfDirectory(QStringList imageFiles, QString shpFile, int null_value = 0, int tile_size = 1024)
{
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动
	const char *pszDriverName = NULL;
	if (QFileInfo(shpFile).completeSuffix().toLower().compare("shp"))
	{
		pszDriverName = "ESRI Shapefile";
	}
	else if (QFileInfo(shpFile).completeSuffix().toLower().compare("kml"))
	{
		pszDriverName = "KML";
	}
	else
	{
		pszDriverName = "ESRI Shapefile";
	}
	OGRSFDriver *poDriver;
	poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
		pszDriverName );
	if( poDriver == NULL )
	{
		return false;//创建文件驱动失败
	}

	if (QFileInfo(shpFile).exists())
	{
		poDriver->DeleteDataSource(shpFile.toLatin1());
	}
	OGRDataSource *poDS;
	poDS = poDriver->CreateDataSource( shpFile.toLatin1(), NULL );
	if( poDS == NULL )
	{
		return false;//创建源失败
	}

	OGRSpatialReference geoSRS;
	geoSRS.SetWellKnownGeogCS( "WGS84" );
	
	// 多边形图形
	OGRLayer *poLayer;
	poLayer = poDS->CreateLayer( "spatial_index", &geoSRS, wkbPolygon, NULL );
	if( poLayer == NULL )
	{
		return false;//创建图层对象失败
	} 

	OGRFieldDefn oField( "id", OFTString );
	oField.SetWidth(32);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		return false;
	}

	//OGRFieldDefn pathField( "path", OFTString );
	//pathField.SetWidth(254);
	//if( poLayer->CreateField( &pathField ) != OGRERR_NONE )
	//{
	//	return false;
	//}

	OGRFieldDefn fileField( "file", OFTString );
	fileField.SetWidth(254);
	if( poLayer->CreateField( &fileField ) != OGRERR_NONE )
	{
		return false;
	}

	int percent = 0;
	//设置光标位置
	HANDLE hOut; 
	hOut = GetStdHandle(STD_OUTPUT_HANDLE); 
	//得到当前光标位置
	COORD pos= GetConsoleCursorPosition(hOut);
	SetConsoleCursorPosition(hOut, pos);
	cout<<percent<<"%";

	for(int i = 0;i < imageFiles.size();i++)
	{
		vector<Point2f> boundary;
		QString boundaryFile = QBeforeLast(imageFiles[i], '.') + ".txt";
		if (QFileInfo(boundaryFile).exists())
		{
			fstream ifs;
			ifs.open(boundaryFile.toLatin1(), ios_base::in);
			double x, y;
			while (ifs>>x>>y)
			{
				boundary.push_back(Point2f(x,y));
			}
		}
		
		if (boundary.size() < 4)
		{
			boundary = getClipedBoundary(imageFiles[i].toLatin1(), null_value, tile_size);
		}

		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		QString strId;
		QTextStream(&strId) <<i+1;
		poFeature->SetField( "id", (const char *)strId.toLatin1() );
		//poFeature->SetField( "path", (const char *)imageFiles[i].toLatin1() );
		poFeature->SetField( "file", (const char *)QAfterLast(imageFiles[i], '\\').toLatin1() );

		OGRPolygon polygon;
		OGRLinearRing *pOGRLinearRing = new OGRLinearRing();
		int nPoint = (int)boundary.size();
		OGRRawPoint* pointsList = new OGRRawPoint[nPoint];
		for(int j = 0;j < nPoint;j++)
		{
			pointsList[j].x = boundary[j].x;
			pointsList[j].y = boundary[j].y;
		}
		//// 闭合
		//pointsList[nPoint].x = boundary[0].x;
		//pointsList[nPoint].y = boundary[0].y;

		pOGRLinearRing->setNumPoints(nPoint);
		pOGRLinearRing->setPoints(nPoint, pointsList);
		pOGRLinearRing->closeRings();
		polygon.addRing(pOGRLinearRing);

		poFeature->SetGeometry( &polygon ); 

		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			return false;
		}
		delete []pointsList;
		OGRFeature::DestroyFeature( poFeature );

		int tmpPercent = (int)((i + 1) / (double)(imageFiles.size()) * 100 + 0.5);
		if(tmpPercent > percent)
		{
			percent = tmpPercent;
			// 在保留光标位置输出，覆盖原输出内容
			SetConsoleCursorPosition(hOut, pos);
			cout<<percent<<"%";
			fflush( stdout );
		}
	}
	// 在保留光标位置输出，覆盖原输出内容
	SetConsoleCursorPosition(hOut, pos);
	cout<<"100%";

	OGRDataSource::DestroyDataSource( poDS );
	return true;
}

void getClipedBoundary2Txt(QString imageFile, int null_value = 0, int tile_size = 1024)
{
	QString outFile = QBeforeLast(imageFile, '.') + ".txt";
	fstream ofs;
	// 占位
	if (QFileInfo(outFile).exists())
	{
		return;
	}
	ofs.open(outFile.toLatin1(), ios_base::out);
	ofs.close();
	vector<Point2f> boundary = getClipedBoundary(imageFile.toLatin1(), null_value, tile_size);
	ofs.open(outFile.toLatin1(), ios_base::out);
	for (int i=0;i < (int)boundary.size();++i)
	{
		if (i > 0)
		{
			ofs<<endl;
		}
		ofs<<boundary[i].x<<"\t"<<boundary[i].y;
	}
	ofs.close();
}

void CreateClipedSpatialIndexOfDirectory(QString strPath, QString shpFile, QStringList filters = QStringList("*.tif"), int null_value = 0, int tile_size = 1024)
{
	QStringList imageFiles;
	QFindFile(strPath, filters, imageFiles);
	CreateClipedSpatialIndexOfDirectory(imageFiles, shpFile, null_value, tile_size);
}

void CreateClipedBoundaryOfDirectory(QString strPath, QStringList filters = QStringList("*.tif"), int null_value = 0, int tile_size = 1024)
{
	QStringList imageFiles;
	QFindFile(strPath, filters, imageFiles);

	int percent = 0;
	//设置光标位置
	HANDLE hOut; 
	hOut = GetStdHandle(STD_OUTPUT_HANDLE); 
	//得到当前光标位置
	COORD pos= GetConsoleCursorPosition(hOut);
	SetConsoleCursorPosition(hOut, pos);
	cout<<percent<<"%";

	for(int i = 0;i < imageFiles.size();i++)
	{
		getClipedBoundary2Txt(imageFiles[i], null_value, tile_size);

		int tmpPercent = (int)((i + 1) / (double)(imageFiles.size()) * 100 + 0.5);
		if(tmpPercent > percent)
		{
			percent = tmpPercent;
			// 在保留光标位置输出，覆盖原输出内容
			SetConsoleCursorPosition(hOut, pos);
			cout<<percent<<"%";
			fflush( stdout );
		}
	}
	// 在保留光标位置输出，覆盖原输出内容
	SetConsoleCursorPosition(hOut, pos);
	cout<<"100%";
}


vector<Point2f> getBoundary(const char* szFilename, int null_value = 0, int tile_size = 1024)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDatasetH hSrcDs = GDALOpen(szFilename, GA_ReadOnly );
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

	vector<Point2f> boundary(4);

	boundary[0] = linesample2Wolrd(ul, adfThisGeoTransform, poCT);
	boundary[1] = linesample2Wolrd(ur, adfThisGeoTransform, poCT);
	boundary[2] = linesample2Wolrd(lr, adfThisGeoTransform, poCT);
	boundary[3] = linesample2Wolrd(ll, adfThisGeoTransform, poCT);
	
	GDALClose( hSrcDataset );
	return boundary;
}

bool CreateSpatialIndex(QStringList imageFiles, QString shpFile, int null_value = 0, int tile_size = 1024)
{
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriver *poDriver;
	poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
		pszDriverName );
	if( poDriver == NULL )
	{
		return false;//创建文件驱动失败
	}

	if (QFileInfo(shpFile).exists())
	{
		poDriver->DeleteDataSource(shpFile.toLatin1());
	}
	OGRDataSource *poDS;
	poDS = poDriver->CreateDataSource( shpFile.toLatin1(), NULL );
	if( poDS == NULL )
	{
		return false;//创建源失败
	}

	OGRSpatialReference geoSRS;
	geoSRS.SetWellKnownGeogCS( "WGS84" );

	// 多边形图形
	OGRLayer *poLayer;
	poLayer = poDS->CreateLayer( "spatial_index", &geoSRS, wkbPolygon, NULL );
	if( poLayer == NULL )
	{
		return false;//创建图层对象失败
	} 

	OGRFieldDefn oField( "id", OFTString );
	oField.SetWidth(32);
	if( poLayer->CreateField( &oField ) != OGRERR_NONE )
	{
		return false;
	}

	//OGRFieldDefn pathField( "path", OFTString );
	//pathField.SetWidth(254);
	//if( poLayer->CreateField( &pathField ) != OGRERR_NONE )
	//{
	//	return false;
	//}

	OGRFieldDefn fileField( "file", OFTString );
	fileField.SetWidth(254);
	if( poLayer->CreateField( &fileField ) != OGRERR_NONE )
	{
		return false;
	}

	int percent = 0;
	//设置光标位置
	HANDLE hOut; 
	hOut = GetStdHandle(STD_OUTPUT_HANDLE); 
	//得到当前光标位置
	COORD pos= GetConsoleCursorPosition(hOut);
	SetConsoleCursorPosition(hOut, pos);
	cout<<percent<<"%";

	for(int i = 0;i < imageFiles.size();i++)
	{
		vector<Point2f> boundary;
		QString boundaryFile = QBeforeLast(imageFiles[i], '.') + ".txt";
		if (QFileInfo(boundaryFile).exists())
		{
			fstream ifs;
			ifs.open(boundaryFile.toLatin1(), ios_base::in);
			double x, y;
			while (ifs>>x>>y)
			{
				boundary.push_back(Point2f(x,y));
			}
		}

		if (boundary.size() < 4)
		{
			boundary = getBoundary(imageFiles[i].toLatin1(), null_value, tile_size);
		}

		OGRFeature *poFeature;
		poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
		QString strId;
		QTextStream(&strId) <<i+1;
		poFeature->SetField( "id", (const char *)strId.toLatin1() );
		//poFeature->SetField( "path", (const char *)imageFiles[i].toLatin1() );
		poFeature->SetField( "file", (const char *)QAfterLast(QDir::toNativeSeparators(imageFiles[i]), '\\').toLatin1() );

		OGRPolygon polygon;
		OGRLinearRing *pOGRLinearRing = new OGRLinearRing();
		int nPoint = (int)boundary.size();
		OGRRawPoint* pointsList = new OGRRawPoint[nPoint];
		for(int j = 0;j < nPoint;j++)
		{
			pointsList[j].x = boundary[j].x;
			pointsList[j].y = boundary[j].y;
		}
		//// 闭合
		//pointsList[nPoint].x = boundary[0].x;
		//pointsList[nPoint].y = boundary[0].y;

		pOGRLinearRing->setNumPoints(nPoint);
		pOGRLinearRing->setPoints(nPoint, pointsList);
		pOGRLinearRing->closeRings();
		polygon.addRing(pOGRLinearRing);

		poFeature->SetGeometry( &polygon ); 

		if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
		{
			return false;
		}
		delete []pointsList;
		OGRFeature::DestroyFeature( poFeature );

		int tmpPercent = (int)((i + 1) / (double)(imageFiles.size()) * 100 + 0.5);
		if(tmpPercent > percent)
		{
			percent = tmpPercent;
			// 在保留光标位置输出，覆盖原输出内容
			SetConsoleCursorPosition(hOut, pos);
			cout<<percent<<"%";
			fflush( stdout );
		}
	}
	// 在保留光标位置输出，覆盖原输出内容
	SetConsoleCursorPosition(hOut, pos);
	cout<<"100%";

	OGRDataSource::DestroyDataSource( poDS );
	return true;
}

void CreateSpatialIndexOfDirectory(QString strPath, QString shpFile, QStringList filters = QStringList("*.tif"), int null_value = 0, int tile_size = 1024)
{
	QStringList imageFiles;
	QFindFile(strPath, filters, imageFiles);
	CreateSpatialIndex(imageFiles, shpFile, null_value, tile_size);
}

}

/************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: GetRasterBoundary  [-f filter] [-o outfile] [-nv nullvalue] [-c] [-h]\n"
		"-f filter: filter.\n"
		"-shp shapefile: output shapefile name. \n    If it is not provided, only the true boundaries would be calculated.\n"
		"-c: enable true boundary calculation.\n"
		"-h: show help information.\n");

	if( pszErrorMsg != NULL )
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

int main( int argc, char** argv )
{
	//CreateSpatialIndexOfDirectory("I:\\testdata\\2000_UTM_std", "I:\\testdata\\2000_UTM_std\\si.shp", 0, 1024);
	//CreateSpatialIndexOfDirectory("D:\\workspace\\Landsat", "D:\\workspace\\Landsat\\si.shp", 0, 1024);
	//CreateSpatialIndexOfDirectory("D:\\workspace\\HJ\\db1", "D:\\workspace\\HJ\\db1\\si.shp", 0, 1024);
	//CreateSpatialIndexOfDirectory("D:\\workspace\\HJ\\db1", "D:\\workspace\\HJ\\db1\\si.shp", QStringList("*.tif"), 0, 1024);
	//CreateBoundaryOfDirectory("D:\\workspace\\HJ\\db1", QStringList("*.tif"), 0, 1024);

	const char *pszFilter = "*.tif";
	const char *pszShapeFile = NULL;
	bool bClip = false;
	int null_value = 0;

	if (argc == 0)
	{
		rasterBoundary::CreateClipedBoundaryOfDirectory("I:\\20140523\\fusion", QStringList("*.tif"), 0, 1024);
	}
	if (argc == 1)
	{
		rasterBoundary::CreateClipedBoundaryOfDirectory("I:\\20140523\\fusion", QStringList("*PAN*.tif"), 0, 1024);
		//CreateClipedBoundaryOfDirectory("I:\\test", QStringList("*.tif"), 0, 1024);
		//CreateSpatialIndexOfDirectory("I:\\test", "index.shp", QStringList(pszFilter), null_value, 1024);
	}

	if (argc > 1)
	{
		/* -------------------------------------------------------------------- */
		/*      Parse arguments.                                                */
		/* -------------------------------------------------------------------- */
		for( int i = 1; i < argc; i++ )
		{
			if( 0 == _stricmp(argv[i],"-filter")
				|| 0 == _stricmp(argv[i],"-f"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszFilter = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-o") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszShapeFile = argv[++i] ;
			}
			else if(0 == _stricmp(argv[i],"-nv") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				null_value = CPLAtofM( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-c") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
				bClip = true;
			}
			else if(0 == _stricmp(argv[i],"-h") )
			{
				Usage(0);
			}
			else
			{
				Usage();
			}
		}

		if (!bClip)
		{
			//CreateBoundaryOfDirectory("", QStringList(pszFilter), null_value, 1024);
			rasterBoundary::CreateSpatialIndexOfDirectory("", pszShapeFile, QStringList(pszFilter), null_value, 1024);
		}
		else if(!pszShapeFile)
		{
			rasterBoundary::CreateClipedBoundaryOfDirectory("", QStringList(pszFilter), null_value, 1024);
		}
		else{
			rasterBoundary::CreateClipedSpatialIndexOfDirectory("", pszShapeFile, QStringList(pszFilter), null_value, 1024);
		}
	}
	return 0;
}