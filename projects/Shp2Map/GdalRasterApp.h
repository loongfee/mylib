#pragma once

// opencv
#include "cv.h"
#include "highgui.h"
#include <opencv2/core/core.hpp>
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#include "opencv2/imgproc/imgproc.hpp"

// QT includes
#include <QMessageBox>
#include <QString>
#include <QObject>

// GDAL includes
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"

#include <iostream>
#include <wtypes.h>
#include <vector>
using namespace std;

class GdalRasterApp
{
public:
	struct RasterBuf
	{
		void * data;
		int iBufTileX;
		int iBufTileY;
		int iBufWidth;
		int iBufHeight;
		int iBufOffsetX;
		int iBufOffsetY;
		int iBandCount;
		GDALDataType eDataType;

		RasterBuf()  //默认构造函数
		{
			data = NULL;
			iBufWidth = -1;
			iBufHeight = -1;

			iBufOffsetX = -1;
			iBufOffsetY = -1;

			iBufTileX = -1;
			iBufTileY = -1;

			iBandCount = 0;

			eDataType = GDT_Byte;
		}
     
		 RasterBuf& operator=(const RasterBuf& s)//重载运算符
		 {
			data = s.data;
			iBufWidth = s.iBufWidth;
			iBufHeight = s.iBufHeight;
			iBufOffsetX = s.iBufOffsetX;
			iBufOffsetY = s.iBufOffsetY;
			iBufTileX = s.iBufTileX;
			iBufTileY = s.iBufTileY;
			iBandCount = s.iBandCount;
			eDataType = s.eDataType;
			return *this;
		 }
		RasterBuf(const RasterBuf& s)//复制构造函数
		{
			data = s.data;
			iBufWidth = s.iBufWidth;
			iBufHeight = s.iBufHeight;
			iBufOffsetX = s.iBufOffsetX;
			iBufOffsetY = s.iBufOffsetY;
			iBufTileX = s.iBufTileX;
			iBufTileY = s.iBufTileY;
			iBandCount = s.iBandCount;
			eDataType = s.eDataType;
		}
	};
public:
	GdalRasterApp(void);
	~GdalRasterApp(void);

	bool open(const char* filename);
	bool close();

	double* getGeoTransform(){return m_adfGeoTransform;};
	char* getGetProjectionRef(){return m_pszTargetSRS;};
	GDALDataset* getDataset(){return m_pDataset;};
	GDALDataType getDataType(){return m_eDataType;};
	int width(){return m_iWidth;};
	int height(){return m_iHeight;};
	int nBand(){return m_iBandCount;};

	int getTileWidth(){return m_iTileWidth;};
	int getTileHeight(){return m_iTileHeight;};
	int getTileCountX(){return m_iTileCountX;};
	int getTileCountY(){return m_iTileCountY;};

	void setTileWidth(int iTileWidth);
	void setTileHeight(int iTileHeight);

	//bool getTileData(int iTileX, int iTileY, RasterBuf *Buf,
	//	int nBandCount,
	//	int *panBandMap = 0,
	//	int nPixelSpace = 0,
	//	int nLineSpace = 0,
	//	int nBandSpace = 0);
	GdalRasterApp::RasterBuf * getTileData(int iTileX, int iTileY, 
								int nBandCount,
								int *panBandMap = 0,
								int nPixelSpace = 0,
								int nLineSpace = 0,
								int nBandSpace = 0);

	bool getTileOffset(int iTileX, int iTileY, int &offsetX, int &offsetY);
	bool readTile(int iTileX, int iTileY,
		int nBandCount,
		int *panBandMap = 0,
		int nPixelSpace = 0,
		int nLineSpace = 0,
		int nBandSpace = 0);

	bool readTile(int iTileX, int iTileY,
							 int nBandCount,
							RasterBuf *Buf,
							int *panBandMap = 0,
							int nPixelSpace = 0,
							int nLineSpace = 0,
							int nBandSpace = 0);
	bool checkTileValidity(int iTileX, int iTileY);

	bool getPixel(int x, int y, void** data);

	GdalRasterApp::RasterBuf* getBufInfo(){return m_CacheList.back();};
	std::vector<GdalRasterApp::RasterBuf*> getCacheInfoList(){return m_CacheList;};

	template<class T1, class T2>
	void linesample2world(T1 linesample, T2 &world)
	{
		world.x = m_adfGeoTransform[0] + linesample.x*m_adfGeoTransform[1] + linesample.y*m_adfGeoTransform[2];
		world.y = m_adfGeoTransform[3] + linesample.x*m_adfGeoTransform[4] + linesample.y*m_adfGeoTransform[5];
	};

	template<class T1, class T2>
	void world2linesample(T1 world, T2 &linesample)
	{
		double m = m_adfGeoTransform[1]*m_adfGeoTransform[5] - m_adfGeoTransform[2]*m_adfGeoTransform[4];
		linesample.x = (m_adfGeoTransform[5]*world.x - m_adfGeoTransform[2]*world.y + m_adfGeoTransform[2]*m_adfGeoTransform[3]
		-m_adfGeoTransform[0]*m_adfGeoTransform[5]) / (m + 10e-10);
		linesample.y = (m_adfGeoTransform[4]*world.x - m_adfGeoTransform[1]*world.y + m_adfGeoTransform[1]*m_adfGeoTransform[3]
		-m_adfGeoTransform[0]*m_adfGeoTransform[4]) / (-m + 10e-10);
	};
protected:
	bool m_IsOpened;

	// raster information
	double m_adfGeoTransform[6];
	char *m_pszTargetSRS;
	GDALDataset* m_pDataset;
	GDALDataType m_eDataType;
	int m_iWidth;
	int m_iHeight;
	int m_iBandCount;

	// Tile Information
	int m_iTileWidth;
	int m_iTileHeight;
	int m_iTileCountX;
	int m_iTileCountY;

	// Buf Information
	//RasterBufInfo m_BufInfo;
	//BYTE* pBufData;
	//BYTE* m_pBufData;
	//int m_iBufTileX;
	//int m_iBufTileY;
	//int m_iBufWidth;
	//int m_iBufHeight;
	//int m_iBufOffsetX;
	//int m_iBufOffsetY;
	//bool updateBufInfo(int iTileX, int iTileY);
	bool updateBufInfo(RasterBuf *bufInfo,int iTileX, int iTileY);
	//void releaseBuf();
	void releaseBuf(RasterBuf *cache);
	void releaseCaches();

	std::vector<RasterBuf*> m_CacheList;
	//std::vector<RasterBufInfo> m_CacheInfoList;

	int m_CacheSize;
};