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
//#include "opencv2/nonfree/nonfree.hpp"
#include "opencv2/imgproc/imgproc_c.h"
//#include "opencv2/legacy/legacy.hpp"
//#include "opencv2/legacy/compat.hpp"
#include "opencv2/imgproc/imgproc.hpp"

// GDAL includes
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"

// ossim
#include <ossim/base/ossimTDpt.h>
#include <ossim/base/ossimDrect.h>
#include <ossim/base/ossimIrect.h>

#include <armadillo>
#include "pca.h"

#include <iostream>
#include <wtypes.h>
#include <vector>
using namespace std;

namespace mylib{
#define MAXVALUE 65535

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

	bool getRect2CvMat(const ossimIrect &rect, cv::Mat& outMat, int band);
	bool getRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, int band, double scale = 1.0, double stretchRatio = 0.01);
	bool getRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, int band, ossimDpt scale = ossimDpt(1.0, 1.0),
		double stretchRatio = 0.01, bool bStretch = true);
	bool getRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, vector<unsigned int> band, ossimDpt scale = ossimDpt(1.0, 1.0),
		double stretchRatio = 0.01, bool bStretch = true);
	bool getCombinedRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, vector<unsigned int> band, vector<double> weights, ossimDpt scale = ossimDpt(1.0, 1.0),
		double stretchRatio = 0.01, bool bStretch = true);
	bool getPrincipalRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, double scale = 1.0);

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
	ossimIrect getBoundary(){return m_Boundary;};

	void setTileWidth(int iTileWidth);
	void setTileHeight(int iTileHeight);
	void setNullValue(double null_value){
		m_bIgnorNull = true;
		m_nullValue = null_value;};
	void disableNullValue() {m_bIgnorNull = false;};

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

	ossimIrect getBoundary()const{ return m_Boundary;};

	GdalRasterApp::RasterBuf* getBufInfo(){return m_CacheList.back();};
	std::vector<GdalRasterApp::RasterBuf*> getCacheInfoList(){return m_CacheList;};

	template<class T1, class T2>
	void linesample2world(T1 linesample, T2 &world)const 
	{
		world.x = m_adfGeoTransform[0] + linesample.x*m_adfGeoTransform[1] + linesample.y*m_adfGeoTransform[2];
		world.y = m_adfGeoTransform[3] + linesample.x*m_adfGeoTransform[4] + linesample.y*m_adfGeoTransform[5];
	};

	template<class T1, class T2>
	void world2linesample(T1 world, T2 &linesample)const 
	{
		double m = m_adfGeoTransform[1]*m_adfGeoTransform[5] - m_adfGeoTransform[2]*m_adfGeoTransform[4];
		linesample.x = (m_adfGeoTransform[5]*world.x - m_adfGeoTransform[2]*world.y + m_adfGeoTransform[2]*m_adfGeoTransform[3]
		-m_adfGeoTransform[0]*m_adfGeoTransform[5]) / (m + 10e-20);
		linesample.y = (m_adfGeoTransform[4]*world.x - m_adfGeoTransform[1]*world.y + m_adfGeoTransform[1]*m_adfGeoTransform[3]
		-m_adfGeoTransform[0]*m_adfGeoTransform[4]) / (-m + 10e-20);
	};

	template<class T>
	void world2lonlat(T &pt)const 
	{
		if( m_poWorld2Lonlat == NULL || !m_poWorld2Lonlat->Transform( 1, &pt.x, &pt.y ) )
		{
			printf( "Transformation failed.\n" );
		}
	}

	template<class T>
	void lonlat2world(T &pt)const
	{
		if( m_poLonlat2World == NULL || !m_poLonlat2World->Transform( 1, &pt.x, &pt.y ) )
		{
			printf( "Transformation failed.\n" );
		}
	}

	template<class T1, class T2>
	void linesample2lonlat(T1 linesample, T2 &lonlat)const 
	{
		linesample2world(linesample, lonlat);
		world2lonlat(lonlat);
	};

	template<class T1, class T2>
	void lonlat2linesample(T1 lonlat, T2 &linesample)const 
	{
		lonlat2world(lonlat);
		world2linesample(lonlat, linesample);
	};

	template<class T>
	cv::Mat gdalData2FloatcvMat(int nw, int nh, void* pRectData)
	{
		cv::Mat outMat = cv::Mat(cv::Size(nw, nh), CV_32FC1);
		//cv::Mat outMat = cv::Mat(cv::Size(nw, nh), CV_8UC1);
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0; j < nw; ++j)
			{
				T data = ((T*)pRectData)[i*nw + j];
				//outMat.data[i*nw + j] = static_cast<float>(data);
				//outMat.data[i*nw + j] = static_cast<GByte>(data);
				outMat.at<float>(i, j) = static_cast<float>(data);
			}
		}
		return outMat;
	}

	template<class T>
	cv::Mat gdalData2DoublecvMat(int nw, int nh, void* pRectData)
	{
		cv::Mat outMat = cv::Mat(cv::Size(nw, nh), CV_64FC1);
		//cv::Mat outMat = cv::Mat(cv::Size(nw, nh), CV_8UC1);
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0; j < nw; ++j)
			{
				T data = ((T*)pRectData)[i*nw + j];
				outMat.data[i*nw + j] = static_cast<double>(data);
				//outMat.data[i*nw + j] = static_cast<GByte>(data);
			}
		}
		return outMat;
	}

	template<class T>
	void stretchRect2Byte(int nw, int nh, void* pRectData, GByte* &pByteData, double stretchRatio = 0.02)
	{
		//if (pByteData)
		//{
		//	delete []pByteData;
		//	pByteData = NULL;
		//}
		//pByteData = new GByte [nw * nh];

		T vmin = 1e10;
		T vmax = -1e10;
		// find min and max values
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0;j < nw;++j)
			{
				T data = ((T*)pRectData)[i*nw + j];
				if (m_bIgnorNull && data == m_nullValue)
				{
					continue;
				}
				if (data < vmin)
				{
					vmin = data;
				}
				if (data > vmax)
				{
					vmax = data;
				}
			}
		}

		//const int nbins = 65536;
		const int nbins = 256;
		double bin_step = (vmax - vmin) / (double)nbins;
		//bin_step = 1.0;

		// histogram statistics
		int histogram[nbins] = {0};  // the rect is expected to be small, thus int would be enough
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0;j < nw;++j)
			{
				T data = ((T*)pRectData)[i*nw + j];
				if (m_bIgnorNull && data == m_nullValue)
				{
					continue;
				}
				// scale to bins
				int v = int((data - vmin)/bin_step+0.5);
				v = v > (nbins-1) ? (nbins-1) : v;
				v = v < 0 ? 0 : v;
				histogram[v]++;
			}
		}

		// find the stretch percent position
		int stretch_count = int(stretchRatio * nh * nw + 0.5);
		int low_pos = 0;
		int low_count = 0;
		for (int i = 0;i < nbins;++i)
		{
			low_count += histogram[i];
			if (low_count > stretch_count)
			{
				low_pos = i;
				break;
			}
		}
		int high_pos = nbins - 1;
		int high_count = 0;
		for (int i = nbins - 1;i > 0;--i)
		{
			high_count += histogram[i];
			if (high_count > stretch_count)
			{
				high_pos = i;
				break;
			}
		}

		// high and low to data
		if (low_pos == high_pos)
		{
			if (low_pos > 0)
			{
				low_pos--;
			}
			else
			{
				high_pos++;
			}
		}
		T vlow = low_pos*bin_step + vmin;
		T vhigh = high_pos*bin_step + vmin;
		double vlen = vhigh - vlow;
		if (vlen < 1e-6)
		{
			vlen += 1e-6;
		}

		double g = 0.8;

		// stretch to byte
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0;j < nw;++j)
			{
				T data = ((T*)pRectData)[i*nw + j];
				if (m_bIgnorNull && data == m_nullValue)
				{
					pByteData[i*nw + j] = m_nullValue;
					continue;
				}
				int oData = (int)(pow((data - vlow) / vlen, g) * 255.0 + 0.5);
				oData = oData > 255 ? 255 : oData;
				oData = oData < 0 ? 0 : oData;
				pByteData[i*nw + j] = (GByte)oData;
			}
		}
	}


	template<class T>
	void stretchRect2Byte(int nw, int nh, void* pRectData, arma::mat &X, int iCol, double stretchRatio = 0.02)
	{
		//if (pByteData)
		//{
		//	delete []pByteData;
		//	pByteData = NULL;
		//}
		//pByteData = new GByte [nw * nh];

		T vmin = -1e10;
		T vmax = 1e10;
		// find min and max values
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0;j < nw;++j)
			{
				T data = ((T*)pRectData)[i*nw + j];
				if (data < vmin)
				{
					vmin = data;
				}
				if (data > vmax)
				{
					vmax = data;
				}
			}
		}

		const int nbins = 65536;
		//const int nbins = 256;
		double bin_step = (vmax - vmin) / (double)nbins;
		//bin_step = 1.0;

		// histogram statistics
		int histogram[nbins] = {0};  // the rect is expected to be small, thus int would be enough
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0;j < nw;++j)
			{
				T data = ((T*)pRectData)[i*nw + j];
				// scale to bins
				int v = int((data - vmin)/bin_step+0.5);
				v = v > (nbins-1) ? (nbins-1) : v;
				v = v < 0 ? 0 : v;
				histogram[v]++;
			}
		}

		// find the stretch percent position
		int stretch_count = int(stretchRatio * nh * nw + 0.5);
		int low_pos = 0;
		int low_count = 0;
		for (int i = 0;i < nbins;++i)
		{
			low_count += histogram[i];
			if (low_count > stretch_count)
			{
				low_pos = i;
				break;
			}
		}
		int high_pos = nbins - 1;
		int high_count = 0;
		for (int i = nbins - 1;i > 0;--i)
		{
			high_count += histogram[i];
			if (high_count > stretch_count)
			{
				high_pos = i;
				break;
			}
		}

		// high and low to data
		if (low_pos == high_pos)
		{
			if (low_pos > 0)
			{
				low_pos--;
			}
			else
			{
				high_pos++;
			}
		}
		T vlow = low_pos*bin_step + vmin;
		T vhigh = high_pos*bin_step + vmin;
		double vlen = vhigh - vlow;
		if (vlen < 1e-6)
		{
			vlen += 1e-6;
		}

		// stretch to byte
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0;j < nw;++j)
			{
				T data = ((T*)pRectData)[i*nw + j];
				int oData = (int)((data - vlow) / vlen * 255.0 + 0.5);
				oData = oData > 255 ? 255 : oData;
				oData = oData < 0 ? 0 : oData;
				X(i*nw + j, iCol) = (GByte)oData;
			}
		}
	}
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
	ossimIrect m_Boundary;
	OGRCoordinateTransformation *m_poWorld2Lonlat;
	OGRCoordinateTransformation *m_poLonlat2World;

	bool m_bIgnorNull;
	double m_nullValue;

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
}