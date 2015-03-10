#include "GdalRasterApp.h"
#include <armadillo>

namespace mylib{
//const int over_lap = 4;
#define HEIGHT_OVER_LAP 0
#define WIDTH_OVER_LAP 0

GdalRasterApp::GdalRasterApp(void)
	:m_poLonlat2World(NULL),
	m_poWorld2Lonlat(NULL),
	m_bIgnorNull(true),
	m_nullValue(0.0)
{
	m_pszTargetSRS = NULL;
	m_pDataset = NULL;
	//m_pBufData = NULL;

	//m_iTileWidth = 512;
	//m_iTileHeight = 512;
	m_iTileWidth = 2048;
	m_iTileHeight = 2048;

	m_iTileCountX = 0;
	m_iTileCountY = 0;

	m_CacheSize = 3;
	m_CacheList.resize(m_CacheSize);
	for (int i = 0; i < m_CacheSize; i++)
	{
		m_CacheList[i] = new RasterBuf;
		//m_CacheList[i]->data = NULL;
	}
}

GdalRasterApp::~GdalRasterApp(void)
{
	releaseCaches();
	//releaseBuf();
}

bool GdalRasterApp::open(const char* filename)
{
	GDALAllRegister();
	m_pDataset	= (GDALDataset*) GDALOpen(filename,  GA_ReadOnly);
	if (! m_pDataset)
	{
		cout<<"打开文件："<<filename<<"失败。"<<endl;
		return false;
	}

	m_iWidth = m_pDataset->GetRasterXSize();
	m_iHeight = m_pDataset->GetRasterYSize();
	m_Boundary = ossimIrect(ossimIpt(0, 0), ossimIpt(m_iWidth-1, m_iHeight-1));
	m_iBandCount = m_pDataset->GetRasterCount();
	m_eDataType = m_pDataset->GetRasterBand(1)->GetRasterDataType();

	m_pDataset->GetGeoTransform( m_adfGeoTransform );			
	m_pszTargetSRS = CPLStrdup(GDALGetProjectionRef(m_pDataset));
	try
	{
		OGRSpatialReference *geoSRS;
		OGRSpatialReference prjSRS;
		prjSRS.importFromWkt(&m_pszTargetSRS);
		geoSRS = prjSRS.CloneGeogCS();
		//geoSRS.SetWellKnownGeogCS( "WGS84" );

		m_poWorld2Lonlat = OGRCreateCoordinateTransformation( &prjSRS, geoSRS );
		m_poLonlat2World = OGRCreateCoordinateTransformation( geoSRS, &prjSRS );
	}
	catch (...)
	{
		std::cerr<<"SRS unknown: "<<filename<<endl;
	}
	m_pszTargetSRS = CPLStrdup(GDALGetProjectionRef(m_pDataset));

	setTileWidth(m_iTileWidth);
	setTileHeight(m_iTileHeight);

	return true;
}

bool GdalRasterApp::close()
{
	releaseCaches();
	if (m_pDataset)
	{
		GDALClose(m_pDataset);
	}
	return true;
}

// Normalizes a given image into a value range between 0 and 255.
cv::Mat norm_0_255(const cv::Mat& src) {
	// Create and return normalized image:
	cv::Mat dst;
	switch (src.channels()) {
	case 1:
		cv::normalize(src, dst, 0, 255, cv::NORM_MINMAX, CV_8UC1);
		break;
	case 3:
		cv::normalize(src, dst, 0, 255, cv::NORM_MINMAX, CV_8UC3);
		break;
	default:
		src.copyTo(dst);
		break;
	}
	return dst;
}

// Converts the images given in src into a row matrix.
cv::Mat asRowMatrix(const vector<cv::Mat>& src, int rtype, double alpha = 1, double beta = 0) {
	// Number of samples:
	size_t n = src.size();
	// Return empty matrix if no matrices given:
	if (n == 0)
		return cv::Mat();
	// dimensionality of (reshaped) samples
	size_t d = src[0].total();
	// Create resulting data matrix:
	cv::Mat data(n, d, rtype);
	// Now copy data:
	for (int i = 0; i < n; i++) {
		//
		if (src[i].empty()) {
			string error_message = cv::format("Image number %d was empty, please check your input data.", i);
			CV_Error(CV_StsBadArg, error_message);
		}
		// Make sure data can be reshaped, throw a meaningful exception if not!
		if (src[i].total() != d) {
			string error_message = cv::format("Wrong number of elements in matrix #%d! Expected %d was %d.", i, d, src[i].total());
			CV_Error(CV_StsBadArg, error_message);
		}
		// Get a hold of the current row:
		cv::Mat xi = data.row(i);
		// Make reshape happy by cloning for non-continuous matrices:
		if (src[i].isContinuous()) {
			src[i].reshape(1, 1).convertTo(xi, rtype, alpha, beta);
		}
		else {
			src[i].clone().reshape(1, 1).convertTo(xi, rtype, alpha, beta);
		}
	}
	return data;
}

bool GdalRasterApp::getRect2CvMat(const ossimIrect &rect, cv::Mat& outMat, int band)
{
	if (!m_Boundary.pointWithin(rect.ul()) || !m_Boundary.pointWithin(rect.lr()))
	{
		return false;
	}
	if (band >= m_iBandCount)
	{
		cerr<<"GdalRasterApp"<<"::getRect2CvMat Warning not enough bands in slave, only "<< m_iBandCount <<endl;
		band = 0;
	}
	GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(band+1);
	outMat = cv::Mat(cv::Size(rect.width(), rect.height()), CV_8UC1);
	//outMat.data = static_cast<uchar*>(imageData->getBuf(0));
//#pragma omp critical
{
	pRasterBand->RasterIO( GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(), 
			outMat.data, rect.width(), rect.height(),  m_eDataType, 0, 0);

}
	return true;
}

bool GdalRasterApp::getPrincipalRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, double scale/* = 1.0*/)
{
	if (!m_Boundary.pointWithin(rect.ul()) || !m_Boundary.pointWithin(rect.lr()))
	{
		return false;
	}
	int nOutWidth = int(rect.width() * scale + 0.5);
	int nOutHeight = int(rect.height() * scale + 0.5);
	int nBands = m_iBandCount;
	if (nBands == 1)
	{
		return getRect2CvMatByte(rect, outMat, 0, scale);
	}
	arma::mat X(nOutWidth * nOutHeight, nBands);
	for (int i = 0;i < nBands;++i)
	{

		GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(i+1);
		int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
		void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
		pRasterBand->RasterIO( GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(), 
			pData, nOutWidth, nOutHeight,  m_eDataType, 0, 0);

		double stretchRatio = 0.01;
		if (m_eDataType == GDT_Byte)
		{
			stretchRect2Byte<GByte>(nOutWidth, nOutHeight, pData, X, i, stretchRatio);
		}
		else if (m_eDataType == GDT_UInt16)
		{
			stretchRect2Byte<GUInt16>(nOutWidth, nOutHeight, pData, X, i, stretchRatio);
		}
		else if (m_eDataType == GDT_Int16)
		{
			stretchRect2Byte<GInt16>(nOutWidth, nOutHeight, pData, X, i, stretchRatio);
		}
		else if (m_eDataType == GDT_UInt32)
		{
			stretchRect2Byte<GUInt32>(nOutWidth, nOutHeight, pData, X, i, stretchRatio);
		}
		else if (m_eDataType == GDT_Int32)
		{
			stretchRect2Byte<GInt32>(nOutWidth, nOutHeight, pData, X, i, stretchRatio);
		}
		else if (m_eDataType == GDT_Float32)
		{
			stretchRect2Byte<float>(nOutWidth, nOutHeight, pData, X, i, stretchRatio);
		}
		else if (m_eDataType == GDT_Float64)
		{
			stretchRect2Byte<double>(nOutWidth, nOutHeight, pData, X, i, stretchRatio);
		}

		delete []pData;
		pData = NULL;
	}

	arma::mat coeff;
	arma::mat score;
	arma::princomp(coeff, score, X);
	outMat = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
	for (int i = 0;i < nOutHeight*nOutWidth;++i)
	{
		outMat.data[i] = score(i, 0);
	}
	return true;
}

bool GdalRasterApp::getRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, int band, double scale/* = 1.0*/,
	double stretchRatio/* = 0.01*/)
{
	ossimIrect clipRect = rect.clipToRect(m_Boundary);
	if (clipRect.width() <= 1 || clipRect.height() <= 1)
	{
		return false;
	}
	//if (!m_Boundary.pointWithin(rect.ul()) || !m_Boundary.pointWithin(rect.lr()))
	//{
	//	return false;
	//}
	if (band >= m_iBandCount)
	{
		//cerr<<"GdalRasterApp"<<"::getRect2CvMat Warning not enough bands in slave, only "<< m_iBandCount <<endl;
		band = 0;
	}
	int nOutWidth = int(rect.width() * scale + 0.5);
	int nOutHeight = int(rect.height() * scale + 0.5);
	GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(band+1);
	int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
	void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
	//outMat.data = static_cast<uchar*>(imageData->getBuf(0));
	//#pragma omp critical
	{
		pRasterBand->RasterIO( GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(), 
			pData, nOutWidth, nOutHeight,  m_eDataType, 0, 0);
	}

	outMat = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
	if (m_eDataType == GDT_Byte)
	{
		stretchRect2Byte<GByte>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_UInt16)
	{
		stretchRect2Byte<GUInt16>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_Int16)
	{
		stretchRect2Byte<GInt16>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_UInt32)
	{
		stretchRect2Byte<GUInt32>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_Int32)
	{
		stretchRect2Byte<GInt32>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_Float32)
	{
		stretchRect2Byte<float>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_Float64)
	{
		stretchRect2Byte<double>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}

	delete []pData;
	pData = NULL;
	return true;
}

bool GdalRasterApp::getRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, int band, ossimDpt scale/* = ossimDpt(1.0, 1.0)*/,
	double stretchRatio/* = 0.01*/, bool bStretch/* = true*/)
{
	ossimIrect clipRect = rect.clipToRect(m_Boundary);
	if (clipRect.width() <= 1 || clipRect.height() <= 1)
	{
		return false;
	}
	//if (!m_Boundary.pointWithin(rect.ul()) || !m_Boundary.pointWithin(rect.lr()))
	//{
	//	return false;
	//}
	if (band >= m_iBandCount)
	{
		//cerr<<"GdalRasterApp"<<"::getRect2CvMat Warning not enough bands in slave, only "<< m_iBandCount <<endl;
		band = 0;
	}
	int nOutWidth = int(rect.width() * scale.x + 0.5);
	int nOutHeight = int(rect.height() * scale.y + 0.5);
	GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(band + 1);
	int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
	void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
	//outMat.data = static_cast<uchar*>(imageData->getBuf(0));
	//#pragma omp critical
	{
		pRasterBand->RasterIO(GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(),
			pData, nOutWidth, nOutHeight, m_eDataType, 0, 0);
	}

	outMat = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
	if (m_eDataType == GDT_Byte)
	{
		if (bStretch)
		{
			stretchRect2Byte<GByte>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else
		{
			memcpy(outMat.data, pData, nOutWidth*nOutHeight);
		}
	}
	else if (m_eDataType == GDT_UInt16)
	{
		stretchRect2Byte<GUInt16>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_Int16)
	{
		stretchRect2Byte<GInt16>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_UInt32)
	{
		stretchRect2Byte<GUInt32>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_Int32)
	{
		stretchRect2Byte<GInt32>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_Float32)
	{
		stretchRect2Byte<float>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}
	else if (m_eDataType == GDT_Float64)
	{
		stretchRect2Byte<double>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
	}

	delete[]pData;
	pData = NULL;
	return true;
}


bool GdalRasterApp::getCombinedRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, vector<unsigned int> bands, vector<double> weights,
	ossimDpt scale/* = ossimDpt(1.0, 1.0)*/,
	double stretchRatio/* = 0.01*/, bool bStretch/* = true*/)
{
	ossimIrect clipRect = rect.clipToRect(m_Boundary);
	if (clipRect.width() <= 1 || clipRect.height() <= 1)
	{
		return false;
	}
	//if (!m_Boundary.pointWithin(rect.ul()) || !m_Boundary.pointWithin(rect.lr()))
	//{
	//	return false;
	//}

	vector<unsigned int> bandList;
	// check bands
	for (size_t i = 0; i < bands.size(); i++)
	{
		int newBand = 0;
		if (bands[i] < m_iBandCount && bands[i] >= 0)
		{
			newBand = bands[i];
		}

		// check exist
		bool band_exist = false;
		for (size_t j = 0; j < bandList.size(); j++)
		{
			if (bandList[j] == newBand)
			{
				band_exist = true;
				break;
			}
		}

		if (!band_exist)
		{
			bandList.push_back(newBand);
		}
	}

	int nbands = (int)bandList.size();
	if (weights.size() != nbands)
	{
		weights.clear();
		for (size_t i = 0; i < bandList.size(); i++)
		{
			weights.push_back(1.0/(double)nbands);
		}
	}
	else
	{
		double weight_sum = 0.0;
		for (size_t i = 0; i < bandList.size(); i++)
		{
			weight_sum += weights[i];
		}
		for (size_t i = 0; i < bandList.size(); i++)
		{
			weights[i] /= (double)weight_sum;
		}
	}

	if (bandList.size() < 1)
	{
		return false;
	}
	else if (bandList.size() == 1)
	{
		int nOutWidth = int(rect.width() * scale.x + 0.5);
		int nOutHeight = int(rect.height() * scale.y + 0.5);
		GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(bandList[0] + 1);
		int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
		void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
		pRasterBand->RasterIO(GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(),
			pData, nOutWidth, nOutHeight, m_eDataType, 0, 0);

		outMat = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
		if (m_eDataType == GDT_Byte)
		{
			if (bStretch)
			{
				stretchRect2Byte<GByte>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
			}
			else
			{
				memcpy(outMat.data, pData, nOutWidth*nOutHeight);
			}
		}
		else if (m_eDataType == GDT_UInt16)
		{
			stretchRect2Byte<GUInt16>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_Int16)
		{
			stretchRect2Byte<GInt16>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_UInt32)
		{
			stretchRect2Byte<GUInt32>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_Int32)
		{
			stretchRect2Byte<GInt32>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_Float32)
		{
			stretchRect2Byte<float>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_Float64)
		{
			stretchRect2Byte<double>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}

		delete[]pData;
		pData = NULL;
	}
	else
	{
		vector<cv::Mat> db;
		int nOutWidth = int(rect.width() * scale.x + 0.5);
		int nOutHeight = int(rect.height() * scale.y + 0.5);
		int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
		//void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
		//cv::Mat combined = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_32FC1);
		cv::Mat combined;
		for (size_t i = 0; i < bandList.size(); i++)
		{
			GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(i + 1);
			int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
			void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
			pRasterBand->RasterIO(GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(),
				pData, nOutWidth, nOutHeight, m_eDataType, 0, 0);

			cv::Mat mt = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_32FC1);

			//cv::Mat outMat = cv::Mat(cv::Size(nw, nh), CV_32FC1);
			//cv::Mat outMat = cv::Mat(cv::Size(nw, nh), CV_8UC1);
			for (int ii = 0; ii < nOutHeight; ii++)
			{
				for (int jj = 0; jj < nOutWidth; ++jj)
				{
					int pos = ii * nOutWidth + jj;
					mt.at<float>(ii, jj) = static_cast<float>(SRCVAL(pData, m_eDataType, pos));
					//mt.data[ii*nOutWidth + jj] = static_cast<GByte>(SRCVAL(pData, m_eDataType, pos));
				}
			}

			//if (m_eDataType == GDT_Byte)
			//{
			//	mt = gdalData2FloatcvMat<GByte>(nOutWidth, nOutHeight, pData);
			//}
			//else if (m_eDataType == GDT_UInt16)
			//{
			//	mt = gdalData2FloatcvMat<GUInt16>(nOutWidth, nOutHeight, pData);
			//}
			//else if (m_eDataType == GDT_Int16)
			//{
			//	mt = gdalData2FloatcvMat<GInt16>(nOutWidth, nOutHeight, pData);
			//}
			//else if (m_eDataType == GDT_UInt32)
			//{
			//	mt = gdalData2FloatcvMat<GUInt32>(nOutWidth, nOutHeight, pData);
			//}
			//else if (m_eDataType == GDT_Int32)
			//{
			//	mt = gdalData2FloatcvMat<GInt32>(nOutWidth, nOutHeight, pData);
			//}
			//else if (m_eDataType == GDT_Float32)
			//{
			//	mt = gdalData2FloatcvMat<float>(nOutWidth, nOutHeight, pData);
			//}
			//else if (m_eDataType == GDT_Float64)
			//{
			//	mt = gdalData2FloatcvMat<double>(nOutWidth, nOutHeight, pData);
			//}
			if (countNonZero(mt) < 1)
			{
				delete[]pData;
				return false;
			}

			//combined = combined + weights[i] * mt;
			if (0 == i)
			{
				combined = weights[i] * mt;
			}
			else
			{
				combined = combined + weights[i] * mt;
			}

			//char buf[256];
			//sprintf(buf, "%d_0.png\0", i + 1);
			//cv::Mat mt_btye = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
			//cv::imwrite(buf, mt);
			//stretchRect2Byte<float>(nOutWidth, nOutHeight, mt.data, mt_btye.data, stretchRatio);
			//cv::imwrite(buf, mt_btye);
			////db.push_back(mt);
			delete[]pData;
		}

		outMat = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
		stretchRect2Byte<float>(nOutWidth, nOutHeight, combined.data, outMat.data, stretchRatio);
		//outMat = norm_0_255(combined);
	}
	return true;
}

bool GdalRasterApp::getRect2CvMatByte(const ossimIrect &rect, cv::Mat& outMat, vector<unsigned int> bands, ossimDpt scale/* = ossimDpt(1.0, 1.0)*/,
	double stretchRatio/* = 0.01*/, bool bStretch/* = true*/)
{
	ossimIrect clipRect = rect.clipToRect(m_Boundary);
	if (clipRect.width() <= 1 || clipRect.height() <= 1)
	{
		return false;
	}
	//if (!m_Boundary.pointWithin(rect.ul()) || !m_Boundary.pointWithin(rect.lr()))
	//{
	//	return false;
	//}

	vector<unsigned int> bandList;
	// check bands
	for (size_t i = 0; i < bands.size(); i++)
	{
		int newBand = 0;
		if (bands[i] < m_iBandCount && bands[i] >= 0)
		{
			newBand = bands[i];
		}

		// check exist
		bool band_exist = false;
		for (size_t j = 0; j < bandList.size(); j++)
		{
			if (bandList[j] == newBand)
			{
				band_exist = true;
				break;
			}
		}

		if (!band_exist)
		{
			bandList.push_back(newBand);
		}
	}

	if (bandList.size() < 1)
	{
		return false;
	}
	else if (bandList.size() == 1)
	{
		int nOutWidth = int(rect.width() * scale.x + 0.5);
		int nOutHeight = int(rect.height() * scale.y + 0.5);
		GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(bandList[0] + 1);
		int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
		void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
		pRasterBand->RasterIO(GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(),
			pData, nOutWidth, nOutHeight, m_eDataType, 0, 0);

		outMat = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
		if (m_eDataType == GDT_Byte)
		{
			if (bStretch)
			{
				stretchRect2Byte<GByte>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
			}
			else
			{
				memcpy(outMat.data, pData, nOutWidth*nOutHeight);
			}
		}
		else if (m_eDataType == GDT_UInt16)
		{
			stretchRect2Byte<GUInt16>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_Int16)
		{
			stretchRect2Byte<GInt16>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_UInt32)
		{
			stretchRect2Byte<GUInt32>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_Int32)
		{
			stretchRect2Byte<GInt32>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_Float32)
		{
			stretchRect2Byte<float>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}
		else if (m_eDataType == GDT_Float64)
		{
			stretchRect2Byte<double>(nOutWidth, nOutHeight, pData, outMat.data, stretchRatio);
		}

		delete[]pData;
		pData = NULL;
	}
	else
	{
		vector<cv::Mat> db;
		int nOutWidth = int(rect.width() * scale.x + 0.5);
		int nOutHeight = int(rect.height() * scale.y + 0.5);
		int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
		//void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
		for (size_t i = 0; i < bandList.size(); i++)
		{
			GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(i + 1);
			int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
			void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
			pRasterBand->RasterIO(GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(),
				pData, nOutWidth, nOutHeight, m_eDataType, 0, 0);

			cv::Mat mt;

			if (m_eDataType == GDT_Byte)
			{
				mt = gdalData2FloatcvMat<GByte>(nOutWidth, nOutHeight, pData);
			}
			else if (m_eDataType == GDT_UInt16)
			{
				mt = gdalData2FloatcvMat<GUInt16>(nOutWidth, nOutHeight, pData);
			}
			else if (m_eDataType == GDT_Int16)
			{
				mt = gdalData2FloatcvMat<GInt16>(nOutWidth, nOutHeight, pData);
			}
			else if (m_eDataType == GDT_UInt32)
			{
				mt = gdalData2FloatcvMat<GUInt32>(nOutWidth, nOutHeight, pData);
			}
			else if (m_eDataType == GDT_Int32)
			{
				mt = gdalData2FloatcvMat<GInt32>(nOutWidth, nOutHeight, pData);
			}
			else if (m_eDataType == GDT_Float32)
			{
				mt = gdalData2FloatcvMat<float>(nOutWidth, nOutHeight, pData);
			}
			else if (m_eDataType == GDT_Float64)
			{
				mt = gdalData2FloatcvMat<double>(nOutWidth, nOutHeight, pData);
			}
			if (countNonZero(mt) < 1)
			{
				delete[]pData;
				return false;
			}
			//if (!getRect2CvMatByte(rect, mt, bandList[i], scale, stretchRatio))
			//{
			//	delete[]pData;
			//	return false;
			//}

			char buf[256];
			sprintf(buf, "%d_0.png\0", i + 1);
			//cv::imwrite(buf, mt);
			db.push_back(mt);
			delete[]pData;
		}
		// Build a matrix with the observations in row:
		cv::Mat data = asRowMatrix(db, CV_32FC1);

		// Number of components to keep for the PCA:
		int num_components = bandList.size();

		// Perform a PCA:
		//cv::PCA pca(data, cv::Mat(), CV_PCA_DATA_AS_ROW, num_components);
		cv::PCA pca(data, cv::Mat(), CV_PCA_DATA_AS_ROW, 0);

		// And copy the PCA results:
		cv::Mat mean = pca.mean.clone();
		cv::Mat eigenvalues = pca.eigenvalues.clone();
		cv::Mat eigenvectors = pca.eigenvectors.clone();
		//for (size_t i = 0; i < bandList.size(); i++)
		//{
		//	cout << eigenvalues.at<float>(i) << endl;
		//	cv::Mat componentMat = norm_0_255(pca.eigenvectors.row(i)).reshape(1, db[i].rows);
		//	char buf[256];
		//	sprintf(buf, "%d.png\0", i + 1);
		//	//cv::imwrite(buf, pca.eigenvectors.row(i).reshape(1, db[i].rows));
		//	cv::imwrite(buf, componentMat);
		//}
		//outMat = norm_0_255(pca.eigenvectors.row(0)).reshape(1, db[0].rows);
		outMat = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
		stretchRect2Byte<float>(nOutWidth, nOutHeight, pca.eigenvectors.row(0).data, outMat.data, stretchRatio);

		//int nOutWidth = int(rect.width() * scale.x + 0.5);
		//int nOutHeight = int(rect.height() * scale.y + 0.5);
		//int nDataTypeSize = GDALGetDataTypeSize(m_eDataType) / 8;
		//void* pData = new GByte[nOutWidth * nOutHeight * nDataTypeSize];
		////stats::pca pca(bandList.size());
		//stats::pca pca(nOutWidth * nOutHeight);
		//for (size_t i = 0; i < bandList.size(); i++)
		//{
		//	GDALRasterBand *pRasterBand = m_pDataset->GetRasterBand(bandList[i] + 1);
		//	pRasterBand->RasterIO(GF_Read, rect.ul().x, rect.ul().y, rect.width(), rect.height(),
		//		pData, nOutWidth, nOutHeight, m_eDataType, 0, 0);
		//	std::vector<double> bandData(nOutWidth * nOutHeight);
		//	for (int ih = 0; ih < nOutHeight; ih++)
		//	{
		//		for (int jw = 0; jw < nOutWidth; ++jw)
		//		{
		//			double data;
		//			if (m_eDataType == GDT_Byte)
		//			{
		//				data = ((GByte*)pData)[ih*nOutWidth + jw];
		//			}
		//			else if (m_eDataType == GDT_UInt16)
		//			{
		//				data = ((GUInt16*)pData)[ih*nOutWidth + jw];
		//			}
		//			else if (m_eDataType == GDT_Int16)
		//			{
		//				data = ((GInt16*)pData)[ih*nOutWidth + jw];
		//			}
		//			else if (m_eDataType == GDT_UInt32)
		//			{
		//				data = ((GUInt32*)pData)[ih*nOutWidth + jw];
		//			}
		//			else if (m_eDataType == GDT_Int32)
		//			{
		//				data = ((GInt32*)pData)[ih*nOutWidth + jw];
		//			}
		//			else if (m_eDataType == GDT_Float32)
		//			{
		//				data = ((float*)pData)[ih*nOutWidth + jw];
		//			}
		//			else if (m_eDataType == GDT_Float64)
		//			{
		//				data = ((double*)pData)[ih*nOutWidth + jw];
		//			}

		//		}
		//	}
		//	pca.add_record(bandData);
		//}
		//pca.solve();
		//vector<double> prinvec1 = pca.get_principal(0);


		//outMat = cv::Mat(cv::Size(nOutWidth, nOutHeight), CV_8UC1);
		//stretchRect2Byte<double>(nOutWidth, nOutHeight, &prinvec1[0], outMat.data, stretchRatio);
		//
		//delete[]pData;
		//pData = NULL;
	}
	return true;
}

void GdalRasterApp::setTileWidth(int iTileWidth)
{
	m_iTileWidth = iTileWidth;
	m_iTileCountX = ceil((m_iWidth - WIDTH_OVER_LAP) / (double)(iTileWidth - WIDTH_OVER_LAP));
	releaseCaches();
}

void GdalRasterApp::setTileHeight(int iTileHeight)
{
	m_iTileHeight = iTileHeight;
	m_iTileCountY = ceil((m_iHeight - HEIGHT_OVER_LAP)/ (double)(iTileHeight - HEIGHT_OVER_LAP));
	releaseCaches();
}

bool GdalRasterApp::checkTileValidity(int iTileX, int iTileY)
{
	if (iTileX >= m_iTileCountX || iTileY >= m_iTileCountY
		|| iTileX < 0 || iTileY <0)
	{
		return false;
	}
	return true;
}

bool GdalRasterApp::getTileOffset(int iTileX, int iTileY, int &offsetX, int &offsetY)
{
	if (!checkTileValidity(iTileX, iTileY))
	{
		cout<<"块号不合法！"<<endl;
		return false;
	}
	offsetX = iTileX*(m_iTileWidth - WIDTH_OVER_LAP);
	offsetY = iTileY*(m_iTileHeight - HEIGHT_OVER_LAP);
	return true;
}

//bool GdalRasterApp::getTileData(int iTileX, int iTileY, 
//								RasterBuf *pBuf,
//								int nBandCount,
//								int *panBandMap/* = 0*/,
//								int nPixelSpace/* = 0*/,
//								int nLineSpace/* = 0*/,
//								int nBandSpace/* = 0*/)
//{
//	// check if there are any empty cache blocks?
//	for (int i = m_CacheSize-1; i >= 0; i--)
//	{
//		if (NULL == m_CacheList[i].data)
//		{			
//			if(!readTile(iTileX, iTileY, nBandCount, &m_CacheList[i]))
//			{
//				return false;
//			}
//			if (i != (m_CacheSize - 1))
//			{
//				RasterBuf tmp = m_CacheList[m_CacheSize-1];
//				m_CacheList[m_CacheSize-1] = m_CacheList[i];
//				m_CacheList[i] = tmp;
//				tmp.data = NULL;
//			}
//			pBuf = &m_CacheList.back();
//			return true;
//		}
//	}
//
//	// check if the required tile exists in cache
//	for (int i = m_CacheSize-1; i >= 0; i--)
//	{
//		if (iTileX == m_CacheList[i].iBufTileX && iTileY == m_CacheList[i].iBufTileY)
//		{
//			// if exist in cache
//			if (i != (m_CacheSize - 1))
//			{
//				RasterBuf tmp = m_CacheList[m_CacheSize-1];
//				m_CacheList[m_CacheSize-1] = m_CacheList[i];
//				m_CacheList[i] = tmp;
//				tmp.data = NULL;
//			}
//			pBuf = &m_CacheList.back();
//			return true;
//		}
//	}
//
//	// if the required tile does not exist in cache
//	RasterBuf p_data = m_CacheList[0];
//	if(!readTile(iTileX, iTileY, nBandCount, &p_data))
//	{
//		return false;
//	}
//	m_CacheList.erase(m_CacheList.begin());
//	m_CacheList.push_back(p_data);
//	pBuf = &m_CacheList.back();
//	return true;
//}

GdalRasterApp::RasterBuf * GdalRasterApp::getTileData(int iTileX, int iTileY, 
								int nBandCount,
								int *panBandMap/* = 0*/,
								int nPixelSpace/* = 0*/,
								int nLineSpace/* = 0*/,
								int nBandSpace/* = 0*/)
{
	// check if there are any empty cache blocks?
	for (int i = m_CacheSize-1; i >= 0; i--)
	{
		if (NULL == m_CacheList[i]->data)
		{			
			if(!readTile(iTileX, iTileY, nBandCount, m_CacheList[i]))
			{
				return NULL;
			}
			if (i != (m_CacheSize - 1))
			{
				RasterBuf* tmp = m_CacheList[m_CacheSize-1];
				m_CacheList[m_CacheSize-1] = m_CacheList[i];
				m_CacheList[i] = tmp;
				tmp = NULL;
			}
			return m_CacheList.back();
		}
	}

	// check if the required tile exists in cache
	for (int i = m_CacheSize-1; i >= 0; i--)
	{
		if (iTileX == m_CacheList[i]->iBufTileX && iTileY == m_CacheList[i]->iBufTileY)
		{
			// if exist in cache
			if (i != (m_CacheSize - 1))
			{
				RasterBuf* tmp = m_CacheList[m_CacheSize-1];
				m_CacheList[m_CacheSize-1] = m_CacheList[i];
				m_CacheList[i] = tmp;
				tmp = NULL;
			}
			return m_CacheList.back();
		}
	}

	// if the required tile does not exist in cache
	RasterBuf *p_data = m_CacheList[0];
	if(!readTile(iTileX, iTileY, nBandCount, p_data))
	{
		return false;
	}
	m_CacheList.erase(m_CacheList.begin());
	m_CacheList.push_back(p_data);
	return m_CacheList.back();
}

//bool GdalRasterApp::readTile(int iTileX, int iTileY,
//							 int nBandCount,
//							 int *panBandMap/* = 0*/,
//							 int nPixelSpace/* = 0*/,
//							 int nLineSpace/* = 0*/,
//							 int nBandSpace/* = 0*/)
//{
//	releaseBuf();
//	if (!updateBufInfo(iTileX, iTileY))
//	{
//		return false;
//	}
//	m_pBufData = new BYTE[m_BufInfo.iBufWidth*m_BufInfo.iBufHeight*nBandCount];
//
//	//int BufSizeX = m_iTileWidth;
//	//int BufSizeY = m_iTileHeight;
//	//// 末尾小块
//	//if (iTileX == m_iTileCountX - 1)
//	//{
//	//	BufSizeX = (m_iTileCountX - 1) % m_iTileWidth + 1;
//	//}
//	//if (iTileY == m_iTileCountY - 1)
//	//{
//	//	BufSizeY = (m_iTileCountY - 1) % m_iTileHeight + 1;
//	//}	
//
//	int offsetX, offsetY;
//	getTileOffset(iTileX, iTileY, offsetX, offsetY);
//	m_pDataset->RasterIO( GF_Read, offsetX, offsetY, m_BufInfo.iBufWidth, m_BufInfo.iBufHeight, m_pBufData, m_BufInfo.iBufWidth, m_BufInfo.iBufHeight, getDataType(),
//		nBandCount, panBandMap, nPixelSpace, nLineSpace, nBandSpace);
//	return true;
//}


bool GdalRasterApp::readTile(int iTileX, int iTileY,
							 int nBandCount,
							 RasterBuf *pcache,
							 int *panBandMap/* = 0*/,
							 int nPixelSpace/* = 0*/,
							 int nLineSpace/* = 0*/,
							 int nBandSpace/* = 0*/)
{
	int nBandDataSize = GDALGetDataTypeSize( m_eDataType ) / 8;
	if (!updateBufInfo(pcache, iTileX, iTileY))
	{
		return false;
	}
	//if (sizeof(pcache->data) != pcache->iBufWidth* pcache->iBufHeight*nBandCount*nBandDataSize)
	//{
	//	delete []pcache->data;
	//	if (GDT_Byte == m_eDataType)
	//	{
	//		pcache->data = new GByte[pcache->iBufWidth*pcache->iBufHeight*nBandCount];
	//	}
	//	else if(GDT_UInt16 == m_eDataType)
	//	{
	//		pcache->data = new GUInt16[pcache->iBufWidth*pcache->iBufHeight*nBandCount];
	//	}
	//	else if(GDT_Int16 == m_eDataType)
	//	{
	//		pcache->data = new GInt16[pcache->iBufWidth*pcache->iBufHeight*nBandCount];
	//	}
	//	else if(GDT_UInt32 == m_eDataType)
	//	{
	//		pcache->data = new GUInt32[pcache->iBufWidth*pcache->iBufHeight*nBandCount];
	//	}
	//	else if(GDT_Int32 == m_eDataType)
	//	{
	//		pcache->data = new GInt32[pcache->iBufWidth*pcache->iBufHeight*nBandCount];
	//	}
	//	else
	//	{
	//		pcache->data = new GByte[pcache->iBufWidth*pcache->iBufHeight*nBandCount];
	//	}
	//}
	if (sizeof(pcache->data) != pcache->iBufWidth* pcache->iBufHeight*nBandCount*nBandDataSize)
	{
		delete []pcache->data;
		pcache->data = new GByte[pcache->iBufWidth*pcache->iBufHeight*nBandCount*nBandDataSize];
	}
	int offsetX, offsetY;
	getTileOffset(iTileX, iTileY, offsetX, offsetY);
	//GByte *temp;
	//temp = new GByte[pcache->iBufWidth*pcache->iBufHeight*nBandCount*nBandDataSize];
	
	//nPixelSpace = sizeof(GByte)*nBandCount;
	//nLineSpace = sizeof(GByte)*nBandCount*pcache->iBufWidth;
	//nBandSpace = sizeof(GByte);
	//m_pDataset->RasterIO( GF_Read, offsetX, offsetY, pcache->iBufWidth, pcache->iBufHeight, (GByte*)pcache->data, pcache->iBufWidth, pcache->iBufHeight, m_eDataType,
	//	nBandCount, panBandMap, nPixelSpace, nLineSpace, nBandSpace);
	nPixelSpace = nBandCount*nBandDataSize;
	nLineSpace = nBandCount*pcache->iBufWidth*nBandDataSize;
	nBandSpace = nBandDataSize;
	m_pDataset->RasterIO( GF_Read, offsetX, offsetY, pcache->iBufWidth, pcache->iBufHeight, pcache->data, pcache->iBufWidth, pcache->iBufHeight, m_eDataType,
		nBandCount, panBandMap, nPixelSpace, nLineSpace, nBandSpace);
	return true;
}

//bool GdalRasterApp::updateBufInfo(int iTileX, int iTileY)
//{
//	return updateBufInfo(m_BufInfo, iTileX, iTileY);
//}

bool GdalRasterApp::updateBufInfo(RasterBuf *pcache,int iTileX, int iTileY)
{
	if (!checkTileValidity(iTileX, iTileY))
	{
		cout<<"块号不合法！"<<endl;
		return false;
	}

	pcache->iBufWidth = m_iTileWidth;
	pcache->iBufHeight = m_iTileHeight;
	// 末尾小块
	if (iTileX == m_iTileCountX - 1)
	{
		pcache->iBufWidth = (m_iWidth - 1) % m_iTileWidth + 1;
	}
	if (iTileY == m_iTileCountY - 1)
	{
		pcache->iBufHeight = (m_iHeight - 1) % m_iTileHeight + 1;
	}

	getTileOffset(iTileX, iTileY, pcache->iBufOffsetX, pcache->iBufOffsetY);

	pcache->iBufTileX = iTileX;
	pcache->iBufTileY = iTileY;

	pcache->iBandCount = m_iBandCount;

	pcache->eDataType = m_eDataType;
	return true;
}

//void GdalRasterApp::releaseBuf()
//{
//	if (m_pBufData)
//	{
//		delete []m_pBufData;
//		m_pBufData = NULL;
//	}
//
//	m_BufInfo.iBufWidth = -1;
//	m_BufInfo.iBufHeight = -1;
//
//	m_BufInfo.iBufOffsetX = -1;
//	m_BufInfo.iBufOffsetY = -1;
//
//	m_BufInfo.iBufTileX = -1;
//	m_BufInfo.iBufTileY = -1;
//}

void GdalRasterApp::releaseBuf(RasterBuf *pcache)
{
	if (pcache->data)
	{
		delete []pcache->data;
		pcache->data = NULL;
	}

	pcache->iBufWidth = -1;
	pcache->iBufHeight = -1;

	pcache->iBufOffsetX = -1;
	pcache->iBufOffsetY = -1;

	pcache->iBufTileX = -1;
	pcache->iBufTileY = -1;

	pcache->iBandCount = 0;
	pcache->eDataType = m_eDataType;
}

void GdalRasterApp::releaseCaches()
{
	for (int i = 0; i < m_CacheSize; i++)
	{
		releaseBuf(m_CacheList[i]);
	}
}

bool GdalRasterApp::getPixel(int x, int y, void** data)
{
	int iTileX = (x-1) / m_iTileWidth;
	int iTileY = (y-1) / m_iTileHeight;
	RasterBuf *pBuf = getTileData(iTileX, iTileY,m_iBandCount);
	if (!pBuf)
	{
		int nBandDataSize = GDALGetDataTypeSize( m_eDataType ) / 8;
		*data = new GByte[m_iBandCount*nBandDataSize];
		memset(*data, 0, m_iBandCount*nBandDataSize);
		return false;
	}

	int xx = x - pBuf->iBufOffsetX;
	int yy = y - pBuf->iBufOffsetY;

	if (xx < 0 || xx >= pBuf->iBufWidth || yy < 0 || yy >= pBuf->iBufHeight)
	{
		int nBandDataSize = GDALGetDataTypeSize( m_eDataType ) / 8;
		*data = new GByte[m_iBandCount*nBandDataSize];
		memset(*data, 0, m_iBandCount*nBandDataSize);
		return false;
	}

	int nBandDataSize = GDALGetDataTypeSize( m_eDataType ) / 8;
	*data = new GByte[m_iBandCount*nBandDataSize];

	int pos = yy * pBuf->iBufWidth + xx;
	memcpy(*data, (GByte*)pBuf->data+pos*m_iBandCount*nBandDataSize, m_iBandCount*nBandDataSize);

	GInt16 r = *((GInt16*)*data);
	GInt16 g = *((GInt16*)*data + 1);
	GInt16 b = *((GInt16*)*data + 2);
	return true;
}
}