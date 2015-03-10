#include "GdalRasterApp.h"

//const int over_lap = 4;
#define HEIGHT_OVER_LAP 0
#define WIDTH_OVER_LAP 0

GdalRasterApp::GdalRasterApp(void)
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
	m_iBandCount = m_pDataset->GetRasterCount();
	m_eDataType = m_pDataset->GetRasterBand(1)->GetRasterDataType();

	m_pDataset->GetGeoTransform( m_adfGeoTransform );			
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

	if (xx < 0 || xx >= pBuf->iBufWidth || y < 0 || y >= pBuf->iBufHeight)
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