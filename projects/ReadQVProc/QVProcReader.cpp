#include "QVProcReader.h"
#include "gdal_priv.h"
#include "cpl_string.h"
#include "strUtil.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <QDir>
#include <QXmlStreamWriter>
#include <QMapIterator>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;


QVProcReader::QVProcReader()
	:theReferenceLine(0),
	time_ready(false)
{

}

bool QVProcReader::read(const char* filename, const char* outPath)
{
	int last_pos = 0;
	//string strPath = SBeforeLast(string(filename), '\\');
	string strPath = string(outPath);
	GDALAllRegister();
	//OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持
	ifstream fs(filename, std::ios_base::binary);
	//fs.open(filename, std::ios_base::binary);
	if(!fs)
	{
		cerr<<"open error!"<<endl;
		return false;
	}
	fs.seekg(0, ios::end);
	LONGLONG file_size = fs.tellg();
	fs.seekg(0, ios::beg);
	
	// 读QUIVIMAGE_HEAD_INFO
	QUIVIMAGE_HEAD_INFO theHeaderInfo;
	readHeadInfo(fs, theHeaderInfo);

	int line_offset = 0;
	//if (0 == strcmp(theHeaderInfo.satellite_id, "HJ-1A"))
	//{
	//	if (0 == strcmp(theHeaderInfo.sensor_id, "CCD-1"))
	//	{
	//		//line_offset = 1835;
	//		line_offset = 0;
	//	}
	//	else if (0 == strcmp(theHeaderInfo.sensor_id, "CCD-2"))
	//	{
	//		//line_offset = 1835;
	//	}
	//}
	//else if (0 == strcmp(theHeaderInfo.satellite_id, "HJ-1B"))
	//{
	//	if (0 == strcmp(theHeaderInfo.sensor_id, "CCD-1"))
	//	{
	//		//line_offset = 2268;
	//		//line_offset = 1835;
	//		line_offset = 0;
	//	}
	//	else if (0 == strcmp(theHeaderInfo.sensor_id, "CCD-2"))
	//	{
	//		//line_offset = 2390;
	//	}
	//}

	fstream ephFile;
	ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::out);
	ephFile.close();

	fstream nadirsFile;
	nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::out);
	nadirsFile.close();

	theHeaderInfo.band_num = theHeaderInfo.gray_image_flag?1:3;
	int lineDataByte = theHeaderInfo.band_num * theHeaderInfo.data_width * theHeaderInfo.sample_bit_count/8;
	theHeaderInfo.line_num = (file_size-180)/(92+lineDataByte);
	theHeaderInfo.sample_num = theHeaderInfo.data_width;

	//theHeaderInfo.band_num = theBands;


	//writeDescXML(theHeaderInfo, (strPath+"\\desc.xml").c_str());

	int percent = 0;
	printf("\r%d%%", percent);
	fflush(stdout);

	//int tileSize = 2*theHeaderInfo.data_width;
	int tileSize = 1*theHeaderInfo.data_width;
	const int OVERLAP = 0;
	int row,col;
	row=(int)(theHeaderInfo.line_num / tileSize);
	col=(int)(theHeaderInfo.sample_num / tileSize);
	if(row * tileSize < theHeaderInfo.line_num ) row =row +1;
	if(col * tileSize < theHeaderInfo.sample_num ) col =col +1;

	string tilFile = strPath + "\\IMAGE.TIL";
	ofstream ftil(tilFile.c_str(),ios::out);
	ftil<<"bandId = \"RGB\""<<endl;
	ftil<<"numTiles = "<<col*row<<";"<<endl;
	ftil<<"tileSizeX = "<<col<<";"<<endl;
	ftil<<"tileSizeY = "<<row<<";"<<endl;
	ftil<<"tileUnits = \"Pixels\";"<<endl;
	ftil<<"tileOverlap = 0;"<<endl;

	GDALDatasetH hDstDS = new GDALDatasetH;
	const char         *pszFormat = "GTiff";
	char               *pszTargetSRS = NULL;
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALDataType eDT;
	if (16 == theHeaderInfo.sample_bit_count)
	{
		eDT = GDALDataType::GDT_Int16;
	}else
	{
		eDT = GDALDataType::GDT_Byte;
	}

	vector<int> panBandMap;
	for (int i=0;i < theHeaderInfo.band_num;++i)
	{
		panBandMap.push_back(i+1);
	}
	vector<QUIVIMAGE_AUX_INFO> auxList;

	GByte *pBuf = NULL;
	int iLine = 0;
	
	int iWidth = theHeaderInfo.sample_num;
	int iHeight = tileSize;
	int lastLine = 0;
	int start_line = 0;
	double last_time = 0.0;
	while (!fs.eof()) {
		// 读取星历参数
		QUIVIMAGE_AUX_INFO aux;
		readQuivAuxInfo(fs, aux);
#if 0
		if (1 == aux.valid_flag)
		{
			// 如果星历参数有效
			//ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
			//ephFile<<iLine+1<<" "
			//	<<aux.satpos.x<<" "
			//	<<aux.satpos.y<<" "
			//	<<aux.satpos.z<<" "
			//	<<aux.satpos.vx<<" "
			//	<<aux.satpos.vy<<" "
			//	<<aux.satpos.vz<<" "
			//	<<aux.satatt.roll<<" "
			//	<<aux.satatt.pitch<<" "
			//	<<aux.satatt.yaw<<endl;
			//ephFile.close();

			if (0 != timeCompare(lastTime, aux.line_time))
			{
				int centerLine = int((lastLine+iLine)*0.5+0.5);

				ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
				//ephFile<<centerLine+1-1792<<" "
				//ephFile<<centerLine+1-1831<<" "
				ephFile<<centerLine<<" "
					<<aux.satpos.x<<" "
					<<aux.satpos.y<<" "
					<<aux.satpos.z<<" "
					<<aux.satpos.vx<<" "
					<<aux.satpos.vy<<" "
					<<aux.satpos.vz<<" "
					<<aux.satatt.roll<<" "
					<<aux.satatt.pitch<<" "
					<<aux.satatt.yaw<<" "
					<<aux.nadir_pos.latitude<<" "
					<<aux.nadir_pos.longitude<<" "
					<<time2String(aux.line_time)<<endl;
				ephFile.close();

				lastTime = aux.line_time;
				start_line = iLine;
			}
			lastLine = iLine;
		}
#else
		if (0 == aux.valid_flag)
		{
			//if (aux.line_count == 0 && !time_ready)
			//{
			//	// 参考行
			//	theReferenceTime = aux.gps_time;
			//	theReferenceLine = iLine;
			//	//theReferenceLine = 0;
			//	time_ready = true;
			//}
			//if (aux.gps_time != last_time && time_ready)
			//{
			//	//if (aux.line_count == 0)
			//	//{
			//	//	// 参考行
			//	//	theReferenceTime = aux.att_time;
			//	//}
			//	// tn = t0 + n*4.369e-3 + delta_t
			//	short line_num = (short)((aux.gps_time - theReferenceTime) / 4.369e-3  + theReferenceLine+ 0.5);// - 1792;
			//	//short line_num = aux.line_count;
			//	ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
			//	//velEci2Ecr(aux.satpos.x, aux.satpos.y, aux.satpos.z, 
			//	//	aux.satpos.vx, aux.satpos.vy, aux.satpos.vz, aux.gps_time);
			//	//posEci2Ecr(aux.satpos.x, aux.satpos.y, aux.satpos.z, aux.gps_time);
			//	ephFile<<line_num<<" "
			//		<<setprecision(25)
			//		<<aux.satpos.x<<" "
			//		<<aux.satpos.y<<" "
			//		<<aux.satpos.z<<" "
			//		<<aux.satpos.vx<<" "
			//		<<aux.satpos.vy<<" "
			//		<<aux.satpos.vz<<" "
			//		<<aux.satatt.roll<<" "
			//		<<aux.satatt.pitch<<" "
			//		<<aux.satatt.yaw<<" "
			//		//<<aux.satatt.vroll<<" "
			//		//<<aux.satatt.vpitch<<" "
			//		//<<aux.satatt.vyaw<<" "
			//		<<aux.gps_time<<endl;
			//	ephFile.close();
			//}
			//if (time_ready)
			//{
			//	last_time = aux.gps_time;
			//}
			if (iLine == 0 || aux.line_count == 0)
			{
				aux.line_num = iLine - aux.line_count - line_offset;
				if (aux.satpos.x == 0 && aux.satpos.y == 0 && aux.satpos.y ==0)
				{
					//iLine++;
					//// 无效行
					//continue;
				}
				else
				{
					auxList.push_back(aux);
				}
				//short line_num = (short)((aux.gps_time - theReferenceTime) / 4.369e-3  + theReferenceLine+ 0.5);// - 1792;
				//short line_num = iLine;
				ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
				//velEci2Ecr(aux.satpos.x, aux.satpos.y, aux.satpos.z, 
				//	aux.satpos.vx, aux.satpos.vy, aux.satpos.vz, aux.gps_time);
				//posEci2Ecr(aux.satpos.x, aux.satpos.y, aux.satpos.z, aux.gps_time);
				ephFile<<aux.line_num<<" "
					<<setprecision(25)
					<<aux.satpos.x<<" "
					<<aux.satpos.y<<" "
					<<aux.satpos.z<<" "
					<<aux.satpos.vx<<" "
					<<aux.satpos.vy<<" "
					<<aux.satpos.vz<<" "
					<<aux.satatt.roll<<" "
					<<aux.satatt.pitch<<" "
					<<aux.satatt.yaw<<" "
					<<aux.satatt.vroll<<" "
					<<aux.satatt.vpitch<<" "
					<<aux.satatt.vyaw<<" "
					//<<aux.satatt.vroll<<" "
					//<<aux.satatt.vpitch<<" "
					//<<aux.satatt.vyaw<<" "
					<<aux.gps_time<<endl;
				ephFile.close();
				double lat, lon, alt;
				ecef2lla(aux.satpos.x, aux.satpos.y, aux.satpos.z, lat, lon, alt);
				nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::app);
				nadirsFile<<aux.line_num<<" "
					<<lat<<" "
					<<lon<<endl;
				nadirsFile.close();
			}
		}
#endif
		int iTile = iLine / tileSize;
		if (iTile*tileSize == iLine)
		{
			iWidth = theHeaderInfo.sample_num;
			iHeight = min(tileSize, theHeaderInfo.line_num-iTile*tileSize);
			const int nBuf = 2048;
			char pszDstFile[nBuf];
			sprintf_s(pszDstFile, nBuf, "%s_%d.tif", SBeforeLast(tilFile, '.').c_str(), iTile+1);
			hDstDS = GDALCreate( hDriver, pszDstFile, iWidth, iHeight, theHeaderInfo.band_num, eDT, NULL );
			pBuf = new GByte[iWidth*iHeight*theHeaderInfo.band_num];

			ftil<<"BEGIN_GROUP = TILE_"<<iTile+1<<endl;
			ftil<<"\tfilename = \""<<SAfterLast(string(pszDstFile), '\\')<<"\""<<endl;
			ftil<<"\tULColOffset = "<<0<<";"<<endl;
			ftil<<"\tULRowOffset = "<<iTile*tileSize<<";"<<endl;
			ftil<<"\tLRColOffset = "<<iWidth-1<<";"<<endl;
			ftil<<"\tLRRowOffset = "<<iTile*tileSize+iHeight-1<<";"<<endl;
			ftil<<"END_GROUP = TILE_"<<iTile+1<<endl;
		}
		//cout<<fs.tellg()<<endl;
		if (file_size - fs.tellg() < lineDataByte)
		{
			((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, theHeaderInfo.band_num, &panBandMap[0],theHeaderInfo.band_num, lineDataByte, 1);
			//int iBand = 1;
			//int band1_offset = -3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth - band1_offset, iHeight, pBuf+theHeaderInfo.band_num*band1_offset, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 2;
			//int band2_offset = -4;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,band2_offset,iWidth, iHeight-band2_offset, pBuf, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			GDALClose( hDstDS );
			CPLFree( pBuf );
			break;
		}

		fs.read((char*)(pBuf+(iLine-iTile*tileSize)*iWidth*theHeaderInfo.band_num), lineDataByte);		
		if (0 == (iLine+1)%tileSize || iLine == theHeaderInfo.line_num)
		{
			((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, theHeaderInfo.band_num, &panBandMap[0],theHeaderInfo.band_num, lineDataByte, 1);
			//int iBand = 1;
			//int band1_offset = -3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth - band1_offset, iHeight, pBuf+theHeaderInfo.band_num*band1_offset, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 2;
			//int band2_offset = -4;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,band2_offset,iWidth, iHeight-band2_offset, pBuf, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			GDALClose( hDstDS );
			CPLFree( pBuf );
		}
		iLine++;

		// 计算当前进度
		streampos pos_ = fs.tellg(); //读取文件指针的位置
		int tmpPercent = (int)(pos_ / (double)(file_size) * 100 + 0.5);
		if(tmpPercent > percent)
		{
			percent = tmpPercent;
			printf("\r%d%%", percent);
			fflush(stdout);
		}
	}
	fs.close();
	ftil<<"END;";
	ftil.close();


	writeDescXML(theHeaderInfo, auxList, (strPath+"\\desc.xml").c_str());

	printf("\r%d%%\n", 100);
	fflush(stdout);
	return true;
}

bool QVProcReader::read_band_seperate(const char* filename, const char* outPath)
{
	int last_pos = 0;
	//string strPath = SBeforeLast(string(filename), '\\');
	string strPath = string(outPath);
	GDALAllRegister();
	//OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持
	ifstream fs(filename, std::ios_base::binary);
	if(!fs)
	{
		cerr<<"open error!"<<endl;
		return false;
	}
	fs.seekg(0, ios::end);
	LONGLONG file_size = fs.tellg();
	fs.seekg(0, ios::beg);

	// 读QUIVIMAGE_HEAD_INFO
	QUIVIMAGE_HEAD_INFO theHeaderInfo;
	readHeadInfo(fs, theHeaderInfo);

	fstream ephFile;
	ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::out);
	ephFile.close();

	fstream nadirsFile;
	nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::out);
	nadirsFile.close();

	theHeaderInfo.band_num = theHeaderInfo.gray_image_flag?1:3;
	int lineDataByte = theHeaderInfo.band_num * theHeaderInfo.data_width * theHeaderInfo.sample_bit_count/8;
	theHeaderInfo.line_num = (file_size-180)/(92+lineDataByte);
	theHeaderInfo.sample_num = theHeaderInfo.data_width;

	//theHeaderInfo.band_num = theBands;

	//writeDescXML(theHeaderInfo, (strPath+"\\desc.xml").c_str());

	int percent = 0;
	printf("\r%d%%", percent);
	fflush(stdout);

	//int tileSize = 2*theHeaderInfo.data_width;
	int tileSize = 1*theHeaderInfo.data_width;
	const int OVERLAP = 0;
	int row,col;
	row=(int)(theHeaderInfo.line_num / tileSize);
	col=(int)(theHeaderInfo.sample_num / tileSize);
	if(row * tileSize < theHeaderInfo.line_num ) row =row +1;
	if(col * tileSize < theHeaderInfo.sample_num ) col =col +1;
	
	vector<GDALDatasetH> hDstDSList;
	const char         *pszFormat = "GTiff";
	char               *pszTargetSRS = NULL;
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALDataType eDT;
	if (16 == theHeaderInfo.sample_bit_count)
	{
		eDT = GDALDataType::GDT_Int16;
	}else
	{
		eDT = GDALDataType::GDT_Byte;
	}
	int dataTypeSize = GDALGetDataTypeSize(eDT) / 8;

	vector<int> panBandMap;
	for (int i=0;i < theHeaderInfo.band_num;++i)
	{
		panBandMap.push_back(i+1);

		const int nBuf = 2048;
		char pszDstFile[nBuf];
		sprintf_s(pszDstFile, nBuf, "%s\\IMAGE_B%d.tif", strPath.c_str(), i+1);
		GDALDatasetH hDstDS = new GDALDatasetH;
		hDstDS = GDALCreate( hDriver, pszDstFile, theHeaderInfo.sample_num, theHeaderInfo.line_num, 1, eDT, NULL );
		hDstDSList.push_back(hDstDS);
	}

	GByte* pBuf = NULL;
	vector<QUIVIMAGE_AUX_INFO> auxList;

	int iLine = 0;
	int iWidth = theHeaderInfo.sample_num;
	int iHeight = tileSize;
	int lastLine = 0;
	int start_line = 0;
	double last_time = 0.0;
	while (!fs.eof()) {
		// 读取星历参数
		QUIVIMAGE_AUX_INFO aux;
		readQuivAuxInfo(fs, aux);
		if (0 == aux.valid_flag)
		{
			if (iLine == 0 || aux.line_count == 0)
			{
				aux.line_num = iLine - aux.line_count;
				if (aux.satpos.x == 0 && aux.satpos.y == 0 && aux.satpos.y ==0)
				{
					//iLine++;
					//// 无效行
					//continue;
				}
				else
				{
					auxList.push_back(aux);
				}
				ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
				ephFile<<aux.line_num<<" "
					<<setprecision(25)
					<<aux.satpos.x<<" "
					<<aux.satpos.y<<" "
					<<aux.satpos.z<<" "
					<<aux.satpos.vx<<" "
					<<aux.satpos.vy<<" "
					<<aux.satpos.vz<<" "
					<<aux.satatt.roll<<" "
					<<aux.satatt.pitch<<" "
					<<aux.satatt.yaw<<" "
					<<aux.satatt.vroll<<" "
					<<aux.satatt.vpitch<<" "
					<<aux.satatt.vyaw<<" "
					<<aux.gps_time<<endl;
				ephFile.close();
				double lat, lon, alt;
				ecef2lla(aux.satpos.x, aux.satpos.y, aux.satpos.z, lat, lon, alt);
				nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::app);
				nadirsFile<<aux.line_num<<" "
					<<lat<<" "
					<<lon<<endl;
				nadirsFile.close();
			}
		}

		int iTile = iLine / tileSize;
		if (iTile*tileSize == iLine)
		{
			iWidth = theHeaderInfo.sample_num;
			iHeight = min(tileSize, theHeaderInfo.line_num-iTile*tileSize);
			pBuf = new GByte[iWidth*iHeight*theHeaderInfo.band_num*dataTypeSize];
		}

		if (file_size - fs.tellg() < lineDataByte || iLine == theHeaderInfo.line_num)
		{
			for (int ib = 0;ib < theHeaderInfo.band_num;++ib)
			{
				int iBand = 1;
				((GDALDataset*)hDstDSList[ib])->RasterIO(GF_Write,0, tileSize*iTile, iWidth, iHeight, pBuf+ib, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
				GDALClose( hDstDSList[ib] );
			}
			CPLFree( pBuf );
			break;
		}

		fs.read((char*)(pBuf+(iLine-iTile*tileSize)*iWidth*theHeaderInfo.band_num), lineDataByte);

		if (0 == (iLine+1)%tileSize && iLine != theHeaderInfo.line_num)
		{
			for (int ib = 0;ib < theHeaderInfo.band_num;++ib)
			{
				int iBand = 1;
				((GDALDataset*)hDstDSList[ib])->RasterIO(GF_Write, 0, tileSize*iTile, iWidth, iHeight, pBuf+ib, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			}
			CPLFree( pBuf );
		}
		iLine++;

		// 计算当前进度
		streampos pos_ = fs.tellg(); //读取文件指针的位置
		int tmpPercent = (int)(pos_ / (double)(file_size) * 100 + 0.5);
		if(tmpPercent > percent)
		{
			percent = tmpPercent;
			printf("\r%d%%", percent);
			fflush(stdout);
		}
	}
	fs.close();


	writeDescXML(theHeaderInfo, auxList, (strPath+"\\desc.xml").c_str());

	printf("\r%d%%", 100);
	fflush(stdout);
	return true;

}

bool QVProcReader::read_by_scene(const char* filename, const char* outPath)
{
	int last_pos = 0;
	//string strPath = SBeforeLast(string(filename), '\\');
	string strPath = string(outPath);
	GDALAllRegister();
	//OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持
	fstream hdfile;
	hdfile.open(strPath + "\\hd.dat", std::ios_base::out);
	hdfile.close();

	ifstream fs;
	fs.open(filename, std::ios_base::binary);
	if(!fs)
	{
		cerr<<"open error!"<<endl;
		return false;
	}
	fs.seekg(0, ios::end);
	LONGLONG file_size = fs.tellg();
	fs.seekg(0, ios::beg);

	// 读QUIVIMAGE_HEAD_INFO
	QUIVIMAGE_HEAD_INFO theHeaderInfo;
	readHeadInfo(fs, theHeaderInfo);
	
	theHeaderInfo.band_num = theHeaderInfo.gray_image_flag?1:3;
	int lineDataByte = theHeaderInfo.band_num * theHeaderInfo.data_width * theHeaderInfo.sample_bit_count/8;
	theHeaderInfo.line_num = (file_size-180)/(92+lineDataByte);
	theHeaderInfo.sample_num = theHeaderInfo.data_width;

	writeDescXML(theHeaderInfo, (strPath+"\\desc.xml").c_str());

	int percent = 0;
	printf("\r%d%%", percent);
	fflush(stdout);

	int sceneSize = 1*theHeaderInfo.data_width;
	const int OVERLAP = 0;
	int nScenes=(int)(theHeaderInfo.line_num / sceneSize);
	if(nScenes * sceneSize < theHeaderInfo.line_num ) nScenes =nScenes +1;

	fstream ephFile;
	ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::out);
	ephFile.close();

	fstream nadirsFile;
	nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::out);
	nadirsFile.close();

	GDALDatasetH hDstDS = new GDALDatasetH;
	const char         *pszFormat = "GTiff";
	char               *pszTargetSRS = NULL;
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALDataType eDT;
	if (16 == theHeaderInfo.sample_bit_count)
	{
		eDT = GDALDataType::GDT_Int16;
	}else
	{
		eDT = GDALDataType::GDT_Byte;
	}

	vector<int> panBandMap;
	for (int i=0;i < theHeaderInfo.band_num;++i)
	{
		panBandMap.push_back(i+1);
	}
	vector<QUIVIMAGE_AUX_INFO> auxList;
	GByte *pBuf = NULL;
	int iLine = 0;
	while (!fs.eof()) {
		// 读取星历参数
		QUIVIMAGE_AUX_INFO aux;
		readQuivAuxInfo(fs, aux);
#if 1
		if (1 == aux.valid_flag)
		{
			// 如果星历参数有效
			//ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
			//ephFile<<iLine+1<<" "
			//	<<aux.satpos.x<<" "
			//	<<aux.satpos.y<<" "
			//	<<aux.satpos.z<<" "
			//	<<aux.satpos.vx<<" "
			//	<<aux.satpos.vy<<" "
			//	<<aux.satpos.vz<<" "
			//	<<aux.satatt.roll<<" "
			//	<<aux.satatt.pitch<<" "
			//	<<aux.satatt.yaw<<endl;
			//ephFile.close();

			//if (0 != timeCompare(lastTime, aux.line_time))
			//{
				//int centerLine = int((lastLine+iLine)*0.5+0.5);

				ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
				//ephFile<<centerLine+1-1792<<" "
				//ephFile<<centerLine+1-1831<<" "
				ephFile << iLine << " "
					<<aux.satpos.x<<" "
					<<aux.satpos.y<<" "
					<<aux.satpos.z<<" "
					<<aux.satpos.vx<<" "
					<<aux.satpos.vy<<" "
					<<aux.satpos.vz<<" "
					<<aux.satatt.roll<<" "
					<<aux.satatt.pitch<<" "
					<<aux.satatt.yaw<<" "
					<<aux.nadir_pos.latitude<<" "
					<<aux.nadir_pos.longitude<<" "
					<<time2String(aux.line_time)<<endl;
				ephFile.close();

			//	lastTime = aux.line_time;
			//	start_line = iLine;
			//}
			//lastLine = iLine;
		}
#else
		if (0 == aux.valid_flag)
		{
			if (iLine == 0 || aux.line_count == 0)
			{
				//if (aux.satpos.x == 0 && aux.satpos.y == 0 && aux.satpos.y ==0)
				//{
				//	// 无效行
				//	continue;
				//}
				aux.line_num = iLine - aux.line_count;
				if (aux.satpos.x == 0 && aux.satpos.y == 0 && aux.satpos.y ==0)
				{
					//iLine++;
					//// 无效行
					//continue;
				}
				else
				{
					auxList.push_back(aux);
				}
				int line_num;
				line_num = iLine - aux.line_count;// - line_offset;
				ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
				ephFile<<line_num<<" "
					<<setprecision(25)
					<<aux.satpos.x<<" "
					<<aux.satpos.y<<" "
					<<aux.satpos.z<<" "
					<<aux.satpos.vx<<" "
					<<aux.satpos.vy<<" "
					<<aux.satpos.vz<<" "
					<<aux.satatt.roll<<" "
					<<aux.satatt.pitch<<" "
					<<aux.satatt.yaw<<" "
					<<aux.satatt.vroll<<" "
					<<aux.satatt.vpitch<<" "
					<<aux.satatt.vyaw<<" "
					<<aux.gps_time<<endl;
				ephFile.close();
				double lat, lon, alt;
				ecef2lla(aux.satpos.x, aux.satpos.y, aux.satpos.z, lat, lon, alt);
				nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::app);
				nadirsFile<<line_num<<" "
					<<lat<<" "
					<<lon<<endl;
				nadirsFile.close();
			}
		}
#endif
		int iScene = iLine / sceneSize;
		int iWidth = theHeaderInfo.sample_num;
		int iHeight = sceneSize;
		int line_offset = iScene * sceneSize;
		// 末尾景
		if (iScene == nScenes - 1)
		{
			iHeight = theHeaderInfo.line_num - line_offset;
		}
		if (iScene*sceneSize == iLine)
		{
			char buf[1024];
			sprintf_s(buf, "%s\\Scene%02d", strPath.c_str(), iScene+1);
			if (!QDir(buf).exists())
			{
				_mkdir(buf);
			}
			string sceneDir(buf);
			string sceneImage = sceneDir + "\\IMAGE.TIF";

			hDstDS = GDALCreate( hDriver, sceneImage.c_str(), iWidth, iHeight, theHeaderInfo.band_num, eDT, NULL );
			pBuf = new GByte[iWidth*iHeight*theHeaderInfo.band_num];
		}
		//cout<<fs.tellg()<<endl;
		if (file_size - fs.tellg() < lineDataByte || iLine == theHeaderInfo.line_num-1)
		{
			((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, theHeaderInfo.band_num, &panBandMap[0],theHeaderInfo.band_num, lineDataByte, 1);
			//int iBand = 1;
			//int band1_offset = 3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,band1_offset,0,iWidth - band1_offset, iHeight, pBuf, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 2;
			//int band2_offset = 4;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight-band2_offset, pBuf+lineDataByte*band2_offset, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			////((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,band2_offset,iWidth, iHeight-band2_offset, pBuf, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			GDALClose( hDstDS );
			delete[]pBuf;
			pBuf = NULL;
			//CPLFree( pBuf );
			break;
		}

		int npos = fs.tellg();
		fs.read((char*)(pBuf+(iLine-iScene*sceneSize)*iWidth*theHeaderInfo.band_num), lineDataByte);	
		if (iLine-iScene*sceneSize == iHeight-1)
		{
			((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, theHeaderInfo.band_num, &panBandMap[0],theHeaderInfo.band_num, lineDataByte, 1);
			//int iBand = 1;
			//int band1_offset = 1;
			////((GDALDataset*)hDstDS)->RasterIO(GF_Write,band1_offset,0,iWidth - band1_offset, iHeight, pBuf, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth - band1_offset, iHeight, pBuf+theHeaderInfo.band_num*band1_offset, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 2;
			//int band2_offset = 2;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight-band2_offset, pBuf+lineDataByte*band2_offset, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			////((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,band2_offset,iWidth, iHeight-band2_offset, pBuf, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			GDALClose(hDstDS);
			delete[]pBuf;
			pBuf = NULL;
			//CPLFree( pBuf );
		}
		iLine++;

		// 计算当前进度
		streampos pos_ = fs.tellg(); //读取文件指针的位置
		int tmpPercent = (int)(pos_ / (double)(file_size) * 100 + 0.5);
		if(tmpPercent > percent)
		{
			percent = tmpPercent;
			printf("\r%d%%", percent);
			fflush(stdout);
		}
	}
	fs.close();
	writeDescXML(theHeaderInfo, auxList, (strPath+"\\desc.xml").c_str());

	// write nScenes aux files respectively
	for (int iScene = 0;iScene < nScenes;++iScene)
	{
		char buf[1024];
		sprintf_s(buf, "%s\\Scene%02d", strPath.c_str(), iScene+1);
		if (!QDir(buf).exists())
		{
			_mkdir(buf);
		}
		string sceneDir(buf);
		
		int iWidth = theHeaderInfo.sample_num;
		int iHeight = sceneSize;
		int line_offset = iScene * sceneSize;
		// 末尾景
		if (iScene == nScenes - 1)
		{
			iHeight = theHeaderInfo.line_num - line_offset;
		}

		QUIVIMAGE_HEAD_INFO sceneHeaderInfo = theHeaderInfo;
		sceneHeaderInfo.line_num = iHeight;
		writeDescXML(sceneHeaderInfo, (sceneDir+"\\desc.xml").c_str());
		int overlap_lines = 3000;

		fstream inEphFile;
		inEphFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::in);
		fstream outEphFile;
		outEphFile.open((sceneDir+"\\eph.txt").c_str(), std::ios_base::out);		
		char ephLine[2048];
		while (inEphFile.getline(ephLine, 2048))
		{
			int scene_line = atoi(SBeforeFirst(string(ephLine), ' ').c_str()) - line_offset;
			if (scene_line > -overlap_lines && scene_line < sceneSize + overlap_lines)
			{
				outEphFile<<scene_line<<" "
					<<SAfterFirst(string(ephLine), ' ')<<endl;
			}
		}
		inEphFile.close();
		outEphFile.close();

		fstream inNadirFile;
		inNadirFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::in);
		fstream outNadirFile;
		outNadirFile.open((sceneDir+"\\nadirs.txt").c_str(), std::ios_base::out);		
		char nadirLine[2048];
		while (inNadirFile.getline(nadirLine, 2048))
		{
			int scene_line = atoi(SBeforeFirst(string(nadirLine), ' ').c_str()) - line_offset;
			if (scene_line > -overlap_lines && scene_line < sceneSize + overlap_lines)
			{
				outNadirFile<<scene_line<<" "
					<<SAfterFirst(string(nadirLine), ' ')<<endl;
			}
		}
		inNadirFile.close();
		outNadirFile.close();

		writeDescXML(theHeaderInfo, auxList, (sceneDir+"\\desc.xml").c_str(), line_offset);
	}

	printf("\r%d%%", 100);
	fflush(stdout);
	return true;
}

bool QVProcReader::read_by_scene_band_seperate(const char* filename, const char* outPath)
{
	string basename = SBeforeLast(SAfterLast(string(filename), '\\'), '.');
	int last_pos = 0;
	//string strPath = SBeforeLast(string(filename), '\\');
	string strPath = string(outPath);
	strPath = strPath + "\\" + basename;
	if(!QDir(strPath.c_str()).exists())
	{
		_mkdir(strPath.c_str());
	}
	GDALAllRegister();
	//OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持
	ifstream fs;
	fs.open(filename, std::ios_base::binary);
	if(!fs)
	{
		cerr<<"open error!"<<endl;
		return false;
	}
	fs.seekg(0, ios::end);
	LONGLONG file_size = fs.tellg();
	fs.seekg(0, ios::beg);

	// 读QUIVIMAGE_HEAD_INFO
	QUIVIMAGE_HEAD_INFO theHeaderInfo;
	readHeadInfo(fs, theHeaderInfo);

	theHeaderInfo.band_num = theHeaderInfo.gray_image_flag?1:3;
	int lineDataByte = theHeaderInfo.band_num * theHeaderInfo.data_width * theHeaderInfo.sample_bit_count/8;
	theHeaderInfo.line_num = (file_size-180)/(92+lineDataByte);
	theHeaderInfo.sample_num = theHeaderInfo.data_width;

	writeDescXML(theHeaderInfo, (strPath+"\\desc.xml").c_str());

	int percent = 0;
	printf("\r%d%%", percent);
	fflush(stdout);

	int sceneSize = 1*theHeaderInfo.data_width;
	const int OVERLAP = 0;
	int nScenes=(int)(theHeaderInfo.line_num / sceneSize);
	if(nScenes * sceneSize < theHeaderInfo.line_num ) nScenes =nScenes +1;

	fstream ephFile;
	ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::out);
	ephFile.close();

	fstream nadirsFile;
	nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::out);
	nadirsFile.close();

	vector<GDALDatasetH> hDstDSList(theHeaderInfo.band_num);
	const char         *pszFormat = "GTiff";
	char               *pszTargetSRS = NULL;
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALDataType eDT;
	if (16 == theHeaderInfo.sample_bit_count)
	{
		eDT = GDALDataType::GDT_Int16;
	}else
	{
		eDT = GDALDataType::GDT_Byte;
	}

	vector<int> panBandMap;
	for (int i=0;i < theHeaderInfo.band_num;++i)
	{
		panBandMap.push_back(i+1);
	}

	vector<QUIVIMAGE_AUX_INFO> auxList;
	GByte *pBuf = NULL;
	int iLine = 0;
	while (!fs.eof()) {
		// 读取星历参数
		QUIVIMAGE_AUX_INFO aux;
		readQuivAuxInfo(fs, aux);
		if (0 == aux.valid_flag)
		{
			if (iLine == 0 || aux.line_count == 0)
			{
				aux.line_num = iLine - aux.line_count;
				if (aux.satpos.x == 0 && aux.satpos.y == 0 && aux.satpos.y ==0)
				{
					//iLine++;
					//// 无效行
					//continue;
				}
				else
				{
					auxList.push_back(aux);
				}
				int line_num;
				line_num = iLine - aux.line_count;// - line_offset;
				ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
				ephFile<<line_num<<" "
					<<setprecision(25)
					<<aux.satpos.x<<" "
					<<aux.satpos.y<<" "
					<<aux.satpos.z<<" "
					<<aux.satpos.vx<<" "
					<<aux.satpos.vy<<" "
					<<aux.satpos.vz<<" "
					<<aux.satatt.roll<<" "
					<<aux.satatt.pitch<<" "
					<<aux.satatt.yaw<<" "
					<<aux.satatt.vroll<<" "
					<<aux.satatt.vpitch<<" "
					<<aux.satatt.vyaw<<" "
					<<aux.gps_time<<endl;
				ephFile.close();
				double lat, lon, alt;
				ecef2lla(aux.satpos.x, aux.satpos.y, aux.satpos.z, lat, lon, alt);
				nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::app);
				nadirsFile<<line_num<<" "
					<<lat<<" "
					<<lon<<endl;
				nadirsFile.close();
			}
		}
		int iScene = iLine / sceneSize;
		int iWidth = theHeaderInfo.sample_num;
		int iHeight = sceneSize;
		int line_offset = iScene * sceneSize;
		// 末尾景
		if (iScene == nScenes - 1)
		{
			iHeight = theHeaderInfo.line_num - line_offset;
		}
		if (iScene*sceneSize == iLine)
		{
			char buf[1024];
			sprintf_s(buf, "%s\\Scene%02d", strPath.c_str(), iScene+1);
			if (!QDir(buf).exists())
			{
				_mkdir(buf);
			}

			for (int i=0;i < theHeaderInfo.band_num;++i)
			{
				panBandMap.push_back(i+1);

				const int nBuf = 2048;
				char pszDstFile[nBuf];
				sprintf_s(pszDstFile, nBuf, "%s\\IMAGE_B%d.tif", buf, i+1);
				hDstDSList[i] = new GDALDatasetH;
				hDstDSList[i] = GDALCreate( hDriver, pszDstFile, iWidth, iHeight, 1, eDT, NULL );
			}
			pBuf = new GByte[iWidth*iHeight*theHeaderInfo.band_num];
		}
		//cout<<fs.tellg()<<endl;
		if (file_size - fs.tellg() < lineDataByte || iLine == theHeaderInfo.line_num-1)
		{
			for (int ib = 0;ib < theHeaderInfo.band_num;++ib)
			{
				int iBand = 1;
				((GDALDataset*)hDstDSList[ib])->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf+ib, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
				GDALClose( hDstDSList[ib] );
			}
			//int iBand = 1;
			//int band1_offset = 3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,band1_offset,0,iWidth - band1_offset, iHeight, pBuf, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 2;
			//int band2_offset = 4;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight-band2_offset, pBuf+lineDataByte*band2_offset, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			////((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,band2_offset,iWidth, iHeight-band2_offset, pBuf, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);

			delete[]pBuf;
			pBuf = NULL;
			break;
		}

		int npos = fs.tellg();
		fs.read((char*)(pBuf+(iLine-iScene*sceneSize)*iWidth*theHeaderInfo.band_num), lineDataByte);	
		if (iLine-iScene*sceneSize == iHeight-1)
		{

			for (int ib = 0;ib < theHeaderInfo.band_num;++ib)
			{
				int iBand = 1;
				((GDALDataset*)hDstDSList[ib])->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf+ib, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
				GDALClose( hDstDSList[ib] );
			}
			//int iBand = 1;
			//int band1_offset = 1;
			////((GDALDataset*)hDstDS)->RasterIO(GF_Write,band1_offset,0,iWidth - band1_offset, iHeight, pBuf, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth - band1_offset, iHeight, pBuf+theHeaderInfo.band_num*band1_offset, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 2;
			//int band2_offset = 2;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight-band2_offset, pBuf+lineDataByte*band2_offset, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			////((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,band2_offset,iWidth, iHeight-band2_offset, pBuf, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
			//iBand = 3;
			//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);

			delete[]pBuf;
			pBuf = NULL;
		}
		iLine++;

		// 计算当前进度
		streampos pos_ = fs.tellg(); //读取文件指针的位置
		int tmpPercent = (int)(pos_ / (double)(file_size) * 100 + 0.5);
		if(tmpPercent > percent)
		{
			percent = tmpPercent;
			printf("\r%d%%", percent);
			fflush(stdout);
		}
	}
	fs.close();
	writeDescXML(theHeaderInfo, auxList, (strPath+"\\desc.xml").c_str());

	// write nScenes aux files respectively
	for (int iScene = 0;iScene < nScenes;++iScene)
	{
		char buf[1024];
		sprintf_s(buf, "%s\\Scene%02d", strPath.c_str(), iScene+1);
		if (!QDir(buf).exists())
		{
			_mkdir(buf);
		}
		string sceneDir(buf);

		int iWidth = theHeaderInfo.sample_num;
		int iHeight = sceneSize;
		int line_offset = iScene * sceneSize;
		// 末尾景
		if (iScene == nScenes - 1)
		{
			iHeight = theHeaderInfo.line_num - line_offset;
		}

		QUIVIMAGE_HEAD_INFO sceneHeaderInfo = theHeaderInfo;
		sceneHeaderInfo.line_num = iHeight;
		writeDescXML(sceneHeaderInfo, (sceneDir+"\\desc.xml").c_str());
		int overlap_lines = 3000;

		fstream inEphFile;
		inEphFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::in);
		fstream outEphFile;
		outEphFile.open((sceneDir+"\\eph.txt").c_str(), std::ios_base::out);		
		char ephLine[2048];
		while (inEphFile.getline(ephLine, 2048))
		{
			int scene_line = atoi(SBeforeFirst(string(ephLine), ' ').c_str()) - line_offset;
			if (scene_line > -overlap_lines && scene_line < sceneSize + overlap_lines)
			{
				outEphFile<<scene_line<<" "
					<<SAfterFirst(string(ephLine), ' ')<<endl;
			}
		}
		inEphFile.close();
		outEphFile.close();

		fstream inNadirFile;
		inNadirFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::in);
		fstream outNadirFile;
		outNadirFile.open((sceneDir+"\\nadirs.txt").c_str(), std::ios_base::out);		
		char nadirLine[2048];
		while (inNadirFile.getline(nadirLine, 2048))
		{
			int scene_line = atoi(SBeforeFirst(string(nadirLine), ' ').c_str()) - line_offset;
			if (scene_line > -overlap_lines && scene_line < sceneSize + overlap_lines)
			{
				outNadirFile<<scene_line<<" "
					<<SAfterFirst(string(nadirLine), ' ')<<endl;
			}
		}
		inNadirFile.close();
		outNadirFile.close();

		writeDescXML(theHeaderInfo, auxList, (sceneDir+"\\desc.xml").c_str(), line_offset);
	}

	printf("\r%d%%", 100);
	fflush(stdout);
	return true;
}

bool QVProcReader::readHeadInfo(ifstream& fs, QUIVIMAGE_HEAD_INFO& header)
{
	//	接收该数据的地面站的标识 16
	fs.read((char*)header.station_id, sizeof(header.station_id));
	//	对应的卫星标识 16
	fs.read((char*)header.satellite_id, sizeof(header.satellite_id));
	//	对应的传感器标识 16
	fs.read((char*)header.sensor_id, sizeof(header.sensor_id));
	//	对应的工作模式信息 16
	fs.read((char*)header.work_mode, sizeof(header.work_mode));
	//	任务标识 20
	fs.read((char*)header.JobTaskID, sizeof(header.JobTaskID));
	//	轨道号 4
	fs.read((char*)&header.orbit_num, sizeof(header.orbit_num));
	//	通道号 4
	fs.read((char*)&header.channel_num, sizeof(header.channel_num));
	//	本任务单传输快视文件总数 4
	fs.read((char*)&header.totalfiles, sizeof(header.totalfiles));
	//	本次传输的文件在所属任务中的文件序号 4
	fs.read((char*)&header.files_num, sizeof(header.files_num));
	//	进行快数处理时的降采样数,上一步的降采样数 4
	fs.read((char*)&header.desample_num, sizeof(header.desample_num));
	//	上一步传进来的宽度,(如果我在做降采样的话),这里会进一步变小 4
	fs.read((char*)&header.data_width, sizeof(header.data_width));
	//	快视处理后的图像量化比特数，只有8Bit和16Bit两种情况 4
	fs.read((char*)&header.sample_bit_count, sizeof(header.sample_bit_count));
	//	当此值为真时，该图像为单通道灰度图像，否则为多通道伪彩色图像 4
	fs.read((char*)&header.gray_image_flag, sizeof(header.gray_image_flag));

	//	数据接收的开始时间 28
	readTimeInfo(fs, header.start_time);
	//  数据接收的结束时间 28 
	readTimeInfo(fs, header.end_time);

	//	对应的地面站在地球上的位置信息 8
	fs.read((char*)&header.station_pos.longitude, sizeof(header.station_pos.longitude));
	fs.read((char*)&header.station_pos.latitude, sizeof(header.station_pos.latitude));
	return true;
}

bool QVProcReader::readTimeInfo(ifstream& fs, _TIMESTAMP& t)
{
	//	year 4
	fs.read((char*)&t.year, sizeof(int));
	//	mon 4
	fs.read((char*)&t.mon, sizeof(int));
	//	day 4
	fs.read((char*)&t.day, sizeof(int));
	//	hour 4
	fs.read((char*)&t.hour, sizeof(int));
	//	min 4
	fs.read((char*)&t.min, sizeof(int));
	//	sec 4
	fs.read((char*)&t.sec, sizeof(int));
	//	nsec 4
	fs.read((char*)&t.nsec, sizeof(int));
	return true;
}

bool QVProcReader::readSAT_POS(ifstream& fs, SAT_POS& pos)
{
	//	x 4
	fs.read((char*)&pos.x, sizeof(float));
	//	y 4
	fs.read((char*)&pos.y, sizeof(float));
	//	z 4
	fs.read((char*)&pos.z, sizeof(float));
	//	vx 4
	fs.read((char*)&pos.vx, sizeof(float));
	//	vy 4
	fs.read((char*)&pos.vy, sizeof(float));
	//	vz 4
	fs.read((char*)&pos.vz, sizeof(float));
	return true;
}

bool QVProcReader::readSAT_ATT(ifstream& fs, SAT_ATT& att)
{
	//	roll 4
	fs.read((char*)&att.roll, sizeof(float));
	//	pitch 4
	fs.read((char*)&att.pitch, sizeof(float));
	//	yaw 4
	fs.read((char*)&att.yaw, sizeof(float));
	const float pi = 3.14159265358979f;
	att.roll *= pi/180.0f;
	att.pitch *= pi/180.0f;
	att.yaw *= pi/180.0f;
	return true;
}

void QVProcReader::BitInverse(char* chrs, int nBit)
{
	char* tmp = new char[nBit];
	int result = 0;
	for (int i = 0;i < nBit;++i)
	{
		tmp[i] = chrs[nBit-1-i];
	}
	memcpy(chrs, tmp, nBit);
}

bool QVProcReader::readQuivAuxInfo(ifstream& fs, QUIVIMAGE_AUX_INFO& aux)
{
	bool bState = true;
	//	该值为=1时有效，其他情况下辅助信息为无效值 4
	fs.read((char*)&aux.valid_flag, sizeof(aux.valid_flag));


	//// 识别码 6
	//char identifier[6];
	//fs.read((char*)identifier, sizeof(identifier));
	//// 路标识
	//char pathId;
	//fs.read((char*)&pathId, 1);

	//// 星务对时码 秒 4
	//unsigned int Star_s;
	//fs.read((char*)&Star_s, sizeof(Star_s));
	//// 星务对时码 毫秒 2
	//unsigned short int Star_ms;
	//fs.read((char*)&Star_ms, sizeof(Star_ms));

	//// 累计秒数 秒 4
	//unsigned int GPS_s;
	//fs.read((char*)&GPS_s, sizeof(GPS_s));
	//// 累计秒数 毫秒 2
	//unsigned short int GPS_ms;
	//fs.read((char*)&GPS_ms, sizeof(GPS_ms));
	//// X轴位置坐标 4字节补码（0.1m） 4
	//signed int x;
	//fs.read((char*)&x, sizeof(x));
	//// Y轴位置坐标 4字节补码（0.1m） 4
	//signed int y;
	//fs.read((char*)&y, sizeof(y));
	//// Z轴位置坐标 4字节补码（0.1m） 4
	//signed int z;
	//fs.read((char*)&z, sizeof(z));
	//// X轴速度坐标 4字节补码（0.01m/s） 4
	//signed int vx;
	//fs.read((char*)&vx, sizeof(vx));
	//// Y轴速度坐标 4字节补码（0.01m/s） 4
	//signed int vy;
	//fs.read((char*)&vy, sizeof(vy));
	//// Z轴速度坐标 4字节补码（0.01m/s） 4
	//signed int vz;
	//fs.read((char*)&vz, sizeof(vz));

	//// GPS对时码 秒 4
	//unsigned int cGPS_s;
	//fs.read((char*)&cGPS_s, sizeof(cGPS_s));
	//// GPS对时码 毫秒 2
	//unsigned short int cGPS_ms;
	//fs.read((char*)&cGPS_ms, sizeof(cGPS_ms));

	//// 姿控系统时间标准双精度型浮点数 秒 8
	//double att_t;
	//fs.read((char*)&att_t, sizeof(att_t));
	//// 滚动角（roll）标准浮点数 弧度 4
	//float roll;
	//fs.read((char*)&roll, sizeof(roll));
	//// 俯仰角（pitch）标准浮点数 弧度 4
	//float pitch;
	//fs.read((char*)&pitch, sizeof(pitch));
	//// 偏航角（yaw）标准浮点数 弧度 4
	//float yaw;
	//fs.read((char*)&yaw, sizeof(yaw));
	//// 滚动角速度（vroll）标准浮点数 弧度/秒 4
	//float vroll;
	//fs.read((char*)&vroll, sizeof(vroll));
	//// 俯仰角速度（pitch）标准浮点数 弧度/秒 4
	//float vpitch;
	//fs.read((char*)&vpitch, sizeof(vpitch));
	//// 偏航角速度（yaw）标准浮点数 弧度/秒 4
	//float vyaw;
	//fs.read((char*)&vyaw, sizeof(vyaw));

	//// CCD 行计数值 2
	//short iLine;
	//fs.read((char*)&iLine, sizeof(iLine));

	//// 暗像元 2
	//short Dark_pixel;
	//fs.read((char*)&Dark_pixel, sizeof(Dark_pixel));

	//// 后 3 字节填零
	//fs.ignore(3);

	//int iCount = 0;
	//if (1 == aux.valid_flag)
	//{
	//	cout<<++iCount<<endl;
	//}


	if (0 == aux.valid_flag)
	{
		// 辅助信息为无效值
		// 识别码 6
		char identifier[6];
		fs.read((char*)identifier, sizeof(identifier));
		//cerr << identifier << endl;
		// 路标识
		char pathId;
		fs.read((char*)&pathId, 1);

		// 星务对时码 秒 4
		unsigned int Star_s;
		//int Star_s;
		int nn = sizeof(Star_s);
		fs.read((char*)&Star_s, sizeof(Star_s));
		BitInverse((char*)&Star_s, sizeof(Star_s));
		// 星务对时码 毫秒 2
		unsigned short int Star_ms;
		fs.read((char*)&Star_ms, sizeof(Star_ms));
		BitInverse((char*)&Star_ms, sizeof(Star_ms));

		// 累计秒数 秒 4
		unsigned int GPS_s;
		fs.read((char*)&GPS_s, sizeof(GPS_s));
		BitInverse((char*)&GPS_s, sizeof(GPS_s));
		// 累计秒数 毫秒 2
		unsigned short int GPS_ms;
		fs.read((char*)&GPS_ms, sizeof(GPS_ms));
		BitInverse((char*)&GPS_ms, sizeof(GPS_ms));
		// X轴位置坐标 4字节补码（0.1m） 4
		int x;
		fs.read((char*)&x, sizeof(x));
		BitInverse((char*)&x, sizeof(x));
		// Y轴位置坐标 4字节补码（0.1m） 4
		int y;
		fs.read((char*)&y, sizeof(y));
		BitInverse((char*)&y, sizeof(y));
		// Z轴位置坐标 4字节补码（0.1m） 4
		int z;
		fs.read((char*)&z, sizeof(z));
		BitInverse((char*)&z, sizeof(z));
		// X轴速度坐标 4字节补码（0.01m/s） 4
		int vx;
		fs.read((char*)&vx, sizeof(vx));
		BitInverse((char*)&vx, sizeof(vx));
		// Y轴速度坐标 4字节补码（0.01m/s） 4
		int vy;
		fs.read((char*)&vy, sizeof(vy));
		BitInverse((char*)&vy, sizeof(vy));
		// Z轴速度坐标 4字节补码（0.01m/s） 4
		int vz;
		fs.read((char*)&vz, sizeof(vz));
		BitInverse((char*)&vz, sizeof(vz));

		// GPS对时码 秒 4
		unsigned int cGPS_s;
		fs.read((char*)&cGPS_s, sizeof(cGPS_s));
		BitInverse((char*)&cGPS_s, sizeof(cGPS_s));
		// GPS对时码 毫秒 2
		unsigned short int cGPS_ms;
		fs.read((char*)&cGPS_ms, sizeof(cGPS_ms));
		BitInverse((char*)&cGPS_ms, sizeof(cGPS_ms));

		// 姿控系统时间标准双精度型浮点数 秒 8
		double att_t;
		float d1, d2;
		fs.read((char*)&d1, 4);
		fs.read((char*)&d2, 4);
		BitInverse((char*)&d1, 4);
		BitInverse((char*)&d2, 4);
		memcpy((char*)&att_t+4, &d1, 4);
		memcpy((char*)&att_t, &d2, 4);
		//fs.read((char*)&att_t, sizeof(att_t));
		//BitInverse((char*)&att_t);
		// 滚动角（roll）标准浮点数 弧度 4
		float roll;
		fs.read((char*)&roll, sizeof(roll));
		BitInverse((char*)&roll, sizeof(roll));
		// 俯仰角（pitch）标准浮点数 弧度 4
		float pitch;
		fs.read((char*)&pitch, sizeof(pitch));
		BitInverse((char*)&pitch, sizeof(pitch));
		// 偏航角（yaw）标准浮点数 弧度 4
		float yaw;
		fs.read((char*)&yaw, sizeof(yaw));
		BitInverse((char*)&yaw, sizeof(yaw));
		// 滚动角速度（vroll）标准浮点数 弧度/秒 4
		float vroll;
		fs.read((char*)&vroll, sizeof(vroll));
		BitInverse((char*)&vroll, sizeof(vroll));
		// 俯仰角速度（pitch）标准浮点数 弧度/秒 4
		float vpitch;
		fs.read((char*)&vpitch, sizeof(vpitch));
		BitInverse((char*)&vpitch, sizeof(vpitch));
		// 偏航角速度（yaw）标准浮点数 弧度/秒 4
		float vyaw;
		fs.read((char*)&vyaw, sizeof(vyaw));
		BitInverse((char*)&vyaw, sizeof(vyaw));

		// CCD 行计数值 2
		//int nnn = sizeof(unsigned short);
		//unsigned short iLine;
		short line_count;
		//char c[2];
		//fs.read(c+1, 1);
		//fs.read(c+0, 1);
		//fs.read(c, 2);
		//memcpy((char*)&line_count, c, 2);
		//cout<<c[0]<<","<<c[1]<<endl;

		fs.read((char*)&line_count, sizeof(line_count));
		BitInverse((char*)&line_count, sizeof(line_count));

		// 暗像元 2
		unsigned short Dark_pixel;
		fs.read((char*)&Dark_pixel, sizeof(Dark_pixel));
		BitInverse((char*)&Dark_pixel, sizeof(Dark_pixel));

		// 后 3 字节填零
		fs.ignore(3);

		bState = false;

		//cout<<iLine<<endl;

		//////////////////////////////////////////////////////////////////////////
		aux.star_time = Star_s + Star_ms * 0.001;
		aux.gps_c_time = cGPS_s + cGPS_ms * 0.001;
		aux.gps_time = GPS_s + GPS_ms * 0.001;
		aux.satpos.x = x * 0.1f;
		aux.satpos.y = y * 0.1f;
		aux.satpos.z = z * 0.1f;
		aux.satpos.vx = vx * 0.01f;
		aux.satpos.vy = vy * 0.01f;
		aux.satpos.vz = vz * 0.01f;

		aux.satatt.roll = roll;
		aux.satatt.pitch = pitch;
		aux.satatt.yaw = yaw;
		aux.satatt.vroll = vroll;
		aux.satatt.vpitch = vpitch;
		aux.satatt.vyaw = vyaw;

		aux.att_time = att_t;
		aux.line_count = line_count;
		
		//cout<<line_count<<endl;
	}
	else if (1 == aux.valid_flag)
	{
		//	对应的星下点的位置 8
		fs.read((char*)&aux.nadir_pos.longitude, sizeof(aux.nadir_pos.longitude));
		fs.read((char*)&aux.nadir_pos.latitude, sizeof(aux.nadir_pos.latitude));
		//	该距离线对应的左侧的位置（精度和纬度） 8
		fs.read((char*)&aux.line_left_pos.longitude, sizeof(aux.line_left_pos.longitude));
		fs.read((char*)&aux.line_left_pos.latitude, sizeof(aux.line_left_pos.latitude));
		//	该距离线对应的右侧的位置（精度和纬度） 8
		fs.read((char*)&aux.line_right_pos.longitude, sizeof(aux.line_right_pos.longitude));
		fs.read((char*)&aux.line_right_pos.latitude, sizeof(aux.line_right_pos.latitude));
		//	卫星的x,y,z坐标以及vx,vy,vz 24
		readSAT_POS(fs, aux.satpos);
		//	卫星的姿态参数,俯仰滚动偏航，r,p,y 12
		readSAT_ATT(fs, aux.satatt);
		//  照射该脉冲的时间信息 28 
		readTimeInfo(fs, aux.line_time);
	}
	return bState;
}

int QVProcReader::timeCompare(const _TIMESTAMP& t1, const _TIMESTAMP& t2)const
{
	int diff = 0;
	if (0 != (diff = t1.year-t2.year))
	{
		return diff;
	}
	if (0!= (diff = t1.mon-t2.mon))
	{
		return diff;
	}
	if (0!= (diff = t1.day-t2.day))
	{
		return diff;
	}
	if (0!= (diff = t1.hour-t2.hour))
	{
		return diff;
	}
	if (0!= (diff = t1.min-t2.min))
	{
		return diff;
	}
	if (0!= (diff = t1.sec-t2.sec))
	{
		return diff;
	}
	if (0!= (diff = t1.nsec-t2.nsec))
	{
		return diff;
	}
	 return diff;
}

string QVProcReader::time2String(const _TIMESTAMP& t)const
{
	const int bufSize = 2048;
	char buf[bufSize];
	//sprintf_s(buf, bufSize, "%d-%d-%dT%d:%d:%lf", t.year, t.mon, t.day, t.hour, t.min, t.sec + 0.001*t.millsec + 1e-6*t.microsec);
	sprintf_s(buf, bufSize, "%4d-%2d-%2d %2d:%2d:%9lf", t.year, t.mon, t.day, t.hour, t.min, t.sec + 0.001*t.millsec + 1e-6*t.microsec);
	return string(buf);
}

void QVProcReader::convertTimeStamp(const string& time_stamp,
	double& ti) const
{
	int    year, month, day, hour, minute;
	double second;

	//***
	// Time stamps are in the format: "yyyy-mm-dd hh:mm:ss.ssssss"
	//***
	int converted = sscanf(time_stamp.c_str(),
		"%4d-%2d-%2d %2d:%2d:%9lf",
		&year, &month, &day,
		&hour, &minute, &second);

	if (converted != 6)
	{
		ti = -NAN;
	}
	else
	{
		ti = (((((year - 2002.0)*12.0 + month - 1.0)*365.0 + day - 1.0)*24.0
			+ hour)*60.0 + minute)*60.0 + second;
	}
}

double QVProcReader::convertTimeStamp(const string& time_stamp) const
{
	double ti;
	convertTimeStamp(time_stamp, ti);
	return ti;
}

void QVProcReader::writeDescXML(const QUIVIMAGE_HEAD_INFO& headerInfo, 
								const vector<QUIVIMAGE_AUX_INFO>& auxList, 
								const char* descXML,
								int line_offset/* = 0*/)const
{
	// get a test document
	pugi::xml_document doc;
	doc.load("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");

	pugi::xml_node DATASET_NODE = doc.append_child("DATASET");
	pugi::xml_node BASE_INFO_NODE = DATASET_NODE.append_child("BASE_INFO");
	BASE_INFO_NODE.append_child("PROCESSING").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("STATION_ID").append_child(pugi::node_pcdata).set_value(headerInfo.station_id);
	BASE_INFO_NODE.append_child("SPACECRAFT_ID").append_child(pugi::node_pcdata).set_value(headerInfo.satellite_id);
	BASE_INFO_NODE.append_child("SENSOR_ID").append_child(pugi::node_pcdata).set_value(headerInfo.sensor_id);
	BASE_INFO_NODE.append_child("RESAMPLE_RATE").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.sample_bit_count).c_str());
	BASE_INFO_NODE.append_child("START_TIME").append_child(pugi::node_pcdata).set_value(time2String(headerInfo.start_time).c_str());
	BASE_INFO_NODE.append_child("STOP_TIME").append_child(pugi::node_pcdata).set_value(time2String(headerInfo.end_time).c_str());
	BASE_INFO_NODE.append_child("COLUMNS").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.sample_num).c_str());
	BASE_INFO_NODE.append_child("LINES").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.line_num).c_str());
	BASE_INFO_NODE.append_child("BANDS").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.band_num).c_str());
	BASE_INFO_NODE.append_child("PYRAMID_LEVELS").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.desample_num).c_str());
	BASE_INFO_NODE.append_child("IMAGE_TILE_SIDE").append_child(pugi::node_pcdata).set_value("512");
	BASE_INFO_NODE.append_child("IFCITY").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("IFPROVINCE").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("IFEPH").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("IFNADIR").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("GRID_TILE_SIDE").append_child(pugi::node_pcdata).set_value("32");
	//end of BASE_INFO

	vector<QUIVIMAGE_AUX_INFO>::const_iterator aux;

	// Ephemeris
	pugi::xml_node Ephemeris_NODE = DATASET_NODE.append_child("Ephemeris");
	for (aux = auxList.begin();
		aux != auxList.end();
		aux++)
	{
		//<Point>
		// <TIME>349234357.6785600185</TIME>
		// <Location>
		//  <X>-3763220.5000000000</X>
		//  <Y>2643178.2500000000</Y>
		//  <Z>5469612.5000000000</Z>
		// </Location>
		// <Velocity>
		//  <X>-3372.7299804688</X>
		//  <Y>4885.5000000000</Y>
		//  <Z>-4670.6000976562</Z>
		// </Velocity>
		//</Point>
		pugi::xml_node Point_NODE = Ephemeris_NODE.append_child("Point");
		Point_NODE.append_child("TIME").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->gps_time).c_str());
		pugi::xml_node Location_NODE = Point_NODE.append_child("Location");
		Location_NODE.append_child("X").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satpos.x).c_str());
		Location_NODE.append_child("Y").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satpos.y).c_str());
		Location_NODE.append_child("Z").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satpos.z).c_str());
		pugi::xml_node Velocity_NODE = Point_NODE.append_child("Velocity");
		Velocity_NODE.append_child("X").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satpos.vx).c_str());
		Velocity_NODE.append_child("Y").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satpos.vy).c_str());
		Velocity_NODE.append_child("Z").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satpos.vz).c_str());
	}


	// Attitudes
	pugi::xml_node Attitudes_NODE = DATASET_NODE.append_child("Attitudes");
	for (aux = auxList.begin();
		aux != auxList.end();
		aux++)
	{
		pugi::xml_node Attitude_NODE = Attitudes_NODE.append_child("Attitude");
		Attitude_NODE.append_child("TIME").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->att_time).c_str());
		pugi::xml_node Angle_NODE = Attitude_NODE.append_child("Angle");
		Angle_NODE.append_child("YAW").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satatt.yaw).c_str());
		Angle_NODE.append_child("PITCH").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satatt.pitch).c_str());
		Angle_NODE.append_child("ROLL").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satatt.roll).c_str());
		pugi::xml_node Angle_Velocity_NODE = Attitude_NODE.append_child("Angle_Velocity");
		Angle_Velocity_NODE.append_child("YAW").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satatt.vyaw).c_str());
		Angle_Velocity_NODE.append_child("PITCH").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satatt.vpitch).c_str());
		Angle_Velocity_NODE.append_child("ROLL").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->satatt.vroll).c_str());
	}

	// Line Times
	pugi::xml_node Line_Times_NODE = DATASET_NODE.append_child("Line_Times");
	for (aux = auxList.begin();
		aux != auxList.end();
		aux++)
	{
		pugi::xml_node Line_Time_NODE = Line_Times_NODE.append_child("Line_Time");
		Line_Time_NODE.append_child("TIME").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.10lf") % aux->star_time).c_str());
		Line_Time_NODE.append_child("LINE").append_child(pugi::node_pcdata).set_value(boost::str(boost::format("%.1lf") % (aux->line_num - line_offset)).c_str());
	}

	doc.save_file(descXML, "\t");
}


void QVProcReader::writeDescXML(const QUIVIMAGE_HEAD_INFO& headerInfo, const char* descXML)const
{
	// get a test document
	pugi::xml_document doc;
	doc.load("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");

	pugi::xml_node DATASET_NODE = doc.append_child("DATASET");
	pugi::xml_node BASE_INFO_NODE = DATASET_NODE.append_child("BASE_INFO");
	BASE_INFO_NODE.append_child("PROCESSING").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("STATION_ID").append_child(pugi::node_pcdata).set_value(headerInfo.station_id);
	BASE_INFO_NODE.append_child("SPACECRAFT_ID").append_child(pugi::node_pcdata).set_value(headerInfo.satellite_id);
	BASE_INFO_NODE.append_child("SENSOR_ID").append_child(pugi::node_pcdata).set_value(headerInfo.sensor_id);
	BASE_INFO_NODE.append_child("RESAMPLE_RATE").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.sample_bit_count).c_str());
	BASE_INFO_NODE.append_child("START_TIME").append_child(pugi::node_pcdata).set_value(time2String(headerInfo.start_time).c_str());
	BASE_INFO_NODE.append_child("STOP_TIME").append_child(pugi::node_pcdata).set_value(time2String(headerInfo.end_time).c_str());
	BASE_INFO_NODE.append_child("COLUMNS").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.sample_num).c_str());
	BASE_INFO_NODE.append_child("LINES").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.line_num).c_str());
	BASE_INFO_NODE.append_child("BANDS").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.band_num).c_str());
	BASE_INFO_NODE.append_child("PYRAMID_LEVELS").append_child(pugi::node_pcdata).set_value(boost::lexical_cast<std::string>(headerInfo.desample_num).c_str());
	BASE_INFO_NODE.append_child("IMAGE_TILE_SIDE").append_child(pugi::node_pcdata).set_value("512");
	BASE_INFO_NODE.append_child("IFCITY").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("IFPROVINCE").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("IFEPH").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("IFNADIR").append_child(pugi::node_pcdata).set_value("0");
	BASE_INFO_NODE.append_child("GRID_TILE_SIDE").append_child(pugi::node_pcdata).set_value("32");
	//end of BASE_INFO

	doc.save_file(descXML, "\t");
}