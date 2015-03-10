#include "HJ1QVProcReader.h"
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

namespace QVProc{
HJ1QVProcReader::HJ1QVProcReader()
	:theReferenceLine(0),
	time_ready(false)
{

}
//
//void posEci2Ecr(double& x, double &y, double &z, double t)
//{
//	//ω=7292115×10-11rads-1±0.150×10-11rads-1
//	double omega = 7.292115e-5;
//	double a = t * omega;
//	double sa = sin(a);
//	double ca = cos(a);
//
//	double xx = ca *x + sa * y;
//	double yy = -sa * x + ca * y;
//
//	x = xx;
//	y = yy;
//}
//
//void velEci2Ecr(double x, double y, double z, 
//				 double &vx, double &vy, double &vz, double t)
//{
//	//ω=7292115×10-11rads-1±0.150×10-11rads-1
//	double omega = 7.292115e-5;
//	double a = t * omega;
//	double sa = sin(a);
//	double ca = cos(a);
//
//	double vxx = ca *vx + sa * vy;
//	double vyy = -sa * vx + ca * vy;
//
//	double xx = -omega*sa*x + omega*ca * y;
//	double yy = -omega*ca * x - omega*sa * y;
//
//	vx = vxx + xx;
//	vy = vyy + yy;
//}
//
//void ecef2lla(double x, double y, double z,
//			  double& lat, double& lon, double& alt)
//{
//	double b, ep, p, th, n;
//
//	/* Constants (WGS ellipsoid) */
//	const double a = 6378137;
//	const double e = 8.1819190842622e-2;
//	const double pi = 3.141592653589793;
//
//	/* Calculation */
//	b = sqrt(pow(a,2)*(1-pow(e,2)));
//	ep = sqrt((pow(a,2)-pow(b,2))/pow(b,2));
//	p = sqrt(pow(x,2)+pow(y,2));
//	th = atan2(a*z, b*p);
//	lon = atan2(y, x);
//	lat = atan2((z+ep*ep*b*pow(sin(th),3)), (p-e*e*a*pow(cos(th),3)));
//	n = a/sqrt(1-e*e*pow(sin(lat),2));
//	alt = p/cos(lat)-n;	
//	lat = (lat*180)/pi;
//	lon = (lon*180)/pi;
//}

bool HJ1QVProcReader::read_by_scene(const char* filename, const char* outPath)
{
	int last_pos = 0;

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

	//writeDescXML(theHeaderInfo, (strPath+"\\desc.xml").c_str());

	int percent = 0;
	printf("\r%d%%", percent);
	fflush(stdout);

	int sceneSize = 1*theHeaderInfo.data_width;
	const int OVERLAP = 0;
	int nScenes=(int)(theHeaderInfo.line_num / sceneSize);
	if(nScenes * sceneSize < theHeaderInfo.line_num ) nScenes =nScenes +1;

	//fstream ephFile;
	//ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::out);
	//ephFile.close();

	//fstream nadirsFile;
	//nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::out);
	//nadirsFile.close();

	//GDALDatasetH hDstDS = new GDALDatasetH;
	//const char         *pszFormat = "GTiff";
	//char               *pszTargetSRS = NULL;
	//GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
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
		streampos pos_current = fs.tellg();
		if ((long long)pos_current + 92 + lineDataByte > file_size)
		{
			break;
		}

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
				//int line_num;
				//line_num = iLine - aux.line_count;// - line_offset;
				//ephFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::app);
				//ephFile<<line_num<<" "
				//	<<setprecision(25)
				//	<<aux.satpos.x<<" "
				//	<<aux.satpos.y<<" "
				//	<<aux.satpos.z<<" "
				//	<<aux.satpos.vx<<" "
				//	<<aux.satpos.vy<<" "
				//	<<aux.satpos.vz<<" "
				//	<<aux.satatt.roll<<" "
				//	<<aux.satatt.pitch<<" "
				//	<<aux.satatt.yaw<<" "
				//	<<aux.satatt.vroll<<" "
				//	<<aux.satatt.vpitch<<" "
				//	<<aux.satatt.vyaw<<" "
				//	<<aux.gps_time<<endl;
				//ephFile.close();
				//double lat, lon, alt;
				//ecef2lla(aux.satpos.x, aux.satpos.y, aux.satpos.z, lat, lon, alt);
				//nadirsFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::app);
				//nadirsFile<<line_num<<" "
				//	<<lat<<" "
				//	<<lon<<endl;
				//nadirsFile.close();
			}
		}
#endif
		// 跳过数据区
		fs.seekg(pos_current + streampos(92 + lineDataByte), ios_base::beg);
		// 计算当前进度
		streampos pos_ = fs.tellg(); //读取文件指针的位置
		int tmpPercent = (int)(pos_ / (double)(file_size)* 100 + 0.5);
		if (tmpPercent > percent)
		{
			percent = tmpPercent;
			printf("\r%d%%", percent);
			fflush(stdout);
		}
		iLine++;

		//int iScene = iLine / sceneSize;
		//int iWidth = theHeaderInfo.sample_num;
		//int iHeight = sceneSize;
		//int line_offset = iScene * sceneSize;
		//// 末尾景
		//if (iScene == nScenes - 1)
		//{
		//	iHeight = theHeaderInfo.line_num - line_offset;
		//}
		//if (iScene*sceneSize == iLine)
		//{
		//	char buf[1024];
		//	sprintf_s(buf, "%s\\Scene%02d", strPath.c_str(), iScene+1);
		//	if (!QDir(buf).exists())
		//	{
		//		_mkdir(buf);
		//	}
		//	string sceneDir(buf);
		//	string sceneImage = sceneDir + "\\IMAGE.TIF";

		//	hDstDS = GDALCreate( hDriver, sceneImage.c_str(), iWidth, iHeight, theHeaderInfo.band_num, eDT, NULL );
		//	pBuf = new GByte[iWidth*iHeight*theHeaderInfo.band_num];
		//}
		////cout<<fs.tellg()<<endl;
		//if (file_size - fs.tellg() < lineDataByte || iLine == theHeaderInfo.line_num-1)
		//{
		//	((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, theHeaderInfo.band_num, &panBandMap[0],theHeaderInfo.band_num, lineDataByte, 1);
		//	//int iBand = 1;
		//	//int band1_offset = 3;
		//	//((GDALDataset*)hDstDS)->RasterIO(GF_Write,band1_offset,0,iWidth - band1_offset, iHeight, pBuf, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
		//	//iBand = 2;
		//	//int band2_offset = 4;
		//	//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight-band2_offset, pBuf+lineDataByte*band2_offset, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
		//	////((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,band2_offset,iWidth, iHeight-band2_offset, pBuf, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
		//	//iBand = 3;
		//	//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
		//	GDALClose( hDstDS );
		//	delete[]pBuf;
		//	pBuf = NULL;
		//	//CPLFree( pBuf );
		//	break;
		//}

		//int npos = fs.tellg();
		//fs.read((char*)(pBuf+(iLine-iScene*sceneSize)*iWidth*theHeaderInfo.band_num), lineDataByte);	
		//if (iLine-iScene*sceneSize == iHeight-1)
		//{
		//	((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, theHeaderInfo.band_num, &panBandMap[0],theHeaderInfo.band_num, lineDataByte, 1);
		//	//int iBand = 1;
		//	//int band1_offset = 1;
		//	////((GDALDataset*)hDstDS)->RasterIO(GF_Write,band1_offset,0,iWidth - band1_offset, iHeight, pBuf, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
		//	//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth - band1_offset, iHeight, pBuf+theHeaderInfo.band_num*band1_offset, iWidth-band1_offset, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
		//	//iBand = 2;
		//	//int band2_offset = 2;
		//	//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight-band2_offset, pBuf+lineDataByte*band2_offset, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
		//	////((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,band2_offset,iWidth, iHeight-band2_offset, pBuf, iWidth, iHeight-band2_offset, eDT, 1, &iBand,theHeaderInfo.band_num, lineDataByte, 1);
		//	//iBand = 3;
		//	//((GDALDataset*)hDstDS)->RasterIO(GF_Write,0,0,iWidth, iHeight, pBuf, iWidth, iHeight, eDT, 1, &iBand, theHeaderInfo.band_num, lineDataByte, 1);
		//	GDALClose(hDstDS);
		//	delete[]pBuf;
		//	pBuf = NULL;
		//	//CPLFree( pBuf );
		//}
		//iLine++;

		//// 计算当前进度
		//streampos pos_ = fs.tellg(); //读取文件指针的位置
		//int tmpPercent = (int)(pos_ / (double)(file_size) * 100 + 0.5);
		//if(tmpPercent > percent)
		//{
		//	percent = tmpPercent;
		//	printf("\r%d%%", percent);
		//}
	}
	fs.close();
	writeDescXML(theHeaderInfo, auxList, outPath);

	//// write nScenes aux files respectively
	//for (int iScene = 0;iScene < nScenes;++iScene)
	//{
	//	char buf[1024];
	//	sprintf_s(buf, "%s\\Scene%02d", strPath.c_str(), iScene+1);
	//	if (!QDir(buf).exists())
	//	{
	//		_mkdir(buf);
	//	}
	//	string sceneDir(buf);
	//	
	//	int iWidth = theHeaderInfo.sample_num;
	//	int iHeight = sceneSize;
	//	int line_offset = iScene * sceneSize;
	//	// 末尾景
	//	if (iScene == nScenes - 1)
	//	{
	//		iHeight = theHeaderInfo.line_num - line_offset;
	//	}

	//	QUIVIMAGE_HEAD_INFO sceneHeaderInfo = theHeaderInfo;
	//	sceneHeaderInfo.line_num = iHeight;
	//	writeDescXML(sceneHeaderInfo, (sceneDir+"\\desc.xml").c_str());
	//	int overlap_lines = 3000;

	//	fstream inEphFile;
	//	inEphFile.open((strPath+"\\eph.txt").c_str(), std::ios_base::in);
	//	fstream outEphFile;
	//	outEphFile.open((sceneDir+"\\eph.txt").c_str(), std::ios_base::out);		
	//	char ephLine[2048];
	//	while (inEphFile.getline(ephLine, 2048))
	//	{
	//		int scene_line = atoi(SBeforeFirst(string(ephLine), ' ').c_str()) - line_offset;
	//		if (scene_line > -overlap_lines && scene_line < sceneSize + overlap_lines)
	//		{
	//			outEphFile<<scene_line<<" "
	//				<<SAfterFirst(string(ephLine), ' ')<<endl;
	//		}
	//	}
	//	inEphFile.close();
	//	outEphFile.close();

	//	fstream inNadirFile;
	//	inNadirFile.open((strPath+"\\nadirs.txt").c_str(), std::ios_base::in);
	//	fstream outNadirFile;
	//	outNadirFile.open((sceneDir+"\\nadirs.txt").c_str(), std::ios_base::out);		
	//	char nadirLine[2048];
	//	while (inNadirFile.getline(nadirLine, 2048))
	//	{
	//		int scene_line = atoi(SBeforeFirst(string(nadirLine), ' ').c_str()) - line_offset;
	//		if (scene_line > -overlap_lines && scene_line < sceneSize + overlap_lines)
	//		{
	//			outNadirFile<<scene_line<<" "
	//				<<SAfterFirst(string(nadirLine), ' ')<<endl;
	//		}
	//	}
	//	inNadirFile.close();
	//	outNadirFile.close();

	//	writeDescXML(theHeaderInfo, auxList, (sceneDir+"\\desc.xml").c_str(), line_offset);
	//}

	printf("\r%d%%\n", 100);
	fflush(stdout);
	return true;
}

bool HJ1QVProcReader::readHeadInfo(ifstream& fs, QUIVIMAGE_HEAD_INFO& header)
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

bool HJ1QVProcReader::readTimeInfo(ifstream& fs, _TIMESTAMP& t)
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

bool HJ1QVProcReader::readSAT_POS(ifstream& fs, SAT_POS& pos)
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

bool HJ1QVProcReader::readSAT_ATT(ifstream& fs, SAT_ATT& att)
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

void HJ1QVProcReader::BitInverse(char* chrs, int nBit)
{
	char* tmp = new char[nBit];
	int result = 0;
	for (int i = 0;i < nBit;++i)
	{
		tmp[i] = chrs[nBit-1-i];
	}
	memcpy(chrs, tmp, nBit);
}

bool HJ1QVProcReader::readQuivAuxInfo(ifstream& fs, QUIVIMAGE_AUX_INFO& aux)
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

int HJ1QVProcReader::timeCompare(const _TIMESTAMP& t1, const _TIMESTAMP& t2)const
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

string HJ1QVProcReader::time2String(const _TIMESTAMP& t)const
{
	const int bufSize = 2048;
	char buf[bufSize];
	sprintf_s(buf, bufSize, "%d-%d-%dT%d:%d:%d:%d", t.year, t.mon, t.day, t.hour, t.min, t.sec, t.nsec);
	return string(buf);
}

void HJ1QVProcReader::writeDescXML(const QUIVIMAGE_HEAD_INFO& headerInfo, 
								const vector<QUIVIMAGE_AUX_INFO>& auxList, 
								const char* descXML,
								int line_offset/* = 0*/)const
{
	// get a test document
	pugi::xml_document doc;
	doc.load("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");

	pugi::xml_node DATASET_NODE =  doc.append_child("DATASET");
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

}