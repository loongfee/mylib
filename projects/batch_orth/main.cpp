#include <stdlib.h>
#include <math.h>
#include <direct.h>


/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"

#include <ossim/projection/ossimMapProjection.h>
#include <ossim/base/ossimPreferences.h>

#include <QDir>

#include <fstream>



#include <strUtil.h>
#include <fileUtil.h>
#include <mprojectdefine.h>
#include <ossim/parallel/ossimMpi.h>
//#include <mpi.h>

using namespace std;
using namespace mylib;

#ifndef _WIN64
#pragma comment(lib, "ossim20.lib")
#pragma comment(lib, "blas_win32_MT.lib")
#pragma comment(lib, "lapack_win32_MT.lib")
#else
#pragma comment(lib, "ossim20x64.lib")
#pragma comment(lib, "blas_win64_MT.lib")
#pragma comment(lib, "lapack_win64_MT.lib")
#endif
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "ossim_plugin.lib")
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")
//#pragma comment(lib, "mpi.lib")


ossimFilename virtual_gcpfile = "virtual_gcp.txt";
ossimFilename virtual_chkfile = "virtual_chk.txt";
char *pszScaleType = "DEFAULT";
char *pszSampleType = "BICUBIC";
char *pszOutFileType = "TIFF";
char *pszLogFile = "log.txt";

const char *pszInputPath = "";
const char *pszFilter = "*.rpb";
const char *pszOutPath = "orth";
const char *pszDemPath = "";
const char *pszPluginPath = "";
const char *pszPreferenceFile = "preference.txt";
const char *pszProjectFile = "projection.txt";
bool bAuto = false;

/************************************************************************/
/* 有利函数系数结构                                                      */
/************************************************************************/
struct RPC_Struct
{
	double lineScale;
	double sampScale;
	double latScale;
	double lonScale;
	double hgtScale;
	double lineOffset;
	double sampOffset;
	double latOffset;
	double lonOffset;
	double hgtOffset;
	double lineNumCoef[20];
	double lineDenCoef[20];
	double sampNumCoef[20];
	double sampDenCoef[20];
};

/************************************************************************/
/* 根据图像坐标和高程计算经纬度                                           */
/*
  输入参数： 
      double x：      图像x坐标
	  double y：      图像y坐标
	  double hgt：    高程（米）。如无精确高程数据，可以直接设为0，也
				       可以用图像范围内平均高程（RPC_Struct.hgtOffset）近
					   似，推荐使用后者。
	  RPC_Struct rs： 有利函数系数结构
  输出参数：
      double* lon：   对应的经度
	  double *lat：   对应的纬度
*/
/************************************************************************/
void image2ground(double x, double y, double hgt,
				  RPC_Struct rs,
				  double* lon, double *lat)
{
	double U = (y - rs.lineOffset) / rs.lineScale;
	double V = (x - rs.sampOffset) / rs.sampScale;
	double nhgt = (hgt - rs.hgtOffset) / rs.hgtScale;
	static const int    MAX_NUM_ITERATIONS  = 10;
	static const double CONVERGENCE_EPSILON = 0.1;

	double epsilonU = CONVERGENCE_EPSILON/rs.lineScale;
	double epsilonV = CONVERGENCE_EPSILON/rs.sampScale;
	int iteration = 0;
	double Pu, Qu, Pv, Qv;
	double dPu_dLat, dQu_dLat, dPv_dLat, dQv_dLat;
	double dPu_dLon, dQu_dLon, dPv_dLon, dQv_dLon;
	double Uc, Vc;
	double deltaU, deltaV;
	double dU_dLat, dU_dLon, dV_dLat, dV_dLon, W;
	double deltaLat, deltaLon;
	double nlat = 0.0;
	double nlon = 0.0;
	do
	{
		double P = nlat;
		double L = nlon;
		double H = nhgt;
		const double* c = rs.lineNumCoef;
		Pu = c[ 0]       + c[ 1]*L     + c[ 2]*P     + c[ 3]*H     +
			c[ 4]*L*P   + c[ 5]*L*H   + c[ 6]*P*H   + c[ 7]*L*L   +
			c[ 8]*P*P   + c[ 9]*H*H   + c[10]*L*P*H + c[11]*L*L*L +
			c[12]*L*P*P + c[13]*L*H*H + c[14]*L*L*P + c[15]*P*P*P +
			c[16]*P*H*H + c[17]*L*L*H + c[18]*P*P*H + c[19]*H*H*H;
		c = rs.lineDenCoef;
		Qu = c[ 0]       + c[ 1]*L     + c[ 2]*P     + c[ 3]*H     +
			c[ 4]*L*P   + c[ 5]*L*H   + c[ 6]*P*H   + c[ 7]*L*L   +
			c[ 8]*P*P   + c[ 9]*H*H   + c[10]*L*P*H + c[11]*L*L*L +
			c[12]*L*P*P + c[13]*L*H*H + c[14]*L*L*P + c[15]*P*P*P +
			c[16]*P*H*H + c[17]*L*L*H + c[18]*P*P*H + c[19]*H*H*H;
		c = rs.sampNumCoef;
		Pv = c[ 0]       + c[ 1]*L     + c[ 2]*P     + c[ 3]*H     +
			c[ 4]*L*P   + c[ 5]*L*H   + c[ 6]*P*H   + c[ 7]*L*L   +
			c[ 8]*P*P   + c[ 9]*H*H   + c[10]*L*P*H + c[11]*L*L*L +
			c[12]*L*P*P + c[13]*L*H*H + c[14]*L*L*P + c[15]*P*P*P +
			c[16]*P*H*H + c[17]*L*L*H + c[18]*P*P*H + c[19]*H*H*H;
		c = rs.sampDenCoef;
		Qv = c[ 0]       + c[ 1]*L     + c[ 2]*P     + c[ 3]*H     +
			c[ 4]*L*P   + c[ 5]*L*H   + c[ 6]*P*H   + c[ 7]*L*L   +
			c[ 8]*P*P   + c[ 9]*H*H   + c[10]*L*P*H + c[11]*L*L*L +
			c[12]*L*P*P + c[13]*L*H*H + c[14]*L*L*P + c[15]*P*P*P +
			c[16]*P*H*H + c[17]*L*L*H + c[18]*P*P*H + c[19]*H*H*H;
		Uc = Pu/Qu;
		Vc = Pv/Qv;

		deltaU = U - Uc;
		deltaV = V - Vc;
		if ((fabs(deltaU) > epsilonU) || (fabs(deltaV) > epsilonV))
		{
			c = rs.lineNumCoef;
			dPu_dLat = c[2] + c[4]*L + c[6]*H + 2*c[8]*P + c[10]*L*H + 2*c[12]*L*P +
				c[14]*L*L + 3*c[15]*P*P + c[16]*H*H + 2*c[18]*P*H;
			c = rs.lineDenCoef;
			dQu_dLat = c[2] + c[4]*L + c[6]*H + 2*c[8]*P + c[10]*L*H + 2*c[12]*L*P +
				c[14]*L*L + 3*c[15]*P*P + c[16]*H*H + 2*c[18]*P*H;
			c = rs.sampNumCoef;
			dPv_dLat = c[2] + c[4]*L + c[6]*H + 2*c[8]*P + c[10]*L*H + 2*c[12]*L*P +
				c[14]*L*L + 3*c[15]*P*P + c[16]*H*H + 2*c[18]*P*H;
			c = rs.sampDenCoef;
			dQv_dLat = c[2] + c[4]*L + c[6]*H + 2*c[8]*P + c[10]*L*H + 2*c[12]*L*P +
				c[14]*L*L + 3*c[15]*P*P + c[16]*H*H + 2*c[18]*P*H;
			c = rs.lineNumCoef;
			dPu_dLon = c[1] + c[4]*P + c[5]*H + 2*c[7]*L + c[10]*P*H + 3*c[11]*L*L +
				c[12]*P*P + c[13]*H*H + 2*c[14]*P*L + 2*c[17]*L*H;
			c = rs.lineDenCoef;
			dQu_dLon = c[1] + c[4]*P + c[5]*H + 2*c[7]*L + c[10]*P*H + 3*c[11]*L*L +
				c[12]*P*P + c[13]*H*H + 2*c[14]*P*L + 2*c[17]*L*H;
			c = rs.sampNumCoef;
			dPv_dLon = c[1] + c[4]*P + c[5]*H + 2*c[7]*L + c[10]*P*H + 3*c[11]*L*L +
				c[12]*P*P + c[13]*H*H + 2*c[14]*P*L + 2*c[17]*L*H;
			c = rs.sampDenCoef;
			dQv_dLon = c[1] + c[4]*P + c[5]*H + 2*c[7]*L + c[10]*P*H + 3*c[11]*L*L +
				c[12]*P*P + c[13]*H*H + 2*c[14]*P*L + 2*c[17]*L*H;
			dU_dLat = (Qu*dPu_dLat - Pu*dQu_dLat)/(Qu*Qu);
			dU_dLon = (Qu*dPu_dLon - Pu*dQu_dLon)/(Qu*Qu);
			dV_dLat = (Qv*dPv_dLat - Pv*dQv_dLat)/(Qv*Qv);
			dV_dLon = (Qv*dPv_dLon - Pv*dQv_dLon)/(Qv*Qv);
			W = dU_dLon*dV_dLat - dU_dLat*dV_dLon;
			deltaLat = (dU_dLon*deltaV - dV_dLon*deltaU) / W;
			deltaLon = (dV_dLat*deltaU - dU_dLat*deltaV) / W;
			nlat += deltaLat;
			nlon += deltaLon;
		}
		iteration++;      
	} while (((fabs(deltaU)>epsilonU) || (fabs(deltaV)>epsilonV))
		&& (iteration < MAX_NUM_ITERATIONS));
	if (iteration == MAX_NUM_ITERATIONS)
	{
		printf("WARNING: 达到最大迭代次数，结果可能不精确。\n");
	}
	*lat = nlat * rs.latScale + rs.latOffset;
	*lon = nlon * rs.lonScale + rs.lonOffset;
}

void image2ground(double x, double y, double hgt,
				  ossimRpcModel::rpcModelStruct rs,
				  double* lon, double *lat)
{
	double U = (y - rs.lineOffset) / rs.lineScale;
	double V = (x - rs.sampOffset) / rs.sampScale;
	double nhgt = (hgt - rs.hgtOffset) / rs.hgtScale;
	static const int    MAX_NUM_ITERATIONS  = 10;
	static const double CONVERGENCE_EPSILON = 0.1;

	double epsilonU = CONVERGENCE_EPSILON/rs.lineScale;
	double epsilonV = CONVERGENCE_EPSILON/rs.sampScale;
	int iteration = 0;
	double Pu, Qu, Pv, Qv;
	double dPu_dLat, dQu_dLat, dPv_dLat, dQv_dLat;
	double dPu_dLon, dQu_dLon, dPv_dLon, dQv_dLon;
	double Uc, Vc;
	double deltaU, deltaV;
	double dU_dLat, dU_dLon, dV_dLat, dV_dLon, W;
	double deltaLat, deltaLon;
	double nlat = 0.0;
	double nlon = 0.0;
	do
	{
		double P = nlat;
		double L = nlon;
		double H = nhgt;
		const double* c = rs.lineNumCoef;
		Pu = c[ 0]       + c[ 1]*L     + c[ 2]*P     + c[ 3]*H     +
			c[ 4]*L*P   + c[ 5]*L*H   + c[ 6]*P*H   + c[ 7]*L*L   +
			c[ 8]*P*P   + c[ 9]*H*H   + c[10]*L*P*H + c[11]*L*L*L +
			c[12]*L*P*P + c[13]*L*H*H + c[14]*L*L*P + c[15]*P*P*P +
			c[16]*P*H*H + c[17]*L*L*H + c[18]*P*P*H + c[19]*H*H*H;
		c = rs.lineDenCoef;
		Qu = c[ 0]       + c[ 1]*L     + c[ 2]*P     + c[ 3]*H     +
			c[ 4]*L*P   + c[ 5]*L*H   + c[ 6]*P*H   + c[ 7]*L*L   +
			c[ 8]*P*P   + c[ 9]*H*H   + c[10]*L*P*H + c[11]*L*L*L +
			c[12]*L*P*P + c[13]*L*H*H + c[14]*L*L*P + c[15]*P*P*P +
			c[16]*P*H*H + c[17]*L*L*H + c[18]*P*P*H + c[19]*H*H*H;
		c = rs.sampNumCoef;
		Pv = c[ 0]       + c[ 1]*L     + c[ 2]*P     + c[ 3]*H     +
			c[ 4]*L*P   + c[ 5]*L*H   + c[ 6]*P*H   + c[ 7]*L*L   +
			c[ 8]*P*P   + c[ 9]*H*H   + c[10]*L*P*H + c[11]*L*L*L +
			c[12]*L*P*P + c[13]*L*H*H + c[14]*L*L*P + c[15]*P*P*P +
			c[16]*P*H*H + c[17]*L*L*H + c[18]*P*P*H + c[19]*H*H*H;
		c = rs.sampDenCoef;
		Qv = c[ 0]       + c[ 1]*L     + c[ 2]*P     + c[ 3]*H     +
			c[ 4]*L*P   + c[ 5]*L*H   + c[ 6]*P*H   + c[ 7]*L*L   +
			c[ 8]*P*P   + c[ 9]*H*H   + c[10]*L*P*H + c[11]*L*L*L +
			c[12]*L*P*P + c[13]*L*H*H + c[14]*L*L*P + c[15]*P*P*P +
			c[16]*P*H*H + c[17]*L*L*H + c[18]*P*P*H + c[19]*H*H*H;
		Uc = Pu/Qu;
		Vc = Pv/Qv;

		deltaU = U - Uc;
		deltaV = V - Vc;
		if ((fabs(deltaU) > epsilonU) || (fabs(deltaV) > epsilonV))
		{
			c = rs.lineNumCoef;
			dPu_dLat = c[2] + c[4]*L + c[6]*H + 2*c[8]*P + c[10]*L*H + 2*c[12]*L*P +
				c[14]*L*L + 3*c[15]*P*P + c[16]*H*H + 2*c[18]*P*H;
			c = rs.lineDenCoef;
			dQu_dLat = c[2] + c[4]*L + c[6]*H + 2*c[8]*P + c[10]*L*H + 2*c[12]*L*P +
				c[14]*L*L + 3*c[15]*P*P + c[16]*H*H + 2*c[18]*P*H;
			c = rs.sampNumCoef;
			dPv_dLat = c[2] + c[4]*L + c[6]*H + 2*c[8]*P + c[10]*L*H + 2*c[12]*L*P +
				c[14]*L*L + 3*c[15]*P*P + c[16]*H*H + 2*c[18]*P*H;
			c = rs.sampDenCoef;
			dQv_dLat = c[2] + c[4]*L + c[6]*H + 2*c[8]*P + c[10]*L*H + 2*c[12]*L*P +
				c[14]*L*L + 3*c[15]*P*P + c[16]*H*H + 2*c[18]*P*H;
			c = rs.lineNumCoef;
			dPu_dLon = c[1] + c[4]*P + c[5]*H + 2*c[7]*L + c[10]*P*H + 3*c[11]*L*L +
				c[12]*P*P + c[13]*H*H + 2*c[14]*P*L + 2*c[17]*L*H;
			c = rs.lineDenCoef;
			dQu_dLon = c[1] + c[4]*P + c[5]*H + 2*c[7]*L + c[10]*P*H + 3*c[11]*L*L +
				c[12]*P*P + c[13]*H*H + 2*c[14]*P*L + 2*c[17]*L*H;
			c = rs.sampNumCoef;
			dPv_dLon = c[1] + c[4]*P + c[5]*H + 2*c[7]*L + c[10]*P*H + 3*c[11]*L*L +
				c[12]*P*P + c[13]*H*H + 2*c[14]*P*L + 2*c[17]*L*H;
			c = rs.sampDenCoef;
			dQv_dLon = c[1] + c[4]*P + c[5]*H + 2*c[7]*L + c[10]*P*H + 3*c[11]*L*L +
				c[12]*P*P + c[13]*H*H + 2*c[14]*P*L + 2*c[17]*L*H;
			dU_dLat = (Qu*dPu_dLat - Pu*dQu_dLat)/(Qu*Qu);
			dU_dLon = (Qu*dPu_dLon - Pu*dQu_dLon)/(Qu*Qu);
			dV_dLat = (Qv*dPv_dLat - Pv*dQv_dLat)/(Qv*Qv);
			dV_dLon = (Qv*dPv_dLon - Pv*dQv_dLon)/(Qv*Qv);
			W = dU_dLon*dV_dLat - dU_dLat*dV_dLon;
			deltaLat = (dU_dLon*deltaV - dV_dLon*deltaU) / W;
			deltaLon = (dV_dLat*deltaU - dU_dLat*deltaV) / W;
			nlat += deltaLat;
			nlon += deltaLon;
		}
      iteration++;      
   } while (((fabs(deltaU)>epsilonU) || (fabs(deltaV)>epsilonV))
            && (iteration < MAX_NUM_ITERATIONS));
   if (iteration == MAX_NUM_ITERATIONS)
   {
	   printf("WARNING: 达到最大迭代次数，结果可能不精确。\n");
   }
   *lat = nlat * rs.latScale + rs.latOffset;
   *lon = nlon * rs.lonScale + rs.lonOffset;
}

bool readRPCFile_GF(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct)
{
	std::fstream fs(rpcFile.c_str());
	char buf[2048];
	vector<char> chList;
	chList.push_back(' ');
	chList.push_back(',');
	chList.push_back(';');
	chList.push_back('=');
	chList.push_back('(');
	chList.push_back(')');
	while (fs.getline(buf, 2048))
	{
		ossimString strLine(buf);
		std::vector<ossimString> strList;
		splitString(strLine, chList, strList);
		if (strList[0].contains("lineOffset"))
		{
			rpcStruct.lineOffset = strList[1].toDouble();
		}
		else if (strList[0].contains("sampOffset"))
		{
			rpcStruct.sampOffset = strList[1].toDouble();
		}
		else if (strList[0].contains("latOffset"))
		{
			rpcStruct.latOffset = strList[1].toDouble();
		}
		else if (strList[0].contains("longOffset"))
		{
			rpcStruct.lonOffset = strList[1].toDouble();
		}
		else if (strList[0].contains("heightOffset"))
		{
			rpcStruct.hgtOffset = strList[1].toDouble();
		}
		else if (strList[0].contains("lineScale"))
		{
			rpcStruct.lineScale = strList[1].toDouble();
		}
		else if (strList[0].contains("sampScale"))
		{
			rpcStruct.sampScale = strList[1].toDouble();
		}
		else if (strList[0].contains("latScale"))
		{
			rpcStruct.latScale = strList[1].toDouble();
		}
		else if (strList[0].contains("longScale"))
		{
			rpcStruct.lonScale = strList[1].toDouble();
		}
		else if (strList[0].contains("heightScale"))
		{
			rpcStruct.hgtScale = strList[1].toDouble();
		}
		else if (strList[0].contains("lineNumCoef"))
		{
			for (int i = 0;i < 20;++i)
			{
				fs.getline(buf, 2048);
				ossimString strLine(buf);
				std::vector<ossimString> strList;
				splitString(strLine, chList, strList);
				rpcStruct.lineNumCoef[i] = strList[0].toDouble();
			}
		}
		else if (strList[0].contains("lineDenCoef"))
		{
			for (int i = 0;i < 20;++i)
			{
				fs.getline(buf, 2048);
				ossimString strLine(buf);
				std::vector<ossimString> strList;
				splitString(strLine, chList, strList);
				rpcStruct.lineDenCoef[i] = strList[0].toDouble();
			}
		}
		else if (strList[0].contains("sampNumCoef"))
		{
			for (int i = 0;i < 20;++i)
			{
				fs.getline(buf, 2048);
				ossimString strLine(buf);
				std::vector<ossimString> strList;
				splitString(strLine, chList, strList);
				rpcStruct.sampNumCoef[i] = strList[0].toDouble();
			}
		}
		else if (strList[0].contains("sampDenCoef"))
		{
			for (int i = 0;i < 20;++i)
			{
				fs.getline(buf, 2048);
				ossimString strLine(buf);
				std::vector<ossimString> strList;
				splitString(strLine, chList, strList);
				rpcStruct.sampDenCoef[i] = strList[0].toDouble();
			}
		}
	}
	fs.close();
	rpcStruct.type = 'B';
	//rpcStruct.type = 'A';
	return true;
}

void create3DGridPoints(const ossimDrect& imageBounds,
						const ossimSensorModel& proj,
						double height,
						ossimTieGptSet*& pTieGptSet,
						ossim_uint32 xSamples,
						ossim_uint32 ySamples,
						bool latlon,
						bool shiftTo0Flag)
{
	if (NULL == pTieGptSet)
	{
		pTieGptSet = new ossimTieGptSet;
	}
	ossim_uint32 x,y;
	ossim_float64 w = imageBounds.width();
	ossim_float64 h = imageBounds.height();
	ossimGpt gpt;
	ossimGpt defaultGround;
	if(ySamples < 1) ySamples = 12;
	if(xSamples < 1) xSamples = 12;
	srand(time(0));
	double xnorm;
	double ynorm;
	ossimDpt ul = imageBounds.ul();
	ossimDpt shiftTo0(-ul.x,
		-ul.y);
	for(y = 0; y < ySamples; ++y)
	{
		for(x = 0; x < xSamples; ++x)
		{
			ossimDpt imagePoint;
			if(ySamples > 1)
			{
				ynorm = (double)y/(double)(ySamples - 1);
			}
			else
			{
				ynorm = 0.0;
			}
			if(xSamples > 1)
			{
				xnorm = (double)x/(double)(xSamples - 1);
			}
			else
			{
				xnorm = 0.0;
			}

			ossimDpt dpt((w-1)*xnorm + ul.x,
				(h-1)*ynorm + ul.y);

			proj.lineSampleToWorld(dpt, gpt);
			//proj.lineSampleHeightToWorld(dpt, height, gpt);
			ossimDpt tmpDpt;
			proj.worldToLineSample(gpt, tmpDpt);
			//ossimGpt gptTemp;
			//proj.lineSampleToWorld(dpt,
			//	gptTemp);
			//ossimDpt dptTemp;
			//proj.worldToLineSample(gptTemp, dptTemp);
			//ossimDpt res = dptTemp - dpt;
			//gpt.changeDatum(defaultGround.datum());
			if(shiftTo0Flag)
			{
				imagePoint = dpt + shiftTo0;
			}
			else
			{
				imagePoint = dpt;
			}
			//gpt.hgt = height;
			imagePoint.y += 27000;

			if (latlon)
			{
				//dpt = ossimDpt(gpt.lat, gpt.lon);
				dpt = ossimDpt(gpt.lon, gpt.lat);
			}
			else
			{
				dpt = proj.m_proj->forward(gpt);
			}
			ossimGpt groundPoint(dpt.x,dpt.y, gpt.hgt);

			ossimString strId;
			char tmpStr[256];
			sprintf(tmpStr, "%d", y * xSamples + x + 1);
			strId = tmpStr;
			ossimTieGpt *aTiePt = new ossimTieGpt(groundPoint, imagePoint, 1.0, strId);
			pTieGptSet->addTiePoint(aTiePt);
		}
	}
}

void create3DGridPoints(const double& minLat,
						const double& maxLat,
						const double& minLon,
						const double& maxLon,
						const ossimSensorModel& proj,
						double height,
						ossimTieGptSet*& pTieGptSet,
						ossim_uint32 xSamples,
						ossim_uint32 ySamples,
						bool latlon)
{
	if (NULL == pTieGptSet)
	{
		pTieGptSet = new ossimTieGptSet;
	}

	ossim_uint32 x,y;
	ossim_float64 w = maxLat - minLat;
	ossim_float64 h = maxLon - minLon;
	ossimGpt gpt;
	ossimGpt defaultGround;
	if(ySamples < 1) ySamples = 12;
	if(xSamples < 1) xSamples = 12;
	srand(time(0));
	double xnorm;
	double ynorm;

	for(y = 0; y < ySamples; ++y)
	{
		for(x = 0; x < xSamples; ++x)
		{
			ossimDpt imagePoint;
			if(ySamples > 1)
			{
				ynorm = (double)y/(double)(ySamples - 1);
			}
			else
			{
				ynorm = 0.0;
			}
			if(xSamples > 1)
			{
				xnorm = (double)x/(double)(xSamples - 1);
			}
			else
			{
				xnorm = 0.0;
			}

			ossimGpt gpt((w-1)*xnorm + minLat,
				(h-1)*ynorm + minLon, height);
			ossimDpt dpt;
			proj.worldToLineSample(gpt, dpt);
			imagePoint = dpt;
			//gpt.hgt = height;

			if (latlon)
			{
				//dpt = ossimDpt(gpt.lat, gpt.lon);
				dpt = ossimDpt(gpt.lon, gpt.lat);
			}
			else
			{
				dpt = proj.m_proj->forward(gpt);
			}
			ossimGpt groundPoint(dpt.x,dpt.y, gpt.hgt);

			ossimString strId;
			char tmpStr[256];
			sprintf(tmpStr, "%d", y * xSamples + x + 1);
			strId = tmpStr;
			ossimTieGpt *aTiePt = new ossimTieGpt(groundPoint, imagePoint, 1.0, strId);
			pTieGptSet->addTiePoint(aTiePt);
		}
	}
}

void create3DGridPoints(ossimFilename rpcFile, 
							 ossimFilename imgFile,
							 ossimFilename projectionFile,
							 ossimFilename demFile)
{
	MyProject prj;
	prj.theMgr = ossimElevManager::instance();
	if(!prj.theMgr->loadElevationPath(ossimFilename(demFile)))
	{
		//cout<<"warning: 加载DEM失败！"<<endl;
		//return;
	}

	// 如果需要读取投影
	if (projectionFile.exists())
	{
		prj.m_MapProjection.addFile(projectionFile);
		prj.m_MapPar = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(prj.m_MapProjection));
		// 测试投影有效性
		if(!prj.m_MapPar)
		{
			// 如果无效，则清空
			cout<<"Warning: 无效的投影文件"<<endl;
		}
	}
	else{
		// 无投影文件
		cout<<"Warning: 未找到相应的投影文件投影文件("<<projectionFile.file()<<")"<<endl;
	}

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demFile));
	prj.m_DemPath = ossimFilename(demFile);

	ossimRpcModel::rpcModelStruct rpcStruct;
	readRPCFile_GF(rpcFile, rpcStruct);


	prj.m_ImgFileNameUnc = imgFile;
	prj.m_ModelType = RPCType;
	if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
	{
		cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
		return;
	}

	ossimRpcModel *rpcModel = new ossimRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	
	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(imgFile);
	if(!handler) return;   //应该弹出警告对话框
	ossimTieGptSet* gptSet = new ossimTieGptSet;
	double max_height = rpcStruct.hgtOffset + rpcStruct.hgtScale;
	double min_height = rpcStruct.hgtOffset - rpcStruct.hgtScale;
	double max_lat = rpcStruct.latOffset + rpcStruct.latScale;
	double min_lat = rpcStruct.latOffset - rpcStruct.latScale;
	double max_lon = rpcStruct.lonOffset + rpcStruct.lonScale;
	double min_lon = rpcStruct.lonOffset - rpcStruct.lonScale;
	int nLevels = 1;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height;// + i*(max_height-min_height)/(nLevels-1);
		//create3DGridPoints(min_lat, max_lat, min_lon, max_lon, *prj.m_sensorModel, hgt, gptSet, 15, 15, false);
		//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 3, handler->getImageRectangle().height(), false, false);
		//create3DGridPoints(ossimIrect(0, 2100, 12000, 2110), *prj.m_sensorModel, hgt, gptSet, 50, 5, false, false);
		//create3DGridPoints(ossimIrect(0, 9100, 12000, 9102), *prj.m_sensorModel, hgt, gptSet, 200, 1, false, false);
		//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 100, 3, false, false);
		create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 100, 2, false, false);
	}

	ossimTieGptSet* chkSet = new ossimTieGptSet;
	nLevels = 6;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height + i*(max_height-min_height)/(nLevels-1);
		//create3DGridPoints(min_lat, max_lat, min_lon, max_lon, *prj.m_sensorModel, hgt, chkSet, 20, 20, false);
		create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, chkSet, 20, 20, false, false);
	}

	handler->close();
	delete handler;

	prj.SavePointToFile(virtual_gcpfile, gptSet, NULL);

	prj.OutputReport("report.txt", prj.m_sensorModel, gptSet, chkSet, true);
}

void GF_RPC_Custom(ossimFilename rpcFile, 
				   ossimFilename imgFile,
				   ossimFilename projectionFile,
				   ossimFilename outFile,
				   ossimFilename demFile)
{
	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	prj.ReadGcpAndProjection(ossimFilename(virtual_gcpfile));
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demFile));//
	prj.m_DemPath=ossimFilename(demFile);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = imgFile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);
	prj.InitiateSensorModel(imgFile);

	int num = static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());
	ossimRpcSolver *solver = new ossimRpcSolver(true, false);
	vector < ossimDpt > imagePoints(num);
	vector < ossimGpt > groundControlPoints(num);
	ossimGpt gpt;
	for(int i = 0;i < num;i++)
	{
		groundControlPoints[i] = prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint();
		ossimGpt gpt = prj.m_MapPar->inverse(ossimDpt(groundControlPoints[i].lat, groundControlPoints[i].lon));
		gpt.hgt = groundControlPoints[i].hgt;
		groundControlPoints[i] = gpt;
		imagePoints[i] = prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint();
	}
	solver->solveCoefficients(imagePoints, groundControlPoints);

	ossimImageGeometry *imageGeom = solver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimRpcModel *rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom, "projection.");
	fstream rpcStruct;
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.m_sensorModel->saveState(prj.geom);
	rpcModel->loadState(prj.geom);
	rpcStruct.open("rpcStruct0.txt", ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	
	prj.Orthograph(outFile);
}

void LogUnknown(QString strFileName)
{
	fstream ofs;
	ofs.open(pszLogFile, ios_base::app);
	ofs<<endl<<"unknown file type: "<<strFileName.toStdString();
	ofs.close();
}

void Rpc_GF(ossimFilename rpcFile, 
			ossimFilename imgFile,
			ossimKeywordlist prjKlw,
			ossimFilename outFile,
			ossimFilename demFile)
{
	//create3DGridPoints(rpcFile, imgFile, projectionFile, demFile);
	//GF_RPC_Custom(rpcFile, imgFile, projectionFile, outFile, demFile);
	//return;
	MyProject prj;
	prj.theMgr = ossimElevManager::instance();
	if(!prj.theMgr->loadElevationPath(ossimFilename(demFile)))
	{
		//cout<<"warning: 加载DEM失败！"<<endl;
		//return;
	}

	prj.m_MapProjection = prjKlw;
	prj.m_MapPar = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(prj.m_MapProjection));
	// 测试投影有效性
	if(!prj.m_MapPar)
	{
		// 如果无效，则清空
		cout<<"Warning: 无效的投影文件"<<endl;
	}

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demFile));
	prj.m_DemPath = ossimFilename(demFile);

	ossimRpcModel::rpcModelStruct rpcStruct;
	readRPCFile_GF(rpcFile, rpcStruct);


	prj.m_ImgFileNameUnc = imgFile;
	prj.m_OutBandList.clear();
	//选择输出波段
	//int nbands = 4;
	//for (int i = 0;i < nbands;++i)
	//{
	//	prj.m_OutBandList.push_back(i);
	//}

	prj.m_ModelType = RPCType;
	if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
	{
		cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
		return;
	}

	ossimRpcModel *rpcModel = new ossimRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	/*ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	if(!handler) return;   //应该弹出警告对话框
	prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	handler->close();
	delete handler;*/
	// for debug
	ossimDpt ls;
	prj.m_sensorModel->worldToLineSample(ossimGpt(47.894170,103.603532,1161.015000), ls);
	ossimGpt gpt;
	prj.m_sensorModel->lineSampleToWorld(ossimDpt(10770.0, 12658.5), gpt);
	ossimDpt dpt = prj.m_MapPar->forward(gpt);
	double lon, lat;
	image2ground(10770.0, 12658.5, 1161.3772336225384, rpcStruct, &lon, &lat);
	ossimDpt dpt1 = prj.m_MapPar->forward(ossimGpt(lat, lon, 1161.3772336225384));
	ossimDpt r = dpt - dpt1;

	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.m_ScaleType = pszScaleType;
	prj.m_SampleType = pszSampleType;
	prj.m_FileOutType = pszOutFileType;


#if OSSIM_HAS_MPI
	ossimMpi::instance()->initialize(NULL, NULL);
#endif
	ossimMpi::instance()->setEnabledFlag(true);

	if (ossimMpi::instance()->getRank() == 0)
	{
		ossimNotify(ossimNotifyLevel_INFO)
			<< "MPI running with "
			<< ossimMpi::instance()->getNumberOfProcessors()
			<< " processors..." << std::endl;
	}
	prj.Orthograph(outFile);
	ossimMpi::instance()->finalize();

}


void Rpc_GF(ossimFilename rpcFile, 
			ossimFilename imgFile,
			ossimFilename projectionFile,
			ossimFilename outFile,
			ossimFilename demFile)
{
	ossimKeywordlist prjKwl;
	// 如果需要读取投影
	if (projectionFile.exists())
	{
		prjKwl.addFile(projectionFile);
	}
	else{
		// 无投影文件
		cout<<"Warning: 未找到相应的投影文件投影文件("<<projectionFile.file()<<")"<<endl;
		return;
	}
	Rpc_GF(rpcFile, imgFile, prjKwl, outFile, demFile);
}

void Batch_GF(QString strPath, QStringList filters = QStringList("*.tif"), QString outPath = "", QString pszProjectFile = "projection.txt", QString demPath = "")
{
	cout<<"********************************************************"<<endl
		<<"-------------------"<<endl
		<<"batch information\n"
		<<"-------------------"<<endl
		<<"outpath:\t"<<outPath.toStdString()<<endl
		<<"filter:\t"<<filters[0].toStdString()<<endl
		<<"dem path:\t"<<demPath.toStdString()<<endl
		<<"auto projection:\t"<<bAuto<<endl
		<<"projection file:\t"<<pszProjectFile.toStdString()<<endl
		<<"preference file:\t"<<pszPreferenceFile<<endl
		<<"plugin path:\t"<<pszPluginPath<<endl
		<<"scale type:\t"<<pszScaleType<<endl
		<<"sample type:\t"<<pszSampleType<<endl
		<<"format:\t"<<pszOutFileType<<endl
		<<"********************************************************"<<endl;
	QStringList rpbFiles;
	QStringList allRpbFiles;
	QFindFile(strPath, filters, allRpbFiles);
	for (int i = 0; i < (int)allRpbFiles.size(); i++)
	{
		//QString imgFile = QBeforeLast(allRpbFiles[i], '.') + ".tiff";
		QString imgFile = QBeforeLast(allRpbFiles[i], '.') + ".tif";
		if (QFileInfo(imgFile).exists())
		{
			rpbFiles.push_back(allRpbFiles[i]);
		}
	}
	if (rpbFiles.size() > 0)
	{
		if (!QDir(outPath).exists())
		{
			_mkdir(outPath.toLatin1());
		}
	}

	for(int i = 0;i < rpbFiles.size();i++)
	{
		QString rpbFile = QDir::toNativeSeparators(rpbFiles[i]);
		//QString imgFile = QBeforeLast(rpbFile, '.') + ".tiff";
		QString imgFile = QBeforeLast(rpbFile, '.') + ".tif";
		if (!bAuto)
		{
			QString extention = "tif";
			if (0 == QString(pszOutFileType).toUpper().compare("TIFF"))
			{
				extention = "tif";
			}
			else if(0 == QString(pszOutFileType).toUpper().compare("PIX"))
			{
				extention = "pix";
			}
			else if(0 == QString(pszOutFileType).toUpper().compare("IMG"))
			{
				extention = "img";
			}
			else
			{
				extention = "tif";
			}
			QString outFile = outPath + "\\" + QBeforeLast(QAfterLast(rpbFile, '\\'), '.') + "." + extention;
			cout<<i+1<<" / "<<rpbFiles.size()<<"..."<<endl;
			if (ossimFilename(outFile.toLatin1()).exists())
			{
				//cout<<"here"<<endl;
				continue;
			}

			Rpc_GF(ossimFilename(rpbFiles[i].toLatin1()),
				ossimFilename(imgFile.toLatin1()),
				ossimFilename(pszProjectFile.toLatin1()),
				ossimFilename(outFile.toLatin1()),
				ossimFilename(demPath.toLatin1()));
		}
		else
		{
			QString rpbFile = QDir::toNativeSeparators(rpbFiles[i]);
			string fileName = QAfterLast(rpbFile, '\\').toLatin1();
			vector<string> results;
			SplitString(fileName, "_", results, false);
			// GF1_WFV1_E100.0_N4.8_20140313_L1A0000182461.rpb
			//  0   1     2     3      4          5
			int UTM_Zone = 47;
			if (results.size() < 6)
			{
				continue;
			}
			if(results[2][0] == 'E')
			{
				double lon = QAfterFirst(results[2].c_str(), 'E').toDouble();
				UTM_Zone = (int)((lon - 0.001) / 6.0);
				UTM_Zone = UTM_Zone + 31;
			}
			else if(results[2][0] == 'W')
			{
				double lon = QAfterFirst(results[2].c_str(), 'W').toDouble();
				UTM_Zone = (int)((lon - 0.001) / 6.0);
				UTM_Zone = UTM_Zone + 30;
			}
			else
			{
				continue;
			}

			double metersPerPixel = 16.0;
			pszOutFileType = "TIFF";
			if (ossimString(results[0]).upcase().contains("GF1"))
			{
				// GF1
				if (ossimString(results[1]).upcase().contains("WFV"))
				{
					pszOutFileType = "IMG";
					metersPerPixel = 16.0;
				}
				else if (ossimString(results[1]).upcase().contains("PMS"))
				{
					if (ossimString(results[5]).upcase().contains("MSS"))
					{
						metersPerPixel = 8.0;
					}
					else if (ossimString(results[5]).upcase().contains("PAN"))
					{
						metersPerPixel = 2.0;
					}
					else
					{
						LogUnknown(rpbFile);
						continue;
					}
				}
				else
				{
					LogUnknown(rpbFile);
					continue;
				}
			}
			else if(ossimString(results[0]).upcase().contains("ZY02C"))
			{
				if (ossimString(results[1]).upcase().contains("HRC"))
				{
					metersPerPixel = 2.36;
				}
				else if (ossimString(results[1]).upcase().contains("PMS"))
				{
					if (ossimString(results[5]).upcase().contains("MUX"))
					{
						metersPerPixel = 10.0;
					}
					else if (ossimString(results[5]).upcase().contains("PAN"))
					{
						metersPerPixel = 5.0;
					}
					else
					{
						LogUnknown(rpbFile);
						continue;
					}
				}
				else
				{
					LogUnknown(rpbFile);
					continue;
				}
			}
			else if(ossimString(results[0]).upcase().contains("ZY3"))
			{

				if (ossimString(results[1]).upcase().contains("NAD"))
				{
					metersPerPixel = 2.1;
				}
				else if (ossimString(results[1]).upcase().contains("MUX"))
				{
					metersPerPixel = 5.8;
				}
				else if (ossimString(results[1]).upcase().contains("TLC"))
				{
					if (ossimString(results[5]).upcase().contains("NAD"))
					{
						metersPerPixel = 2.1;
					}
					else if (ossimString(results[5]).upcase().contains("MUX"))
					{
						metersPerPixel = 5.8;
					}
				}
				else
				{
					LogUnknown(rpbFile);
					continue;
				}
			}
			else
			{
				LogUnknown(rpbFile);
				continue;
			}


			QString extention = "tif";
			if (0 == QString(pszOutFileType).toUpper().compare("TIFF"))
			{
				extention = "tif";
			}
			else if(0 == QString(pszOutFileType).toUpper().compare("PIX"))
			{
				extention = "pix";
			}
			else if(0 == QString(pszOutFileType).toUpper().compare("IMG"))
			{
				extention = "img";
			}
			else
			{
				extention = "tif";
			}
			QString outFile = outPath + "\\" + QBeforeLast(QAfterLast(rpbFile, '\\'), '.') + "." + extention;
			cout<<i+1<<" / "<<rpbFiles.size()<<"..."<<endl;
			if (ossimFilename(outFile.toLatin1()).exists())
			{
				//cout<<"here"<<endl;
				continue;
			}

			ossimUtmProjection* utm = new ossimUtmProjection;
			utm->setZone(UTM_Zone);
			utm->setMetersPerPixel(ossimDpt(metersPerPixel, metersPerPixel));
			ossimKeywordlist prjKwl;
			utm->saveState(prjKwl);
			Rpc_GF(ossimFilename(rpbFiles[i].toLatin1()),
				ossimFilename(imgFile.toLatin1()),
				prjKwl,
				ossimFilename(outFile.toLatin1()),
				ossimFilename(demPath.toLatin1()));
		}
	}
}

/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)
{
	printf( 
		"Usage: GetRasterBoundary  [-o outpath] [-f filter] [-dem dempath]\n"
		"\t[-prj projectionfile] [-format outfiletype]\n"
		"\t[-pref preference file] [-plugin plugin path]\n"
		"\t[-scale scaletype] [-sample sampletype]\n"
		"  -o outpath\t: output directory\n"
		"  -f filter\t: filter\n"
		"  -dem dempatht: dem path\n"
		"  -prj projectionfile: projection file\n"
		"  -pref preference file: preference file\n"
		"  -plugin plugin path: plugin path\n"
		"  -auto: auto select projection based on image filename\n"
		"  -scale scaletype: scaletype, e.g. UCHAR, UINT16, UINT8, FLOAT32, FLOAT64...\n"
		"  -sample sampletype: sampletype, e.g. NEAREST_NEIGHBOR, BILINEAR, BICUBIC\n"
		"  -format outfiletype: scaletype, e.g. TIFF, PIX, IMG\n");

	if( pszErrorMsg != NULL )
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

int main( int argc, char** argv )
{
	ossimString cwd = getcwd(NULL, 0);
	int null_value = 0;
	
	if (argc == 0)
	{
		Batch_GF("E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320", QStringList(pszFilter), 
			"E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320\\orth", pszProjectFile, "");
	}
	if (argc == 1)
	{
		pszPreferenceFile = "D:\\opensource\\ossim\\preference.txt";
		ossimInit::instance()->loadPlugins(ossimFilename(pszPluginPath) + "\\ossimgdal_plugin.dll");
		ossimPreferences::instance()->loadPreferences(pszPreferenceFile);
		ossimInit::instance()->initialize();
		//pszOutFileType = "IMG";
		//pszScaleType = "UCHAR";
		//bAuto = true;
		//Batch_GF("G:\\SJ9A", QStringList("*MUX.rpb"), "G:\\SJ9A\\orth", "G:\\SJ9A\\projection10.txt", "");
		//Batch_GF(pszInputPath, QStringList(pszFilter), pszOutPath, pszProjectFile, pszDemPath);
		//Batch_GF("I:\\GF1_test", QStringList("*MSS1.rpb"), "I:\\GF1_test\\orth", "I:\\GF1_test\\projection.txt", "");
		//create3DGridPoints("E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320\\IMAGE0.rpb", 
		//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320\\IMAGE0.tif",
		//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320\\projection.txt",
		//	"D:\\workspace\\dem\\srtm90");
		//Batch_GF("E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320", QStringList(pszFilter), 
		//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320\\orth", "E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320\\projection.txt", "");
		//pszDemPath = "D:\\workspace\\dem\\srtm90";
		pszDemPath = "D:\\workspace\\dem\\aster30";
		//Batch_GF("E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02", QStringList(pszFilter), 
		//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\orth", "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\projection.txt", pszDemPath);
		Batch_GF("E:\\20140509\\pan", QStringList(pszFilter), 
			"E:\\20140509\\pan\\orth", "E:\\20140509\\pan\\projection.txt", pszDemPath);
	}

	if (argc > 1)
	{
		/* -------------------------------------------------------------------- */
		/*      Parse arguments.                                                */
		/* -------------------------------------------------------------------- */
		for( int i = 1; i < argc; i++ )
		{
			if( 0 == _stricmp(argv[i],"-o") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutPath = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-f") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszFilter = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-dem") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszDemPath = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-prj") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszProjectFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-plugin") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszPluginPath = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-pref") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszPreferenceFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-auto") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
				bAuto = true;
			}
			else if( 0 == _stricmp(argv[i],"-scale") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszScaleType = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-sample") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszSampleType = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-format") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutFileType = argv[++i] ;
			}
			else
			{
				Usage();
			}
		}

		ossimInit::instance()->loadPlugins(ossimFilename(pszPluginPath) + "\\ossimgdal_plugin.dll");
		ossimPreferences::instance()->loadPreferences(pszPreferenceFile);
		Batch_GF(pszInputPath, QStringList(pszFilter), pszOutPath, pszProjectFile, pszDemPath);
	}
	return 1;
}