#include "AlosApp.h"

namespace mylib{

bool AlosModel::setProjection()
{
	m_MapProjection.clear();
	ossimString tmp;
	//double a,b;
	char* prefix="";
	//中央纬度
	double centerLat = 0.0;	
	//半球
	char theHemisphere = m_UTMZone[m_UTMZone.size() - 1];
	//投影带
	ossim_uint32 theZone = m_UTMZone.beforePos(m_UTMZone.size() - 1).toInt();
	//中央经线
	double centerLon = (double)((theZone - 30) * 6 - 3);

	m_MapProjection.add(prefix,ossimKeywordNames::ORIGIN_LATITUDE_KW, centerLat,true);
	m_MapProjection.add(prefix,ossimKeywordNames::CENTRAL_MERIDIAN_KW,centerLon,true);

	//使用椭球模型
	m_EllipsoidModel = "WE";
	const 	ossimEllipsoid * theelli=ossimEllipsoidFactory::instance()->create(m_EllipsoidModel.c_str());
	theelli->saveState(m_MapProjection, prefix);

	ossimDpt	theMetersPerPixel;
	//X方向分辨率
	theMetersPerPixel.lon = m_PixelSize;
	//Y方向分辨率
	theMetersPerPixel.lat = m_PixelSize;

	//设置分辨率 单位为米
	m_MapProjection.add(prefix,	ossimKeywordNames::PIXEL_SCALE_XY_KW,theMetersPerPixel.toString().c_str(),true);
	m_MapProjection.add(prefix,	ossimKeywordNames::PIXEL_SCALE_UNITS_KW,	ossimUnitTypeLut::instance()->getEntryString(OSSIM_METERS),true); 

	ossimDpt	theFalseEastingNorthing;
	//东向偏移值
	theFalseEastingNorthing.lon = 500000.0;
	//北向偏移值
	theFalseEastingNorthing.lat = 0.0;

	m_MapProjection.add(prefix,ossimKeywordNames::FALSE_EASTING_NORTHING_KW,theFalseEastingNorthing.toString().c_str(), true);
	m_MapProjection.add(prefix,	ossimKeywordNames::FALSE_EASTING_NORTHING_UNITS_KW,	ossimUnitTypeLut::instance()->getEntryString(OSSIM_METERS), true);
	bool theElevationLookupFlag=true;
	m_MapProjection.add(prefix,	ossimKeywordNames::ELEVATION_LOOKUP_FLAG_KW,	ossimString::toString(theElevationLookupFlag), true);

	//设置UTM投影
	m_MapProjection.add(prefix,ossimKeywordNames::ZONE_KW,theZone,true);
	m_MapProjection.add(prefix,ossimKeywordNames::HEMISPHERE_KW,theHemisphere,true);
	m_MapProjection.add(prefix,	"type","ossimUtmProjection",true);



	m_MapPar=PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(m_MapProjection));
	if(m_MapPar==NULL)
	{
		cout<<"您的投影信息设置有误！"<<endl;
		return false;
	}
	return true;
}

bool AlosModel::readRPCFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct)
{
	fstream fin;
	fin.open(rpcFile.c_str(), ios_base::in);
	string strtmp0;
	std::getline(fin, strtmp0);
	QString strtmp(strtmp0.c_str());

	/*
	―――――――――――――――――――――――――――――――――――――
	FIELD                         SIZE         VALUE RANGE         UNITS
	―――――――――――――――――――――――――――――――――――――
	LINE_OFF                       6      000000 to 999999        pixels
	SAMP_OFF                       5       00000 to 99999         pixels
	LAT_OFF                        8         ±90.0000            degrees
	LONG_OFF                       9        ±180.0000            degrees
	HEIGHT_OFF                     5          ±9999              meters
	LINE_SCALE                     6      000001 to 999999        pixels
	SAMP_SCALE                     5       00001 to 99999         pixels
	LAT_SCALE                      8         ±90.0000            degrees
	LONG_SCALE                     9        ±180.0000            degrees
	HEIGHT_SCALE                   5          ±9999              meters
	LINE_NUM_COEFF1               12       ±9.999999E±9
	･･                         ･･            ･･
	LINE_NUM_COEFF20              12       ±9.999999E±9
	LINE_DEN_COEFF1               12       ±9.999999E±9
	･･                         ･･            ･･
	LINE_DEN_COEFF20              12       ±9.999999E±9
	SAMP_NUM_COEFF1               12       ±9.999999E±9
	･･                         ･･            ･･
	SAMP_NUM_COEFF20              12       ±9.999999E±9
	SAMP_DEN_COEFF1               12       ±9.999999E±9
	･･                         ･･            ･･
	SAMP_DEN_COEFF20              12       ±9.999999E±9
	―――――――――――――――――――――――――――――――――――――
	*/

	//ossim_float64 LINE_OFF;
	//ossim_float64 SAMP_OFF;
	//ossim_float64 LAT_OFF;
	//ossim_float64 LONG_OFF;
	//ossim_float64 HEIGHT_OFF;
	//ossim_float64 LINE_SCALE;
	//ossim_float64 SAMP_SCALE;
	//ossim_float64 LAT_SCALE;
	//ossim_float64 LONG_SCALE;
	//ossim_float64 HEIGHT_SCALE;
	//std::vector<double> LINE_NUM_COEFF;
	//std::vector<double> LINE_DEN_COEFF;
	//std::vector<double> SAMP_NUM_COEFF;
	//std::vector<double> SAMP_DEN_COEFF;
	int nCoeff = 20;

	int pos = 0;
	
	int word_width = 6;
	rpcStruct.lineOffset = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 5;
	rpcStruct.sampOffset = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 8;
	rpcStruct.latOffset = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 9;
	rpcStruct.lonOffset = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 5;
	rpcStruct.hgtOffset = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 6;
	rpcStruct.lineScale = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 5;
	rpcStruct.sampScale = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 8;
	rpcStruct.latScale = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 9;
	rpcStruct.lonScale = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;

	word_width = 5;
	rpcStruct.hgtScale = strtmp.mid(pos, word_width).toDouble();
	pos+=word_width;



	int i;
	for(i = 0;i < nCoeff;i++)
	{
		word_width = 12;
		rpcStruct.lineNumCoef[i] = strtmp.mid(pos, word_width).toDouble();
		pos+=word_width;
	}
	for(i = 0;i < nCoeff;i++)
	{
		word_width = 12;
		rpcStruct.lineDenCoef[i] = strtmp.mid(pos, word_width).toDouble();
		pos+=word_width;
	}
	for(i = 0;i < nCoeff;i++)
	{
		word_width = 12;
		rpcStruct.sampNumCoef[i] = strtmp.mid(pos, word_width).toDouble();
		pos+=word_width;
	}
	for(i = 0;i < nCoeff;i++)
	{
		word_width = 12;
		rpcStruct.sampDenCoef[i] = strtmp.mid(pos, word_width).toDouble();
		pos+=word_width;
	}
	rpcStruct.type = ossimRpcModel::PolynomialType::B;

	return true;
}


AlosPRISM::AlosPRISM(const ossimFilename& AlosDir)
{
	if(!init(AlosDir))
	{
		cout<<"Alos PRISM数据目录不合法！"<<endl;
		return;
	}
	setProjection();
	readRPCFile(m_FileRPC, m_rpcStruct);
	m_bInitState = true;
}
bool AlosPRISM::init(const ossimFilename& AlosDir)
{//Alos PRISM
	QStringList strListTmp;
	QFindFile(AlosDir.c_str(), QStringList("HDR-ALPSMW*_O1B2R_UW.txt"), strListTmp);
	if(strListTmp.size() != 1)
		return false;
	m_FileHDR = strListTmp[0].toLatin1();
	strListTmp.clear();


	QFindFile(AlosDir.c_str(), QStringList("IMG-ALPSMW*_O1B2R_UW.tif"), strListTmp);
	if(strListTmp.size() != 1)
		return false;
	m_FileTIF = strListTmp[0].toLatin1();
	strListTmp.clear();

	QFindFile(AlosDir.c_str(), QStringList("RPC-ALPSMW*_O1B2R_UW.txt"), strListTmp);
	if(strListTmp.size() != 1)
		return false;
	m_FileRPC = strListTmp[0].toLatin1();
	strListTmp.clear();

	fstream ff;
	ff.open(m_FileHDR.c_str(), ios_base::in);

	if(!ff)
		return false;

	ossimString strTmp;
	ossimDpt v[4] = {ossimDpt(0.0,0.0)};
	ossimDpt vUtm[4] = {ossimDpt(0.0,0.0)};
	while(getline(ff,strTmp))
	{
		int pos = (int)strTmp.find_first_of('=');
		ossimString strName = strTmp.substr(0,pos);
		ossimString strValue = strTmp.substr(pos + 1,strTmp.length() - pos - 1);
		int nstart = (int)strValue.find_first_of('"') + 1;
		int nend = (int)strValue.find_last_of('"') - 1;
		strValue = strValue.substr(nstart, nend - nstart + 1);

		if(strName.trim(" ") == "SceneID")
			m_SceneID = strValue.trim(" ");

		if(strName.trim(" ") == "ProductID")
			m_ProductID = strValue.trim(" ");

		if(strName.trim(" ") == "Projection")
			m_Projection = strValue.trim(" ");

		if(strName.trim(" ") == "Datum")
			m_Datum = strValue.trim(" ");

		if(strName.trim(" ") == "EllipsoidModel")
			m_EllipsoidModel = strValue.trim(" ");

		if(strName.trim(" ") == "UTMZone")
			m_UTMZone = strValue.trim(" ");

		if(strName.trim(" ") == "PixelSize")
			m_PixelSize = strValue.toDouble();

		if(strName.trim(" ") == "Columns")
			m_Columns = strValue.toDouble();

		if(strName.trim(" ") == "Lines")
			m_Lines = strValue.toDouble();

		if(strName.trim(" ") == "SceneCenterTime")
			m_SceneCenterTime = strValue.trim(" ");

		if(strName.trim(" ") == "SceneLeftTopEasting")
			vUtm[0].lon = strValue.toDouble();
		if(strName.trim(" ") == "SceneLeftTopNorthing")
			vUtm[0].lat = strValue.toDouble();

		if(strName.trim(" ") == "SceneRightTopEasting")
			vUtm[1].lon = strValue.toDouble();
		if(strName.trim(" ") == "SceneRightTopNorthing")
			vUtm[1].lat = strValue.toDouble();

		if(strName.trim(" ") == "SceneRightBottomEasting")
			vUtm[2].lon = strValue.toDouble();
		if(strName.trim(" ") == "SceneRightBottomNorthing")
			vUtm[2].lat = strValue.toDouble();

		if(strName.trim(" ") == "SceneLeftBottomEasting")
			vUtm[3].lon = strValue.toDouble();
		if(strName.trim(" ") == "SceneLeftBottomNorthing")
			vUtm[3].lat = strValue.toDouble();


		if(strName.trim(" ") == "SceneLeftTopLatitude")
			v[0].lon = strValue.toDouble() - 1.0;
		if(strName.trim(" ") == "SceneLeftTopLongitude")
			v[0].lat = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneRightTopLatitude")
			v[1].lon = strValue.toDouble() - 1.0;
		if(strName.trim(" ") == "SceneRightTopLongitude")
			v[1].lat = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneLeftBottomLatitude")
			v[2].lon = strValue.toDouble() - 1.0;
		if(strName.trim(" ") == "SceneLeftBottomLongitude")
			v[2].lat = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneRightBottomLatitude")
			v[3].lon = strValue.toDouble() - 1.0;
		if(strName.trim(" ") == "SceneRightBottomLongitude")
			v[3].lat = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneCenterLatitude")
			m_theRefImgPt.x = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneCenterLongitude")
			m_theRefImgPt.y = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneCenterNorthing")
			m_theRefGndPt.lat = strValue.toDouble();

		if(strName.trim(" ") == "SceneCenterEasting")
			m_theRefGndPt.lon = strValue.toDouble();
	}
	m_theBoundGndPolygon = ossimPolygon(4, vUtm);
	m_theImageClipRect = ossimPolygon(4,v);
	m_OutBandList.clear();
	m_OutBandList.push_back(0);
	return true;
}

AlosAVNIR2::AlosAVNIR2(const ossimFilename& AlosDir)
{
	if(!init(AlosDir))
	{
		cout<<"Alos AVNIR2数据目录不合法！"<<endl;
		return;
	}
	setProjection();
	readRPCFile(m_FileRPC, m_rpcStruct);
	m_bInitState = true;
}

bool AlosAVNIR2::init(const ossimFilename& AlosDir)
{//Alos AVNIR2
	QStringList strListTmp;
	QFindFile(AlosDir.c_str(), strListTmp, QStringList("HDR-ALAV2*_O1B2R_U.txt"));
	if(strListTmp.size() != 1)
		return false;
	m_FileHDR = strListTmp[0].toLatin1();
	strListTmp.clear();

	//dir.GetAllFiles(AlosDir, &strListTmp, wxT("IMG-01-ALAV2*_O1B2R_UW.tif"));
	//if(strListTmp.size() != 1)
	//	return false;
	//m_FileTIF[0] = ossimFilename(strListTmp[0]);
	//strListTmp.clear();

	//dir.GetAllFiles(AlosDir, &strListTmp, wxT("IMG-02-ALAV2*_O1B2R_UW.tif"));
	//if(strListTmp.size() != 1)
	//	return false;
	//m_FileTIF[1] = ossimFilename(strListTmp[0]);
	//strListTmp.clear();

	//dir.GetAllFiles(AlosDir, &strListTmp, wxT("IMG-03-ALAV2*_O1B2R_UW.tif"));
	//if(strListTmp.size() != 1)
	//	return false;
	//m_FileTIF[2] = ossimFilename(strListTmp[0]);
	//strListTmp.clear();

	//dir.GetAllFiles(AlosDir, &strListTmp, wxT("IMG-04-ALAV2*_O1B2R_UW.tif"));
	//if(strListTmp.size() != 1)
	//	return false;
	//m_FileTIF[3] = ossimFilename(strListTmp[0]);
	//strListTmp.clear();

	QFindFile(AlosDir.c_str(), strListTmp, QStringList("IMG-ALAV2*_O1B2R_U.tif"));
	if(strListTmp.size() != 1)
		return false;
	m_FileTIF = strListTmp[0].toLatin1();
	strListTmp.clear();

	QFindFile(AlosDir.c_str(), strListTmp, QStringList("RPC-ALAV2*_O1B2R_U.txt"));
	if(strListTmp.size() != 1)
		return false;
	m_FileRPC = strListTmp[0].toLatin1();
	strListTmp.clear();

	fstream ff;
	ff.open(m_FileHDR.c_str(), ios_base::in);

	if(!ff)
		return false;

	ossimString strTmp;
	ossimDpt v[4] = {ossimDpt(0.0,0.0)};
	ossimDpt vUtm[4] = {ossimDpt(0.0,0.0)};
	while(getline(ff,strTmp))
	{
		int pos = (int)strTmp.find_first_of('=');
		ossimString strName = strTmp.substr(0,pos);
		ossimString strValue = strTmp.substr(pos + 1,strTmp.length() - pos - 1);
		int nstart = (int)strValue.find_first_of('"') + 1;
		int nend = (int)strValue.find_last_of('"') - 1;
		strValue = strValue.substr(nstart, nend - nstart + 1);

		if(strName.trim(" ") == "SceneID")
			m_SceneID = strValue.trim(" ");

		if(strName.trim(" ") == "ProductID")
			m_ProductID = strValue.trim(" ");

		if(strName.trim(" ") == "Projection")
			m_Projection = strValue.trim(" ");

		if(strName.trim(" ") == "Datum")
			m_Datum = strValue.trim(" ");

		if(strName.trim(" ") == "EllipsoidModel")
			m_EllipsoidModel = strValue.trim(" ");

		if(strName.trim(" ") == "UTMZone")
			m_UTMZone = strValue.trim(" ");

		if(strName.trim(" ") == "PixelSize")
			m_PixelSize = strValue.toDouble();

		if(strName.trim(" ") == "Columns")
			m_Columns = strValue.toDouble();

		if(strName.trim(" ") == "Lines")
			m_Lines = strValue.toDouble();

		if(strName.trim(" ") == "SceneCenterTime")
			m_SceneCenterTime = strValue.trim(" ");

		if(strName.trim(" ") == "SceneLeftTopEasting")
			vUtm[0].lon = strValue.toDouble();
		if(strName.trim(" ") == "SceneLeftTopNorthing")
			vUtm[0].lat = strValue.toDouble();

		if(strName.trim(" ") == "SceneRightTopEasting")
			vUtm[1].lon = strValue.toDouble();
		if(strName.trim(" ") == "SceneRightTopNorthing")
			vUtm[1].lat = strValue.toDouble();

		if(strName.trim(" ") == "SceneRightBottomEasting")
			vUtm[2].lon = strValue.toDouble();
		if(strName.trim(" ") == "SceneRightBottomNorthing")
			vUtm[2].lat = strValue.toDouble();

		if(strName.trim(" ") == "SceneLeftBottomEasting")
			vUtm[3].lon = strValue.toDouble();
		if(strName.trim(" ") == "SceneLeftBottomNorthing")
			vUtm[3].lat = strValue.toDouble();


		if(strName.trim(" ") == "SceneLeftTopLatitude")
			v[0].lon = strValue.toDouble() - 1.0;
		if(strName.trim(" ") == "SceneLeftTopLongitude")
			v[0].lat = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneRightTopLatitude")
			v[1].lon = strValue.toDouble() - 1.0;
		if(strName.trim(" ") == "SceneRightTopLongitude")
			v[1].lat = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneLeftBottomLatitude")
			v[2].lon = strValue.toDouble() - 1.0;
		if(strName.trim(" ") == "SceneLeftBottomLongitude")
			v[2].lat = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneRightBottomLatitude")
			v[3].lon = strValue.toDouble() - 1.0;
		if(strName.trim(" ") == "SceneRightBottomLongitude")
			v[3].lat = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneCenterLatitude")
			m_theRefImgPt.x = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneCenterLongitude")
			m_theRefImgPt.y = strValue.toDouble() - 1.0;

		if(strName.trim(" ") == "SceneCenterNorthing")
			m_theRefGndPt.lat = strValue.toDouble();

		if(strName.trim(" ") == "SceneCenterEasting")
			m_theRefGndPt.lon = strValue.toDouble();
	}
	m_theBoundGndPolygon = ossimPolygon(4, vUtm);
	m_theImageClipRect = ossimPolygon(4,v);
	m_OutBandList.clear();
	m_OutBandList.push_back(0);
	m_OutBandList.push_back(1);
	m_OutBandList.push_back(2);
	m_OutBandList.push_back(3);
	return true;
}
}