#include <iostream>
#include <iterator>
#include <fstream>
using namespace std;

#include "..\func.h"
#include "..\mprojectdefine.h"
#include "..\AlosBatch.h"

//void func1()
//{
//	//wxString sourcePath = wxT("Z:\\share\\data_share\\landsat\\gcp_work\\tflong\\data");
//	rspfFilename sourcePath = "D:\\data";
//	//wxString gcpPath = wxT("Z:\\share\\data_share\\landsat\\gcp_result\\20100512");
//	//wxString gcpPath = wxT("Z:\\share\\data_share\\landsat\\gcp_result");
//	rspfFilename gcpPath = "D:\\Result";
//	//wxString destinationPath = wxT("Z:\\share\\data_share\\landsat\\gcp_work\\tflong\\my_gcp_result");
//	rspfFilename destinationPath = "D:\\Result\\gcp_result";
//	rspfFilename elevationpath = "I:\\share\\dzzw\\demnew\\rspf_dem";
//	wxArrayString sourceFiles;
//	wxArrayString gcpFiles;
//	wxArrayString destinationDirs;
//	getFiles(sourcePath, gcpPath, destinationPath, sourceFiles, gcpFiles, destinationDirs);
//
//	AutoOptimize(sourceFiles, gcpFiles, destinationDirs, elevationpath);
//}
void linefeature_test1()
{
	//rspfFilename elevationpath = "D:\\workspace\\dem";
	//rspfFilename gcpfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\points.txt";
	//rspfFilename linefile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\lines.txt";
	//rspfFilename sourcefile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\header.dat";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\report.txt";
	//rspfFilename outfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\rect.tif";

	

	//rspfFilename elevationpath = "D:\\workspace\\dem";
	//rspfFilename gcpfile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\points.txt";
	//rspfFilename linefile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\lines.txt";
	//rspfFilename sourcefile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\header.dat";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\report.txt";
	//rspfFilename outfile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\rect.tif";


	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	rspfFilename gcpfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003816\\points.txt";
	rspfFilename linefile = "D:\\workspace\\testdata\\LineFeature\\LD2010003816\\lines.txt";
	rspfFilename sourcefile = "D:\\workspace\\testdata\\LineFeature\\LD2010003816\\header.dat";
	rspfFilename reportfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003816\\report_points.txt";
	rspfFilename outfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003816\\rect.tif";


	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	//prj.ReadLineAndProjection(linefile, tieLineList);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(tieLineList);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(6);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);
	prj.InitiateSensorModel(sourcefile);
	rspfSensorModel* sensorModel = prj.m_sensorModel;

	if( NULL == sensorModel )
		return;


	int i;
	vector<rspfTieFeature> tieFeatureList;
	for(i = 0;i < (int)prj.m_CtrlGptSet->size();i++)
	{
		rspfTieFeature tieFeature;
		tieFeature.m_featureType = rspfTieFeature::rspfTiePointPoint;
		tieFeature.m_TiePoints.push_back(*prj.m_CtrlGptSet->getTiePoints()[i]);
		tieFeatureList.push_back(tieFeature);
	}
	for(i = 0;i < (int)tieLineList.size();i++)
	{
		rspfTieFeature tieFeature;
		tieFeature.m_featureType = rspfTieFeature::rspfFeatureType::rspfLineFeature;
		tieFeature.m_TiePoints.push_back(tieLineList[i].first);
		tieFeature.m_TiePoints.push_back(tieLineList[i].second);
		tieFeatureList.push_back(tieFeature);
	}

	prj.UpdateSensorModel(tieFeatureList, sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, sensorModel, prj.geom);
	//prj.UpdateSensorModel(tieLineList, sensorModel, prj.geom);

	//bool optimizeFinished = false;
	//int nIteration = 0;
	//int Max_Iteration = 5;
	//int nControl = 10;
	//int nCheck = 0;
	//double metre_threshold = 40.0;
	//NEWMAT::ColumnVector ctrlResidue;

	//for(; !optimizeFinished; optimizeFinished = prj.CheckSenserModel(sensorModel, tieLineList,
	//	metre_threshold, ctrlResidue))
	//{
	//	if(nIteration++ > Max_Iteration)
	//	{
	//		break;
	//	}

	//	prj.UpdateSensorModel(tieLineList, sensorModel, prj.geom);
	//}

	sensorModel->loadState(prj.geom);

	
	NEWMAT::ColumnVector residue0 = prj.getResidue(sensorModel, tieFeatureList);
	double ki0=residue0.SumSquare();
	cout<<residue0<<endl;
	cout<<endl<<ki0<<endl<<endl;

	prj.OutputReport(reportfile, sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	//prj.Orthograph(outfile);
}

void func3()
{
	rspfInit::instance()->initialize();
	rspfElevManager* theMgr;
	theMgr = rspfElevManager::instance();

	rspfDpt fen(2.0,2.0);

	rspfUtmProjection* UtmMerc = new rspfUtmProjection;
	rspfDpt offset(500000.0, 0.0);
	UtmMerc->setUlEastingNorthing(offset);
	UtmMerc->setZone(49);
	UtmMerc->setMetersPerPixel(fen);//setDatum


	rspfKeywordlist kwl;
	ofstream outf;
	outf.open("D:\\gpt_formated.txt");

	UtmMerc->saveState(kwl);
	outf<<"MAP_PROJECTION_BEGIN"<<endl;
	kwl.print(outf);
	outf<<"MAP_PROJECTION_END"<<endl;

	outf.setf(ios::fixed, ios::floatfield);
	outf.precision(6);


	fstream fs;
	fs.open("D:\\gpt.txt",ios_base::in);
	rspfString strID;
	double x,y;
	double lat,lon,hgt;
	while(fs>>strID>>x>>y>>lat>>lon>>hgt)
	{
		rspfGpt gpt(lat,lon,hgt);

		gpt.hgt=theMgr->getHeightAboveEllipsoid(gpt);
		if (rspf::isnan(gpt.hgt)) gpt.hgt=0;
		rspfDpt dpt = UtmMerc->forward(gpt);

		outf<<strID
			<<"\t"<<x
			<<"\t"<<y
			<<"\t"<<dpt.samp
			<<"\t"<<dpt.line
			<<"\t"<<gpt.hgt<<"\n";
	}
	fs.close();
	outf.close();
}
void func4()
{

	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename gcpfile = "D:\\workspace\\spot5_gpts.txt";
	rspfFilename linefile = "D:\\workspace\\spot5_lines.txt";
	rspfFilename sourcefile = "D:\\workspace\\source\\281268-20040602-2.5\\scene01\\imagery.tif";
	rspfFilename reportfile = "D:\\workspace\\spot5_report.txt";

	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));
	prj.ReadLineAndProjection(linefile, tieLineList);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(tieLineList);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(6);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);
	prj.InitiateSensorModel(sourcefile);
	rspfSensorModel* sensorModel = prj.m_sensorModel;

	if( NULL == sensorModel )
		return;

	int i;
	vector<rspfTieFeature> tieFeatureList;
	for(i = 0;i < (int)prj.m_CtrlGptSet->size();i++)
	{
		rspfTieFeature tieFeature;
		tieFeature.m_featureType = rspfTieFeature::rspfFeatureType::rspfGptFeature;
		tieFeature.m_TiePoints.push_back(*prj.m_CtrlGptSet->getTiePoints()[i]);
		tieFeatureList.push_back(tieFeature);
	}
	for(i = 0;i < (int)tieLineList.size();i++)
	{
		rspfTieFeature tieFeature;
		tieFeature.m_featureType = rspfTieFeature::rspfFeatureType::rspfLineFeature;
		tieFeature.m_TiePoints.push_back(tieLineList[i].first);
		tieFeature.m_TiePoints.push_back(tieLineList[i].second);
		tieFeatureList.push_back(tieFeature);
	}

	prj.UpdateSensorModel(tieFeatureList, sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, sensorModel, prj.geom);

	prj.OutputReport(reportfile, sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
}

void func5()
{

	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename gcpfile = "D:\\gpt_formated.txt";
	rspfFilename sourcefile = "D:\\workspace\\source\\th_cat_100530021840553_1\\imagery.tif";
	rspfFilename reportfile = "D:\\workspace\\theos_report.txt";

	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(tieLineList);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(6);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);
	prj.InitiateSensorModel(sourcefile);
	rspfSensorModel* sensorModel = prj.m_sensorModel;

	if( NULL == sensorModel )
		return;

	prj.UpdateSensorModel(*prj.m_CtrlGptSet, sensorModel, prj.geom);

	prj.OutputReport(reportfile, sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
}

void TheosAffine()
{
	wxString fileName = "Z:\\share\\theos\\TH_CAT_100724090339354_1#1P_S285338_53ARC327_W3464\\TH_CAT_100724090339354_1\\IMAGERY.TIF";
	//wxString fileName = "D:\\workspace\\testdata\\beijing\\001\\imagery.tif";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	//rspfFilename elevationpath = "Z:\\share\\alos_data\\beijing\\bj_rspfdem";
	wxString gcpfile = "Z:\\share\\theos\\TH_CAT_100724090339354_1#1P_S285338_53ARC327_W3464\\TH_CAT_100724090339354_1\\gcp.txt";
	//wxString gcpfile = "D:\\workspace\\testdata\\beijing\\001\\gcp.txt";
	rspfFilename outfile = "Z:\\share\\theos\\TH_CAT_100724090339354_1#1P_S285338_53ARC327_W3464\\TH_CAT_100724090339354_1\\rect.tif";
	//rspfFilename outfile = "D:\\rect.tif";
	//rspfFilename outfile = "D:\\workspace\\testdata\\beijing\\001\\rect.tif";
	rspfFilename reportfile = "Z:\\share\\theos\\TH_CAT_100724090339354_1#1P_S285338_53ARC327_W3464\\TH_CAT_100724090339354_1\\report.txt";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\beijing\\001\\report.txt";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = AffineType;
	prj.InitiateSensorModel(rspfFilename(fileName));

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.Orthograph(outfile);
}

bool readRPCFile(wxString rpcFile, rspfRpcModel::rpcModelStruct& rpcStruct)
{
	ifstream fin(rpcFile.c_str());
	string strtmp;
	getline(fin,strtmp);

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

	//rspf_float64 LINE_OFF;
	//rspf_float64 SAMP_OFF;
	//rspf_float64 LAT_OFF;
	//rspf_float64 LONG_OFF;
	//rspf_float64 HEIGHT_OFF;
	//rspf_float64 LINE_SCALE;
	//rspf_float64 SAMP_SCALE;
	//rspf_float64 LAT_SCALE;
	//rspf_float64 LONG_SCALE;
	//rspf_float64 HEIGHT_SCALE;
	//std::vector<double> LINE_NUM_COEFF;
	//std::vector<double> LINE_DEN_COEFF;
	//std::vector<double> SAMP_NUM_COEFF;
	//std::vector<double> SAMP_DEN_COEFF;
	int nCoeff = 20;

	int pos = 0;

	wxString strtmp1;
	strtmp1 = strtmp.substr(pos, 6);
	pos+=6;
	strtmp1.ToDouble(&rpcStruct.lineOffset);

	strtmp1 = strtmp.substr(pos, 5);
	pos+=5;
	strtmp1.ToDouble(&rpcStruct.sampOffset);

	strtmp1 = strtmp.substr(pos, 8);
	pos+=8;
	strtmp1.ToDouble(&rpcStruct.latOffset);

	strtmp1 = strtmp.substr(pos, 9);
	pos+=9;
	strtmp1.ToDouble(&rpcStruct.lonOffset);

	strtmp1 = strtmp.substr(pos, 5);
	pos+=5;
	strtmp1.ToDouble(&rpcStruct.hgtOffset);

	strtmp1 = strtmp.substr(pos, 6);
	pos+=6;
	strtmp1.ToDouble(&rpcStruct.lineScale);

	strtmp1 = strtmp.substr(pos, 5);
	pos+=5;
	strtmp1.ToDouble(&rpcStruct.sampScale);

	strtmp1 = strtmp.substr(pos, 8);
	pos+=8;
	strtmp1.ToDouble(&rpcStruct.latScale);

	strtmp1 = strtmp.substr(pos, 9);
	pos+=9;
	strtmp1.ToDouble(&rpcStruct.lonScale);

	strtmp1 = strtmp.substr(pos, 5);
	pos+=5;
	strtmp1.ToDouble(&rpcStruct.hgtScale);

	int i;
	for(i = 0;i < nCoeff;i++)
	{
		strtmp1 = strtmp.substr(pos, 12);
		pos+=12;
		strtmp1.ToDouble(&rpcStruct.lineNumCoef[i]);
	}
	for(i = 0;i < nCoeff;i++)
	{
		strtmp1 = strtmp.substr(pos, 12);
		pos+=12;
		strtmp1.ToDouble(&rpcStruct.lineDenCoef[i]);
	}
	for(i = 0;i < nCoeff;i++)
	{
		strtmp1 = strtmp.substr(pos, 12);
		pos+=12;
		strtmp1.ToDouble(&rpcStruct.sampNumCoef[i]);
	}
	for(i = 0;i < nCoeff;i++)
	{
		strtmp1 = strtmp.substr(pos, 12);
		pos+=12;
		strtmp1.ToDouble(&rpcStruct.sampDenCoef[i]);
	}
	rpcStruct.type = rspfRpcModel::PolynomialType::B;

	return true;
}
bool readWorldviewRPCFile(wxString rpcFile, rspfRpcModel::rpcModelStruct& rpcStruct)
{
	rspfWorldviewRpcHeader worldviewRpcHeader;
	worldviewRpcHeader.open(rspfFilename(rpcFile));

	rpcStruct.lineOffset = worldviewRpcHeader.theLineOffset;
	rpcStruct.sampOffset = worldviewRpcHeader.theSampOffset;
	rpcStruct.latOffset = worldviewRpcHeader.theLatOffset;
	rpcStruct.lonOffset = worldviewRpcHeader.theLonOffset;
	rpcStruct.hgtOffset = worldviewRpcHeader.theHeightOffset;
	rpcStruct.lineScale = worldviewRpcHeader.theLineScale;
	rpcStruct.sampScale = worldviewRpcHeader.theSampScale;
	rpcStruct.latScale = worldviewRpcHeader.theLatScale;
	rpcStruct.lonScale = worldviewRpcHeader.theLonScale;
	rpcStruct.hgtScale = worldviewRpcHeader.theHeightScale;


	int nCoeff = 20;
	int i;
	for(i = 0;i < nCoeff;i++)
		rpcStruct.lineDenCoef[i] = worldviewRpcHeader.theLineDenCoeff[i];
	for(i = 0;i < nCoeff;i++)
		rpcStruct.lineNumCoef[i] = worldviewRpcHeader.theLineNumCoeff[i];
	for(i = 0;i < nCoeff;i++)
		rpcStruct.sampDenCoef[i] = worldviewRpcHeader.theSampDenCoeff[i];
	for(i = 0;i < nCoeff;i++)
		rpcStruct.sampNumCoef[i] = worldviewRpcHeader.theSampNumCoeff[i];

	if(worldviewRpcHeader.isAPolynomial())
		rpcStruct.type = rspfRpcModel::PolynomialType::A;
	if(worldviewRpcHeader.isBPolynomial())
		rpcStruct.type = rspfRpcModel::PolynomialType::B;

	return true;
}
bool AlosRpc()
{
	//wxString fileName = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\IMG-01-ALPSMW137562600-O1B1___W";
	wxString fileName = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\merge.tif";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	//wxString fileName = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\01.tif";
	wxString rpcFile = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\RPC-ALPSMW137562600-O1B1___W";
	//wxString rpcFile = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\RPC-01-ALPSMW137562600-O1B1___W";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	wxString gcpfile = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\gcp_partial.txt";
	rspfFilename outfile = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\rect.tif";
	rspfFilename reportfile = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\report.txt";

	//const char *pszFormat = "GTiff";
	//GDALDriver *poDriver;
	//GDALDriverH	 hDriver;
	//char **papszMetadata = NULL;
	//char *pszTargetSRS = NULL;
	//double adfThisGeoTransform[6];

	////GDALDataset *poSrcDS;
	//GDALDatasetH hDstDS;
	//GDALDataset *poDstDS;
	//GDALDataType eDT;
	//CPLErr eErr = CE_None, eErr1 = CE_None;
	//unsigned int  bands;
	//GDALAllRegister();

	//GDALDatasetH hDataset = GDALOpen( fileName, GA_ReadOnly );
	//GDALDataset *poDataset_IMG = (GDALDataset *)hDataset;
	//if (poDataset_IMG == NULL ) {
	//	wxString str_error = fileName + wxT(" 没有顺利打开！"); 
	//	cout<<str_error;
	//	GDALDestroyDriverManager();
	//	return false;
	//}

	//bands = GDALGetRasterCount( poDataset_IMG );

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);

	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;




	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	//prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.Orthograph(outfile);


	//rspfRpcSolver solver(true, false);
	/*vector < rspfDpt > imagePoints;
	vector < rspfGpt > groundControlPoints;
	m_MapPar =NULL;
	m_MapPar = PTR_CAST(rspfMapProjection,
		rspfMapProjectionFactory::instance()->createProjection(m_MapProjection));
	if (!m_MapPar) return false;
	int num = static_cast<int>(m_OptCtrlGptSet->getTiePoints().size());
	rspfGpt tGpt;
	for(int i = 0;i < num;i++)
	{

		tGpt = m_MapPar->inverse_do(rspfDpt(m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().lat,
			m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().lon),m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint());
		tGpt.hgt = m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().hgt;
		groundControlPoints.push_back(tGpt);
		imagePoints.push_back(m_OptCtrlGptSet->getTiePoints()[i]->getImagePoint());
	}*/
	//solver.solveInitialCoefficients()
	//solver.solveCoefficients(imagePoints,groundControlPoints);

	//m_sensorModel = solver.createRpcModel();
	//m_sensorModel->m_proj = m_MapPar;
	//m_sensorModel->saveState(geom);

	return true;
}

bool readIkonosRpcFile(wxString rpcFile, rspfRpcModel::rpcModelStruct& rpcStruct)
{
	FILE *pf;
	fopen_s(&pf, rpcFile, "r+");

	char strDescription[255];
	char strUnit[255];
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.lineOffset, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.sampOffset, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.latOffset, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.lonOffset, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.hgtOffset, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.lineScale, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.sampScale, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.latScale, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.lonScale, &strUnit);
	fscanf_s(pf, "%s %lf %s",&strDescription, &rpcStruct.hgtScale, &strUnit);

	int nCoeff = 20;
	int i;
	for(i = 0;i < nCoeff;++i)
		fscanf_s(pf, "%s %lf",&strDescription, &rpcStruct.lineNumCoef[i]);
	for(i = 0;i < nCoeff;++i)
		fscanf_s(pf, "%s %lf",&strDescription, &rpcStruct.lineDenCoef[i]);
	for(i = 0;i < nCoeff;++i)
		fscanf_s(pf, "%s %lf",&strDescription, &rpcStruct.sampNumCoef[i]);
	for(i = 0;i < nCoeff;++i)
		fscanf_s(pf, "%s %lf",&strDescription, &rpcStruct.sampDenCoef[i]);

	rpcStruct.type = rspfRpcModel::PolynomialType::B;
	return true;
}
void IkonosTest()
{
	wxString fileName = "D:\\workspace\\testdata\\Ikonos\\band1.tif";
	wxString rpcFile = "D:\\workspace\\testdata\\Ikonos\\po_10655_blu_0010000_rpc.txt";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	wxString gcpfile = "D:\\workspace\\testdata\\Ikonos\\cpoy0010000.txt";
	rspfFilename outfile = "D:\\workspace\\testdata\\Ikonos\\rect.tif";
	rspfFilename reportfile = "D:\\workspace\\testdata\\Ikonos\\report.txt";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readIkonosRpcFile(rpcFile, rpcStruct);

	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	//prj.GetElevations(rspfFilename(fileName), prj.m_CtrlGptSet);
	//prj.GetElevations(rspfFilename(fileName), prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;




	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.Orthograph(outfile);


	//rspfRpcSolver solver(true, false);
	/*vector < rspfDpt > imagePoints;
	vector < rspfGpt > groundControlPoints;
	m_MapPar =NULL;
	m_MapPar = PTR_CAST(rspfMapProjection,
	rspfMapProjectionFactory::instance()->createProjection(m_MapProjection));
	if (!m_MapPar) return false;
	int num = static_cast<int>(m_OptCtrlGptSet->getTiePoints().size());
	rspfGpt tGpt;
	for(int i = 0;i < num;i++)
	{

	tGpt = m_MapPar->inverse_do(rspfDpt(m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().lat,
	m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().lon),m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint());
	tGpt.hgt = m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().hgt;
	groundControlPoints.push_back(tGpt);
	imagePoints.push_back(m_OptCtrlGptSet->getTiePoints()[i]->getImagePoint());
	}*/
	//solver.solveInitialCoefficients()
	//solver.solveCoefficients(imagePoints,groundControlPoints);

	//m_sensorModel = solver.createRpcModel();
	//m_sensorModel->m_proj = m_MapPar;
	//m_sensorModel->saveState(geom);
}

void AlosAffine()
{
	wxString fileName = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\merge.tif";
	//wxString fileName = "D:\\workspace\\testdata\\beijing\\001\\imagery.tif";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	//rspfFilename elevationpath = "Z:\\share\\alos_data\\beijing\\bj_rspfdem";
	wxString gcpfile = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\gcp_partial.txt";
	//wxString gcpfile = "D:\\workspace\\testdata\\beijing\\001\\gcp.txt";
	rspfFilename outfile = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\rect.tif";
	//rspfFilename outfile = "D:\\workspace\\testdata\\beijing\\001\\rect.tif";
	rspfFilename reportfile = "D:\\workspace\\testdata\\PRISM_1082600_20080824CEOS1B1(RPC)\\report.txt";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\beijing\\001\\report.txt";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = AffineType;
	prj.InitiateSensorModel(rspfFilename(fileName));

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	//prj.Orthograph(outfile);
}

void AlosNoGcp()
{
	wxString fileName = "D:\\workspace\\testdata\\Alos\\1b2\\d1000771-006\\IMG-ALPSMW230332695_O1B2R_UW.tif";
	//wxString fileName = "D:\\workspace\\testdata\\beijing\\001\\imagery.tif";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	//rspfFilename elevationpath = "Z:\\share\\alos_data\\beijing\\bj_rspfdem";
	wxString gcpfile = "D:\\workspace\\testdata\\Alos\\1b2\\gpt20.txt";
	//wxString gcpfile = "D:\\workspace\\testdata\\beijing\\001\\gcp.txt";
	rspfFilename outfile = "D:\\workspace\\testdata\\Alos\\1b2\\rect_20GCPs.tif";
	//rspfFilename outfile = "D:\\workspace\\testdata\\beijing\\001\\rect.tif";
	rspfFilename reportfile = "D:\\workspace\\testdata\\Alos\\1b2\\report.txt";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\beijing\\001\\report.txt";
	wxString rpcFile = "D:\\workspace\\testdata\\Alos\\1b2\\d1000771-006\\RPC-ALPSMW230332695_O1B2R_UW.txt";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	//prj.Orthograph(outfile);
}
void Alos1b1()
{
	wxString fileName = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\merge1.tif";
	//wxString fileName = "D:\\workspace\\testdata\\beijing\\001\\imagery.tif";
	rspfFilename elevationpath = "";//"D:\\workspace\\dem";
	rspfFilename EGMfile = "";//"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	//rspfFilename elevationpath = "Z:\\share\\alos_data\\beijing\\bj_rspfdem";
	wxString gcpfile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\gcp.txt";
	//wxString gcpfile = "D:\\workspace\\testdata\\beijing\\001\\gcp.txt";
	rspfFilename outfile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\rect_NODEM.tif";
	//rspfFilename outfile = "D:\\workspace\\testdata\\beijing\\001\\rect.tif";
	rspfFilename reportfile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\report_NODEM.txt";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\beijing\\001\\report.txt";
	wxString rpcFile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\RPC-ALPSMW230332690-O1B1___W";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.Orthograph(outfile);
}

void AlosNoGcp_MS()
{
	wxString fileName[4] = {"H:\\testdata\\d1001111-006\\IMG-01-ALAV2A220853080_O1B2R_U.tif",
	"H:\\testdata\\d1001111-006\\IMG-01-ALAV2A220853080_O1B2R_U.tif",
	"H:\\testdata\\d1001111-006\\IMG-03-ALAV2A220853080_O1B2R_U.tif",
	"H:\\testdata\\d1001111-006\\IMG-04-ALAV2A220853080_O1B2R_U.tif"};
	//wxString fileName = "D:\\workspace\\testdata\\beijing\\001\\imagery.tif";
	//rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	//rspfFilename elevationpath = "Z:\\share\\alos_data\\beijing\\bj_rspfdem";
	wxString gcpfile = "H:\\testdata\\d1001111-006\\gcp.txt";
	//wxString gcpfile = "D:\\workspace\\testdata\\beijing\\001\\gcp.txt";
	rspfFilename outfile[4] = {"H:\\testdata\\d1001111-006\\rect_01.tif",
	"H:\\testdata\\d1001111-006\\rect_02.tif",
	"H:\\testdata\\d1001111-006\\rect_03.tif",
	"H:\\testdata\\d1001111-006\\rect_04.tif"};
	//rspfFilename outfile = "D:\\workspace\\testdata\\beijing\\001\\rect.tif";
	rspfFilename reportfile = "H:\\testdata\\d1001111-006\\report.txt";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\beijing\\001\\report.txt";
	wxString rpcFile = "H:\\testdata\\d1001111-006\\RPC-ALAV2A220853080_O1B2R_U.txt";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName[0]);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName[0]));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	for(i = 0;i < 4;i++)
	{
		prj.m_ImgFileNameUnc = rspfFilename(fileName[i]);
		prj.Orthograph(outfile[i]);
	}
}

void AlosLine()
{
	wxString fileName = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\merge1.tif";
	//wxString fileName = "D:\\workspace\\testdata\\beijing\\001\\imagery.tif";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	//rspfFilename elevationpath = "Z:\\share\\alos_data\\beijing\\bj_rspfdem";
	wxString gcpfile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\gcp_.txt";
	wxString linefile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\lines4.txt";
	//wxString gcpfile = "D:\\workspace\\testdata\\beijing\\001\\gcp.txt";
	rspfFilename outfile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\rect_line_4+2.tif";
	//rspfFilename outfile = "D:\\workspace\\testdata\\beijing\\001\\rect.tif";
	rspfFilename reportfile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\report_line_0+4.txt";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\beijing\\001\\report.txt";
	wxString rpcFile = "D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\RPC-ALPSMW230332690-O1B1___W";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;

	prj.ReadGcpAndProjection(rspfFilename(gcpfile));
	prj.ReadLineAndProjection(rspfFilename(linefile), tieLineList);

	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);


	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);
	prj.GetElevations(tieLineList);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	int i;
	vector<rspfTieFeature> tieFeatureList;
	for(i = 0;i < (int)prj.m_CtrlGptSet->size();i++)
	{
		rspfTieFeature tieFeature;
		tieFeature.m_featureType = rspfTieFeature::rspfFeatureType::rspfGptFeature;
		tieFeature.m_TiePoints.push_back(*prj.m_CtrlGptSet->getTiePoints()[i]);
		tieFeatureList.push_back(tieFeature);
	}
	for(i = 0;i < (int)tieLineList.size();i++)
	{
		rspfTieFeature tieFeature;
		tieFeature.m_featureType = rspfTieFeature::rspfFeatureType::rspfLineFeature;
		tieFeature.m_TiePoints.push_back(tieLineList[i].first);
		tieFeature.m_TiePoints.push_back(tieLineList[i].second);
		tieFeatureList.push_back(tieFeature);
	}


	prj.UpdateSensorModel(tieFeatureList, prj.m_sensorModel, prj.geom);

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	//prj.Orthograph(outfile);
}

void block_test_l5()
{
	rspfBlockAdjustment ba;

	MyProject prj1,prj2,prj3,prj4,prj5,prj6;
	const int num = 6;
	const int numTotal = 6;
	MyProject prj[numTotal];
	rspfFilename sourceFile[numTotal] = {	"D:\\workspace\\testdata\\BlockAdjustment\\LD2010003650\\header.dat",
		"D:\\workspace\\testdata\\BlockAdjustment\\LD2010003651\\header.dat",
		"D:\\workspace\\testdata\\BlockAdjustment\\LD2010003813\\header.dat",
		"D:\\workspace\\testdata\\BlockAdjustment\\LD2010003814\\header.dat",
		"D:\\workspace\\testdata\\BlockAdjustment\\LD2010003815\\header.dat",
		"D:\\workspace\\testdata\\BlockAdjustment\\LD2010003816\\header.dat"};
	rspfFilename demPath[numTotal] = {	"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem"};
	rspfFilename chkFile[numTotal] = {	"D:\\workspace\\testdata\\BlockAdjustment\\L5-TM-123-031-20080428-LD2010003650.txt",
		"D:\\workspace\\testdata\\BlockAdjustment\\L5-TM-123-032-20080428-LD2010003651.txt",
		"D:\\workspace\\testdata\\BlockAdjustment\\L5-TM-123-033-20080428-LD2010003813.txt",
		"D:\\workspace\\testdata\\BlockAdjustment\\L5-TM-123-034-20080428-LD2010003814.txt",
		"D:\\workspace\\testdata\\BlockAdjustment\\L5-TM-123-035-20080428-LD2010003815.txt",
		"D:\\workspace\\testdata\\BlockAdjustment\\L5-TM-123-036-20080428-LD2010003816.txt"};
	rspfFilename outFile[numTotal] = {	"D:\\workspace\\testdata\\BlockAdjustment\\rect_123_031.TIF",
		"D:\\workspace\\testdata\\BlockAdjustment\\rect_123_032.TIF",
		"D:\\workspace\\testdata\\BlockAdjustment\\rect_123_033.TIF",
		"D:\\workspace\\testdata\\BlockAdjustment\\rect_123_034.TIF",
		"D:\\workspace\\testdata\\BlockAdjustment\\rect_123_035.TIF",
		"D:\\workspace\\testdata\\BlockAdjustment\\rect_123_036.TIF"};
	int i;
	for (i = 0;i < num;++i)
	{
		InitializePrj(prj[i], sourceFile[i], demPath[i], chkFile[i]);
		ba.addSensorModel(prj[i].m_sensorModel);
	}

	vector< rspfBlockTieGpt > gptList;
	// 加载（已知）控制点
	AddGpts("D:\\workspace\\testdata\\BlockAdjustment\\gpt.txt", gptList);
	// 加载同名点
	AddCpts("D:\\workspace\\testdata\\BlockAdjustment\\cpt.txt", gptList);
	// 进行区域网平差
	ba.adjustment(gptList);
	for(i = 0;i < num;i++)
	{
		prj[i].m_sensorModel->saveState(prj[i].geom);
	}

	// 用优化后的模型输出检查点精度
	fstream fs;
	fs.open("D:\\workspace\\testdata\\BlockAdjustment\\residue.txt",ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(6);
	for (i = 0;i < num;++i)
	{
		AppendResidue2File(prj[i].m_sensorModel, prj[i].m_CtrlGptSet, fs);
		fs<<endl;
	}
	fs.close();

	for(i = 0;i < num;i++)
	{
		prj[i].m_OutBandList.clear();
		prj[i].m_OutBandList.push_back(6);
		prj[i].m_OutBandList.push_back(3);
		prj[i].m_OutBandList.push_back(1);
		prj[i].Orthograph(outFile[i]);
	}
}

void block_test2()
{
	rspfBlockAdjustment ba;

	MyProject prj1,prj2,prj3,prj4,prj5,prj6;
	int num = 4;
	MyProject prj[4];
	rspfFilename sourceFile[4] = {	"D:\\testdata\\037testdata\\LD2010001518\\header.dat",
		"D:\\testdata\\037testdata\\LD2010001813\\L71121037_03720030511_HRF.FST",
		"D:\\testdata\\037testdata\\LD2010002712\\header.dat",
		"D:\\testdata\\037testdata\\LD2010001401\\L71123037_03720030525_HRF.FST"};
	rspfFilename demPath[4] = {	"I:\\share\\dzzw\\demnew\\rspf_dem",
		"I:\\share\\dzzw\\demnew\\rspf_dem",
		"I:\\share\\dzzw\\demnew\\rspf_dem",
		"I:\\share\\dzzw\\demnew\\rspf_dem"};
	rspfFilename chkFile[4] = {	"D:\\testdata\\037testdata\\L5-TM-120-037-20090426-LD2010001518.txt",
		"D:\\testdata\\037testdata\\L7-ETM+-121-037-20030511-LD2010001813.txt",
		"D:\\testdata\\037testdata\\L5-TM-122-037-20091017-LD2010002712.txt",
		"D:\\testdata\\037testdata\\L7-ETM+-123-037-20030525-LD2010001401.txt"};

	prj[0].m_ModelType = RPCType;
	for (int i = 0;i < num;++i)
	{
		//prj[i].m_ModelType = ModelType::PolynomialType;
		InitializePrj(prj[i], sourceFile[i], demPath[i], chkFile[i]);
		ba.addSensorModel(prj[i].m_sensorModel);
	}

	vector< rspfBlockTieGpt > gptList;
	// 加载（已知）控制点
	AddGpts("D:\\testdata\\037testdata\\gpt.txt", gptList);
	// 加载同名点
	AddCpts("D:\\testdata\\037testdata\\cpt.txt", gptList);
	// 进行区域网平差
	ba.adjustment(gptList);

	// 用优化后的模型输出检查点精度
	fstream fs;
	fs.open("D:\\testdata\\037testdata\\residue.txt",ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(6);
	for (int i = 0;i < num;++i)
	{
		AppendResidue2File(prj[i].m_sensorModel, prj[i].m_CtrlGptSet, fs);
	}
	fs.close();
}

void block_test3()
{
	MyProject prj;
	prj.m_ModelType = RPCType;
	InitializePrj(prj, "D:\\testdata\\037testdata\\LD2010001518\\header.dat",
		"I:\\share\\dzzw\\demnew\\rspf_dem",
		"D:\\testdata\\037testdata\\L5-TM-120-037-20090426-LD2010001518.txt");

	prj.UpdateSensorModel(*prj.m_OptCtrlGptSet, prj.m_sensorModel, prj.geom);
	cout<<prj.geom;

	fstream fs;
	fs.open("D:\\testdata\\037testdata\\residue.txt",ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(6);
	AppendResidue2File(prj.m_sensorModel, prj.m_CtrlGptSet, fs);
	fs.close();
}
void block_alos()
{
	rspfBlockAdjustment ba;

	MyProject prj1,prj2,prj3,prj4,prj5,prj6;
	const int num = 5;
	const int numTotal = 5;
	MyProject prj[numTotal];
	string foldName = "results_without001002003005";
	rspfFilename sourceFile[numTotal] = {	
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\merge1.tif",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-002\\merge1.tif",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-003\\merge1.tif",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-004\\merge1.tif",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-005\\merge1.tif"};
	rspfFilename demPath[numTotal] = {	
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem"};
	rspfFilename chkFile[numTotal] = {	
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\gcp_std.txt",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-002\\gcp_std.txt",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-003\\gcp_std.txt",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-004\\gcp_std.txt",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-005\\gcp_std.txt"};
	rspfFilename outFile[numTotal] = {	
		"D:\\workspace\\testdata\\Alos\\Block_data\\" + foldName + "\\block_rect001_.TIF",
		"D:\\workspace\\testdata\\Alos\\Block_data\\" + foldName + "\\block_rect002_.TIF",
		"D:\\workspace\\testdata\\Alos\\Block_data\\" + foldName + "\\block_rect003_.TIF",
		"D:\\workspace\\testdata\\Alos\\Block_data\\" + foldName + "\\block_rect004_.TIF",
		"D:\\workspace\\testdata\\Alos\\Block_data\\" + foldName + "\\block_rect005_.TIF"};
	rspfFilename egmPath[numTotal] = {
		"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC",
		"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC",
		"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC",
		"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC",
		"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC"};
	rspfFilename rpcFile[numTotal] = {
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-001\\RPC-ALPSMW230332690-O1B1___W",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-002\\RPC-ALPSMW230332695-O1B1___W",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-003\\RPC-ALPSMW230332700-O1B1___W",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-004\\RPC-ALPSMW232812690-O1B1___W",
		"D:\\workspace\\testdata\\Alos\\Block_data\\d1000771-005\\RPC-ALPSMW232812695-O1B1___W"};
	int i;
	for (i = 0;i < num;++i)
	{
		rspfRpcProjection rpcProjection;
		rspfRpcModel::rpcModelStruct rpcStruct;
		readRPCFile(rpcFile[i], rpcStruct);

		rspfKeywordlist MapProjection;
		prj[i].theMgr = rspfElevManager::instance();
		prj[i].theMgr->loadElevationPath(rspfFilename(demPath[i]));
		prj[i].theMgr->openMGH(egmPath[i]);
		prj[i].m_DemPath=rspfFilename(demPath[i]);
		prj[i].ReadGcpAndProjection(rspfFilename(chkFile[i]));

		prj[i].GetElevations(prj[i].m_CtrlGptSet);
		prj[i].GetElevations(prj[i].m_ChkGptSet);

		rspfMapProjection *MapPar = prj[i].m_MapPar;
		MapProjection = prj[i].m_MapProjection;

		prj[i].m_ImgFileNameUnc = rspfFilename(sourceFile[i]);
		prj[i].m_OutBandList.clear();
		prj[i].m_OutBandList.push_back(0);

		prj[i].m_ModelType = RPCType;
		prj[i].InitiateSensorModel(rspfFilename(sourceFile[i]));
		rspfRpcModel *rpcModel = new rspfRpcModel;
		rpcModel->setAttributes(rpcStruct);
		prj[i].m_sensorModel = rpcModel;
		prj[i].m_sensorModel->m_proj = prj[i].m_MapPar;

		ba.addSensorModel(prj[i].m_sensorModel);
	}

	vector< rspfBlockTieGpt > gptList;
	// 加载（已知）控制点
	AddGpts(rspfFilename("D:\\workspace\\testdata\\Alos\\Block_data\\" + foldName + "\\gpt_std.txt"), gptList);
	prj[0].GetElevations(gptList);
	// 加载同名点
	AddCpts(rspfFilename("D:\\workspace\\testdata\\Alos\\Block_data\\" + foldName + "\\tpt.txt"), gptList);
	// 进行区域网平差
	ba.adjustment(gptList);
	for(i = 0;i < num;i++)
	{
		prj[i].m_sensorModel->saveState(prj[i].geom);
	}

	wxString GptFile = "D:\\workspace\\testdata\\Alos\\Block_data\\" + foldName + "\\GptList.txt";
	// 保存控制点和连接点
	saveGptList(GptFile, prj, gptList);

	// 用优化后的模型输出检查点精度
	fstream fs;
	fs.open("D:\\workspace\\testdata\\Alos\\Block_data\\" + wxString(foldName) + "\\residue_.txt",ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(6);
	for (i = 0;i < num;++i)
	{
		AppendResidue2File(prj[i].m_sensorModel, prj[i].m_CtrlGptSet, fs);
		fs<<endl;
	}
	fs.close();

	for(i = 2;i < num;i++)
	{
		prj[i].m_OutBandList.clear();
		//prj[i].m_OutBandList.push_back(6);
		//prj[i].m_OutBandList.push_back(3);
		//prj[i].m_OutBandList.push_back(1);
		prj[i].m_OutBandList.push_back(0);
		prj[i].Orthograph(outFile[i]);
	}
}
void block_testSpot()
{
	//// sopt5 block adjustment start
	//rspfBlockAdjustment ba;

	//// 依次添加SensorModel
	//rspfSensorModel* sensorModel0 = getSensorModel("I:\\share\\testdata\\277270-20040601-10\\scene01\\imagery.tif");
	//rspfSensorModel* sensorModel1 = getSensorModel("I:\\share\\testdata\\278270-20040728-10\\scene01\\imagery.tif");
	//rspfSensorModel* sensorModel2 = getSensorModel("I:\\share\\testdata\\279270-20040523-10\\scene01\\imagery.tif", 115.0);	
	//ba.addSensorModel(sensorModel0);
	//ba.addSensorModel(sensorModel1);
	//ba.addSensorModel(sensorModel2);

	//vector< rspfBlockTieGpt > gptList;
	//// 加载（已知）控制点
	//AddGpts("I:\\share\\testdata\\gpt.txt", gptList);
	//// 加载同名点
	//AddCpts("I:\\share\\testdata\\cpt.txt", gptList);
	//// 进行区域网平差
	//ba.adjustment(gptList);

	//// 用优化后的模型输出检查点精度
	//fstream fs;
	//fs.open("I:\\share\\testdata\\residue.txt",ios_base::out);
	//fs.setf(ios::fixed, ios::floatfield);
	//fs.precision(6);
	//AppendResidueFromGcpFile(sensorModel0, "I:\\share\\testdata\\277270.txt", fs);
	//AppendResidueFromGcpFile(sensorModel1, "I:\\share\\testdata\\278270.txt", fs);
	//AppendResidueFromGcpFile(sensorModel2, "I:\\share\\testdata\\279270.txt", fs);	
	//fs.close();
	//// sopt5 block adjustment end
}

void block_testL5()
{
	//// TM block adjustment start
	//rspfBlockAdjustment ba;

	//// 依次添加SensorModel
	//rspfSensorModel* sensorModel0 = getSensorModel("I:\\share\\balance_data\\12340-t20061101\\header.dat");
	//rspfSensorModel* sensorModel1 = getSensorModel("I:\\share\\balance_data\\12440-t20041118\\header.dat");
	//rspfSensorModel* sensorModel2 = getSensorModel("I:\\share\\balance_data\\12540-t20070814\\header.dat");
	//ba.addSensorModel(sensorModel0);
	//ba.addSensorModel(sensorModel1);
	//ba.addSensorModel(sensorModel2);

	//vector< rspfBlockTieGpt > gptList;
	//// 加载（已知）控制点
	//AddGpts("I:\\share\\balance_data\\gpt.txt", gptList);
	//// 加载同名点
	//AddCpts("I:\\share\\balance_data\\cpt.txt", gptList);
	//// 进行区域网平差
	//ba.adjustment(gptList);

	//// 用优化后的模型输出检查点精度
	//fstream fs;
	//fs.open("I:\\share\\balance_data\\residue.txt",ios_base::out);
	//fs.setf(ios::fixed, ios::floatfield);
	//fs.precision(6);
	//AppendResidueFromGcpFile(sensorModel0, "I:\\share\\balance_data\\12340chk.txt", fs);
	//AppendResidueFromGcpFile(sensorModel1, "I:\\share\\balance_data\\12440chk.txt", fs);
	//AppendResidueFromGcpFile(sensorModel2, "I:\\share\\balance_data\\12540chk.txt", fs);
	//fs.close();
	//// TM block adjustment end
}


rspfRefPtr<rspfIFStream> theInputStream;

bool getMGH(double lon, double lat, double& mgh)
{
	if(!theInputStream.valid())
	{
		return rspf::nan();
	}
	if(theInputStream->fail())
	{
		theInputStream->clear();
		theInputStream->seekg(0);
	}
	rspfEndian endian;

	double xRes = 4.0;
	double yRes = 4.0;
	double theWidth = 360 * xRes;
	double theHeight = 180 * yRes + 1.0;

	double xi = lon * xRes;
	double yi = lat * yRes;
	rspf_sint64 x0 = static_cast<rspf_sint64>(xi);
	rspf_sint64 y0 = static_cast<rspf_sint64>(yi);

	double xt0 = xi - x0;
	double yt0 = yi - y0;
	double xt1 = 1-xt0;
	double yt1 = 1-yt0;

	double w00 = xt1*yt1;
	double w01 = xt0*yt1;
	double w10 = xt1*yt0;
	double w11 = xt0*yt0;


	if ( xi < 0 || yi < 0 ||
		x0 > (theWidth  - 1.0) ||
		y0 > (theHeight  - 1.0) )
	{
		return rspf::nan();
	}

	if(x0 == (theWidth  - 1.0))
	{
		--x0;
	}
	if(y0 == (theHeight  - 1.0))
	{
		--y0;
	}

	rspf_sint16 p[4];

	rspf_uint64 bytesPerLine  = theWidth * sizeof(rspf_sint16);

	std::streampos offset = y0*bytesPerLine + x0*sizeof(rspf_sint16) + 1;

	theInputStream->seekg(offset, ios::beg);
	theInputStream->read((char*)p, sizeof(rspf_sint16));

	// Get the second post.
	theInputStream->read((char*)(p+1), sizeof(rspf_sint16));

	//   offset += (bytesPerLine-2*sizeof(rspf_uint16));

	theInputStream->ignore(bytesPerLine-2*sizeof(rspf_sint16));
	// Get the third post.
	theInputStream->read((char*)(p+2), sizeof(rspf_sint16));

	// Get the fourth post.
	theInputStream->read((char*)(p+3), sizeof(rspf_sint16));

	if(theInputStream->fail())
	{
		theInputStream->clear();
		return rspf::nan();
	}

	double p00 = p[0];
	double p01 = p[1];
	double p10 = p[2];
	double p11 = p[3];

	if (p00<-10000.0)  return rspf::nan();
	if (p01<-10000.0)  return rspf::nan();
	if (p10<-10000.0)  return rspf::nan();
	if (p11<-10000.0)  return rspf::nan();

	if (p00>10000.0)  return rspf::nan();
	if (p01>10000.0)  return rspf::nan();
	if (p10>10000.0)  return rspf::nan();
	if (p11>10000.0)  return rspf::nan();

	//if (p00 == info.theNullHeightValue)
	//	w00 = 0.0;
	//if (p01 == info.theNullHeightValue)
	//	w01 = 0.0;
	//if (p10 == info.theNullHeightValue)
	//	w10 = 0.0;
	//if (p11 == info.theNullHeightValue)
	//	w11 = 0.0;

	double sum_weights = w00 + w01 + w10 + w11;

	if (sum_weights)
	{
		return (p00*w00 + p01*w01 + p10*w10 + p11*w11) / sum_weights;
	}

	return rspf::nan();
}



double getMGH1(ifstream& fs, const rspfGpt& gpt)
{
	//if(!theMGHInputStream.valid())
	//{
	//	return rspf::nan();
	//}
	//if(theMGHInputStream->fail())
	//{
	//	theMGHInputStream->clear();
	//	theMGHInputStream->seekg(0);
	//}

	double xRes = 4.0;
	double yRes = 4.0;
	double theWidth = 360 * xRes;
	double theHeight = 180 * yRes + 1.0;
	double theScale = 100.0;

	double xi;
	if(gpt.lon >= 0)
		xi = gpt.lon * xRes;
	else
		xi = (gpt.lon + 360.0) * xRes;
	double yi = (90.0 - gpt.lat) * yRes;

	rspf_sint64 x0 = static_cast<rspf_sint64>(xi);
	rspf_sint64 y0 = static_cast<rspf_sint64>(yi);

	double xt0 = xi - x0;
	double yt0 = yi - y0;
	double xt1 = 1-xt0;
	double yt1 = 1-yt0;

	double w00 = xt1*yt1;
	double w01 = xt0*yt1;
	double w10 = xt1*yt0;
	double w11 = xt0*yt0;


	if ( xi < 0 || yi < 0 ||
		x0 > (theWidth  - 1.0) ||
		y0 > (theHeight  - 1.0) )
	{
		return rspf::nan();
	}

	if(x0 == (theWidth  - 1.0))
	{
		--x0;
	}
	if(y0 == (theHeight  - 1.0))
	{
		--y0;
	}

	rspf_sint16 p[4];
	rspf_uint8 tmp[2];

	rspf_uint64 bytesPerLine  = theWidth * sizeof(rspf_sint16);

	std::streampos offset = y0*bytesPerLine + x0*sizeof(rspf_sint16);

	fs.seekg(offset, ios::beg);
	fs.read((char*)(tmp+1), sizeof(rspf_uint8));
	fs.read((char*)(tmp), sizeof(rspf_uint8));
	*p = *(rspf_sint16*)(tmp);

	// Get the second post.
	fs.read((char*)(tmp+1), sizeof(rspf_uint8));
	fs.read((char*)(tmp), sizeof(rspf_uint8));
	*(p+1) = *(rspf_sint16*)(tmp);

	//   offset += (bytesPerLine-2*sizeof(T));

	fs.ignore(bytesPerLine-2*sizeof(rspf_sint16));
	// Get the third post.
	fs.read((char*)(tmp+1), sizeof(rspf_uint8));
	fs.read((char*)(tmp), sizeof(rspf_uint8));
	*(p+2) = *(rspf_sint16*)(tmp);

	// Get the fourth post

	fs.read((char*)(tmp+1), sizeof(rspf_uint8));
	fs.read((char*)(tmp), sizeof(rspf_uint8));
	*(p+3) = *(rspf_sint16*)(tmp);

	if(fs.fail())
	{
		fs.clear();
		return rspf::nan();
	}

	double p00 = p[0];
	double p01 = p[1];
	double p10 = p[2];
	double p11 = p[3];

	if (p00<-10000.0)  return rspf::nan();
	if (p01<-10000.0)  return rspf::nan();
	if (p10<-10000.0)  return rspf::nan();
	if (p11<-10000.0)  return rspf::nan();

	if (p00>10000.0)  return rspf::nan();
	if (p01>10000.0)  return rspf::nan();
	if (p10>10000.0)  return rspf::nan();
	if (p11>10000.0)  return rspf::nan();

	//if (p00 == info.theNullHeightValue)
	//	w00 = 0.0;
	//if (p01 == info.theNullHeightValue)
	//	w01 = 0.0;
	//if (p10 == info.theNullHeightValue)
	//	w10 = 0.0;
	//if (p11 == info.theNullHeightValue)
	//	w11 = 0.0;

	double sum_weights = w00 + w01 + w10 + w11;

	if (sum_weights)
	{
		return (p00*w00 + p01*w01 + p10*w10 + p11*w11) / sum_weights / theScale;
	}

	return rspf::nan();
}

void test0000()
{
	cout<<"size of float: "<<sizeof(float)<<endl;
	cout<<"size of char: "<<sizeof(char)<<endl;

	//rspfFilename filename("D:\\workspace\\testdata\\Alos\\WW15MGH.DAC");
	//theMGHInputStream = rspfStreamFactoryRegistry::instance()->createNewIFStream(filename,
	//	ios::in | ios::binary);
	ifstream  mghfile;
	mghfile.open("D:\\workspace\\testdata\\Alos\\WW15MGH.DAC", ios_base::binary);
	//mghfile.open("D:\\workspace\\testdata\\Alos\\WW15MGH.DAC", ios_base::in);

	fstream infile;
	fstream outfile;
	infile.open("D:\\workspace\\testdata\\Alos\\inputData.txt", ios_base::in);
	outfile.open("D:\\workspace\\testdata\\Alos\\outputData.txt",ios_base::out);
	double lat, lon;
	int i,j;
	for(i = 0;i < 721;i++)
		for(j = 0;j < 1440;j++)
		{
			lat = (360 - i) / 4.0;
			if(j <= 720)
				lon = j / 4.0;
			else
				lon = (j - 1440) / 4.0;
			rspfGpt gpt(lat, lon, 0.0);
			gpt.hgt = getMGH1(mghfile, gpt);
			outfile<<gpt.lat<<"\t"<<gpt.lon<<"\t"<<gpt.hgt<<endl;
		}
	infile.close();
	outfile.close();
	mghfile.close();
}
void WorldviewRpc()
{
	wxString fileName = "D:\\workspace\\testdata\\Worldview\\052391246010_01_P001_PSH\\Imagery.tif";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	wxString rpcFile = "D:\\workspace\\testdata\\Worldview\\052391246010_01_P001_PSH\\10MAR31062653-S2AS-052391246010_01_P001.RPB";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	wxString gcpfile = "D:\\workspace\\testdata\\Worldview\\052391246010_01_P001_PSH\\gcp.txt";
	rspfFilename outfile = "D:\\workspace\\testdata\\Worldview\\052391246010_01_P001_PSH\\rect.tif";
	rspfFilename reportfile = "D:\\workspace\\testdata\\Worldview\\052391246010_01_P001_PSH\\report.txt";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readWorldviewRPCFile(rpcFile, rpcStruct);

	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;




	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.Orthograph(outfile);


	//rspfRpcSolver solver(true, false);
	/*vector < rspfDpt > imagePoints;
	vector < rspfGpt > groundControlPoints;
	m_MapPar =NULL;
	m_MapPar = PTR_CAST(rspfMapProjection,
	rspfMapProjectionFactory::instance()->createProjection(m_MapProjection));
	if (!m_MapPar) return false;
	int num = static_cast<int>(m_OptCtrlGptSet->getTiePoints().size());
	rspfGpt tGpt;
	for(int i = 0;i < num;i++)
	{

	tGpt = m_MapPar->inverse_do(rspfDpt(m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().lat,
	m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().lon),m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint());
	tGpt.hgt = m_OptCtrlGptSet->getTiePoints()[i]->getGroundPoint().hgt;
	groundControlPoints.push_back(tGpt);
	imagePoints.push_back(m_OptCtrlGptSet->getTiePoints()[i]->getImagePoint());
	}*/
	//solver.solveInitialCoefficients()
	//solver.solveCoefficients(imagePoints,groundControlPoints);

	//m_sensorModel = solver.createRpcModel();
	//m_sensorModel->m_proj = m_MapPar;
	//m_sensorModel->saveState(geom);

	return;
}
void Theos_Physical()
{

	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename gcpfile = "D:\\workspace\\testdata\\gpt_formated.txt";
	rspfFilename sourcefile = "D:\\workspace\\testdata\\Theos\\th_cat_100530021840553_1\\imagery.tif";
	rspfFilename reportfile = "D:\\workspace\\testdata\\Theos\\th_cat_100530021840553_1\\theos_report.txt";

	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(tieLineList);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(6);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);
	prj.InitiateSensorModel(sourcefile);
	rspfSensorModel* sensorModel = prj.m_sensorModel;

	if( NULL == sensorModel )
		return;

	prj.UpdateSensorModel(*prj.m_CtrlGptSet, sensorModel, prj.geom);

	prj.OutputReport(reportfile, sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
}
void Alos_Rpc()
{
	wxString fileName = "H:\\AlosData\\d1001111-001\\IMG-ALPSMW220853075_O1B2R_UW.tif";
	//wxString fileName = "D:\\workspace\\testdata\\beijing\\001\\imagery.tif";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	//rspfFilename elevationpath = "Z:\\share\\alos_data\\beijing\\bj_rspfdem";
	wxString gcpfile = "H:\\AlosData\\d1001111-001\\gcp.txt";
	//wxString gcpfile = "D:\\workspace\\testdata\\beijing\\001\\gcp.txt";
	rspfFilename outfile = "H:\\AlosData\\d1001111-001\\rect.tif";
	//rspfFilename outfile = "D:\\workspace\\testdata\\beijing\\001\\rect.tif";
	rspfFilename reportfile = "H:\\AlosData\\d1001111-001\\report.txt";
	//rspfFilename reportfile = "D:\\workspace\\testdata\\beijing\\001\\report.txt";
	wxString rpcFile = "H:\\AlosData\\d1001111-001\\RPC-ALPSMW220853075_O1B2R_UW.txt";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	//prj.Orthograph(outfile);
}
void Alos_1B2_Pan_Rpc()
{
	wxString imagePath = "H:\\AlosData\\d1001111-005";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	wxDir dir;
	wxArrayString strListTmp;
	dir.GetAllFiles(imagePath, &strListTmp, wxT("IMG-ALPSMW*_O1B2R_UW.tif"));
	if(strListTmp.size() < 1)
		return;
	rspfFilename fileName = rspfFilename(strListTmp[0]);
	strListTmp.clear();

	dir.GetAllFiles(imagePath, &strListTmp, wxT("RPC-ALPSMW*_O1B2R_UW.txt"));
	if(strListTmp.size() < 1)
		return;
	rspfFilename rpcFile = rspfFilename(strListTmp[0]);
	strListTmp.clear();

	wxString gcpfile = imagePath + "\\gcp.txt";
	wxString reportfile = imagePath + "\\report.txt";
	wxString outfile = imagePath + "\\rect.tif";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(rspfFilename(reportfile), prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	//prj.Orthograph(rspfFilename(outfile));
}
void Alos_1B2_MS_Rpc()
{
	wxString imagePath = "H:\\AlosData\\d1001111-006";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	wxString fileName[4];
	wxDir dir;
	wxArrayString strListTmp;
	dir.GetAllFiles(imagePath, &strListTmp, wxT("IMG-01-ALAV2A*_O1B2R_U.tif"));
	if(strListTmp.size() < 1)
		return;
	fileName[0] = strListTmp[0];

	strListTmp.clear();
	dir.GetAllFiles(imagePath, &strListTmp, wxT("IMG-02-ALAV2A*_O1B2R_U.tif"));
	if(strListTmp.size() < 1)
		return;
	fileName[1] = strListTmp[0];

	strListTmp.clear();
	dir.GetAllFiles(imagePath, &strListTmp, wxT("IMG-03-ALAV2A*_O1B2R_U.tif"));
	if(strListTmp.size() < 1)
		return;
	fileName[2] = strListTmp[0];

	strListTmp.clear();
	dir.GetAllFiles(imagePath, &strListTmp, wxT("IMG-04-ALAV2A*_O1B2R_U.tif"));
	if(strListTmp.size() < 1)
		return;
	fileName[3] = strListTmp[0];

	strListTmp.clear();
	dir.GetAllFiles(imagePath, &strListTmp, wxT("RPC-ALAV2A*_O1B2R_U.txt"));
	if(strListTmp.size() < 1)
		return;
	rspfFilename rpcFile = rspfFilename(strListTmp[0]);


	wxString gcpfile = imagePath + "\\gcp.txt";
	wxString reportfile = imagePath + "\\report.txt";
	wxString outfile[4];
	outfile[0] = imagePath + "\\rect01.tif";
	outfile[1] = imagePath + "\\rect02.tif";
	outfile[2] = imagePath + "\\rect03.tif";
	outfile[3] = imagePath + "\\rect04.tif";

	rspfRpcProjection rpcProjection;
	rspfRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);


	rspfTieGptSet tieGptSet;
	rspfKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.m_ChkGptSet = new rspfTieGptSet;
	vector < rspfTieLine > tieLineList;
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));//
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath=rspfFilename(elevationpath);
	prj.ReadGcpAndProjection(rspfFilename(gcpfile));

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	rspfMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = rspfFilename(fileName[0]);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(rspfFilename(fileName[0]));
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(rpcStruct);
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;
	int num = static_cast<int>(tieGptSet.size());

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	cout<<prj.geom;

	prj.OutputReport(rspfFilename(reportfile), prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	for(i = 0;i < 4;i++)
	{
		prj.m_ImgFileNameUnc = rspfFilename(fileName[i]);
		prj.Orthograph(rspfFilename(outfile[i]));
	}
}



bool GetL5Files(wxString sourcePath, wxArrayString &sourceFiles)
{
	if( sourcePath == wxEmptyString )
		return false;

	wxDir dir;
	dir.GetAllFiles(sourcePath, &sourceFiles, wxT("header.dat"), wxDIR_DEFAULT);
	return true;
}

bool GetL7Files(wxString sourcePath, wxArrayString &sourceFiles)
{
	if( sourcePath == wxEmptyString )
		return false;

	wxDir dir;
	dir.GetAllFiles(sourcePath, &sourceFiles, wxT("*_HRF.FST"), wxDIR_DEFAULT);
	return true;
}



bool setProjectionInfo(MyProject& prj, int centerLongitude)
{

	rspfTempFilename pp;
	pp.generateRandomFile();

	ofstream os(pp.c_str(),ios::out|ios_base::app);
	os<<"MAP_PROJECTION_END"<<endl;
	os<<"MAP_PROJECTION_BEGIN"<<endl;
	os<<"central_meridian:  "<<centerLongitude<<endl;
	os<<"datum:  WGE"<<endl;
	os<<"elevation_lookup_flag:  1"<<endl;
	os<<"ellipse_code:  WE"<<endl;
	os<<"ellipse_name:  WGS 84"<<endl;
	os<<"false_easting_northing:  ( 500000.000000000000000, 0.000000000000000 )"<<endl;
	os<<"false_easting_northing_units:  meters"<<endl;
	os<<"major_axis:  6378137.000000000000000"<<endl;
	os<<"minor_axis:  6356752.314199999900000"<<endl;
	os<<"origin_latitude:  0.000000000000000"<<endl;
	os<<"pixel_scale_units:  meters"<<endl;
	os<<"pixel_scale_xy:  ( 30.000000000000000, 30.000000000000000 )"<<endl;
	os<<"scale_factor:  1.000000000000000"<<endl;
	os<<"tie_point_units:  meters"<<endl;
	os<<"tie_point_xy:  ( 237500.000000000000000, 3463380.000000000000000 )"<<endl;
	os<<"type:  rspfTransMercatorProjection"<<endl;
	os<<"MAP_PROJECTION_END"<<endl;

	prj.m_MapProjection.clear();
	prj.m_MapPar=NULL;
	prj.m_MapProjection.addFile(pp.c_str());
	prj.m_MapPar=PTR_CAST(rspfMapProjection,
		rspfMapProjectionFactory::instance()->createProjection(prj.m_MapProjection));
	return true;
}
void PointNorth()
{
	wxString inputPath = "D:\\workspace\\source\\"; 
	wxString outputPath = "D:\\workspace\\result\\"; 

	wxArrayString sourceFiles;
	sourceFiles.clear();
	GetL5Files(inputPath, sourceFiles);
	GetL7Files(inputPath, sourceFiles);
	int nsize = sourceFiles.GetCount();
	for(int i = 0;i < nsize;i++)
	{
		MyProject prj;
		prj.m_ImgFileNameUnc = rspfFilename(sourceFiles.Item( i ));

		rspfImageHandler *handler   = rspfImageHandlerRegistry::instance()->open(prj.m_ImgFileNameUnc);
		if(!handler) return;

		rspfKeywordlist kwlnew,kk;
		rspfDpt imagesize;

		rspfLandSatModel *lmodel = new rspfLandSatModel(prj.m_ImgFileNameUnc);
		setProjectionInfo(prj, lmodel->getCenterLongitude());
		prj.m_sensorModel = lmodel;
		prj.m_sensorModel->saveState(prj.geom);

		prj.m_OutBandList.clear();
		prj.m_OutBandList.push_back(6);
		prj.m_OutBandList.push_back(3);
		prj.m_OutBandList.push_back(1);
		wxString outName = outputPath;
		outName<<i+1;
		outName += ".tif";
		prj.Orthograph(rspfFilename(outName));
	}
}

void searchPathRow()
{
	wxString inputPath = "Z:\\share\\data_share\\landsat\\source"; 
	wxString outputPath = "D:\\workspace\\result";
	wxArrayString sourceFiles;
	sourceFiles.clear(); 

	GetL7Files(inputPath, sourceFiles);
	int nsize = sourceFiles.GetCount();
	ofstream outfile;
	outfile.open("D:\\workspace\\result\\L7searchResult.txt");

	for(int i = 0;i < nsize;i++)
	{
		MyProject prj;
		prj.m_ImgFileNameUnc = rspfFilename(sourceFiles.Item( i ));

		//rspfImageHandler *handler   = rspfImageHandlerRegistry::instance()->open(prj.m_ImgFileNameUnc);
		//if(!handler) return;

		rspfLandSatModel *lmodel = new rspfLandSatModel(prj.m_ImgFileNameUnc);
		int pathNum = lmodel->getPathNumber();
		int rowNum = lmodel->getRowNumber();

		if(pathNum >= 155 && pathNum <= 126 && rowNum >= 44 && rowNum <= 58)
			outfile<<prj.m_ImgFileNameUnc<<"\t"<<pathNum<<"\t"<<rowNum<<endl;
	}
	outfile.close();
}
bool Alos_PRISM_Polynomial()
{
	rspfFilename AlosDir = "H:\\1204\\data\\d1001856-088_ALPSMW223912670_O1B2R_UW_RPC";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";

	rspfFilename outfile = AlosDir + "\\rect000.tif";
	rspfFilename reportfile = AlosDir+"\\report.txt";

	MyProject prj;
	AlosPRISM alosUti(AlosDir);
	if(!alosUti.getInitState())
	{
		cout<<"warning: Alos PRISM数据\""<<AlosDir<<"\"初始化失败！"<<endl;
		return false;
	}
	prj.theMgr = rspfElevManager::instance();
	if(!prj.theMgr->loadElevationPath(rspfFilename(elevationpath)))
	{
		cout<<"warning: 加载DEM失败！"<<endl;
		return false;
	}
	if(!prj.theMgr->openMGH(EGMfile))
	{
		cout<<"warning: 加载大地水准面数据失败！"<<endl;
		return false;
	}
	prj.ReadGcpAndProjection(rspfFilename(AlosDir+"\\gcp.txt"));
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath = rspfFilename(elevationpath);
	prj.GetElevations(prj.m_CtrlGptSet);

	prj.m_MapProjection = alosUti.getKeywordlist();
	prj.m_MapPar = alosUti.getMapProjection();

	prj.m_ImgFileNameUnc = rspfFilename(alosUti.m_FileTIF);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = PolynomialType;
	prj.m_PolyDegree = "1 x y x2 xy y2 x3 y3 xy2 x2y z xz yz";
	if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
	{
		cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
		return false;
	}
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.m_sensorModel->saveState(prj.geom);
	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	cout<<prj.geom;
	prj.Orthograph(outfile);
	return true;
}
bool Alos_PRISM_Rpc()
{
	rspfFilename AlosDir = "H:\\1204\\data\\d1001856-088_ALPSMW223912670_O1B2R_UW_RPC";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";

	rspfFilename outfile = AlosDir + "\\rect000.tif";

	MyProject prj;
	AlosPRISM alosUti(AlosDir);
	if(!alosUti.getInitState())
	{
		cout<<"warning: Alos PRISM数据\""<<AlosDir<<"\"初始化失败！"<<endl;
		return false;
	}
	prj.theMgr = rspfElevManager::instance();
	if(!prj.theMgr->loadElevationPath(rspfFilename(elevationpath)))
	{
		cout<<"warning: 加载DEM失败！"<<endl;
		return false;
	}
	if(!prj.theMgr->openMGH(EGMfile))
	{
		cout<<"warning: 加载大地水准面数据失败！"<<endl;
		return false;
	}
	prj.ReadGcpAndProjection(rspfFilename(AlosDir+"\\gcp.txt"));
	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(elevationpath));
	prj.theMgr->openMGH(EGMfile);
	prj.m_DemPath = rspfFilename(elevationpath);
	prj.GetElevations(prj.m_CtrlGptSet);

	prj.m_MapProjection = alosUti.getKeywordlist();
	prj.m_MapPar = alosUti.getMapProjection();

	prj.m_ImgFileNameUnc = rspfFilename(alosUti.m_FileTIF);
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
	{
		cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
		return false;
	}
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(alosUti.getRpcModelStruct());
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);
	prj.OutputReport(AlosDir+"\\report.txt", prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	cout<<prj.geom;
	//prj.Orthograph(outfile);
	return true;
}
bool Alos_AVNIR2_Rpc()
{
	rspfFilename AlosDir = "H:\\testdata\\blocktest\\alos\\d1001884-010_ALAV2A254402700_O1B2R_U_RPC";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";

	rspfFilename gcpFile = AlosDir + "\\gcp.txt";
	rspfFilename reportFile = AlosDir + "\\report.txt";
	rspfFilename outFile = AlosDir + "\\rect.tif";

	MyProject prj;
	AlosAVNIR2 alosUti(AlosDir);
	prj.ReadGcpAndProjection(gcpFile);
	if(!alosUti.getInitState())
	{
		cout<<"warning: Alos AVNIR2数据\""<<AlosDir<<"\"初始化失败！"<<endl;
		return false;
	}
	prj.theMgr = rspfElevManager::instance();
	if(!prj.theMgr->loadElevationPath(rspfFilename(elevationpath)))
	{
		cout<<"warning: 加载DEM失败！"<<endl;
		return false;
	}
	if(!prj.theMgr->openMGH(EGMfile))
	{
		cout<<"warning: 加载大地水准面数据失败！"<<endl;
		return false;
	}

	prj.m_ImgFileNameUnc = rspfFilename(alosUti.m_FileTIF);
	prj.m_ImgFileNameout = outFile;
	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(2);
	prj.m_OutBandList.push_back(1);
	prj.m_OutBandList.push_back(0);

	prj.m_ModelType = RPCType;
	if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
	{
		cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
		return false;
	}
	rspfRpcModel *rpcModel = new rspfRpcModel;
	rpcModel->setAttributes(alosUti.getRpcModelStruct());
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	//prj.m_CtrlGptSet = new rspfTieGptSet;
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	prj.m_sensorModel->saveState(prj.geom);

	cout<<prj.geom;

	prj.OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	//prj.Orthograph(outFile);
	return true;
}


bool Alos_Batch()
{
	wxString AlosDir = "H:\\Alos数据";
	wxString OutputDir = "H:\\结果";
	rspfFilename elevationpath = "D:\\workspace\\dem";
	rspfFilename EGMfile = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	AlosBatch ab(AlosDir, OutputDir);
	ab.setElevationPath(elevationpath);
	ab.setEgmFile(EGMfile);
	ab.Run();
	return true;
}


void block_alos_landsat()
{
	rspfBlockAdjustment ba;

	//MyProject prj1,prj2,prj3,prj4,prj5,prj6;
	const int num = 9;
	const int numTotal = 9;
	MyProject prj[numTotal];
	string foldName = "alos_landsat";
	rspfFilename sourceFile[numTotal] = {
		"H:\\testdata\\blocktest\\l7\\118029\\LD2010008272\\L71118029_02920020924_HRF.FST",
		"H:\\testdata\\blocktest\\l7\\119029\\LD2010008685\\L71119029_02920020915_HRF.FST",
		"H:\\testdata\\blocktest\\alos\\d1001884-011_ALAV2A254402690_O1B2R_U_RPC",
		"H:\\testdata\\blocktest\\alos\\d1001884-010_ALAV2A254402700_O1B2R_U_RPC",
		"H:\\testdata\\blocktest\\alos\\d1001884-009_ALAV2A254402710_O1B2R_U_RPC",
		"H:\\testdata\\blocktest\\pan\\d1001883-011_alpsmw254402715_o1b2r_uw_rpc",
		"H:\\testdata\\blocktest\\alos\\d1001884-022_ALAV2A249442710_O1B2R_U_RPC",
		"H:\\testdata\\blocktest\\pan\\d1001883-033_alpsmw251922710_o1b2r_uw_rpc",
		"H:\\testdata\\blocktest\\pan\\d1001883-032_alpsmw251922715_o1b2r_uw_rpc"};
	rspfFilename demPath[numTotal] = {			
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem"};
	rspfFilename chkFile[numTotal] = {	
		"H:\\testdata\\blocktest\\l7\\118029\\LD2010008272\\gcp.txt",
		"H:\\testdata\\blocktest\\l7\\119029\\LD2010008685\\gcp.txt",
		"H:\\testdata\\blocktest\\alos\\d1001884-011_ALAV2A254402690_O1B2R_U_RPC\\gcp.txt",
		"H:\\testdata\\blocktest\\alos\\d1001884-010_ALAV2A254402700_O1B2R_U_RPC\\gcp.txt",
		"H:\\testdata\\blocktest\\alos\\d1001884-009_ALAV2A254402710_O1B2R_U_RPC\\gcp.txt",
		"H:\\testdata\\blocktest\\pan\\d1001883-011_alpsmw254402715_o1b2r_uw_rpc\\gcp.txt",
		"H:\\testdata\\blocktest\\alos\\d1001884-022_ALAV2A249442710_O1B2R_U_RPC\\gcp.txt",
		"H:\\testdata\\blocktest\\pan\\d1001883-033_alpsmw251922710_o1b2r_uw_rpc\\gcp.txt",
		"H:\\testdata\\blocktest\\pan\\d1001883-032_alpsmw251922715_o1b2r_uw_rpc\\gcp.txt"};
	rspfFilename outFile[numTotal] = {
		"H:\\testdata\\blocktest\\" + foldName + "\\118029.TIF",
		"H:\\testdata\\blocktest\\" + foldName + "\\119029.TIF",
		"H:\\testdata\\blocktest\\" + foldName + "\\2690.TIF",
		"H:\\testdata\\blocktest\\" + foldName + "\\2700.TIF",
		"H:\\testdata\\blocktest\\" + foldName + "\\2710.TIF",
		"H:\\testdata\\blocktest\\" + foldName + "\\254402715.TIF",
		"H:\\testdata\\blocktest\\" + foldName + "\\249442710.TIF",
		"H:\\testdata\\blocktest\\" + foldName + "\\251922710.TIF",
		"H:\\testdata\\blocktest\\" + foldName + "\\251922715.TIF"};
	/*rspfFilename egmPath[numTotal] = {
		"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC",
		"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC",
		"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC"};*/
		rspfFilename egmPath[numTotal] = {
			"",
			"",
			""};
		ModelType modelType[numTotal] = {
			OrbitType,
			OrbitType,
			RPCType,
			RPCType,
			RPCType,
			RPCType,
			RPCType,
			RPCType,
			RPCType
		};

	int i;
	for (i = 0;i < num;++i)
	{
		prj[i].m_ImgFileNameUnc = sourceFile[i];
		prj[i].m_ImgFileNameout = outFile[i];
		prj[i].m_ModelType = modelType[i];
		prj[i].ReadGcpAndProjection(rspfFilename(chkFile[i]));
		prj[i].theMgr = rspfElevManager::instance();
		prj[i].theMgr->loadElevationPath(rspfFilename(demPath[i]));
		prj[i].theMgr->openMGH(egmPath[i]);
		prj[i].m_DemPath=rspfFilename(demPath[i]);

		prj[i].GetElevations(prj[i].m_CtrlGptSet);
		prj[i].GetElevations(prj[i].m_ChkGptSet);


		prj[i].InitiateSensorModel();
		prj[i].m_sensorModel->m_proj = prj[i].m_MapPar;

		ba.addSensorModel(prj[i].m_sensorModel);
	}

	vector< rspfBlockTieGpt > gptList;
	// 加载（已知）控制点
	AddGpts(rspfFilename("H:\\testdata\\blocktest\\" + foldName + "\\gpt.txt"), gptList);
	GetBlockElevations(prj, gptList);
	//prj[0].GetElevations(rspfFilename(sourceFile[0]), gptList);
	// 加载同名点
	AddCpts(rspfFilename("H:\\testdata\\blocktest\\" + foldName + "\\tpt.txt"), gptList);
	// 进行区域网平差
	ba.adjustment(gptList, rspfBlockAdjustment::RobustMode::NONE);
	for(i = 0;i < num;i++)
	{
		prj[i].m_sensorModel->saveState(prj[i].geom);
	}

	wxString GptFile = "H:\\testdata\\blocktest\\" + foldName + "\\GptList.txt";
	// 保存控制点和连接点
	saveGptList(GptFile, prj, gptList);

	// 用优化后的模型输出检查点精度
	fstream fs;
	fs.open("H:\\testdata\\blocktest\\" + wxString(foldName) + "\\residue.txt",ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(6);
	for (i = 0;i < num;++i)
	{
		AppendResidue2File(prj[i].m_sensorModel, prj[i].m_CtrlGptSet, fs);
		fs<<endl;
	}
	fs.close();

	for(i = 0;i < num;i++)
	{
		prj[i].m_OutBandList.clear();
		//prj[i].m_OutBandList.push_back(6);
		//prj[i].m_OutBandList.push_back(3);
		//prj[i].m_OutBandList.push_back(1);
		prj[i].m_OutBandList.push_back(0);
		if(!wxFile::Exists(outFile[i]))
			prj[i].Orthograph(outFile[i]);
	}
}
void orthl7()
{
	MyProject prj;
	rspfFilename sourceFile = "H:\\testdata\\blocktest\\l7\\118029\\LD2010008272\\L71118029_02920020924_HRF.FST";
	rspfFilename demPath = "D:\\workspace\\dem";
	rspfFilename chkFile = "H:\\testdata\\blocktest\\l7\\118029\\LD2010008272\\gcp.txt";
	rspfFilename outFile = "H:\\testdata\\blocktest\\l7\\118029\\LD2010008272\\rect.tif";
	rspfFilename egmPath = "D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	//rspfFilename egmPath = "";//"D:\\workspace\\testdata\\Alos\\WW15MGH.DAC";
	rspfFilename reportFile = "H:\\testdata\\blocktest\\l7\\118029\\LD2010008272\\report.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(rspfFilename(chkFile));


	prj.theMgr = rspfElevManager::instance();
	prj.theMgr->loadElevationPath(rspfFilename(demPath));
	prj.theMgr->openMGH(egmPath);
	prj.m_DemPath=rspfFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	rspfDpt imagepoint,cimagepoint;
	rspfGpt goundpoint,tGpt;
	int i;

	NEWMAT::ColumnVector m_CtrlResidue = prj.CalcResidue(prj.m_sensorModel, *prj.m_CtrlGptSet);

	vector<rspfRefPtr<rspfTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<rspfRefPtr<rspfTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse_do(rspfDpt(goundpoint.lat, goundpoint.lon),goundpoint);
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	//double *targetVariance;
	prj.m_sensorModel->Huberk = 200.0;
	//prj.m_sensorModel->robustoptimizeFit(*prj.m_CtrlGptSet,targetVariance, rspfSensorModel::RobustMode::HUBER, m_CtrlResidue);
	prj.m_sensorModel->saveState(prj.geom);

	for(i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	{
		rspfDpt dpt = prj.m_sensorModel->m_proj->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
		rspfGpt gpt(dpt.x,dpt.y);
		prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(rspfGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	}

	prj.OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	prj.Orthograph(outFile);
}

void landsat_block()
{
	rspfBlockAdjustment ba;
	wxString foldName = "noGcp";
	wxString strDir = "D:\\workspace\\landsat_block";
	const int num = 8;
	const int numTotal = 8;
	MyProject prj[numTotal];
	rspfFilename sourceFile[numTotal] = {	strDir+"\\123-033\\header.dat",
		strDir+"\\123-034\\header.dat",
		strDir+"\\124-032\\header.dat",
		strDir+"\\124-033\\header.dat",
		strDir+"\\124-034\\header.dat",
		strDir+"\\124-035\\header.dat",
		strDir+"\\125-033\\header.dat",
		strDir+"\\125-034\\header.dat"};
	rspfFilename demPath[numTotal] = {	"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem",
		"D:\\workspace\\dem"};
	rspfFilename chkFile[numTotal] = {	strDir+"\\123-033\\gcp.txt",
		strDir+"\\123-034\\gcp.txt",
		strDir+"\\124-032\\gcp.txt",
		strDir+"\\124-033\\gcp.txt",
		strDir+"\\124-034\\gcp.txt",
		strDir+"\\124-035\\gcp.txt",
		strDir+"\\125-033\\gcp.txt",
		strDir+"\\125-034\\gcp.txt"};
	rspfFilename outFile[numTotal] = {	strDir+"\\"+foldName+"\\123-033.tif",
		strDir+"\\"+foldName+"\\123-034.tif",
		strDir+"\\"+foldName+"\\124-032.tif",
		strDir+"\\"+foldName+"\\124-033.tif",
		strDir+"\\"+foldName+"\\124-034.tif",
		strDir+"\\"+foldName+"\\124-035.tif",
		strDir+"\\"+foldName+"\\125-033.tif",
		strDir+"\\"+foldName+"\\125-034.tif"};
	int i;
	for (i = 0;i < num;++i)
	{
		InitializePrj(prj[i], sourceFile[i], demPath[i], chkFile[i]);
		ba.addSensorModel(prj[i].m_sensorModel);
	}

	vector< rspfBlockTieGpt > gptList;
	// 加载（已知）控制点
	AddGpts(rspfFilename(strDir+"\\"+foldName+"\\gpt.txt"), gptList);
	// 加载同名点
	AddCpts(rspfFilename(strDir+"\\"+foldName+"\\tpt.txt"), gptList);
	// 进行区域网平差
	ba.adjustment(gptList);
	for(i = 0;i < num;i++)
	{
		prj[i].m_sensorModel->saveState(prj[i].geom);
	}

	wxString GptFile = strDir + "\\" + foldName + "\\GptList.txt";
	// 保存控制点和连接点
	saveGptList(GptFile, prj, gptList);

	// 用优化后的模型输出检查点精度
	fstream fs;
	fs.open(strDir+"\\"+foldName+"\\residue.txt",ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(6);
	for (i = 0;i < num;++i)
	{
		AppendResidue2File(prj[i].m_sensorModel, prj[i].m_CtrlGptSet, fs);
		fs<<endl;
	}
	fs.close();

	for(i = 0;i < num;i++)
	{
		prj[i].m_OutBandList.clear();
		//prj[i].m_OutBandList.push_back(6);
		prj[i].m_OutBandList.push_back(3);
		//prj[i].m_OutBandList.push_back(1);
		prj[i].Orthograph(outFile[i]);
	}
}
