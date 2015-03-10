#include "func.h"

#include <ossim_plugin\radi\radiRpcSolver.h>

#ifndef _WIN64
#pragma comment(lib, "ossim20.lib")
#pragma comment(lib, "blas_win32_MT.lib")
#pragma comment(lib, "lapack_win32_MT.lib")
#else
#pragma comment(lib, "ossim20x64.lib")
#pragma comment(lib, "blas_win64_MT.lib")
#pragma comment(lib, "lapack_win64_MT.lib")
#endif

#pragma comment(lib, "ossim_plugin.lib")
#pragma comment(lib, "OpenThreads.lib")
#pragma comment(lib, "gdal_i.lib")
//#pragma comment(lib, "mylib.lib")

namespace mylib{
//int OptimizeGcp(wxString sourcefile, wxString gcpfile, wxString destinationfile, wxString elevationpath)
//{	
//	double metre_threshold = 60.0;
//	bool optimizeFinished = false;
//	int nIteration = 0;
//	int Max_Iteration = 30;
//	int nControl = 20;
//	int nCheck = 0;
//
//	ossimKeywordlist geom;
//	ossimKeywordlist MapProjection;
//
//	ossimTieGptSet *totalCtrlGptSet = new ossimTieGptSet;
//	ossimTieGptSet *totalChkGptSet = new ossimTieGptSet;
//	ossimTieGptSet *ctrlGptSet = new ossimTieGptSet;
//	ossimTieGptSet *chkGptSet = new ossimTieGptSet;
//	ossimTieGptSet *errGptSet = new ossimTieGptSet;
//
//	NEWMAT::ColumnVector ctrlResidue;
//	NEWMAT::ColumnVector chkResidue;
//	NEWMAT::ColumnVector errResidue;
//
//	if( gcpfile == wxEmptyString)
//	{
//		delete totalCtrlGptSet;
//		delete totalChkGptSet;
//		delete ctrlGptSet;
//		delete chkGptSet;
//		delete errGptSet;
//		return ERR_GCP;
//	}
//	MyProject prj;
//	prj.theMgr = ossimElevManager::instance();
//	prj.theMgr->loadElevationPath(wxString2ossimString(elevationpath));//
//	prj.m_DemPath=wxString2ossimString(elevationpath);
//	prj.ReadGcpAndProjection(wxString2ossimString(gcpfile));
//
//	prj.GetElevations(totalCtrlGptSet);
//
//	ossimMapProjection *MapPar = prj.m_MapPar;
//	MapProjection = prj.m_MapProjection;
//
//	prj.InitiateSensorModel(wxString2ossimString(sourcefile));
//	ossimSensorModel* sensorModel = prj.m_sensorModel;
//	if( NULL == sensorModel )
//		return ERR_SOURCE;
//
//	for(; !optimizeFinished; optimizeFinished = prj.CheckSenserModel(sensorModel, totalCtrlGptSet, ctrlGptSet, chkGptSet, errGptSet,
//		metre_threshold, ctrlResidue, chkResidue, errResidue))
//	{
//		if(nIteration++ > Max_Iteration)
//		{
//			delete totalCtrlGptSet;
//			delete ctrlGptSet;
//			delete chkGptSet;
//			delete errGptSet;
//			return ERR_OPTIMIZE;
//		}
//
//		//控制点优化
//		if(!prj.DistributeOptimize(totalCtrlGptSet, ctrlGptSet, chkGptSet, nControl, nCheck))
//		{
//			delete totalCtrlGptSet;
//			delete ctrlGptSet;
//			delete chkGptSet;
//			delete errGptSet;
//			return ERR_DESTINATION;
//		}
//
//		prj.UpdateSensorModel(*ctrlGptSet, sensorModel, geom);
//	}
//	wxString filenametosave = destinationfile.BeforeLast('-') + wxT(".txt");
//	prj.SaveOptPointToFile(wxString2ossimString(filenametosave),ctrlGptSet);
//
//	wxString reportfile = destinationfile.BeforeLast('-');
//	reportfile += wxT("_RMS.txt");
//	if(!prj.OutputReport(wxString2ossimString(reportfile), sensorModel, ctrlGptSet, chkGptSet))
//	{
//		delete totalCtrlGptSet;
//		delete totalChkGptSet;
//		delete ctrlGptSet;
//		delete chkGptSet;
//		delete errGptSet;
//		return ERR_REPORT;
//	}
//
//	delete totalCtrlGptSet;
//	delete totalChkGptSet;
//	delete ctrlGptSet;
//	delete chkGptSet;
//	delete errGptSet;
//	return SUCCESS;
//}
//
//void AutoOptimize(wxArrayString sourceFiles, wxArrayString gcpFiles, wxArrayString destinationFiles, wxString elevationpath)
//{
//	int m_currentIndex = 0;
//	for(m_currentIndex = 0; m_currentIndex < (int)sourceFiles.size(); m_currentIndex++)
//	{
//		wxString destinationPath = destinationFiles[m_currentIndex].substr(0,destinationFiles[m_currentIndex].find_last_of('\\'));
//		int hResult = OptimizeGcp(sourceFiles[m_currentIndex], gcpFiles[m_currentIndex], destinationFiles[m_currentIndex], elevationpath);
//		switch( hResult )
//		{
//		case ERR_SOURCE:
//			{
//				wxString logfilename = destinationPath + wxT("\\errorllog.txt");
//				wxTextFile logfile;
//				logfile.Create(logfilename);
//				if( logfile.Open(logfilename) )
//				{
//					time_t tm;
//					tm = time(NULL);
//					char tmp[64];
//					strftime(tmp,sizeof(tmp),"%Y/%m/%d %X",localtime(&tm));
//					const wxString templogstring = sourceFiles.Item(m_currentIndex) + wxT("......读取待校正图像失败......") + wxString::FromUTF8(tmp);
//					logfile.AddLine( templogstring );
//				}
//				logfile.Write(wxTextFileType_None);
//				logfile.Close();
//				break;
//			}
//		case ERR_GCP:
//			{
//				wxString logfilename = destinationPath + wxT("\\errorllog.txt");
//				wxTextFile logfile;
//				logfile.Create(logfilename);
//				if( logfile.Open(logfilename) )
//				{
//					time_t tm;
//					tm = time(NULL);
//					char tmp[64];
//					strftime(tmp,sizeof(tmp),"%Y/%m/%d %X",localtime(&tm));
//					const wxString templogstring = sourceFiles.Item(m_currentIndex) + wxT("......读取控制点文件失败......") + wxString::FromUTF8(tmp);
//					logfile.AddLine( templogstring );
//				}
//				logfile.Write(wxTextFileType_None);
//				logfile.Close();
//				break;
//			}
//		case ERR_OPTIMIZE:
//			{
//				wxString logfilename = destinationPath + wxT("\\errorllog.txt");
//				wxTextFile logfile;
//				logfile.Create(logfilename);
//				if( logfile.Open(logfilename) )
//				{
//					time_t tm;
//					tm = time(NULL);
//					char tmp[64];
//					strftime(tmp,sizeof(tmp),"%Y/%m/%d %X",localtime(&tm));
//					const wxString templogstring = gcpFiles.Item(m_currentIndex) + wxT("......控制点优化失败......") + wxString::FromUTF8(tmp);
//					logfile.AddLine( templogstring );
//				}
//				logfile.Write(wxTextFileType_None);
//				logfile.Close();
//				break;
//			}
//		case ERR_DESTINATION:
//			{
//				wxString logfilename = destinationPath + wxT("\\errorllog.txt");
//				wxTextFile logfile;
//				logfile.Create(logfilename);
//				if( logfile.Open(logfilename) )
//				{
//					time_t tm;
//					tm = time(NULL);
//					char tmp[64];
//					strftime(tmp,sizeof(tmp),"%Y/%m/%d %X",localtime(&tm));
//					const wxString templogstring = destinationFiles.Item(m_currentIndex) + wxT("......校正图像失败......") + wxString::FromUTF8(tmp);
//					logfile.AddLine( templogstring );
//				}
//				logfile.Write(wxTextFileType_None);
//				logfile.Close();
//				break;
//			}
//		case ERR_REPORT:
//			{
//				wxString logfilename = destinationPath + wxT("\\errorllog.txt");
//				wxTextFile logfile;
//				logfile.Create(logfilename);
//				if( logfile.Open(logfilename) )
//				{
//					time_t tm;
//					tm = time(NULL);
//					char tmp[64];
//					strftime(tmp,sizeof(tmp),"%Y/%m/%d %X",localtime(&tm));
//					const wxString templogstring = sourceFiles.Item(m_currentIndex) + wxT("......生成误差报告失败......") + wxString::FromUTF8(tmp);
//					logfile.AddLine( templogstring );
//				}
//				logfile.Write(wxTextFileType_None);
//				logfile.Close();
//				break;
//			}
//		case SUCCESS:
//			{
//				wxString logfilename = destinationPath + wxT("\\successfullog.txt");
//				wxTextFile logfile;
//				logfile.Create(logfilename);
//				if( logfile.Open(logfilename) )
//				{
//					time_t tm;
//					tm = time(NULL);
//					char tmp[64];
//					strftime(tmp,sizeof(tmp),"%Y/%m/%d %X",localtime(&tm));
//					const wxString templogstring = sourceFiles.Item(m_currentIndex) + wxT("......控制点优化完成......") + wxString::FromUTF8(tmp);
//					logfile.AddLine( templogstring );
//				}
//				logfile.Write(wxTextFileType_None);
//				logfile.Close();
//				break;
//			}
//		}
//	}
//}


ossimTieGptSet *ImportGptFile(ossimFilename gcpFile)
{
	ossimTieGptSet *gptSet = new ossimTieGptSet;
	ossimTieGpt*		aTiePt;
	ossimElevManager*	theMgr;
	theMgr = ossimElevManager::instance();

	//读取控制点文件，并装载入gptSet
	fstream ifGct;
	ifGct.open(gcpFile.c_str(), ios_base::in);
	while(1)
	{
		int ptName;
		ossimGpt gptCtr;
		ossimDpt dptCtr;
		//double wgtCtr;

		if(!(ifGct>>ptName>>dptCtr.x>>dptCtr.y>>gptCtr.lat>>gptCtr.lon>>gptCtr.hgt))
			break;

		//ossimGpt tGpt;
		//tGpt = m_transMerc->inverse(ossimDpt(gptCtr.lat,gptCtr.lon));
		//tGpt.hgt = theMgr->getHeightAboveEllipsoid(tGpt);
		aTiePt=new ossimTieGpt(gptCtr,dptCtr,0);
		aTiePt->GcpNumberID = ossimString::toString(ptName);
		aTiePt->score = 1.0;

		gptSet->addTiePoint(aTiePt);
	}
	ifGct.close();

	return gptSet;
}

void AddGpts(ossimFilename gptFile, vector<ossimBlockTieGpt> &gptList)
{
	//读取控制点文件，并装载入gptList
	fstream ifGct;
	ifGct.open(gptFile.c_str(), ios_base::in);
	while(1)
	{
		ossimString ptName;
		ossimGpt gptCtr;
		ossimDpt dptCtr;
		int index;
		//double wgtCtr;
		ossimTieGpt*		aTiePt;

		if(!(ifGct>>ptName>>index>>dptCtr.x>>dptCtr.y>>gptCtr.lat>>gptCtr.lon>>gptCtr.hgt))
			break;

		aTiePt=new ossimTieGpt(gptCtr,dptCtr,0);
		aTiePt->GcpNumberID = ptName;

		ossimBlockTieGpt blockTieGpt(*aTiePt);
		blockTieGpt.m_ID = ptName;
		blockTieGpt.m_ImageIndices.push_back(index);
		blockTieGpt.m_DptList.push_back(dptCtr);
		blockTieGpt.m_nType = 1;

		gptList.push_back(blockTieGpt);
	}
	ifGct.close();
}
void AddCpts(ossimFilename cptFile, vector<ossimBlockTieGpt> &gptList)
{
	//读取控制点文件，并装载入gptList
	fstream ifGct;
	ifGct.open(cptFile.c_str(), ios_base::in);
	while(1)
	{
		ossimString ptName;
		ossimGpt gptCtr;
		int index1;
		int index2;
		ossimDpt dpt1;
		ossimDpt dpt2;
		//double wgtCtr;
		ossimBlockTieGpt blockTieGpt;

		if(!(ifGct>>ptName>>index1>>dpt1.x>>dpt1.y>>index2>>dpt2.x>>dpt2.y))
			break;

		blockTieGpt.m_ID = ptName;
		blockTieGpt.m_ImageIndices.push_back(index1);
		blockTieGpt.m_DptList.push_back(dpt1);
		blockTieGpt.m_ImageIndices.push_back(index2);
		blockTieGpt.m_DptList.push_back(dpt2);
		blockTieGpt.m_nType = 0;

		gptList.push_back(blockTieGpt);
	}
	ifGct.close();
}
void ToLatLon(ossimSensorModel* sensorModel, ossimTieGptSet* gptSet)
{
	ossimDpt imagepoint,cimagepoint;
	ossimGpt goundpoint,tGpt;
	int num = gptSet->size();
	for(int i = 0;i < num;++i)
	{
		imagepoint = gptSet->getTiePoints()[i]->getImagePoint();
		goundpoint = gptSet->getTiePoints()[i]->getGroundPoint();
		tGpt = sensorModel->m_proj->inverse(ossimDpt(goundpoint.lat,goundpoint.lon));
		tGpt.hgt = gptSet->getTiePoints()[i]->hgt;
		gptSet->refTiePoints()[i]->setGroundPoint(tGpt);
	}
}
void ToLatLon(ossimMapProjection* proj, ossimTieGptSet* gptSet)
{
	ossimDpt imagepoint, cimagepoint;
	ossimGpt goundpoint, tGpt;
	int num = gptSet->size();
	for (int i = 0; i < num; ++i)
	{
		imagepoint = gptSet->getTiePoints()[i]->getImagePoint();
		goundpoint = gptSet->getTiePoints()[i]->getGroundPoint();
		tGpt = proj->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
		tGpt.hgt = gptSet->getTiePoints()[i]->hgt;
		gptSet->refTiePoints()[i]->setGroundPoint(tGpt);
	}
}

ossimTieGptSet* ReadGptFromFile(ossimFilename gcpFile)
{
	ossimTieGptSet *gptSet = new ossimTieGptSet;
	ossimTieGpt*		aTiePt;

	//读取控制点文件，并装载入gptSet
	fstream ifGct;
	ifGct.open(gcpFile.c_str(), ios_base::in);
	while(1)
	{
		char ptName[256];
		ossimGpt gptCtr;
		ossimDpt dptCtr;
		//double wgtCtr;

		if(!(ifGct>>ptName>>dptCtr.x>>dptCtr.y>>gptCtr.lat>>gptCtr.lon>>gptCtr.hgt))
			break;

		aTiePt=new ossimTieGpt(gptCtr,dptCtr,0);
		aTiePt->GcpNumberID = ossimString(ptName);
		aTiePt->score = 1.0;

		gptSet->addTiePoint(aTiePt);
	}
	ifGct.close();

	return gptSet;
}

void AppendResidue2File(ossimSensorModel* sensorModel, ossimTieGptSet* GptSet, fstream &fs)
{
	//ToLatLon(sensorModel, GptSet);
	//residue = sensorModel->getResidue(*GptSet, true);
	int num = GptSet->size();
	ossimDpt imagepoint;
	ossimGpt goundpoint,tGpt;
	ossimDpt residue1,residue2;
	int i;
	NEWMAT::ColumnVector residue;
	residue.ReSize(num*3);
	for(i = 0;i < num;++i)
	{
		imagepoint = GptSet->getTiePoints()[i]->getImagePoint();
		tGpt.hgt = GptSet->getTiePoints()[i]->hgt;
		sensorModel->lineSampleToWorld(imagepoint,tGpt);
		residue1 = sensorModel->m_proj->forward(tGpt);
		ossimGpt gpt = GptSet->getTiePoints()[i]->getGroundPoint();
		residue2 = ossimDpt(gpt.lat,gpt.lon);

		//ossimDpt imDerp;

		//ossimDpt resIm;
		//ossimTieGpt tiePoint = *GptSet->getTiePoints()[i];
		//resIm = tiePoint.tie - sensorModel->forward(tiePoint);

		fs<<GptSet->getTiePoints()[i]->GcpNumberID<<"\t"
			<<GptSet->getTiePoints()[i]->getImagePoint().x<<"\t"
			<<GptSet->getTiePoints()[i]->getImagePoint().y<<"\t"
			//<<GptSet->getTiePoints()[i]->getGroundPoint().lat<<"\t"
			//<<GptSet->getTiePoints()[i]->getGroundPoint().lon<<"\t"
			<<residue1.lon<<"\t"
			<<residue1.lat<<"\t"
			<<residue1.x - residue2.x<<"\t"
			<<residue1.y - residue2.y<<endl;

		//residue.element(i * 3 + 0) = residue2.x - residue1.x;
		//residue.element(i * 3 + 1) = residue2.y - residue1.y;
		//residue.element(i * 3 + 2) = 0;
	}

	/*for(i = 0;i < static_cast<int>(GptSet->size());++i)
	{
		fs<<GptSet->getTiePoints()[i]->GcpNumberID<<"\t"
			<<GptSet->getTiePoints()[i]->getImagePoint().x<<"\t"
			<<GptSet->getTiePoints()[i]->getImagePoint().y<<"\t"
			<<residue.element(3 * i + 0)<<"\t"
			<<residue.element(3 * i + 1)<<endl;
	}*/
}
void InitializePrj(MyProject& prj, ossimFilename sourceFile, ossimFilename demPath, ossimFilename projectionFile)
{
	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_DemPath = demPath;
	prj.ReadGcpAndProjection(projectionFile, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	*prj.m_OptCtrlGptSet = *prj.m_CtrlGptSet;
	*prj.m_OptChkGptSet = *prj.m_ChkGptSet;
	prj.InitiateSensorModel();
}

bool GetBlockElevations(MyProject prj[], vector< ossimBlockTieGpt >& gptList)
{
	ossimGpt tGpt;
	for(unsigned int i = 0;i < gptList.size();i++)
	{
		if(0 == gptList[i].m_nType)
			continue;
		tGpt = prj[gptList[i].m_ImageIndices[0]].m_MapPar->inverse(ossimDpt(gptList[i].getGroundPoint().lat, gptList[i].getGroundPoint().lon));
		ossimGpt gpt = gptList[i].getGroundPoint();
		gpt.hgt = prj[gptList[i].m_ImageIndices[0]].theMgr->getHeightAboveEllipsoid(tGpt);
		gptList[i].setGroundPoint(gpt);
	}

	return true;
}

bool saveGptList(const char* GptFile, MyProject prj[], vector< ossimBlockTieGpt > gptList)
{
	fstream fs;
	fs.open(GptFile, ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(2);

	int nsize = gptList.size();

	for(int i = 0;i < nsize;i++)
	{
		int imageIndex = gptList[i].m_ImageIndices[0];
		ossimDpt image_point = gptList[i].m_DptList[0];
		ossimGpt lonlat;
		ossimDpt world_point;
		prj[imageIndex].m_sensorModel->lineSampleToWorld(image_point,lonlat);
		world_point = prj[imageIndex].m_sensorModel->m_proj->forward(lonlat);
		
		fs<<i+1<<"\t"<<gptList[i].m_nType<<"\t"<<world_point.lon<<"\t"<<world_point.lat<<endl;
	}
	fs.close();
	return true;
}


bool AppendTieAreas(vector<ossimTieFeature>& tieFeatureList, const vector<ossimDArea> &dptAreaList, const vector<ossimGArea> &gptAreaList)
{
	if(dptAreaList.size() != gptAreaList.size()) return false;
	for(unsigned int i = 0;i < dptAreaList.size();i++)
	{
		ossimTieFeature tieFeature;
		tieFeature.setImageFeature(dptAreaList[i]);
		tieFeature.setGroundFeature(gptAreaList[i]);
		tieFeatureList.push_back(tieFeature);
	}
	return true;
}

bool AppendTieLines(vector<ossimTieFeature>& tieFeatureList, const vector<ossimDLine> &dptLineList, const vector<ossimGLine> &gptLineList)
{
	if(dptLineList.size() != gptLineList.size()) return false;
	for(unsigned int i = 0;i < dptLineList.size();i++)
	{
		ossimTieFeature tieFeature;
		tieFeature.setImageFeature(dptLineList[i]);
		tieFeature.setGroundFeature(gptLineList[i]);
		tieFeatureList.push_back(tieFeature);
	}
	return true;
}

bool AppendTieFeatures(vector<ossimTieFeature>& tieFeatureList, const vector<ossimDFeature> &dptFeatureList, const vector<ossimGFeature> &gptFeatureList)
{
	for(unsigned int i = 0;i < dptFeatureList.size();i++)
	{
		if("NULL" == dptFeatureList[i].strId) continue;
		for(unsigned int j = 0;j <gptFeatureList.size();j++)
		{
			if(dptFeatureList[i].strId == gptFeatureList[j].strId)
			{
				ossimTieFeature tieFeature;
				tieFeature.setImageFeature(dptFeatureList[i]);
				tieFeature.setGroundFeature(gptFeatureList[j]);
				tieFeature.setId(dptFeatureList[i].strId);
				tieFeatureList.push_back(tieFeature);
				break;
				//gptFeatureList.erase(gptFeatureList.begin() + j);
			}
		}
	}
	return true;
}

bool AppendTieFeaturesFromTieGptSet(vector<ossimTieFeature>& tieFeatureList, ossimTieGptSet* pGptSet)
{
	int num = (int)pGptSet->getTiePoints().size();
	for (int i = 0; i < num; i++)
	{
		ossimDFeature dFeature;
		ossimGFeature gFeature;
		dFeature.m_featureType = ossimFeatureType::ossimPointType;
		dFeature.m_Points.push_back(pGptSet->getTiePoints()[i]->getImagePoint());
		gFeature.m_featureType = ossimFeatureType::ossimPointType;
		gFeature.m_Points.push_back(pGptSet->getTiePoints()[i]->getGroundPoint());
		
		ossimTieFeature tieFeature;
		tieFeature.setImageFeature(dFeature);
		tieFeature.setGroundFeature(gFeature);
		char buf[1024];
		int id = pGptSet->getTiePoints()[i]->GcpNumberID.toInt();
		sprintf_s(buf, "1%03d", id);
		tieFeature.setId(buf);
		tieFeatureList.push_back(tieFeature);
	}
	return true;
}

template<class T>
bool ReadDFeatures(string strFilename, vector <T> & featureList)
{
	char strtmp[1024];
	std::vector<string> strList;
	string str;

	fstream os;
	os.open(strFilename.c_str(), ios_base::out);

	while(os.getline(strtmp,1024))
	{
		str = strtmp;
		if(str.empty()) continue;
		strList.clear();
		SplitString(str, "	", strList, false);
		if (strList.size()==1) {
			str = strList[0];
			strList.clear();
			SplitString(str, " ", strList, false);
		}
		if (strList.size() != 3) continue;

		T tmpFeature;
		tmpFeature.strId = strList[0];
		int num = atoi(strList[2].c_str());

		if(0 == strList[1].compare("Point"))
		{
			tmpFeature.m_featureType = ossimFeatureType::ossimPointType;
		}
		else if(0 == strList[1].compare("FreeLine"))
		{
			tmpFeature.m_featureType = ossimFeatureType::ossimFreeLineType;
		}
		else if(0 == strList[1].compare("Polygon"))
		{
			tmpFeature.m_featureType = ossimFeatureType::ossimPolygonType;
		}
		else continue;

		for(int i = 0;i < num;i++)
		{
			if(!os.getline(strtmp, 1024)) return false;
			string strPoint = strtmp;
			if(strPoint.empty()) continue;
			vector<string> strCorList;
			SplitString(strPoint, "	", strCorList, false);
			if(strCorList.size()==1){
				strPoint = strCorList[0];
				strCorList.clear();
				SplitString(strPoint, "	", strCorList, false);
			}
			if(strCorList.size() != 2) continue;
			ossimPoint pt(atof(strCorList[0].c_str()), atof(strCorList[1].c_str()));
			tmpFeature.m_Points.push_back(pt);
		}
		featureList.push_back(tmpFeature);
	}
	return true;
}

ossimProjection* newUtmView(const ossimGpt& centerGround,
							const ossimDpt& metersPerPixel)
{
	ossimUtmProjection* utm = new ossimUtmProjection;

	// we will make it a square pixel in meters
	double averageGsd = (metersPerPixel.x + metersPerPixel.y)*.5;
	utm->setZone(centerGround);
	utm->setMetersPerPixel(ossimDpt(metersPerPixel));


	return utm;
}

ossimRpcModel *createRpcModelFromPoints(ossimTieGptSet* gptSet, const ossimDpt& imageShift)
{
	ossimRpcSolver *rpcSolver = new ossimRpcSolver();

	int nsize = gptSet->size();
	std::vector<ossimDpt> imagePoints;
	std::vector<ossimGpt> groundPoints;
	imagePoints.resize(nsize);
	groundPoints.resize(nsize);
	for(int i = 0;i < nsize;i++)
	{
		imagePoints[i] = gptSet->getTiePoints()[i]->getImagePoint();
		groundPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
	}
	rpcSolver->solveCoefficients(imagePoints, groundPoints, imageShift);

	ossimImageGeometry *imageGeom = rpcSolver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimRpcModel *rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom);
	return rpcModel;
}


ossimRpcModel *createRpcModelFromProjection(const ossimDrect& imageBounds,
							 ossimProjection* imageProj,
							 ossim_uint32 xSamples,
							 ossim_uint32 ySamples,
							 bool shiftTo0Flag)
{
	ossimRpcSolver *rpcSolver = new ossimRpcSolver();
	rpcSolver->solveCoefficients(imageBounds,
							imageProj,
							xSamples,
							ySamples,
							shiftTo0Flag);
	ossimImageGeometry *imageGeom = rpcSolver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimRpcModel *rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom);
	return rpcModel;
}

ossimTieGptSet* createTieGptSet(const ossimDrect& imageBounds,
									  const ossimSensorModel& proj,
									  ossim_uint32 xSamples,
									  ossim_uint32 ySamples,
									  bool latlon,
									  bool shiftTo0Flag)
{
	ossimTieGptSet* pTieGptSet = new ossimTieGptSet;
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

			ossimDpt dpt(w*xnorm + ul.x,
				h*ynorm + ul.y);

			proj.lineSampleToWorld(dpt,
				gpt);
			gpt.changeDatum(defaultGround.datum());
			if(shiftTo0Flag)
			{
				imagePoint = dpt + shiftTo0;
			}
			else
			{
				imagePoint = dpt;
			}
			double h = ossimElevManager::instance()->getHeightAboveMSL(gpt);
			if(ossim::isnan(h) == false)
			{
				gpt.height(h);
			}
			if(gpt.isHgtNan())
			{
				gpt.height(0.0);
			}
			if (latlon)
			{
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
	return pTieGptSet;
}


void createTieGptSet(const ossimDrect& imageBounds,
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

			proj.lineSampleHeightToWorld(dpt, height, gpt);
			//ossimGpt gptTemp;
			//proj.lineSampleToWorld(dpt,
			//	gptTemp);
			//ossimDpt dptTemp;
			//proj.worldToLineSample(gptTemp, dptTemp);
			//ossimDpt res = dptTemp - dpt;
			gpt.changeDatum(defaultGround.datum());
			if(shiftTo0Flag)
			{
				imagePoint = dpt + shiftTo0;
			}
			else
			{
				imagePoint = dpt;
			}
			//gpt.hgt = height;

			if (latlon)
			{
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

			ossimDpt dpt((int)((w - 1)*xnorm + ul.x),
				int((h-1)*ynorm + ul.y));

			proj.lineSampleHeightToWorld(dpt, height, gpt);
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

			if (latlon)
			{
				dpt = ossimDpt(gpt.lat, gpt.lon);
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

//void create3DGridPoints(const ossimDrect& imageBounds,
//						const ossimSensorModel& proj,
//						double height,
//						ossimTieGptSet*& pTieGptSet,
//						ossim_uint32 xSamples,
//						ossim_uint32 ySamples,
//						bool latlon,
//						bool shiftTo0Flag)
//{
//	if (NULL == pTieGptSet)
//	{
//		pTieGptSet = new ossimTieGptSet;
//	}
//	ossim_uint32 x,y;
//	ossim_float64 w = imageBounds.width();
//	ossim_float64 h = imageBounds.height();
//	ossimGpt gpt;
//	ossimGpt defaultGround;
//	if(ySamples < 1) ySamples = 12;
//	if(xSamples < 1) xSamples = 12;
//	srand(time(0));
//	double xnorm;
//	double ynorm;
//	ossimDpt ul = imageBounds.ul();
//	ossimDpt shiftTo0(-ul.x,
//		-ul.y);
//	ossimGpt gpt_ul;
//	ossimGpt gpt_lr;
//	proj.lineSampleToWorld(ossimDpt(0, 0), gpt_ul);
//	proj.lineSampleToWorld(ossimDpt(w-1, h-1), gpt_lr);
//	double xStep = (gpt_lr.lat-gpt_ul.lat)/(xSamples-1);
//	double yStep = (gpt_lr.lon-gpt_ul.lon)/(ySamples-1);
//	for(y = 0; y < ySamples; ++y)
//	{
//		for(x = 0; x < xSamples; ++x)
//		{
//			ossimDpt imagePoint;
//			if(ySamples > 1)
//			{
//				ynorm = (double)y/(double)(ySamples - 1);
//			}
//			else
//			{
//				ynorm = 0.0;
//			}
//			if(xSamples > 1)
//			{
//				xnorm = (double)x/(double)(xSamples - 1);
//			}
//			else
//			{
//				xnorm = 0.0;
//			}
//
//
//			ossimGpt gpt(gpt_ul.lat+x*xStep, gpt_ul.lon+y*yStep, height);
//			ossimDpt dpt = proj.forward(gpt);
//			//proj.worldToLineSample(gpt, dpt);
//
//			//proj.lineSampleHeightToWorld(dpt, height, gpt);
//			//ossimGpt gptTemp;
//			//proj.lineSampleToWorld(dpt,
//			//	gptTemp);
//			//ossimDpt dptTemp;
//			//proj.worldToLineSample(gptTemp, dptTemp);
//			//ossimDpt res = dptTemp - dpt;
//			gpt.changeDatum(defaultGround.datum());
//			if(shiftTo0Flag)
//			{
//				imagePoint = dpt + shiftTo0;
//			}
//			else
//			{
//				imagePoint = dpt;
//			}
//			//gpt.hgt = height;
//
//			if (latlon)
//			{
//				dpt = ossimDpt(gpt.lon, gpt.lat);
//			}
//			else
//			{
//				dpt = proj.m_proj->forward(gpt);
//			}
//			ossimGpt groundPoint(dpt.x,dpt.y, gpt.hgt);
//
//			ossimString strId;
//			char tmpStr[256];
//			sprintf(tmpStr, "%d", y * xSamples + x + 1);
//			strId = tmpStr;
//			ossimTieGpt *aTiePt = new ossimTieGpt(groundPoint, imagePoint, 1.0, strId);
//			pTieGptSet->addTiePoint(aTiePt);
//		}
//	}
//}

//bool ReadFeatures(string strFilename, vector < ossimDFeature > & featureList, ossimDFeature::ossimDFeatureType featureType/* = ossimFeature::ossimUnknown*/)
//{
//	char strtmp[1024];
//	std::vector<string> strList;
//	string str;
//
//	ifstream os(strFilename.c_str());
//
//	vector<string> delimiterList;
//	delimiterList.push_back(" ");
//	delimiterList.push_back("\t");
//	while(os.getline(strtmp,1024))
//	{
//		str = strtmp;
//		if(str.empty()) continue;
//		strList.clear();
//		splitString(str, delimiterList, strList);
//		//SplitString(str, "	", strList, false);
//		//if (strList.size()==1) {
//		//	str = strList[0];
//		//	strList.clear();
//		//	SplitString(str, " ", strList, false);
//		//}
//		if (strList.size() != 3) continue;
//
//		ossimDFeature tmpFeature;
//		tmpFeature.strId = strList[0];
//		int num = atoi(strList[2].c_str());
//
//		if(0 == strList[1].compare("point"))
//		{
//			if(ossimDFeature::ossimDUnknown != featureType && ossimDFeature::ossimDPointType != featureType ) continue;
//			tmpFeature.m_featureType = ossimDFeature::ossimDPointType;
//		}
//		else if(0 == strList[1].compare("line"))
//		{
//			if(ossimDFeature::ossimDUnknown != featureType && ossimDFeature::ossimDLineType != featureType ) continue;
//			tmpFeature.m_featureType = ossimDFeature::ossimDLineType;
//		}
//		else if(0 == strList[1].compare("area"))
//		{
//			if(ossimDFeature::ossimDUnknown != featureType && ossimDFeature::ossimDAreaType != featureType ) continue;
//			tmpFeature.m_featureType = ossimDFeature::ossimDAreaType;
//		}
//		else continue;
//
//		for(int i = 0;i < num;i++)
//		{
//			if(!os.getline(strtmp, 1024)) return false;
//			string strPoint = strtmp;
//			if(strPoint.empty()) continue;
//			vector<string> strCorList;
//			splitString(strPoint, delimiterList, strCorList);
//			if(strCorList.size() != 2) continue;
//			ossimDpt pt(atof(strCorList[0].c_str()), atof(strCorList[1].c_str()));
//			tmpFeature.m_Points.push_back(pt);
//		}
//		featureList.push_back(tmpFeature);
//	}
//	return true;
//}
//
//bool ReadFeatures(string strFilename, vector < ossimGFeature > & featureList, ossimGFeature::ossimGFeatureType featureType/* = ossimGFeature::ossimGUnknown*/)
//{
//	char strtmp[1024];
//	std::vector<string> strList;
//	string str;
//
//	ifstream os(strFilename.c_str());
//	vector<string> delimiterList;
//	delimiterList.push_back(" ");
//	delimiterList.push_back("\t");
//
//	while(os.getline(strtmp,1024))
//	{
//		str = strtmp;
//		if(str.empty()) continue;
//		strList.clear();
//		splitString(str, delimiterList, strList);
//		if (strList.size() != 3) continue;
//
//		ossimGFeature tmpFeature;
//		tmpFeature.strId = strList[0];
//		int num = atoi(strList[2].c_str());
//
//		if(0 == strList[1].compare("point"))
//		{
//			if(ossimGFeature::ossimGUnknown != featureType && ossimGFeature::ossimGPointType != featureType ) continue;
//			tmpFeature.m_featureType = ossimGFeature::ossimGPointType;
//		}
//		else if(0 == strList[1].compare("line"))
//		{
//			if(ossimGFeature::ossimGUnknown != featureType && ossimGFeature::ossimGLineType != featureType ) continue;
//			tmpFeature.m_featureType = ossimGFeature::ossimGLineType;
//		}
//		else if(0 == strList[1].compare("area"))
//		{
//			if(ossimGFeature::ossimGUnknown != featureType && ossimGFeature::ossimGAreaType != featureType ) continue;
//			tmpFeature.m_featureType = ossimGFeature::ossimGAreaType;
//		}
//		else continue;
//
//		for(int i = 0;i < num;i++)
//		{
//			if(!os.getline(strtmp, 1024)) return false;
//			string strPoint = strtmp;
//			if(strPoint.empty()) continue;
//			vector<string> strCorList;
//			splitString(strPoint, delimiterList, strCorList);
//			if(strCorList.size()==1){
//				strPoint = strCorList[0];
//				strCorList.clear();
//				splitString(strPoint, delimiterList, strCorList);
//			}
//			if(strCorList.size() == 3)
//			{
//				ossimGpt pt(atof(strCorList[0].c_str()), atof(strCorList[1].c_str()), atof(strCorList[2].c_str()));
//				tmpFeature.m_Points.push_back(pt);
//			}
//			else if(strCorList.size() == 2)
//			{
//				ossimGpt pt(atof(strCorList[0].c_str()), atof(strCorList[1].c_str()), 0.0f);
//				tmpFeature.m_Points.push_back(pt);
//			}
//		}
//		featureList.push_back(tmpFeature);
//	}
//	return true;
//}

bool FeatureReport(ossimFilename reportfile, ossimSensorModel* sensorModel, vector<ossimTieFeature> tieFeatureList)
{
	int i;
	fstream fs;
	fs.open(reportfile.c_str(),ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(2);

	fs<<"残差报告\n\n"<<"残差单位：像素\n";

	NEWMAT::ColumnVector residue = CalcResidue(sensorModel, tieFeatureList);

	fs<<"清单："<<endl;

	fs<<"标示\t"<<"误差\t"<<"残差X\t"<<"残差Y\t"<<"类型\t"<<endl;

	fs<<"用于计算模型的控制特征数："<<tieFeatureList.size()<<endl;

	//residue = CalcResidue(*m_CtrTieGptSet);
	ossimString strType[10];
	strType[ossimTieFeature::ossimTiePointPoint] = "pp";
	strType[ossimTieFeature::ossimTieLineLine] = "ll";
	strType[ossimTieFeature::ossimTieAreaArea] = "aa";
	strType[ossimTieFeature::ossimTiePointLine] = "pl";
	strType[ossimTieFeature::ossimTiePointArea] = "pa";
	strType[ossimTieFeature::ossimTieLinePoint] = "lp";
	strType[ossimTieFeature::ossimTieLineArea] = "la";
	strType[ossimTieFeature::ossimTieAreaPoint] = "ap";
	strType[ossimTieFeature::ossimTieAreaLine] = "al";
	strType[ossimTieFeature::ossimTieUnknown] = "uk";

	if(residue.Nrows() < 2)
	{
		fs.close();
		return true;
	}

	double minX = fabs(residue.element(0));
	double maxX = fabs(residue.element(0));
	double minY = fabs(residue.element(1));
	double maxY = fabs(residue.element(1));
	double minXY = sqrt(minX * minX + minY * minY);
	double maxXY = sqrt(maxX * maxX + maxY * maxY);
	double sumX_2 = 0.0;
	double sumY_2 = 0.0;
	for(i = 0;i < static_cast<int>(tieFeatureList.size());i++)
	{
		fs<<tieFeatureList[i].getId()<<"\t"
			<<sqrt(residue.element(i * 2 + 0)*residue.element(i * 2 + 0) + residue.element(i * 2 + 1)*residue.element(i * 2 + 1))<<"\t"
			<<residue.element(i * 2 + 0)<<"\t"
			<<residue.element(i * 2 + 1)<<"\t"
			<<strType[tieFeatureList[i].getTieFeatureType()]<<"\n";
		double x = residue.element(i * 2 + 0);
		double y = residue.element(i * 2 + 1);
		double xy = sqrt(x * x + y * y);
		minX = (minX > fabs(x)) ? fabs(x) : minX;
		minY = (minY > fabs(y)) ? fabs(y) : minY;
		minXY = (minXY > xy) ? xy : minXY;
		maxX = (maxX < fabs(x)) ? fabs(x) : maxX;
		maxY = (maxY < fabs(y)) ? fabs(y) : maxY;
		maxXY = (maxXY < xy) ? xy : maxXY;
		sumX_2 += x * x;
		sumY_2 += y * y;
	}

	double rmsX = sqrt(sumX_2 / tieFeatureList.size());
	double rmsY = sqrt(sumY_2 / tieFeatureList.size());
	double rmsXY = sqrt(rmsX * rmsX + rmsY * rmsY);
	fs<<"统计："<<endl;
	fs<<"控制特征个数："<<tieFeatureList.size()<<endl;
	fs<<"minX\t"<<"minY\t"<<"min平面\t"<<"maxX\t"<<"maxY\t"<<"max平面\t"<<"rmsX\t"<<"rmsY\t"<<"rms平面\t"<<endl;
	fs<<minX<<"\t"<<minY<<"\t"<<minXY<<"\t"
		<<maxX<<"\t"<<maxY<<"\t"<<maxXY<<"\t"
		<<rmsX<<"\t"<<rmsY<<"\t"<<rmsXY<<"\n";
	fs.close();

	return true;
}

NEWMAT::ColumnVector CalcResidue(ossimSensorModel* sensorModel,vector<ossimTieFeature> tieFeatureList)
{
	int nFeature = tieFeatureList.size();
	for(int iFeature = 0;iFeature < nFeature;iFeature++)
	{
		int nPoint = tieFeatureList[iFeature].getGroundFeature().m_Points.size();
		for(int i = 0;i < nPoint;++i)
		{
			ossimGpt goundpoint = tieFeatureList[iFeature].getGroundFeature().m_Points[i];
			ossimGpt tGpt = sensorModel->m_proj->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
			tGpt.hgt = goundpoint.hgt;
			tieFeatureList[iFeature].refGroundFeature().m_Points[i] = tGpt;
		}
	}

	NEWMAT::ColumnVector residue =  sensorModel->getResidue(tieFeatureList);

	for(int iFeature = 0;iFeature < nFeature;iFeature++)
	{
		int nPoint = tieFeatureList[iFeature].getGroundFeature().m_Points.size();
		for(int i = 0;i < nPoint;++i)
		{
			ossimDpt imagepoint = sensorModel->m_proj->forward(tieFeatureList[iFeature].getGroundFeature().m_Points[i]);
			tieFeatureList[iFeature].refGroundFeature().m_Points[i] = ossimGpt(imagepoint.x,imagepoint.y,tieFeatureList[iFeature].getGroundFeature().m_Points[i].hgt);
		}
	}
	return residue;
}

ossimMapProjection* CreateProjection(ProjectionParameters proParam, ossimKeywordlist& MapProjection)
{
	//ossimKeywordlist MapProjection;
	char* prefix="";
	// 椭球模型
	ossimString EllipsoidModel;
	if(proParam.DatumName.contains("西安80"))
	{
		EllipsoidModel = "XA";
	}
	else if(proParam.DatumName.contains("北京54"))
	{
		EllipsoidModel = "KA";
	}
	else{
		EllipsoidModel = "WE";
	}

	//中央经线
	MapProjection.add(prefix,ossimKeywordNames::ORIGIN_LATITUDE_KW, proParam.TrueOriginLatitude,true);
	MapProjection.add(prefix,ossimKeywordNames::CENTRAL_MERIDIAN_KW,proParam.TrueOriginLongitude,true);

	//使用椭球模型
	const 	ossimEllipsoid * theelli = ossimEllipsoidFactory::instance()->create(EllipsoidModel.c_str());
	theelli->saveState(MapProjection, prefix);

	ossimDpt	theMetersPerPixel;
	//X方向分辨率
	theMetersPerPixel.lon = proParam.PixelSizeX;
	//Y方向分辨率
	theMetersPerPixel.lat = proParam.PixelSizeY;

	//设置分辨率 单位为米
	MapProjection.add(prefix,	ossimKeywordNames::PIXEL_SCALE_XY_KW, theMetersPerPixel.toString().c_str(),true);
	MapProjection.add(prefix,	ossimKeywordNames::PIXEL_SCALE_UNITS_KW,	ossimUnitTypeLut::instance()->getEntryString(OSSIM_METERS),true); 

	ossimDpt	theFalseEastingNorthing;
	//东向偏移值
	theFalseEastingNorthing.lon = proParam.EastingFalse;
	//北向偏移值
	theFalseEastingNorthing.lat = proParam.NorthingFalse;
	MapProjection.add(prefix,ossimKeywordNames::FALSE_EASTING_NORTHING_KW,theFalseEastingNorthing.toString().c_str(), true);
	MapProjection.add(prefix,	ossimKeywordNames::FALSE_EASTING_NORTHING_UNITS_KW,	ossimUnitTypeLut::instance()->getEntryString(OSSIM_METERS), true);
	bool theElevationLookupFlag=true;
	MapProjection.add(prefix,	ossimKeywordNames::ELEVATION_LOOKUP_FLAG_KW,	ossimString::toString(theElevationLookupFlag), true);

	//设置TM投影
	MapProjection.add(prefix,	ossimKeywordNames::SCALE_FACTOR_KW, 1.0,true);
	MapProjection.add(prefix,	"type","ossimTransMercatorProjection",true);

	ossimMapProjection *MapPar = PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(MapProjection));
	return MapPar;
}


bool readRPBFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct)
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

bool readRPCFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct)
{
	std::fstream fs;
	fs.open(rpcFile.c_str(), ios_base::in);
	char buf[2048];
	while (fs.getline(buf, 2048))
	{
		ossimString strLine(buf);
		std::vector<ossimString> splitString = strLine.split(":", true);
		if (splitString[0].contains("LINE_OFF"))
		{
			rpcStruct.lineOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("SAMP_OFF"))
		{
			rpcStruct.sampOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LAT_OFF"))
		{
			rpcStruct.latOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LONG_OFF"))
		{
			rpcStruct.lonOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("HEIGHT_OFF"))
		{
			rpcStruct.hgtOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LINE_SCALE"))
		{
			rpcStruct.lineScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("SAMP_SCALE"))
		{
			rpcStruct.sampScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LAT_SCALE"))
		{
			rpcStruct.latScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LONG_SCALE"))
		{
			rpcStruct.lonScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("HEIGHT_SCALE"))
		{
			rpcStruct.hgtScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LINE_NUM_COEFF_"))
		{
			int i = splitString[0].after("LINE_NUM_COEFF_").toInt();
			rpcStruct.lineNumCoef[i-1] = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LINE_DEN_COEFF_"))
		{
			int i = splitString[0].after("LINE_DEN_COEFF_").toInt();
			rpcStruct.lineDenCoef[i-1] = splitString[1].toDouble();
		}
		else if (splitString[0].contains("SAMP_NUM_COEFF_"))
		{
			int i = splitString[0].after("SAMP_NUM_COEFF_").toInt();
			rpcStruct.sampNumCoef[i-1] = splitString[1].toDouble();
		}
		else if (splitString[0].contains("SAMP_DEN_COEFF_"))
		{
			int i = splitString[0].after("SAMP_DEN_COEFF_").toInt();
			rpcStruct.sampDenCoef[i-1] = splitString[1].toDouble();
		}
	}
	fs.close();
	//rpcStruct.type = 'B';
	rpcStruct.type = 'A';
	return true;
}

bool readRPBFile(ossimFilename rpcFile, ossimplugins::radiRpcModel::rpcModelStruct& rpcStruct)
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

bool readRPCFile(ossimFilename rpcFile, ossimplugins::radiRpcModel::rpcModelStruct& rpcStruct)
{
	std::fstream fs;
	fs.open(rpcFile.c_str(), ios_base::in);
	char buf[2048];
	while (fs.getline(buf, 2048))
	{
		ossimString strLine(buf);
		std::vector<ossimString> splitString = strLine.split(":", true);
		if (splitString[0].contains("LINE_OFF"))
		{
			rpcStruct.lineOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("SAMP_OFF"))
		{
			rpcStruct.sampOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LAT_OFF"))
		{
			rpcStruct.latOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LONG_OFF"))
		{
			rpcStruct.lonOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("HEIGHT_OFF"))
		{
			rpcStruct.hgtOffset = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LINE_SCALE"))
		{
			rpcStruct.lineScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("SAMP_SCALE"))
		{
			rpcStruct.sampScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LAT_SCALE"))
		{
			rpcStruct.latScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LONG_SCALE"))
		{
			rpcStruct.lonScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("HEIGHT_SCALE"))
		{
			rpcStruct.hgtScale = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LINE_NUM_COEFF_"))
		{
			int i = splitString[0].after("LINE_NUM_COEFF_").toInt();
			rpcStruct.lineNumCoef[i-1] = splitString[1].toDouble();
		}
		else if (splitString[0].contains("LINE_DEN_COEFF_"))
		{
			int i = splitString[0].after("LINE_DEN_COEFF_").toInt();
			rpcStruct.lineDenCoef[i-1] = splitString[1].toDouble();
		}
		else if (splitString[0].contains("SAMP_NUM_COEFF_"))
		{
			int i = splitString[0].after("SAMP_NUM_COEFF_").toInt();
			rpcStruct.sampNumCoef[i-1] = splitString[1].toDouble();
		}
		else if (splitString[0].contains("SAMP_DEN_COEFF_"))
		{
			int i = splitString[0].after("SAMP_DEN_COEFF_").toInt();
			rpcStruct.sampDenCoef[i-1] = splitString[1].toDouble();
		}
	}
	fs.close();
	//rpcStruct.type = 'B';
	rpcStruct.type = 'A';
	return true;
}

NEWMAT::ColumnVector CalcResidue(ossimSensorModel* sensorModel,ossimTieGptSet gptSet, bool bPixel/*=false*/, bool bLonLat/* = false*/)
{
	ossimDpt imagepoint,cimagepoint;
	ossimGpt goundpoint,tGpt;
	ossimDpt residue1,residue2;
	int i;

	int num = static_cast<int>(gptSet.size());
	NEWMAT::ColumnVector residue;
	//residue.ReSize(num*3);
	residue.ReSize(num*2);

	vector<ossimRefPtr<ossimTieGpt> >& theTPV = gptSet.refTiePoints();

	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;


	if(bPixel)
	{
		for(i = 0;i < num;++i)
		{
			ossimGpt gpt = gptSet.getTiePoints()[i]->getGroundPoint();
			if (!bLonLat)
			{
				gpt = sensorModel->m_proj->inverse(ossimDpt(gpt.lat, gpt.lon));
				gpt.hgt = gptSet.getTiePoints()[i]->getGroundPoint().hgt;
			}
			ossimDpt imagepoint = gptSet.getTiePoints()[i]->getImagePoint() - sensorModel->forward(gpt);

			residue.element(i * 2 + 0) = imagepoint.x;
			residue.element(i * 2 + 1) = imagepoint.y;
		}
	}
	else
	{
		for(i = 0;i < num;++i)
		{
			imagepoint = gptSet.getTiePoints()[i]->getImagePoint();
			//sensorModel->lineSampleToWorld(imagepoint,tGpt);
			double hgt = gptSet.getTiePoints()[i]->getGroundPoint().hgt;
			sensorModel->lineSampleHeightToWorld(imagepoint, hgt, tGpt);
			//ossimDpt checkDpt;
			//sensorModel->worldToLineSample(tGpt, checkDpt);
			//checkDpt = checkDpt - imagepoint;
			if(sensorModel->m_proj) tGpt.datum(sensorModel->m_proj->getDatum());	//loong

			if (bLonLat)
			{
				residue1 = sensorModel->m_proj->forward(tGpt);
				residue2 = sensorModel->m_proj->forward(gptSet.getTiePoints()[i]->getGroundPoint());
			}
			else
			{
				//residue1 = ossimDpt(tGpt.lat, tGpt.lon);
				residue1 = sensorModel->m_proj->forward(tGpt);
				ossimGpt gpt = gptSet.getTiePoints()[i]->getGroundPoint();
				residue2 = ossimDpt(gpt.lat,gpt.lon);
			}

			residue.element(i * 2 + 0) = residue2.x - residue1.x;
			residue.element(i * 2 + 1) = residue2.y - residue1.y;
			//residue.element(i * 3 + 2) = gpt.hgt - tGpt.hgt;
		}
	}

	return residue;
}

bool OutputReport(ossimFilename reportfile, ossimSensorModel* sensorModel, ossimTieGptSet* ctrlSet, ossimTieGptSet* chkSet, bool bPixel/*=false*/, bool bLonLat/* = false*/)
{
	int i;
	fstream fs;
	fs.open(reportfile.c_str(), ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(4);

	if(bPixel)
	{
		fs<<"残差报告\n\n"<<"残差单位：像素\n";
	}
	else
	{
		fs<<"残差报告\n\n"<<"残差单位：米\n";
	}

	NEWMAT::ColumnVector residue1 = CalcResidue(sensorModel, *ctrlSet, bPixel, bLonLat);
	NEWMAT::ColumnVector residue2 = CalcResidue(sensorModel, *chkSet, bPixel, bLonLat);

	fs<<"控制点点数："<<ctrlSet->size();
	if(static_cast<int>(ctrlSet->size()) > 0)
	{
		double ResidueX = 0.0;
		double ResidueY = 0.0;
		double ResidueZ = 0.0;

		for(i = 0;i < static_cast<int>(ctrlSet->size());i++)
		{
			ResidueX += residue1.element(i * 2 + 0) * residue1.element(i * 2 + 0);
			ResidueY += residue1.element(i * 2 + 1) * residue1.element(i * 2 + 1);
			//ResidueZ += residue1.element(i * 3 + 2) * residue1.element(i * 3 + 2);
		}
		ResidueX = sqrt(ResidueX / static_cast<int>(ctrlSet->size()));
		ResidueY = sqrt(ResidueY / static_cast<int>(ctrlSet->size()));
		ResidueZ = sqrt(ResidueZ / static_cast<int>(ctrlSet->size()));
		fs<<"\t\t"<<"均方差：X "<<ResidueX<<"\tY "<<ResidueY<<"\tZ "<<ResidueZ<<endl;
	}
	else
		fs<<endl;

	fs<<"检查点点数："<<chkSet->size();
	if(static_cast<int>(chkSet->size()) > 0)
	{
		double ResidueX = 0.0;
		double ResidueY = 0.0;
		double ResidueZ = 0.0;

		for(i = 0;i < static_cast<int>(chkSet->size());i++)
		{
			ResidueX += residue2.element(i * 2 + 0) * residue2.element(i * 2 + 0);
			ResidueY += residue2.element(i * 2 + 1) * residue2.element(i * 2 + 1);
			//ResidueZ += residue2.element(i * 3 + 2) * residue2.element(i * 3 + 2);
		}
		ResidueX = sqrt(ResidueX / static_cast<int>(chkSet->size()));
		ResidueY = sqrt(ResidueY / static_cast<int>(chkSet->size()));
		ResidueZ = sqrt(ResidueZ / static_cast<int>(chkSet->size()));
		fs<<"\t\t"<<"均方差：X "<<ResidueX<<"\tY "<<ResidueY<<"\tZ "<<ResidueZ<<endl;
	}
	else
		fs<<endl;

	fs<<"清单："<<endl;
	fs<<"用于计算模型的控制点点数："<<ctrlSet->size()<<endl;
	fs<<setw(10)<<"标示"<<"\t"
		<<setw(10)<<"误差"<<"\t"
		<<setw(10)<<"残差X"<<"\t"
		<<setw(10)<<"残差Y"<<"\t"
		<<setw(8)<<"类型"<<"\t"
		<<setw(10)<<"高程"<<"\t"
		<<setw(15)<<"实际X"<<"\t"
		<<setw(15)<<"实际Y"<<"\t"
		<<setw(15)<<"对照X"<<"\t"
		<<setw(15)<<"对照Y"<<"\t"<<endl;


	//residue = CalcResidue(*m_CtrTieGptSet);
	for(i = 0;i < static_cast<int>(ctrlSet->size());i++)
	{
		fs<<setw(10)<<ctrlSet->getTiePoints()[i]->GcpNumberID<<"\t"
			<<setw(10)<<sqrt(residue1.element(i * 2 + 0)*residue1.element(i * 2 + 0) + residue1.element(i * 2 + 1)*residue1.element(i * 2 + 1))<<"\t"
			<<setw(10)<<residue1.element(i * 2 + 0)<<"\t"
			<<setw(10)<<residue1.element(i * 2 + 1)<<"\t"
			<<setw(8)<<"GCP"<<"\t";
		if (bPixel)
		{
			fs<<setw(10)<<ctrlSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
				<<setw(15)<<ctrlSet->getTiePoints()[i]->getImagePoint().x<<"\t"
				<<setw(15)<<ctrlSet->getTiePoints()[i]->getImagePoint().y<<"\t"
				<<setw(15)<<ctrlSet->getTiePoints()[i]->getImagePoint().x + residue1.element(i * 2 + 0)<<"\t"
				<<setw(15)<<ctrlSet->getTiePoints()[i]->getImagePoint().y + residue1.element(i * 2 + 1)<<"\n";
		}
		else
		{
			if (bLonLat)
			{
				ossimDpt refDpt = sensorModel->m_proj->forward(ctrlSet->getTiePoints()[i]->getGroundPoint());
				fs<<setw(10)<<ctrlSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
					<<setw(15)<<refDpt.x<<"\t"
					<<setw(15)<<refDpt.y<<"\t"
					<<setw(15)<<refDpt.x + residue1.element(i * 2 + 0)<<"\t"
					<<setw(15)<<refDpt.y + residue1.element(i * 2 + 1)<<"\n";
			}
			else{
				fs<<setw(10)<<ctrlSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
					<<setw(15)<<ctrlSet->getTiePoints()[i]->getGroundPoint().lat<<"\t"
					<<setw(15)<<ctrlSet->getTiePoints()[i]->getGroundPoint().lon<<"\t"
					<<setw(15)<<ctrlSet->getTiePoints()[i]->getGroundPoint().lat + residue1.element(i * 2 + 0)<<"\t"
					<<setw(15)<<ctrlSet->getTiePoints()[i]->getGroundPoint().lon + residue1.element(i * 2 + 1)<<"\n";
			}
		}
	}

	//fs<<endl;

	//residue = CalcResidue(*m_ChkTieGptSet);
	for(i = 0;i < static_cast<int>(chkSet->size());i++)
	{
		fs<<setw(10)<<chkSet->getTiePoints()[i]->GcpNumberID<<"\t"
			<<setw(10)<<sqrt(residue2.element(i * 2 + 0)*residue2.element(i * 2 + 0) + residue2.element(i * 2 + 1)*residue2.element(i * 2 + 1))<<"\t"
			<<setw(10)<<residue2.element(i * 2 + 0)<<"\t"
			<<setw(10)<<residue2.element(i * 2 + 1)<<"\t"
			//<<residue2.element(i * 3 + 2)<<"\t"
			<<setw(8)<<"CHECK"<<"\t";
		if (bPixel)
		{
			fs<<setw(10)<<chkSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
				<<setw(15)<<chkSet->getTiePoints()[i]->getImagePoint().x<<"\t"
				<<setw(15)<<chkSet->getTiePoints()[i]->getImagePoint().y<<"\t"
				<<setw(15)<<chkSet->getTiePoints()[i]->getImagePoint().x + residue2.element(i * 2 + 0)<<"\t"
				<<setw(15)<<chkSet->getTiePoints()[i]->getImagePoint().y + residue2.element(i * 2 + 1)<<"\n";
		}
		else
		{
			if (bLonLat)
			{
				ossimDpt refDpt = sensorModel->m_proj->forward(chkSet->getTiePoints()[i]->getGroundPoint());
				fs<<setw(10)<<chkSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
					<<setw(15)<<refDpt.x<<"\t"
					<<setw(15)<<refDpt.y<<"\t"
					<<setw(15)<<refDpt.x + residue2.element(i * 2 + 0)<<"\t"
					<<setw(15)<<refDpt.y + residue2.element(i * 2 + 1)<<"\n";
			}
			else{
				fs<<setw(10)<<chkSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
					<<setw(15)<<chkSet->getTiePoints()[i]->getGroundPoint().lat<<"\t"
					<<setw(15)<<chkSet->getTiePoints()[i]->getGroundPoint().lon<<"\t"
					<<setw(15)<<chkSet->getTiePoints()[i]->getGroundPoint().lat + residue2.element(i * 2 + 0)<<"\t"
					<<setw(15)<<chkSet->getTiePoints()[i]->getGroundPoint().lon + residue2.element(i * 2 + 1)<<"\n";
			}
		}
	}

	fs.close();

	return true;
}

ossimRpcModel* createRpcModelFromPoints(ossimFilename gcpFile)
{
	ossimTieGptSet* gptSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	mylib::readGcpFile(gcpFile, gptSet, chkSet, &prjKwl);
	mylib::get_elevation(gptSet, prjKwl, 0.0);

	ossimMapProjection* MapPar = PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(prjKwl));

	int num = static_cast<int>(gptSet->getTiePoints().size());
	ossimplugins::radiRpcSolver *solver = new ossimplugins::radiRpcSolver(true, false);
	vector < ossimDpt > imagePoints;
	vector < ossimGpt > groundControlPoints;
	ossimGpt gpt;
	if (MapPar && !MapPar->isGeographic())
	{
		for(int i = 0;i < (int)gptSet->getTiePoints().size();i++)
		{
			ossimGpt gpt = MapPar->inverse(ossimDpt(gptSet->getTiePoints()[i]->getGroundPoint().lat, 
				gptSet->getTiePoints()[i]->getGroundPoint().lon));
			gpt.hgt = gptSet->getTiePoints()[i]->getGroundPoint().hgt;
			groundControlPoints.push_back(gpt);
			//groundControlPoints.push_back(prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint());
			imagePoints.push_back(gptSet->getTiePoints()[i]->getImagePoint());
		}
	}
	else
	{
		for(int i = 0;i < (int)gptSet->getTiePoints().size();i++)
		{
			groundControlPoints.push_back(gptSet->getTiePoints()[i]->getGroundPoint());
			imagePoints.push_back(gptSet->getTiePoints()[i]->getImagePoint());
		}

	}

	//solver->solveCoefficients(imagePoints, groundControlPoints, true);
	//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LS, 1e-5);
	//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::RIDGE, 1e-5);
	solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LASSO, 1e-5);

	ossimImageGeometry *imageGeom = solver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimRpcModel* rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom, "projection.");
	return rpcModel;
}

void Hj1CreateRpcs(ossimSensorModel* sensorModel, ossimIrect theImageClipRect, ossimFilename workfold)
{
	ossimFilename projectionFile = workfold + "\\projection.txt";
	ossimFilename reportFile = workfold + "\\report0.txt";
	ossimFilename rpcFile = workfold + "\\rpcStruct.txt";

	// statistic the max and min Height
	ossim_uint32 xSamples = 20;
	ossim_uint32 ySamples = 20;

	double max_Height = -1.0e10;
	double min_Height = 1.0e10;
	ossim_uint32 x,y;
	ossim_float64 w = theImageClipRect.width();
	ossim_float64 h = theImageClipRect.height();
	ossimGpt gpt;
	ossimGpt defaultGround;
	if(ySamples < 1) ySamples = 12;
	if(xSamples < 1) xSamples = 12;

	double xnorm;
	double ynorm;
	ossimDpt ul = theImageClipRect.ul();
	ossimDpt shiftTo0(-ul.x, -ul.y);
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

			sensorModel->lineSampleToWorld(dpt, gpt);
			if (!gpt.hasNans())
			{
				if (gpt.height() > max_Height)
				{
					max_Height = gpt.height();
				}
				if (gpt.height() < min_Height)
				{
					min_Height = gpt.height();
				}
			}
		}
	}

	min_Height -= 500.0;
	max_Height += 500.0;
	ossimTieGptSet *gptSet = NULL;
	int nLevels = 5;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_Height + i*(max_Height-min_Height)/(nLevels-1);
		create3DGridPoints(theImageClipRect, *sensorModel, hgt, gptSet, 15, 15, true, false);
	}
	ossimTieGptSet *chkSet = NULL;
	nLevels = 7;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_Height + i*(max_Height-min_Height)/(nLevels-1);
		create3DGridPoints(theImageClipRect, *sensorModel, hgt, chkSet, 21, 21, true, false);
	}

	ossimplugins::radiRpcSolver *solver = new ossimplugins::radiRpcSolver(true, false);
	int num = (int)gptSet->getTiePoints().size();
	vector < ossimDpt > imagePoints(num);
	vector < ossimGpt > groundControlPoints(num);
	for(int i = 0;i < num;i++)
	{
		groundControlPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
		imagePoints[i] = gptSet->getTiePoints()[i]->getImagePoint();
	}
	//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LS, 1e-8);
	solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::RIDGE, 1e-4);
	//solver->solveCoefficients(imagePoints, groundControlPoints, true, ossimplugins::radiRpcSolver::LASSO, 1e-5);

	ossimImageGeometry *imageGeom = solver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimRpcModel* rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom, "projection.");


	fstream rpcStructFile;
	rpcStructFile.open(rpcFile.c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStructFile);
	rpcStructFile.close();
	//OutputReport1(reportFile0, prj.m_sensorModel, gptSet, chkSet);
	mylib::OutputReport(reportFile, rpcModel, gptSet, chkSet, false, false);
}
}