#include <iostream>
#include <iterator>
#include <fstream>
using namespace std;


#include <func.h>
#include <mprojectdefine.h>
#include "..\QuickbirdRpcModel.h"

#include <strUtil.h>
#include <fileUtil.h>
using namespace mylib;

void AlosCreateVirtualGpcs(ossimFilename alosFold, ossimFilename virtualGcpFile)
{
	//MyProject prj;
	//AlosPRISM alosUti(alosFold);
	////AlosAVNIR2 alosUti(alosFold);
	//if(!alosUti.getInitState())
	//{
	//	cout<<"warning: Alos PRISM数据\""<<alosFold<<"\"初始化失败！"<<endl;
	//	return;
	//}

	//prj.m_MapProjection = alosUti.getKeywordlist();
	//prj.m_MapPar = alosUti.getMapProjection();

	////double hgtOffset = alosUti.getRpcModelStruct().hgtOffset;
	////double hgtScale = alosUti.getRpcModelStruct().hgtScale;
	////double minH = (hgtScale-hgtOffset)*0.5;
	////double maxH = (hgtScale+hgtOffset)*0.5;
	//prj.m_ImgFileNameUnc = ossimFilename(alosUti.m_FileTIF);
	//prj.m_OutBandList.clear();
	//prj.m_OutBandList.push_back(0);

	//prj.m_ModelType = RPCType;
	//if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
	//{
	//	cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
	//	return;
	//}
	//ossimRpcModel *rpcModel = new ossimRpcModel;
	//rpcModel->setAttributes(alosUti.getRpcModelStruct());
	//prj.m_sensorModel = rpcModel;
	//prj.m_sensorModel->m_proj = prj.m_MapPar;

	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(alosUti.m_FileTIF);
	////double minH = 0.0;
	////double maxH = 3000.0;
	//if(!handler) return;   //应该弹出警告对话框
	//
	////ossimTieGptSet* gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 11, 11);
	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 21, 21);


	//double minH = 1000000.0;
	//double maxH = -1000000.0;
	//for (int i = 0;i < (int)chkSet->getTiePoints().size();++i)
	//{
	//	double h = chkSet->getTiePoints()[i]->hgt;
	//	minH = (h<minH)?h:minH;
	//	maxH = (h>maxH)?h:maxH;
	//}
	//int nElevationLayers = 5;
	//double hStep = (maxH-minH)/(nElevationLayers-1);
	//ossimTieGptSet* gptSet = NULL;
	//for (int i = 0;i < nElevationLayers;++i)
	//{
	//	double h = minH + i * hStep;
	//	createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, h, gptSet, 11, 11);
	//}
	//

	//handler->close();
	//delete handler;

	//ossimFilename ctrFile = alosFold + "\\ctr.txt";
	//ossimFilename chkFile = alosFold + "\\chk.txt";
	//ossimFilename rpcFile = alosFold + "\\rpc.txt";

	//fstream rpcStruct;
	//prj.m_sensorModel->saveState(prj.geom);
	//rpcModel->loadState(prj.geom);
	//rpcStruct.open(rpcFile.c_str(), ios_base::out);
	//rpcModel->saveRpcModelStruct(rpcStruct);
	//rpcStruct.close();

	//prj.SavePointToFile(ctrFile, gptSet, NULL, false);
	//prj.SavePointToFile(chkFile, chkSet, NULL, false);
}


//void batch_AlosCreateVirtualGpcs(QString inputDir)
//{
//	QStringList allFiles;
//	QStringList strFilter;
//	strFilter<<"IMG-ALPSMW*_O1B2R_U*.tif"<<"IMG-alav2*_O1B2R_U*.tif";
//	QFindFile(inputDir, allFiles, strFilter);
//	int nTotal = (int)allFiles.size();
//	for (int i = 0; i < nTotal; i++)
//	{
//		AlosCreateVirtualGpcs(ossimFilename(QBeforeLast(allFiles[i], '\\').toLatin1()), "");
//	}
//}



//void AlosRpc()
//{
//	//ossimFilename workfold = "I:\\testdata\\alos\\hlj\\d1000771-006";
//	ossimFilename workfold = "I:\\testdata\\alos_batch\\multiple\\d1001884-021_alav2a249442720_o1b2r_u_rpc";
//	ossimFilename sourcefile = workfold + "\\IMG-ALPSMW230332695_O1B2R_UW.tif";
//	//ossimFilename elevationpath = "D:\\workspace\\dem";
//	ossimFilename gcpfile = workfold + "\\gcp.txt";
//	ossimFilename ctrfile = workfold + "\\controlp_.txt";
//	ossimFilename chkfile = workfold + "\\chk_.txt";
//	ossimFilename reportfile = workfold + "\\RpcReport.txt";
//
//	MyProject prj;
//	//AlosPRISM alosUti(workfold);
//	AlosAVNIR2 alosUti(workfold);
//	if(!alosUti.getInitState())
//	{
//		cout<<"warning: Alos PRISM数据\""<<workfold<<"\"初始化失败！"<<endl;
//		return;
//	}
//	//prj.theMgr = ossimElevManager::instance();
//	//if(!prj.theMgr->loadElevationPath(ossimFilename(elevationpath)))
//	//{
//	//	cout<<"warning: 加载DEM失败！"<<endl;
//	//	return;
//	//}
//	prj.m_MapProjection = alosUti.getKeywordlist();
//	prj.m_MapPar = alosUti.getMapProjection();
//	//prj.ReadGcpAndProjection(gcpfile);
//	prj.ReadGcpAndProjection(ctrfile, prj.m_CtrlGptSet, prj.m_ChkGptSet);
//	prj.ReadGcpAndProjection(chkfile, prj.m_ChkGptSet, prj.m_ChkGptSet);
//	prj.theMgr = ossimElevManager::instance();
//	//prj.theMgr->loadElevationPath(ossimFilename(elevationpath));
//	//prj.m_DemPath = ossimFilename(elevationpath);
//	prj.GetElevations(prj.m_CtrlGptSet);
//	prj.GetElevations(prj.m_ChkGptSet);
//
//	//prj.SavePointToFile(workfold+"\\gcp_hgt.txt", prj.m_CtrlGptSet, prj.m_ChkGptSet, true);
//	prj.SavePoint2Shape(workfold+"\\ctr_point.shp", prj.m_CtrlGptSet);
//	prj.SavePoint2Shape(workfold+"\\chk_point.shp", prj.m_ChkGptSet);
//
//
//	prj.m_ImgFileNameUnc = ossimFilename(alosUti.m_FileTIF);
//	prj.m_OutBandList.clear();
//	prj.m_OutBandList.push_back(0);
//
//	prj.m_ModelType = RPCType;
//	if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
//	{
//		cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
//		return;
//	}
//	ossimRpcModel *rpcModel = new ossimRpcModel;
//	rpcModel->setAttributes(alosUti.getRpcModelStruct());
//	prj.m_sensorModel = rpcModel;
//	prj.m_sensorModel->m_proj = prj.m_MapPar;
//
//	/*ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
//	if(!handler) return;   //应该弹出警告对话框
//	prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
//	prj.m_sensorModel->m_proj = prj.m_MapPar;
//	handler->close();
//	delete handler;*/
//
//	prj.m_sensorModel->m_bHgtOptimize = true;
//	//prj.m_sensorModel->m_bHgtOptimize = false;
//	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
//	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, true);
//}
//
//
//void AlosRpc2()
//{
//	//ossimFilename workfold = "H:\\testdata\\AlosData\\d1001111-001";
//	ossimFilename workfold = "D:\\workspace\alos\\d1001111-001";
//	ossimFilename sourcefile = workfold + "\\IMG-ALPSMW220853075_O1B2R_UW1.tif";
//	ossimFilename elevationpath = "D:\\workspace\\dem";
//	ossimFilename gcpfile = workfold + "\\gcp.txt";
//	ossimFilename reportfile = workfold + "\\RpcReport.txt";
//
//	MyProject prj;
//	AlosPRISM alosUti(workfold);
//	if(!alosUti.getInitState())
//	{
//		cout<<"warning: Alos PRISM数据\""<<workfold<<"\"初始化失败！"<<endl;
//		return;
//	}
//	prj.theMgr = ossimElevManager::instance();
//	if(!prj.theMgr->loadElevationPath(ossimFilename(elevationpath)))
//	{
//		cout<<"warning: 加载DEM失败！"<<endl;
//		return;
//	}
//
//	prj.m_MapProjection = alosUti.getKeywordlist();
//	prj.m_MapPar = alosUti.getMapProjection();
//
//	prj.ReadGcpAndProjection(gcpfile);
//	prj.theMgr = ossimElevManager::instance();
//	prj.theMgr->loadElevationPath(ossimFilename(elevationpath));
//	prj.m_DemPath = ossimFilename(elevationpath);
//	prj.GetElevations(prj.m_CtrlGptSet);
//	prj.GetElevations(prj.m_ChkGptSet);
//
//	prj.m_ImgFileNameUnc = ossimFilename(alosUti.m_FileTIF);
//	prj.m_OutBandList.clear();
//	prj.m_OutBandList.push_back(0);
//
//	prj.m_ModelType = RPCType;
//	if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
//	{
//		cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
//		return;
//	}
//	ossimRpcModel *rpcModel = new ossimRpcModel;
//	rpcModel->setAttributes(alosUti.getRpcModelStruct());
//	prj.m_sensorModel = rpcModel;
//	prj.m_sensorModel->m_proj = prj.m_MapPar;
//
//	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
//	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
//
//
//	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
//	if(!handler) return;   //应该弹出警告对话框
//	ossimTieGptSet* gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 10, 10);
//	ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 20, 20);
//	handler->close();
//	delete handler;
//
//	int num = static_cast<int>(gptSet->getTiePoints().size());
//	ossimRpcSolver *solver = new ossimRpcSolver(true, false);
//	vector < ossimDpt > imagePoints(num);
//	vector < ossimGpt > groundControlPoints(num);
//	for(int i = 0;i < num;i++)
//	{
//		groundControlPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
//		ossimGpt gpt = prj.m_MapPar->inverse(ossimDpt(groundControlPoints[i].lat, groundControlPoints[i].lon));
//		gpt.hgt = groundControlPoints[i].hgt;
//		groundControlPoints[i] = gpt;
//		imagePoints[i] = gptSet->getTiePoints()[i]->getImagePoint();
//	}
//
//	solver->solveCoefficients(imagePoints, groundControlPoints);
//
//	ossimImageGeometry *imageGeom = solver->createRpcModel();
//	rpcModel = (ossimRpcModel*)imageGeom->getProjection();
//
//	prj.m_sensorModel = rpcModel;
//	prj.m_sensorModel->m_proj = prj.m_MapPar;
//	fstream rpcStruct;
//	prj.m_sensorModel->saveState(prj.geom);
//	rpcModel->loadState(prj.geom);
//	rpcStruct.open((workfold+"\\rpcStruct0.txt").c_str(), ios_base::out);
//	rpcModel->saveRpcModelStruct(rpcStruct);
//	rpcStruct.close();
//	prj.SavePointToFile(workfold+"\\virtualGpt.txt", gptSet, chkSet);
//	prj.OutputReport(workfold+"\\report0.txt", prj.m_sensorModel, gptSet, chkSet, false);
//	return;
//
//	prj.m_sensorModel->m_bHgtOptimize = true;
//	//prj.m_sensorModel->m_bHgtOptimize = false;
//	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
//	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, true);
//}
//
//
//void AlosRpc_createTieGptSet()
//{
//	//ossimFilename workfold = "H:\\testdata\\AlosData\\d1001111-001";
//	//ossimFilename workfold = "D:\\workspace\\alos\\d1001111-001";
//	ossimFilename workfold = "D:\\workspace\\alos\\hlj\\d1001884-021_alav2a249442720_o1b2r_u_rpc";
//	ossimFilename sourcefile = workfold + "\\img-alav2a249442720_o1b2r_u.tif";
//	workfold = "D:\\workspace\\alos\\hlj\\d1000771-006";
//	sourcefile = workfold + "\\IMG-ALPSMW230332695_O1B2R_UW.tif";
//	//ossimFilename elevationpath = "D:\\workspace\\dem-hgt";
//	//ossimFilename elevationpath = "I:\\testdata\\alos\\beijing\\dem-mosaic";
//	ossimFilename elevationpath = "I:\\testdata\\ASTER GDEM V2";
//	ossimFilename gcpfile = workfold + "\\gcp.txt";
//	ossimFilename virtual_gcpfile = workfold + "\\virtual_gcp.txt";
//	ossimFilename virtual_chkfile = workfold + "\\virtual_chk.txt";
//	ossimFilename reportfile = workfold + "\\RpcReport.txt";
//
//	MyProject prj;
//	AlosPRISM alosUti(workfold);
//	//AlosAVNIR2 alosUti(workfold);
//	if(!alosUti.getInitState())
//	{
//		cout<<"warning: Alos PRISM数据\""<<workfold<<"\"初始化失败！"<<endl;
//		return;
//	}
//
//	ossimRpcModel *rpcModel = new ossimRpcModel;
//	rpcModel->setAttributes(alosUti.getRpcModelStruct());
//	prj.m_sensorModel = rpcModel;
//	prj.m_sensorModel->m_proj = prj.m_MapPar;
//
//	ossimImageHandler *handler = ossimImageHandlerRegistry::instance()->open(sourcefile);
//	if (!handler) return;   //应该弹出警告对话框
//	ossimTieGptSet* gptSet = new ossimTieGptSet;
//	double max_height = 300.0;
//	double min_height = 0.0;
//
//	int nLevels = 5;
//	for (int i = 0; i < nLevels; ++i)
//	{
//		double hgt = min_height + i*(max_height - min_height) / (nLevels - 1);
//		create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 15, 15, true, true);
//	}
//	//int nLevels = 1;
//	//for (int i=0;i < nLevels;++i)
//	//{
//	//	double hgt = 0.0;// min_height + i*(max_height - min_height) / (nLevels - 1);
//	//	create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 1, 23999, false, false);
//	//}
//
//	ossimTieGptSet* chkSet = new ossimTieGptSet;
//	nLevels = 10;
//	for (int i = 0; i < nLevels; ++i)
//	{
//		double hgt = min_height + i*(max_height - min_height) / (nLevels - 1);
//		create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, chkSet, 30, 30, true, true);
//	}
//	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 0.0, gptSet, 10, 10, false);
//	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 1000.0, gptSet, 10, 10, false);
//	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 2000.0, gptSet, 10, 10, false);
//	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 3000.0, gptSet, 10, 10, false);
//	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 4000.0, gptSet, 10, 10, false);
//	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 5000.0, gptSet, 10, 10, false);
//	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 20, 20, true);
//	//gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 11, 11, false);
//	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 21, 21, false);
//	handler->close();
//	delete handler;
//
//	prj.SavePointToFile(virtual_gcpfile, gptSet, NULL);
//	prj.SavePointToFile(virtual_chkfile, chkSet, NULL);
//
//	////prj.theMgr = ossimElevManager::instance();
//	////if(!prj.theMgr->loadElevationPath(ossimFilename(elevationpath)))
//	////{
//	////	cout<<"warning: 加载DEM失败！"<<endl;
//	////	return;
//	////}
//
//	//prj.m_MapProjection = alosUti.getKeywordlist();
//	//prj.m_MapPar = alosUti.getMapProjection();
//
//	////prj.ReadGcpAndProjection(gcpfile);
//	////prj.theMgr = ossimElevManager::instance();
//	////prj.theMgr->loadElevationPath(ossimFilename(elevationpath));
//	////prj.m_DemPath = ossimFilename(elevationpath);
//	////prj.GetElevations(prj.m_CtrlGptSet);
//	////prj.GetElevations(prj.m_ChkGptSet);
//
//	////prj.m_ImgFileNameUnc = ossimFilename(alosUti.m_FileTIF);
//	////prj.m_OutBandList.clear();
//	////prj.m_OutBandList.push_back(0);
//
//	//prj.m_ModelType = RPCType;
//	//if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
//	//{
//	//	cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
//	//	return;
//	//}
//	//ossimRpcModel *rpcModel = new ossimRpcModel;
//	//rpcModel->setAttributes(alosUti.getRpcModelStruct());
//	//prj.m_sensorModel = rpcModel;
//	//prj.m_sensorModel->m_proj = prj.m_MapPar;
//
//	////prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
//	////prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
//
//	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
//	//if(!handler) return;   //应该弹出警告对话框
//	//ossimTieGptSet* gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 10, 10, true);
//	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 20, 20, true);
//	////ossimTieGptSet* gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 10, 10, false);
//	////ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 20, 20, false);
//	//handler->close();
//	//delete handler;
//
//	//prj.SavePointToFile(virtual_gcpfile, gptSet, chkSet);
//}


void Landsat_create3DGridPoints()
{

	ossimFilename workfold = "D:\\workspace\\Landsat\\LS5_TM_20090721_050335_050401_146030_FASTB_L2";
	//ossimFilename workfold = "D:\\testdata\\LD2010003816";
	ossimFilename sourcefile = workfold + "\\header.dat";
	ossimFilename elevationpath = "D:\\workspace\\dem-hgt";
	ossimFilename gcpfile = workfold + "\\gcps.txt";
	ossimFilename virtual_gcpfile = workfold + "\\virtual_gcp.txt";
	ossimFilename virtual_chkfile = workfold + "\\virtual_chk.txt";

	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(elevationpath));//
	prj.m_DemPath=ossimFilename(elevationpath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(6);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);
	prj.InitiateSensorModel(sourcefile);

	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	if(!handler) return;   //应该弹出警告对话框
	ossimTieGptSet* gptSet = new ossimTieGptSet;
	double max_height = 4000.0;
	double min_height = 0.0;
	int nLevels = 5;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height + i*(max_height-min_height)/(nLevels-1);
		create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 11, 11, false);
	}

	ossimTieGptSet* chkSet = new ossimTieGptSet;
	nLevels = 6;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height + i*(max_height-min_height)/(nLevels-1);
		create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, chkSet, 21, 21, false);
	}
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 0.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 1000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 2000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 3000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 4000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 5000.0, gptSet, 10, 10, false);
	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 20, 20, true);
	//gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 11, 11, false);
	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 21, 21, false);
	handler->close();
	delete handler;

	prj.SavePointToFile(virtual_gcpfile, gptSet, NULL);
	prj.SavePointToFile(virtual_chkfile, chkSet, NULL);
}


void Spot_create3DGridPoints()
{
	ossimFilename workfold = "D:\\workspace\\Spot5\\beijing";
	workfold = "D:\\workspace\\Spot5\\shisanling";
	workfold = "D:\\workspace\\Spot5\\tianjin";
	//workfold = "D:\\workspace\\Spot5\\279269-20040523-2.5";
	//ossimFilename workfold = "I:\\testdata\\spot5\\wuhan";
	//ossimFilename workfold = "I:\\wuhan\\spot\\2005";
	//ossimFilename workfold = "E:\\beijing";
	MyProject prj;
	ossimFilename sourcefile = workfold + "\\281268-20040602-2.5\\scene01\\imagery.tif";
	sourcefile = workfold + "\\279269-20040523-2.5\\scene01\\imagery.tif";
	sourcefile = workfold + "\\282270-20040518-2.5\\scene01\\imagery.tif";
	//sourcefile = workfold + "\\282270-20040518-10\\scene01\\imagery.tif";
	//sourcefile = workfold + "\\scene01\\imagery.tif";
	//ossimFilename sourcefile = workfold + "\\10\\SCENE01\\imagery.tif";
	//ossimFilename sourcefile = workfold + "\\2.5\\SCENE01\\imagery.tif";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";
	ossimFilename gcpfile = workfold + "\\gcp_auto.txt";
	//gcpfile = workfold + "\\gcp_auto_pci.txt";

	//ossimFilename demPath = "I:\\spot5\\beijing";
	ossimFilename demPath = "D:\\workspace\\dem\\aster30";
	//ossimFilename demPath = "E:\\dem_bj";
	//ossimFilename chkFile = workfold + "\\281268_2m5_gcp09033_hgt.txt";
	ossimFilename chkFile = workfold + "\\gcp_auto_pci.txt";
	//ossimFilename chkFile = workfold + "\\gcp_all_hgt.txt";

	//ossimFilename gcp_hgt_File = workfold + "\\wgs84_hgt.txt";
	//ossimFilename chkFile = workfold + "\\gcp.txt";
	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report_orbit.txt";
	ossimFilename virtual_gcpfile = workfold + "\\virtual_gcp.txt";
	ossimFilename virtual_chkfile = workfold + "\\virtual_chk.txt";
	ossimFilename reportfile = workfold + "\\report.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	//prj.SavePointToFile(gcp_hgt_File, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);


	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	if(!handler) return;   //应该弹出警告对话框
	ossimTieGptSet* gptSet = new ossimTieGptSet;
	double max_height = 1100.0;
	double min_height = 0.0;

	ossimIrect imageRect(handler->getImageRectangle().ul(),
		handler->getImageRectangle().lr());
	//ossimIrect imageRect(handler->getImageRectangle().lr()*0.5, 
	//	handler->getImageRectangle().lr());
	//ossimIrect imageRect(handler->getImageRectangle().ul(),
	//	handler->getImageRectangle().lr()*0.5);
	int nLevels = 5;
	for (int i = 0; i < nLevels; ++i)
	{
		double hgt = min_height + i*(max_height - min_height) / (nLevels - 1);
		create3DGridPoints(imageRect, *prj.m_sensorModel, hgt, gptSet, 15, 15, true, false);
	}
	//int nLevels = 1;
	//for (int i=0;i < nLevels;++i)
	//{
	//	double hgt = 0.0;// min_height + i*(max_height - min_height) / (nLevels - 1);
	//	create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 1, 23999, false, false);
	//}

	ossimTieGptSet* chkSet = new ossimTieGptSet;
	nLevels = 10;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height + i*(max_height-min_height)/(nLevels-1);
		create3DGridPoints(imageRect, *prj.m_sensorModel, hgt, chkSet, 30, 30, true, false);
	}
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 0.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 1000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 2000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 3000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 4000.0, gptSet, 10, 10, false);
	//create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, 5000.0, gptSet, 10, 10, false);
	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 20, 20, true);
	//gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 11, 11, false);
	//ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel, 21, 21, false);
	handler->close();
	delete handler;

	prj.SavePointToFile(virtual_gcpfile, gptSet, NULL);
	prj.SavePointToFile(virtual_chkfile, chkSet, NULL);
}


void Theos_create3DGridPoints()
{
	ossimFilename workfold = "D:\\workspace\\Theos\\112577-116719-TH2014002255(1A1116)";
	//ossimFilename workfold = "E:\\beijing";
	MyProject prj;
	ossimFilename sourcefile = workfold + "\\THEOS1_20140221065605507_14002255-0_01_1_1\\imagery.tif";
	//ossimFilename sourcefile = workfold + "\\10\\SCENE01\\imagery.tif";
	//ossimFilename sourcefile = workfold + "\\2.5\\SCENE01\\imagery.tif";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";
	//ossimFilename gcpfile = workfold + "\\wgs84_hgt.txt";
	ossimFilename gcpfile = workfold + "\\gcps.txt";

	//ossimFilename demPath = "I:\\spot5\\beijing";
	ossimFilename demPath = "D:\\workspace\\dem\\aster30";
	//ossimFilename demPath = "E:\\dem_bj";
	//ossimFilename chkFile = workfold + "\\281268_2m5_gcp09033_hgt.txt";
	ossimFilename chkFile = workfold + "\\wgs84_gcp.txt";
	//ossimFilename chkFile = workfold + "\\gcp_all_hgt.txt";

	//ossimFilename gcp_hgt_File = workfold + "\\wgs84_hgt.txt";
	//ossimFilename chkFile = workfold + "\\gcp.txt";
	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report_orbit.txt";
	ossimFilename virtual_gcpfile = workfold + "\\virtual_gcp.txt";
	ossimFilename virtual_chkfile = workfold + "\\virtual_chk.txt";
	ossimFilename reportfile = workfold + "\\report.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	//prj.SavePointToFile(gcp_hgt_File, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);


	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	if(!handler) return;   //应该弹出警告对话框
	ossimTieGptSet* gptSet = new ossimTieGptSet;
	double max_height = 1100.0;
	double min_height = 0.0;
	int nLevels = 5;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height + i*(max_height-min_height)/(nLevels-1);
		create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, gptSet, 15, 15, true, false);
	}

	ossimTieGptSet* chkSet = new ossimTieGptSet;
	nLevels = 1;
	for (int i=0;i < nLevels;++i)
	{
		double hgt = min_height;// + i*(max_height-min_height)/(nLevels-1);
		create3DGridPoints(handler->getImageRectangle(), *prj.m_sensorModel, hgt, chkSet, 12000, 1, true, false);
	}
	handler->close();
	delete handler;

	prj.SavePointToFile(virtual_gcpfile, gptSet, NULL);
	prj.SavePointToFile(virtual_chkfile, chkSet, NULL);
}


//void QbRpc1()
//{
//	ossimFilename workfold = "E:\\testdata\\QB\\Src\\052407462010_01_P003_PAN";
//	ossimFilename sourcefile = workfold + "\\10SEP14023541-P2AS-052407462010_01_P003.TIL";
//	ossimFilename elevationpath = "D:\\workspace\\dem_qb";
//	ossimFilename gcpfile = workfold + "\\all_30.txt";
//	ossimFilename reportfile = workfold + "\\RpcReport.txt";
//	QuickbirdRpcModel *qbRpcModel = new QuickbirdRpcModel(sourcefile);
//	qbRpcModel->ReadGcpAndProjection(gcpfile, qbRpcModel->m_ctrSet, qbRpcModel->m_chkSet);
//	qbRpcModel->m_sensorModel->m_proj = qbRpcModel->m_MapPar;
//	qbRpcModel->m_theElevManager = ossimElevManager::instance();
//
//	if(!qbRpcModel->m_theElevManager->loadElevationPath(ossimFilename(elevationpath)))
//	{
//		cout<<"warning: 加载DEM失败！"<<endl;
//		return;
//	}
//	qbRpcModel->GetElevations(qbRpcModel->m_ctrSet);
//	qbRpcModel->GetElevations(qbRpcModel->m_chkSet);
//
//	ossimImageGeometry *imageGeom;
//	fstream rpcStruct;
//
//	int num = static_cast<int>(qbRpcModel->m_ctrSet->getTiePoints().size());
//	vector < ossimDpt > imagePoints(num);
//	vector < ossimGpt > groundControlPoints(num);
//
//	for(int i = 0;i < num;i++)
//	{
//		groundControlPoints[i] = qbRpcModel->m_ctrSet->getTiePoints()[i]->getGroundPoint();
//		ossimGpt gpt = qbRpcModel->m_MapPar->inverse(ossimDpt(groundControlPoints[i].lat, groundControlPoints[i].lon));
//		gpt.hgt = groundControlPoints[i].hgt;
//		groundControlPoints[i] = gpt;
//		imagePoints[i] = qbRpcModel->m_ctrSet->getTiePoints()[i]->getImagePoint();
//	}
//	ossimRpcSolver *solver = new ossimRpcSolver(true, false);
//	solver->solveCoefficients(imagePoints, groundControlPoints);
//	imageGeom = solver->createRpcModel();
//	ossimRpcModel* rpcModel = (ossimRpcModel*)imageGeom->getProjection();
//	qbRpcModel->m_sensorModel = rpcModel;
//	qbRpcModel->m_sensorModel->m_proj = qbRpcModel->m_MapPar;
//	qbRpcModel->m_sensorModel->saveState(qbRpcModel->m_geom);
//	rpcModel->loadState(qbRpcModel->m_geom);
//	rpcStruct.open((workfold+"\\rpcStruct.txt").c_str(), ios_base::out);
//	rpcModel->saveRpcModelStruct(rpcStruct);
//	rpcStruct.close();
//	//qbRpcModel->OutputReport(workfold+"\\report.txt", qbRpcModel->m_sensorModel, qbRpcModel->m_ctrSet, qbRpcModel->m_chkSet);
//	MyProject prj;
//	prj.OutputReport(workfold+"\\report.txt", qbRpcModel->m_sensorModel, qbRpcModel->m_ctrSet, qbRpcModel->m_chkSet, false);
//
//
//	for(int i = 0;i < num;i++)
//	{
//		//groundControlPoints[i] = qbRpcModel->m_ctrSet->getTiePoints()[i]->getGroundPoint();
//		//ossimGpt gpt = qbRpcModel->m_MapPar->inverse(ossimDpt(groundControlPoints[i].lat, groundControlPoints[i].lon),groundControlPoints[i]);
//		//gpt.hgt = groundControlPoints[i].hgt;
//		//groundControlPoints[i] = gpt;
//		imagePoints[i] = qbRpcModel->m_ctrSet->getTiePoints()[i]->getImagePoint();
//		groundControlPoints[i] = qbRpcModel->m_ctrSet->getTiePoints()[i]->getGroundPoint();
//	}
//
//	ossimRpcSolver *solver_Xyz2Rc = new ossimRpcSolver(true, false);
//	solver_Xyz2Rc->solveCoefficients(imagePoints, groundControlPoints);
//	imageGeom = solver_Xyz2Rc->createRpcXyz2RcModel(qbRpcModel->m_MapPar);
//	ossimRpcXyz2RcModel* rpcModel_Xyz2Rc = (ossimRpcXyz2RcModel*)imageGeom->getProjection();
//	qbRpcModel->m_sensorModel = rpcModel_Xyz2Rc;
//	qbRpcModel->m_sensorModel->m_proj = qbRpcModel->m_MapPar;
//	qbRpcModel->m_sensorModel->saveState(qbRpcModel->m_geom);
//	rpcModel_Xyz2Rc->loadState(qbRpcModel->m_geom);
//	rpcStruct.open((workfold+"\\rpcStruct_Xyz2Rc.txt").c_str(), ios_base::out);
//	rpcModel_Xyz2Rc->saveRpcModelStruct(rpcStruct);
//	rpcStruct.close();
//	//qbRpcModel->OutputReport(workfold+"\\report_Xyz2Rc.txt", qbRpcModel->m_sensorModel, qbRpcModel->m_ctrSet, qbRpcModel->m_chkSet);
//	prj.OutputReport(workfold+"\\report_Xyz2Rc.txt", qbRpcModel->m_sensorModel, qbRpcModel->m_ctrSet, qbRpcModel->m_chkSet, true);
//
//
//
//	ossimRpcSolver *solver_Rcz2Xy = new ossimRpcSolver(true, false);
//	//solver_Rcz2Xy->solveCoefficients_Rcz2Xy(imagePoints, groundControlPoints);
//	imageGeom = solver_Rcz2Xy->createRpcRcz2XyModel(qbRpcModel->m_MapPar);
//	ossimRpcRcz2XyModel* rpcModel_Rcz2Xy = (ossimRpcRcz2XyModel*)imageGeom->getProjection();
//	qbRpcModel->m_sensorModel = rpcModel_Rcz2Xy;
//	qbRpcModel->m_sensorModel->m_proj = qbRpcModel->m_MapPar;
//	qbRpcModel->m_sensorModel->saveState(qbRpcModel->m_geom);
//	rpcModel_Rcz2Xy->loadState(qbRpcModel->m_geom);
//	rpcStruct.open((workfold+"\\rpcStruct_Rcz2Xy.txt").c_str(), ios_base::out);
//	rpcModel_Rcz2Xy->saveRpcModelStruct(rpcStruct);
//	rpcStruct.close();
//	qbRpcModel->OutputReport(workfold+"\\report_Rcz2Xy.txt", qbRpcModel->m_sensorModel, qbRpcModel->m_ctrSet, qbRpcModel->m_chkSet);
//
//
//	return;
//}

void QbRpc2()
{
	//ossimFilename workfold = "H:\\testdata\\QB\\Src\\052407462010_01_P003_PAN";
	//ossimFilename sourcefile = workfold + "\\10SEP14023541-P2AS-052407462010_01_P003.TIL";
	//ossimFilename elevationpath = "D:\\workspace\\dem_qb";
	ossimFilename workfold = "G:\\project_qb\\source\\5\\2\\052407492010_02_P001";
	ossimFilename sourcefile = workfold + "\\10OCT29043027-P2AS-052407492010_02_P001.TIL";
	ossimFilename elevationpath = "G:\\QB\\DEM";
	ossimFilename gcpfile = workfold + "\\gcp.txt";
	ossimFilename reportfile = workfold + "\\RpcReport.txt";
	QuickbirdRpcModel *qbRpcModel = new QuickbirdRpcModel(sourcefile);
	qbRpcModel->ReadGcpAndProjection(gcpfile, qbRpcModel->m_ctrSet, qbRpcModel->m_chkSet);
	qbRpcModel->m_sensorModel->m_proj = qbRpcModel->m_MapPar;
	qbRpcModel->m_theElevManager = ossimElevManager::instance();

	if(!qbRpcModel->m_theElevManager->loadElevationPath(ossimFilename(elevationpath)))
	{
		cout<<"warning: 加载DEM失败！"<<endl;
		return;
	}
	qbRpcModel->GetElevations(qbRpcModel->m_ctrSet);
	qbRpcModel->GetElevations(qbRpcModel->m_chkSet);
	qbRpcModel->UpdateSensorModel(*qbRpcModel->m_ctrSet, qbRpcModel->m_sensorModel, qbRpcModel->m_geom);
	qbRpcModel->OutputReport(workfold+"\\report0.txt", qbRpcModel->m_sensorModel, qbRpcModel->m_ctrSet, qbRpcModel->m_chkSet);

	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	if(!handler) return;   //应该弹出警告对话框
	ossimTieGptSet* gptSet = createTieGptSet(handler->getImageRectangle(),*qbRpcModel->m_sensorModel, 10, 10);
	ossimTieGptSet* chkSet = createTieGptSet(handler->getImageRectangle(),*qbRpcModel->m_sensorModel, 20, 20);
	handler->close();
	delete handler;

	//int num = static_cast<int>(gptSet->getTiePoints().size());
	//ossimRpcSolver *solver = new ossimRpcSolver(true, false);
	//vector < ossimDpt > imagePoints(num);
	//vector < ossimGpt > groundControlPoints(num);
	//for(int i = 0;i < num;i++)
	//{
	//	//groundControlPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
	//	//ossimGpt gpt = qbRpcModel->m_MapPar->inverse(ossimDpt(groundControlPoints[i].lat, groundControlPoints[i].lon),groundControlPoints[i]);
	//	//gpt.hgt = groundControlPoints[i].hgt;
	//	//groundControlPoints[i] = gpt;
	//	imagePoints[i] = gptSet->getTiePoints()[i]->getImagePoint();
	//	groundControlPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
	//}

	//solver->solveCoefficients(imagePoints, groundControlPoints);
	//ossimImageGeometry *imageGeom = solver->createRpcXyz2RcModel(qbRpcModel->m_MapPar);
	//ossimRpcXyz2RcModel* rpcModel = (ossimRpcXyz2RcModel*)imageGeom->getProjection();

	//rpcModel->m_optimizeType = ossimRpcXyz2RcModel::imageAffine;
	//qbRpcModel->m_sensorModel = rpcModel;
	//qbRpcModel->m_sensorModel->m_proj = qbRpcModel->m_MapPar;
	fstream rpcStruct;
	qbRpcModel->m_sensorModel->saveState(qbRpcModel->m_geom);
	ossimRpcModel *rpcModel = new ossimRpcModel;
	rpcModel->loadState(qbRpcModel->m_geom);
	rpcStruct.open((workfold+"\\rpcStruct0.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	qbRpcModel->SavePointToFile(workfold+"\\virtualGpt.txt", gptSet, chkSet);
	MyProject prj;
	prj.OutputReport(workfold+"\\report1.txt", qbRpcModel->m_sensorModel, gptSet, chkSet, true);
	//qbRpcModel->OutputReport(workfold+"\\report1.txt", qbRpcModel->m_sensorModel, gptSet, chkSet, true);
	return;
}


void Hj1Rpc()
{
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossim_plugin.dll");
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimgdal_plugin.dll");
	ossimInit::instance()->initialize();
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401050432_201401050443";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_SYC_201401050300_201401050309";
	ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401070348_201401070354";
	//ossimFilename workfold = "E:\\beijing";
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\IMAGE.TIF";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";

	//ossimFilename demPath = "I:\\spot5\\beijing";
	ossimFilename demPath = "D:\\workspace\\dem-hgt";
	//ossimFilename demPath = "E:\\dem_bj";
	//ossimFilename chkFile = workfold + "\\281268_2m5_gcp09033_hgt.txt";
	ossimFilename gcpfile = workfold + "\\gcp.txt";
	//ossimFilename chkFile = workfold + "\\gcp_all_hgt.txt";

	ossimFilename gcp_hgt_File = workfold + "\\gcp_hgt.txt";
	ossimFilename outFile = workfold + "\\rect_rpc.tif";
	ossimFilename reportFile = workfold + "\\report_orbit.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));

	prj.theMgr = ossimElevManager::instance();
	//prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.SavePointToFile(gcp_hgt_File, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.InitiateSensorModel();
	prj.m_sensorModel->setImageRect(ossimDrect(ossimDpt(0.0,0.0), ossimDpt(11999, 11999)));
	prj.m_sensorModel->m_proj = prj.m_MapPar;


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
	ossimRpcModel* rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom, "projection.");

	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	//// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;fstream rpcStruct;
	fstream rpcStruct;
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.m_sensorModel->saveState(prj.geom);
	rpcModel->loadState(prj.geom);
	rpcStruct.open((workfold+"\\rpcStruct0.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	prj.OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, true);


	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);
	prj.Orthograph(outFile);
}


void Rcp_Generic()
{
	//ossimFilename workfold = "D:\\workspace\\ZY\\ZY3_01a_hsnfavm_001126_20121011_111248_0007_SASMAC_CHN_sec_rel_004_1210127866";
	//ossimFilename imageFile = workfold + "\\ZY3_01a_hsnfavm_001126_20121011_111248_0007_SASMAC_CHN_sec_rel_004_1210127866.tif";
	//ossimFilename rpcFile = workfold + "\\ZY3_01a_hsnfavm_001126_20121011_111248_0007_SASMAC_CHN_sec_rel_004_1210127866_rpc.txt";
	ossimFilename workfold = "I:\\testdata\\zy\\ZY3_01a_hsnnavp_006179_20120323_112052_0007_SASMAC_CHN_sec_rel_001_1203303678";
	ossimFilename imageFile = workfold + "\\ZY3_01a_hsnnavp_006179_20120323_112052_0007_SASMAC_CHN_sec_rel_001_1203303678.tif";
	ossimFilename rpcFile = workfold + "\\ZY3_01a_hsnnavp_006179_20120323_112052_0007_SASMAC_CHN_sec_rel_001_1203303678_rpc.txt";
	ossimFilename outFile = workfold + "\\rect.tif";
	//ossimFilename elevationpath = "D:\\workspace\\ossimdem";
	ossimFilename elevationpath = "D:\\workspace\\dem-hgt";
	//ossimFilename elevationpath = "I:\\zy\\ASTGTM_N18E109V";

	ossimFilename gcpfile = workfold + "\\features\\virtualGcps.txt";
	ossimFilename reportfile = workfold + "\\reports\\RpcReport.txt";

	MyProject prj;
	prj.theMgr = ossimElevManager::instance();
	if(!prj.theMgr->loadElevationPath(ossimFilename(elevationpath)))
	{
		cout<<"warning: 加载DEM失败！"<<endl;
		return;
	}
	prj.ReadGcpAndProjection(gcpfile);
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(elevationpath));
	prj.m_DemPath = ossimFilename(elevationpath);
	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimRpcModel::rpcModelStruct rpcStruct;
	readRPCFile(rpcFile, rpcStruct);

	//ossimIrect bounds = prhandler->getBoundingRect();
	//ossimGpt centerGround;
	//inputProjection->lineSampleToWorld(bounds.midPoint(), centerGround);
	//ossimProjection* newUtmView(const ossimGpt& centerGround, const ossimDpt& metersPerPixel);

	ossimKeywordlist kwl;
	ProjectionParameters proParam;
	proParam.DatumName = "WGS84";
	proParam.TrueOriginLongitude = 111.0;
	proParam.TrueOriginLatitude = 0.0;
	proParam.EastingFalse = 500000;
	proParam.NorthingFalse = 0;
	//proParam.PixelSizeX = 2.1;
	//proParam.PixelSizeY = 2.1;
	//proParam.PixelSizeX = 5.831721842186000870;
	//proParam.PixelSizeY = 5.831721842186000870;
	proParam.PixelSizeX = 2.100596582524536071;
	proParam.PixelSizeY = 2.100596582524536071;
	ossimMapProjection* projection = CreateProjection(proParam, kwl);
	prj.m_MapProjection = kwl;
	prj.m_MapPar = projection;

	prj.m_ImgFileNameUnc = imageFile;
	prj.m_OutBandList.clear();
	//选择输出波段
	int nbands = 4;
	for (int i = 0;i < nbands;++i)
	{
		prj.m_OutBandList.push_back(i);
	}

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

	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	prj.Orthograph(outFile);
}



void SpotRpcTest1()
{
	//ossimFilename workfold = "E:\\beijing";
	ossimFilename workfold = "D:\\workspace\\SPOT5\\beijing";
	ossimFilename sourceFile = workfold + "\\281268-20040602-2.5\\scene01\\imagery.tif";
	//ossimFilename demPath = "E:\\dem_bj";
	//ossimFilename demPath = "D:\\workspace\\dem";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";
	ossimFilename demPath = "D:\\workspace\\dem\\aster30";
	ossimFilename gcpfile = workfold + "\\281268_2m5_gcp09033_hgt.txt";
	//ossimFilename gcpfile = workfold + "\\gcp_80_12.txt";
	//ossimFilename chkFile = workfold + "\\gcp.txt";
	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile0 = workfold + "\\report0.txt";
	ossimFilename reportFile1 = workfold + "\\report1.txt";



	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));//
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);
	prj.InitiateSensorModel(sourceFile);

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
	rpcStruct.open((workfold+"\\rpcStruct0.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	//// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;
	prj.OutputReport(reportFile0, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);

	//double tmp;
	//for(int i = 0;i < num;i++)
	//{
	//	groundControlPoints.push_back(prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint());
	//	imagePoints.push_back(prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint());
	//}
	//prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, &tmp);		//用优化控制点进行模型优化
	//prj.m_sensorModel->updateModel();
	//prj.m_sensorModel->saveState(prj.geom);
	//prj.OutputReport(reportfileall, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	//for(int i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	//{
	//	ossimDpt dpt = prj.m_MapPar->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
	//	ossimGpt gpt(dpt.x,dpt.y);
	//	prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(ossimGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	//}
	prj.m_sensorModel->m_bHgtOptimize = false;

	//// 全参数优化
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.m_sensorModel->saveState(prj.geom);
	//rpcModel->loadState(prj.geom);

	//rpcStruct.open((workfold+"\\rpcStruct1.txt").c_str(), ios_base::out);
	//rpcModel->saveRpcModelStruct(rpcStruct);
	//rpcStruct.close();
	//prj.OutputReport(reportFile1, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);
	// 像方仿射变换
	prj.Orthograph(outFile);

}

void SpotRpcTest2()
{
	//ossimFilename workfold = "E:\\beijing";
	ossimFilename workfold = "D:\\workspace\\testdata\\SPOT5\\beijing";
	ossimFilename sourceFile = workfold + "\\281268-20040602-2.5\\scene01\\imagery.tif";
	//ossimFilename demPath = "E:\\dem_bj";
	ossimFilename demPath = "D:\\workspace\\dem_bj";
	ossimFilename gcpfile = workfold + "\\gcp_80_12.txt";
	//ossimFilename gcpfile = workfold + "\\281268_2m5_gcp09033_hgt.txt";
	//ossimFilename gcpfile = workfold + "\\gcp.txt";
	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile0 = workfold + "\\report0.txt";
	ossimFilename reportFile1 = workfold + "\\report1.txt";



	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));//
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);
	prj.InitiateSensorModel(sourceFile);
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);

	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourceFile);
	if(!handler) return;   //应该弹出警告对话框
	ossimTieGptSet* gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel);
	prj.SavePointToFile(workfold+"\\virtualGpt.txt", gptSet, gptSet);
	handler->close();
	delete handler;

	int num = static_cast<int>(gptSet->getTiePoints().size());
	ossimRpcSolver *solver = new ossimRpcSolver(true, false);
	vector < ossimDpt > imagePoints(num);
	vector < ossimGpt > groundControlPoints(num);
	for(int i = 0;i < num;i++)
	{
		groundControlPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
		ossimGpt gpt = prj.m_MapPar->inverse(ossimDpt(groundControlPoints[i].lat, groundControlPoints[i].lon));
		gpt.hgt = groundControlPoints[i].hgt;
		groundControlPoints[i] = gpt;
		imagePoints[i] = gptSet->getTiePoints()[i]->getImagePoint();
	}

	solver->solveCoefficients(imagePoints, groundControlPoints);
	ossimImageGeometry *imageGeom = solver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimRpcModel* rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom);

	fstream rpcStruct;
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.m_sensorModel->saveState(prj.geom);
	rpcModel->loadState(prj.geom);
	rpcStruct.open((workfold+"\\rpcStruct0.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	prj.OutputReport1(reportFile0, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	prj.OutputReport(reportFile0, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);

	//// 优化
	//rpcModel->m_optimizeType = ossimRpcModel::allParam;
	//prj.m_sensorModel = rpcModel;
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//prj.UpdateSensorModel(*gptSet, prj.m_sensorModel, prj.geom);
	//prj.m_sensorModel->saveState(prj.geom);
	//rpcModel->loadState(prj.geom);

	// 像方仿射变换
	//rpcModel->m_optimizeType = ossimRpcModel::imageAffine;
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);

	prj.m_sensorModel->saveState(prj.geom);
	rpcModel->loadState(prj.geom);
	rpcStruct.open((workfold+"\\rpcStruct1.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	prj.OutputReport1(reportFile0, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	prj.OutputReport(reportFile1, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);

}

//void AlosRpcTest()
//{
//	ossimFilename workfold = "H:\\testdata\\AlosData\\d1001111-001";
//	ossimFilename sourcefile = workfold + "\\IMG-ALPSMW220853075_O1B2R_UW.tif";
//	ossimFilename elevationpath = "D:\\workspace\\dem";
//	ossimFilename gcpfile = workfold + "\\features\\virtualGcps.txt";
//	ossimFilename reportfile = workfold + "\\reports\\RpcReport.txt";
//
//	MyProject prj;
//	AlosPRISM alosUti(workfold);
//	if(!alosUti.getInitState())
//	{
//		cout<<"warning: Alos PRISM数据\""<<workfold<<"\"初始化失败！"<<endl;
//		return;
//	}
//	prj.theMgr = ossimElevManager::instance();
//	if(!prj.theMgr->loadElevationPath(ossimFilename(elevationpath)))
//	{
//		cout<<"warning: 加载DEM失败！"<<endl;
//		return;
//	}
//	prj.ReadGcpAndProjection(gcpfile);
//	prj.theMgr = ossimElevManager::instance();
//	prj.theMgr->loadElevationPath(ossimFilename(elevationpath));
//	prj.m_DemPath = ossimFilename(elevationpath);
//	prj.GetElevations(prj.m_CtrlGptSet);
//	prj.GetElevations(prj.m_ChkGptSet);
//
//	prj.m_MapProjection = alosUti.getKeywordlist();
//	prj.m_MapPar = alosUti.getMapProjection();
//
//	prj.m_ImgFileNameUnc = ossimFilename(alosUti.m_FileTIF);
//	prj.m_OutBandList.clear();
//	prj.m_OutBandList.push_back(0);
//
//	prj.m_ModelType = RPCType;
//	if(!prj.InitiateSensorModel(prj.m_ImgFileNameUnc))
//	{
//		cout<<"warning: 未找到适合\""<<prj.m_ImgFileNameUnc<<"\"的模型！"<<endl;
//		return;
//	}
//	ossimRpcModel *rpcModel = new ossimRpcModel;
//	rpcModel->setAttributes(alosUti.getRpcModelStruct());
//	prj.m_sensorModel = rpcModel;
//	prj.m_sensorModel->m_proj = prj.m_MapPar;
//
//	/*ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
//	if(!handler) return;   //应该弹出警告对话框
//	prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
//	prj.m_sensorModel->m_proj = prj.m_MapPar;
//	handler->close();
//	delete handler;*/
//
//	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
//	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
//}

void LandsatRpcTest1()
{
	ossimFilename workfold = "D:\\workspace\\testdata\\Landsat\\LD2010003816";
	//ossimFilename workfold = "D:\\testdata\\LD2010003816";
	ossimFilename sourcefile = workfold + "\\header.dat";
	ossimFilename elevationpath = "D:\\workspace\\dem";
	ossimFilename gcpfile = workfold + "\\gcp60+60.txt";
	ossimFilename reportfile = workfold + "\\reports\\report.txt";
	ossimFilename reportfileall = workfold + "\\reports\\reportall.txt";
	ossimFilename outfile = workfold + "\\rect_points.tif";

	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(elevationpath));//
	prj.m_DemPath=ossimFilename(elevationpath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(6);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);
	prj.InitiateSensorModel(sourcefile);

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
	ossimRpcModel* rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom, "projection.");

	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	//// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;fstream rpcStruct;
	fstream rpcStruct;
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.m_sensorModel->saveState(prj.geom);
	rpcModel->loadState(prj.geom);
	rpcStruct.open((workfold+"\\rpcStruct0.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);

	//double tmp;
	//for(int i = 0;i < num;i++)
	//{
	//	groundControlPoints.push_back(prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint());
	//	imagePoints.push_back(prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint());
	//}
	//prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, &tmp);		//用优化控制点进行模型优化
	//prj.m_sensorModel->updateModel();
	//prj.m_sensorModel->saveState(prj.geom);
	//prj.OutputReport(reportfileall, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	//for(int i = 0;i < static_cast<int>(prj.m_CtrlGptSet->getTiePoints().size());i++)
	//{
	//	ossimDpt dpt = prj.m_MapPar->forward(*prj.m_CtrlGptSet->getTiePoints()[i]);
	//	ossimGpt gpt(dpt.x,dpt.y);
	//	prj.m_CtrlGptSet->refTiePoints()[i]->setGroundPoint(ossimGpt(dpt.x,dpt.y,prj.m_CtrlGptSet->getTiePoints()[i]->hgt));
	//}

	// 全参数优化
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.m_sensorModel->saveState(prj.geom);
	rpcModel->loadState(prj.geom);
	// 像方仿射变换
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.OutputReport(reportfileall, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, true);
	//prj.Orthograph(outfile);

}

void LandsatRpcTest2()
{
	ossimFilename workfold = "D:\\workspace\\testdata\\BlockAdjustment\\LD2010003650";
	//ossimFilename workfold = "D:\\workspace\\testdata\\Landsat\\LD2010003816";
	//ossimFilename workfold = "D:\\testdata\\LD2010003816";
	ossimFilename sourcefile = workfold + "\\header.dat";
	ossimFilename elevationpath = "D:\\workspace\\dem";
	ossimFilename tmProjection = "D:\\workspace\\testdata\\Landsat\\TM_Projection.txt";
	ossimFilename gcpfile = workfold + "\\gcp70.txt";
	ossimFilename reportfile = workfold + "\\reports\\report.txt";
	ossimFilename reportfile1 = workfold + "\\reports\\report1.txt";
	ossimFilename outfile = workfold + "\\rect_points.tif";

	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	prj.ReadGcpAndProjection(ossimFilename(gcpfile));
	if(!prj.m_MapPar)
	{
		prj.CreateL5Projection(sourcefile);
	}
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(elevationpath));//
	prj.m_DemPath=ossimFilename(elevationpath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(6);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);
	prj.InitiateSensorModel(sourcefile);
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);

	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	if(!handler) return;   //应该弹出警告对话框
	ossimTieGptSet* gptSet = createTieGptSet(handler->getImageRectangle(),*prj.m_sensorModel);
	handler->close();
	delete handler;

	int num = static_cast<int>(gptSet->getTiePoints().size());
	ossimRpcSolver *solver = new ossimRpcSolver(true, false);
	vector < ossimDpt > imagePoints(num);
	vector < ossimGpt > groundControlPoints(num);
	for(int i = 0;i < num;i++)
	{
		groundControlPoints[i] = gptSet->getTiePoints()[i]->getGroundPoint();
		ossimGpt gpt = prj.m_MapPar->inverse(ossimDpt(groundControlPoints[i].lat, groundControlPoints[i].lon));
		gpt.hgt = groundControlPoints[i].hgt;
		groundControlPoints[i] = gpt;
		imagePoints[i] = gptSet->getTiePoints()[i]->getImagePoint();
	}

	solver->solveCoefficients(imagePoints, groundControlPoints);

	ossimImageGeometry *imageGeom = solver->createRpcModel();
	ossimKeywordlist geom;
	imageGeom->saveState(geom);
	ossimRpcModel* rpcModel = new ossimRpcModel();
	rpcModel->loadState(geom);


	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	fstream rpcStruct;
	prj.m_sensorModel->saveState(prj.geom);
	rpcModel->loadState(prj.geom);
	rpcStruct.open((workfold+"\\rpcStruct0.txt").c_str(), ios_base::out);
	rpcModel->saveRpcModelStruct(rpcStruct);
	rpcStruct.close();
	prj.SavePointToFile(workfold+"\\virtualGpt.txt", gptSet, prj.m_ChkGptSet);
	prj.OutputReport(workfold+"\\report0.txt", prj.m_sensorModel, gptSet, prj.m_ChkGptSet, false);
	return;

	// 像方仿射变换
	prj.m_sensorModel = rpcModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;
	prj.m_sensorModel->m_bHgtOptimize = false;
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet, false);
	//prj.Orthograph(outfile);

}