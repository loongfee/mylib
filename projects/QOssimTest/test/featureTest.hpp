#include <iostream>
#include <iterator>
#include <fstream>
using namespace std;

#include <func.h>
#include <mprojectdefine.h>
#include "..\QuickbirdRpcModel.h"
//#include "..\AlosBatch.h"


#include <ogrsf_frmts.h>
#include <gdal.h>
#include <ogr_api.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>


#include <ossim_plugin\ossimEnvisatAsarModel.h>
#include <ossim_plugin\ossimRadarSatModel.h>
#include <ossim_plugin\ossimGeometricSarSensorModel.h>
#include <ossim_plugin\ossimRadarSat2RPCModel.h>
using namespace ossimplugins;
using namespace mylib;

template<class T>
bool shape2points(string shpFile, vector<T>& gptList)
{
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动

	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");//得到shp文件的处理器

	OGRDataSource* poDS = poDriver->Open( shpFile.c_str(), NULL );//打开文件

	//int iLayerCount = poDS->GetLayerCount();

	OGRLayer* poLayer = poDS->GetLayer(0);//获取shp图层

	//读取几何和属性值
	OGRFeature *poFeature;
	poLayer->ResetReading();//确保是从该层的开头开始
	gptList.clear();
	while ((poFeature=poLayer->GetNextFeature())!=NULL)
	{
		OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
		int iField;
		for( iField = 0; iField < poFDefn->GetFieldCount(); iField++ )
		{
			OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
			printf( "%s=", poFieldDefn->GetNameRef() );
			if( poFieldDefn->GetType() == OFTInteger )
			{
				printf( "%d\n", poFeature->GetFieldAsInteger( iField ) );
			}
			else if( poFieldDefn->GetType() == OFTReal )
			{
				printf( "%.3f\n", poFeature->GetFieldAsDouble(iField) );
			}
			else if( poFieldDefn->GetType() == OFTString )
			{
				printf( "%s\n", poFeature->GetFieldAsString(iField) );
			}
			else
			{
				printf( "%s\n", poFeature->GetFieldAsString(iField) );
			}
		}

		OGRGeometry *poGeometry;
		poGeometry = poFeature->GetGeometryRef();
		OGRwkbGeometryType geoType = poGeometry->getGeometryType();
		if( poGeometry != NULL 
			&& wkbFlatten(poGeometry->getGeometryType()) == wkbPoint )//点文件
		{
			OGRPoint *poPoint = (OGRPoint *) poGeometry;
			printf( "%.3f,%.3f\n", poPoint->getX(), poPoint->getY() );
			T gpt(poPoint->getX(), poPoint->getY());
			gptList.push_back(gpt);
		}
		else
		{
			printf( "Error: the shape file is not a point layer.\n" );
			return false;
		}
		printf("\n");
	}
	OGRDataSource::DestroyDataSource( poDS );
	OGRCleanupAll();//资源清理
	return true;
}


#include "GdalRasterApp.h"
ossimTieGptSet* shape2gcps(string srcShpFile, string refShpFile)
{
	vector<ossimDpt> dptList;
	vector<ossimGpt> gptList;
	shape2points<ossimDpt>(srcShpFile, dptList);
	shape2points<ossimGpt>(refShpFile, gptList);
	
	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	if (dptList.size() != gptList.size())
	{
		return gcpSet;
	}

	string imageFile = "D:\\workspace\\jiading\\radarsat\\imagery_HH_LS_8_Lee.tif";
	//string imageFile = "D:\\workspace\\Landsat\\features\\edge\\source1.tif";
	GdalRasterApp gdalApp;
	gdalApp.open(imageFile.c_str());
	int Height = gdalApp.getDataset()->GetRasterYSize();
	for (int i = 0; i < (int)dptList.size(); i++)
	{
		dptList[i].y = Height - dptList[i].y;
	}

	for (int i = 0; i < (int)dptList.size(); i++)
	{
		ossimTieGpt *aTiePt = new ossimTieGpt(gptList[i], dptList[i], 0);
		char strId[64];
		sprintf_s(strId, "%d", i+1);
		aTiePt->GcpNumberID = ossimString(strId);
		gcpSet->addTiePoint(aTiePt);
	}
	return gcpSet;
}

void feature_radarsat2()
{	
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossim_plugin.dll");
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimgdal_plugin.dll");
	ossimInit::instance()->initialize();

	ossimFilename workfold = "D:\\workspace\\Radarsat2";
	ossimFilename sourceFile = workfold + "\\imagery_HH.tif";
	//ossimFilename sourceFile = workfold + "\\imagery_HH_LS_8.tif";
	ossimFilename demPath = "D:\\workspace\\dem-hgt";
	ossimFilename gcpfile = workfold + "\\gcps.txt";
	ossimFilename projection_file = workfold + "\\projection.txt";
	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile0 = workfold + "\\report.txt";
	
	ossimFilename image_point_shape = "D:\\workspace\\jiading\\features\\point_source.shp";
	ossimFilename ground_point_shape = "D:\\workspace\\jiading\\features\\point_reference.shp";
	
	ossimFilename srcFeatureFile = "D:\\workspace\\jiading\\features\\source.txt";
	ossimFilename refFeatureFile = "D:\\workspace\\jiading\\features\\reference.txt";
	ossimFilename reportfile = workfold + "\\reports\\report.txt";
	ossimFilename featureReport = workfold + "\\reports\\feature_report.txt";
	ossimFilename outfile = workfold + "\\rect_points.tif";

	ossimFilename point_shape_input = "D:\\workspace\\jiading\\shapes\\points_input.shp";
	ossimFilename line_shape_input = "D:\\workspace\\jiading\\shapes\\lines_input.shp";
	ossimFilename area_shape_input = "D:\\workspace\\jiading\\shapes\\areas_input.shp";
	ossimFilename point_shape_output = "D:\\workspace\\jiading\\shapes\\points_output.shp";
	ossimFilename line_shape_output = "D:\\workspace\\jiading\\shapes\\lines_output.shp";
	ossimFilename area_shape_output = "D:\\workspace\\jiading\\shapes\\areas_output.shp";

	// 加载高程数据
	ossimElevManager::instance()->loadElevationPath(demPath);

	// 创建RadarSat2 RPC模型
	ossimRadarSat2RPCModel* pRadarSat2RPCModel = new ossimRadarSat2RPCModel();
	pRadarSat2RPCModel->open(sourceFile);
	ossimRpcModel* pRpcModel = new ossimRpcModel;
	ossimRpcModel::rpcModelStruct rpcStruct;
	pRadarSat2RPCModel->getRpcParameters(rpcStruct);
	pRpcModel->setAttributes(rpcStruct);
		

	MyProject prj;
	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(sourceFile);
	prj.m_sensorModel = pRpcModel;
	// 指定输出投影
	prj.m_MapProjection.addFile(projection_file);
	ossimRefPtr<ossimProjection> proj;
	proj = ossimProjectionFactoryRegistry::instance()->createProjection(prj.m_MapProjection);
	prj.m_MapPar = PTR_CAST(ossimMapProjection, proj.get());
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	prj.m_CtrlGptSet = shape2gcps(image_point_shape, ground_point_shape);
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	//prj.ReadGcpAndProjection(ossimFilename(gcpfile));
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));//
	prj.m_DemPath=ossimFilename(demPath);


	// 分别读取点线面特征

	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDUnknown);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(srcFeatureFile.c_str(), imageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(refFeatureFile.c_str(), groundFeatureList, ossimFeatureType::ossimUnknown);
	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDPointType);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGPointType);
	//prj.ReadImagePoints(srcFeatureFile, imageFeatureList);
	//prj.ReadGroundPoints(refFeatureFile, groundFeatureList);
	//prj.ReadImageLines(imageLineFile, imageFeatureList);
	//prj.ReadGroundLines(groundLineFile, groundFeatureList);
	//prj.ReadImageAreas(imageAreaFile, imageFeatureList);
	//prj.ReadGroundAreas(groundAreaFile, groundFeatureList);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);
	
	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	if( NULL == prj.m_sensorModel )
		return;
	vector<ossimTieFeature> tieFeatureList;
	//根据序号连接对应的特征
	AppendTieFeatures(tieFeatureList, imageFeatureList, groundFeatureList);

	// 加入控制点
	//AppendTieFeaturesFromTieGptSet(tieFeatureList, prj.m_CtrlGptSet);

	//AppendTieAreas(tieFeatureList, dptAreaList, gptAreaList);
	prj.GetElevations(tieFeatureList);
	

	prj.UpdateSensorModel(tieFeatureList, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(tieLineList, prj.m_sensorModel, prj.geom);

	prj.m_sensorModel->loadState(prj.geom);


	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	FeatureReport(featureReport, prj.m_sensorModel, tieFeatureList);
	
	// 加入控制点
	AppendTieFeaturesFromTieGptSet(tieFeatureList, prj.m_CtrlGptSet);

	prj.SaveFeaturetoShape(point_shape_input, tieFeatureList, ossimFeatureType::ossimPointType);
	prj.SaveFeaturetoShape(line_shape_input, tieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(area_shape_input, tieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.UpdateFeatures(prj.m_sensorModel, tieFeatureList);
	prj.SaveFeaturetoShape(area_shape_output, tieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.SaveFeaturetoShape(line_shape_output, tieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(point_shape_output, tieFeatureList, ossimFeatureType::ossimPointType);

	prj.Orthograph(outfile);
}

void feature_landsat()
{
	string workfold = "D:\\workspace\\Landsat\\";
	workfold = "D:\\workspace\\Landsat\\xinjiang";
	//ossimFilename workfold = "C:\\LD2010003816";
	string sourcefile = workfold + "\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat";
	sourcefile = workfold + "\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat";
	//ossimFilename elevationpath = "D:\\workspace\\Landsat\\dem_xj";
	string elevationpath = "D:\\workspace\\dem\\aster30";

	string srcFeatureFile = workfold + "\\features\\source.txt";
	string refFeatureFile = workfold + "\\features\\reference.txt";
	//srcFeatureFile = workfold + "\\source.txt";
	//refFeatureFile = workfold + "\\reference.txt";

	string gcpfile = workfold + "\\gcps\\gcps_hgt_26.txt";
	string chkPointFile = workfold + "\\features\\chkPoints.txt";
	string imagePointFile = workfold + "\\features\\imagePoints.txt";
	string groundPointFile = workfold + "\\features\\groundPoints.txt";
	string imageLineFile = workfold + "\\features\\2\\imageLines.txt";
	string groundLineFile = workfold + "\\features\\2\\groundLines.txt";
	//ossimFilename imageAreaFile = workfold + "\\features\\imageAreas.txt";
	//ossimFilename groundAreaFile = workfold + "\\features\\groundAreas.txt";
	string imageAreaFile = workfold + "\\features\\source.txt";
	string groundAreaFile = workfold + "\\features\\reference.txt";
	string reportfile = workfold + "\\reports\\report.txt";
	string featureReport = workfold + "\\reports\\feature_report.txt";
	string outfile = workfold + "\\rect_points.tif";

	imageLineFile = workfold + "\\source.txt";
	groundLineFile = workfold + "\\reference.txt";

	string point_shape_input = workfold + "\\shapes\\points_input.shp";
	string line_shape_input = workfold + "\\shapes\\lines_input.shp";
	string area_shape_input = workfold + "\\shapes\\areas_input.shp";
	string point_shape_output = workfold + "\\shapes\\points_output.shp";
	string line_shape_output = workfold + "\\shapes\\lines_output.shp";
	string area_shape_output = workfold + "\\shapes\\areas_output.shp";


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

	// 分别读取点线面特征

	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDUnknown);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(srcFeatureFile.c_str(), imageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(refFeatureFile.c_str(), groundFeatureList, ossimFeatureType::ossimUnknown);
	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDPointType);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGPointType);
	//prj.ReadImagePoints(srcFeatureFile, imageFeatureList);
	//prj.ReadGroundPoints(refFeatureFile, groundFeatureList);
	//prj.ReadImageLines(imageLineFile, imageFeatureList);
	//prj.ReadGroundLines(groundLineFile, groundFeatureList);
	//prj.ReadImageAreas(imageAreaFile, imageFeatureList);
	//prj.ReadGroundAreas(groundAreaFile, groundFeatureList);

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
	//ossimSensorModel* sensorModel = prj.m_sensorModel;

	if( NULL == prj.m_sensorModel )
		return;
	vector<ossimTieFeature> tieFeatureList;
	//根据序号连接对应的特征
	AppendTieFeatures(tieFeatureList, imageFeatureList, groundFeatureList);

	// 加入控制点
	//AppendTieFeaturesFromTieGptSet(tieFeatureList, prj.m_CtrlGptSet);

	//AppendTieAreas(tieFeatureList, dptAreaList, gptAreaList);
	prj.GetElevations(tieFeatureList);
	
	// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;


	prj.UpdateSensorModel(tieFeatureList, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(tieLineList, prj.m_sensorModel, prj.geom);

	prj.m_sensorModel->loadState(prj.geom);

	//NEWMAT::ColumnVector residue0 = prj.getResidue(prj.m_sensorModel, tieFeatureList);
	//double ki0=residue0.SumSquare();
	//cout<<residue0<<endl;
	//cout<<endl<<ki0<<endl<<endl;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	FeatureReport(featureReport, prj.m_sensorModel, tieFeatureList);

	//
	vector<ossimDFeature> CheckimageFeatureList;
	vector<ossimGFeature> CheckgroundFeatureList;
	//ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	//ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(imageLineFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundLineFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	vector<ossimTieFeature> ChecktieFeatureList;
	AppendTieFeatures(ChecktieFeatureList, CheckimageFeatureList, CheckgroundFeatureList);

	prj.SaveFeaturetoShape(point_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPointType);
	prj.SaveFeaturetoShape(line_shape_input, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(area_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.UpdateFeatures(prj.m_sensorModel, ChecktieFeatureList);
	prj.SaveFeaturetoShape(area_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.SaveFeaturetoShape(line_shape_output, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(point_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPointType);

	prj.Orthograph(outfile);
}


void feature_gf1()
{
	string workfold = "D:\\workspace\\GF\\GF1_WFV3_E75.9_N42.2_20140831_L1A0000322879";
	//ossimFilename workfold = "C:\\LD2010003816";
	string sourcefile = workfold + "\\GF1_WFV3_E75.9_N42.2_20140831_L1A0000322879\\GF1_WFV3_E75.9_N42.2_20140831_L1A0000322879.tiff";
	//ossimFilename elevationpath = "D:\\workspace\\Landsat\\dem_xj";
	string elevationpath = "D:\\workspace\\dem-hgt";

	string srcFeatureFile = workfold + "\\features\\source.txt";
	string refFeatureFile = workfold + "\\features\\reference.txt";

	string projection_file = workfold + "\\projection.txt";
	string gcpfile = workfold + "\\gcps\\gcps_hgt_26.txt";
	string chkPointFile = workfold + "\\features\\chkPoints.txt";
	string imagePointFile = workfold + "\\features\\imagePoints.txt";
	string groundPointFile = workfold + "\\features\\groundPoints.txt";
	string imageLineFile = workfold + "\\features\\2\\imageLines.txt";
	string groundLineFile = workfold + "\\features\\2\\groundLines.txt";
	//ossimFilename imageAreaFile = workfold + "\\features\\imageAreas.txt";
	//ossimFilename groundAreaFile = workfold + "\\features\\groundAreas.txt";
	string imageAreaFile = workfold + "\\features\\source.txt";
	string groundAreaFile = workfold + "\\features\\reference.txt";
	string reportfile = workfold + "\\reports\\report.txt";
	string featureReport = workfold + "\\reports\\feature_report.txt";
	string outfile = workfold + "\\rect_points.tif";

	string point_shape_input = workfold + "\\shapes\\points_input.shp";
	string line_shape_input = workfold + "\\shapes\\lines_input.shp";
	string area_shape_input = workfold + "\\shapes\\areas_input.shp";
	string point_shape_output = workfold + "\\shapes\\points_output.shp";
	string line_shape_output = workfold + "\\shapes\\lines_output.shp";
	string area_shape_output = workfold + "\\shapes\\areas_output.shp";

	// 加载高程数据
	ossimElevManager::instance()->loadElevationPath(elevationpath);

	// 创建RadarSat2 RPC模型
	ossimplugins::radiRpcModel* pRpcModel = new ossimplugins::radiRpcModel();
	if (!pRpcModel->parseRpcFile(sourcefile))
	{
		return;
	}

	ossimKeywordlist MapProjection;
	MyProject prj;
	prj.m_ModelType = RPCType;
	prj.InitiateSensorModel(sourcefile);
	prj.m_sensorModel = pRpcModel;
	// 指定输出投影
	prj.m_MapProjection.addFile(projection_file);
	ossimRefPtr<ossimProjection> proj;
	proj = ossimProjectionFactoryRegistry::instance()->createProjection(prj.m_MapProjection);
	prj.m_MapPar = PTR_CAST(ossimMapProjection, proj.get());
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	prj.m_CtrlGptSet = new ossimTieGptSet;
	prj.m_ChkGptSet = new ossimTieGptSet;
	//vector < ossimTieLine > tieLineList;
	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	//prj.ReadGcpAndProjection(ossimFilename(gcpfile));
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(elevationpath));//
	prj.m_DemPath = ossimFilename(elevationpath);

	// 分别读取点线面特征

	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDUnknown);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(srcFeatureFile.c_str(), imageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(refFeatureFile.c_str(), groundFeatureList, ossimFeatureType::ossimUnknown);
	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDPointType);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGPointType);
	//prj.ReadImagePoints(srcFeatureFile, imageFeatureList);
	//prj.ReadGroundPoints(refFeatureFile, groundFeatureList);
	//prj.ReadImageLines(imageLineFile, imageFeatureList);
	//prj.ReadGroundLines(groundLineFile, groundFeatureList);
	//prj.ReadImageAreas(imageAreaFile, imageFeatureList);
	//prj.ReadGroundAreas(groundAreaFile, groundFeatureList);

	//prj.GetElevations(prj.m_CtrlGptSet);
	//prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(2);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(1);

	if (NULL == prj.m_sensorModel)
		return;
	vector<ossimTieFeature> tieFeatureList;
	//根据序号连接对应的特征
	AppendTieFeatures(tieFeatureList, imageFeatureList, groundFeatureList);

	// 加入控制点
	//AppendTieFeaturesFromTieGptSet(tieFeatureList, prj.m_CtrlGptSet);

	//AppendTieAreas(tieFeatureList, dptAreaList, gptAreaList);
	prj.GetElevations(tieFeatureList);

	// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;


	prj.UpdateSensorModel(tieFeatureList, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(tieLineList, prj.m_sensorModel, prj.geom);

	prj.m_sensorModel->loadState(prj.geom);

	//NEWMAT::ColumnVector residue0 = prj.getResidue(prj.m_sensorModel, tieFeatureList);
	//double ki0=residue0.SumSquare();
	//cout<<residue0<<endl;
	//cout<<endl<<ki0<<endl<<endl;

	//prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	FeatureReport(featureReport, prj.m_sensorModel, tieFeatureList);

	//
	vector<ossimDFeature> CheckimageFeatureList;
	vector<ossimGFeature> CheckgroundFeatureList;
	ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(imageLineFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundLineFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	vector<ossimTieFeature> ChecktieFeatureList;
	AppendTieFeatures(ChecktieFeatureList, CheckimageFeatureList, CheckgroundFeatureList);

	prj.SaveFeaturetoShape(point_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPointType);
	prj.SaveFeaturetoShape(line_shape_input, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(area_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.UpdateFeatures(prj.m_sensorModel, ChecktieFeatureList);
	prj.SaveFeaturetoShape(area_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.SaveFeaturetoShape(line_shape_output, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(point_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPointType);

	prj.Orthograph(outfile);
}

void feature_landsat1()
{
	string workfold = "D:\\workspace\\Landsat\\feature_test";
	workfold = "D:\\workspace\\Landsat\\xinjiang";
	//ossimFilename workfold = "C:\\LD2010003816";
	string sourcefile = workfold + "\\LS5_TM_20100717_050010_050035_145033_FASTB_L2\\header.dat";
	sourcefile = workfold + "\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat";
	//ossimFilename elevationpath = "D:\\workspace\\Landsat\\dem_xj";
	string elevationpath = "D:\\workspace\\dem\\srtm90";

	string srcFeatureFile = workfold + "\\features\\source.txt";
	string refFeatureFile = workfold + "\\features\\reference.txt";

	string gcpfile = workfold + "\\gcps\\chk100.txt";
	string chkPointFile = workfold + "\\features\\chkPoints.txt";
	string imagePointFile = workfold + "\\features\\imagePoints.txt";
	string groundPointFile = workfold + "\\features\\groundPoints.txt";
	string imageLineFile = workfold + "\\features\\2\\imageLines.txt";
	string groundLineFile = workfold + "\\features\\2\\groundLines.txt";
	//ossimFilename imageAreaFile = workfold + "\\features\\imageAreas.txt";
	//ossimFilename groundAreaFile = workfold + "\\features\\groundAreas.txt";
	string imageAreaFile = workfold + "\\features\\source.txt";
	string groundAreaFile = workfold + "\\features\\reference.txt";
	string reportfile = workfold + "\\reports\\report.txt";
	string featureReport = workfold + "\\reports\\feature_report.txt";
	string outfile = workfold + "\\rect_points.tif";

	string point_shape_input = workfold + "\\shapes\\points_input.shp";
	string line_shape_input = workfold + "\\shapes\\lines_input.shp";
	string area_shape_input = workfold + "\\shapes\\areas_input.shp";
	string point_shape_output = workfold + "\\shapes\\points_output.shp";
	string line_shape_output = workfold + "\\shapes\\lines_output.shp";
	string area_shape_output = workfold + "\\shapes\\areas_output.shp";


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

	// 分别读取点线面特征

	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDUnknown);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(srcFeatureFile.c_str(), imageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(refFeatureFile.c_str(), groundFeatureList, ossimFeatureType::ossimUnknown);
	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDPointType);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGPointType);
	//prj.ReadImagePoints(srcFeatureFile, imageFeatureList);
	//prj.ReadGroundPoints(refFeatureFile, groundFeatureList);
	//prj.ReadImageLines(imageLineFile, imageFeatureList);
	//prj.ReadGroundLines(groundLineFile, groundFeatureList);
	//prj.ReadImageAreas(imageAreaFile, imageFeatureList);
	//prj.ReadGroundAreas(groundAreaFile, groundFeatureList);

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
	//ossimSensorModel* sensorModel = prj.m_sensorModel;

	if( NULL == prj.m_sensorModel )
		return;
	vector<ossimTieFeature> tieFeatureList;
	//根据序号连接对应的特征
	AppendTieFeatures(tieFeatureList, imageFeatureList, groundFeatureList);

	// 加入控制点
	AppendTieFeaturesFromTieGptSet(tieFeatureList, prj.m_CtrlGptSet);

	//AppendTieAreas(tieFeatureList, dptAreaList, gptAreaList);
	prj.GetElevations(tieFeatureList);

	// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;


	prj.UpdateSensorModel(tieFeatureList, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(tieLineList, prj.m_sensorModel, prj.geom);

	prj.m_sensorModel->loadState(prj.geom);

	//NEWMAT::ColumnVector residue0 = prj.getResidue(prj.m_sensorModel, tieFeatureList);
	//double ki0=residue0.SumSquare();
	//cout<<residue0<<endl;
	//cout<<endl<<ki0<<endl<<endl;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	FeatureReport(featureReport, prj.m_sensorModel, tieFeatureList);

	//
	vector<ossimDFeature> CheckimageFeatureList;
	vector<ossimGFeature> CheckgroundFeatureList;
	ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(imageLineFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundLineFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	vector<ossimTieFeature> ChecktieFeatureList;
	AppendTieFeatures(ChecktieFeatureList, CheckimageFeatureList, CheckgroundFeatureList);
	prj.GetElevations(ChecktieFeatureList);

	prj.SaveFeaturetoShape(point_shape_input, tieFeatureList, ossimFeatureType::ossimPointType);
	prj.SaveFeaturetoShape(line_shape_input, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(area_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.UpdateFeatures(prj.m_sensorModel, ChecktieFeatureList);
	prj.SaveFeaturetoShape(area_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.SaveFeaturetoShape(line_shape_output, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(point_shape_output, tieFeatureList, ossimFeatureType::ossimPointType);

	//prj.Orthograph(outfile);
}

void straightline_landsat()
{
	//ossimFilename elevationpath = "D:\\workspace\\dem";
	//ossimFilename gcpfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\points.txt";
	//ossimFilename linefile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\lines.txt";
	//ossimFilename sourcefile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\header.dat";
	//ossimFilename reportfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\report.txt";
	//ossimFilename outfile = "D:\\workspace\\testdata\\LineFeature\\LD2010003865\\rect.tif";



	//ossimFilename elevationpath = "D:\\workspace\\dem";
	//ossimFilename gcpfile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\points.txt";
	//ossimFilename linefile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\lines.txt";
	//ossimFilename sourcefile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\header.dat";
	//ossimFilename reportfile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\report.txt";
	//ossimFilename outfile = "D:\\workspace\\testdata\\LineFeature\\LD2010001519\\rect.tif";

	ossimFilename workfold = "D:\\workspace\\LD2010003816";
	//ossimFilename workfold = "C:\\LD2010003816";
	ossimFilename sourcefile = workfold + "\\header.dat";
	ossimFilename elevationpath = "D:\\workspace\\dem-hgt";
	ossimFilename gcpfile = workfold + "\\features\\points_bak.txt";
	ossimFilename chkPointFile = workfold + "\\features\\chkPoints.txt";
	ossimFilename imagePointFile = workfold + "\\features\\imagePoints.txt";
	ossimFilename groundPointFile = workfold + "\\features\\groundPoints.txt";
	ossimFilename imageLineFile = workfold + "\\features\\2\\imageLines.txt";
	ossimFilename groundLineFile = workfold + "\\features\\2\\groundLines.txt";
	//ossimFilename imageAreaFile = workfold + "\\features\\imageAreas.txt";
	//ossimFilename groundAreaFile = workfold + "\\features\\groundAreas.txt";
	ossimFilename imageAreaFile = workfold + "\\features\\12area\\imageAreas.txt";
	ossimFilename groundAreaFile = workfold + "\\features\\12area\\groundAreas.txt";
	ossimFilename reportfile = workfold + "\\reports\\report.txt";
	ossimFilename featureReport = workfold + "\\reports\\feature_report.txt";
	ossimFilename outfile = workfold + "\\rect_points.tif";

	ossimFilename point_shape_input = workfold + "\\shapes\\points_input.shp";
	ossimFilename line_shape_input = workfold + "\\shapes\\lines_input.shp";
	ossimFilename area_shape_input = workfold + "\\shapes\\areas_input.shp";
	ossimFilename point_shape_output = workfold + "\\shapes\\points_output.shp";
	ossimFilename line_shape_output = workfold + "\\shapes\\lines_output.shp";
	ossimFilename area_shape_output = workfold + "\\shapes\\areas_output.shp";


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

	// 分别读取点线面特征

	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDUnknown);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(imageLineFile.c_str(), imageFeatureList, ossimFeatureType::ossimStraightLineType);
	ReadFeatures<ossimGFeature, ossimGpt>(groundLineFile.c_str(), groundFeatureList, ossimFeatureType::ossimStraightLineType);
	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDPointType);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGPointType);
	prj.ReadImagePoints(imagePointFile, imageFeatureList);
	prj.ReadGroundPoints(groundPointFile, groundFeatureList);
	//prj.ReadImageLines(imageLineFile, imageFeatureList);
	//prj.ReadGroundLines(groundLineFile, groundFeatureList);
	//prj.ReadImageAreas(imageAreaFile, imageFeatureList);
	//prj.ReadGroundAreas(groundAreaFile, groundFeatureList);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	ossimMapProjection *MapPar = prj.m_MapPar;
	MapProjection = prj.m_MapProjection;

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(4);
	prj.m_OutBandList.push_back(3);
	prj.m_OutBandList.push_back(2);
	prj.InitiateSensorModel(sourcefile);
	//ossimSensorModel* sensorModel = prj.m_sensorModel;

	if( NULL == prj.m_sensorModel )
		return;
	vector<ossimTieFeature> tieFeatureList;
	//根据序号连接对应的特征
	AppendTieFeatures(tieFeatureList, imageFeatureList, groundFeatureList);
	/*
	for(unsigned int i = 0;i < (int)prj.m_CtrlGptSet->size();i++)
	{
		ossimTieFeature tieFeature;
		tieFeature.setId(*prj.m_CtrlGptSet->getTiePoints()[i]->GcpNumberID);
		tieFeature.setImageFeature(ossimDptFeature(prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint()));
		tieFeature.setGroundFeature(ossimGptFeature(prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint()));
		tieFeatureList.push_back(tieFeature);
	}
	for(i = 0;i < (int)prj.m_ChkGptSet->size();i++)
	{
		ossimTieFeature tieFeature;
		tieFeature.setImageFeature(ossimDptFeature(prj.m_ChkGptSet->getTiePoints()[i]->getImagePoint()));
		tieFeature.setGroundFeature(ossimGptFeature(prj.m_ChkGptSet->getTiePoints()[i]->getGroundPoint()));
		tieFeatureList.push_back(tieFeature);
	}
	
	for(i = 0;i < (int)tieLineList.size();i++)
	{
		ossimTieFeature tieFeature;
		tieFeature.m_featureType = ossimTieFeature::ossimFeatureType::ossimLineFeature;
		tieFeature.m_TiePoints.push_back(tieLineList[i].first);
		tieFeature.m_TiePoints.push_back(tieLineList[i].second);
		tieFeatureList.push_back(tieFeature);
	}*/

	//AppendTieAreas(tieFeatureList, dptAreaList, gptAreaList);
	prj.GetElevations(tieFeatureList);



	// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;


	//prj.UpdateSensorModel(tieFeatureList, prj.m_sensorModel, prj.geom);
	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(tieLineList, prj.m_sensorModel, prj.geom);

	prj.m_sensorModel->loadState(prj.geom);

	//NEWMAT::ColumnVector residue0 = prj.getResidue(prj.m_sensorModel, tieFeatureList);
	//double ki0=residue0.SumSquare();
	//cout<<residue0<<endl;
	//cout<<endl<<ki0<<endl<<endl;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	FeatureReport(featureReport, prj.m_sensorModel, tieFeatureList);

	//
	vector<ossimDFeature> CheckimageFeatureList;
	vector<ossimGFeature> CheckgroundFeatureList;
	ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimDFeature, ossimDpt>(imageLineFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundLineFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	//ReadFeatures(imageAreaFile.c_str(), CheckimageFeatureList, ossimDFeature::ossimDAreaType);
	//ReadFeatures(groundAreaFile.c_str(), CheckgroundFeatureList, ossimGFeature::ossimGAreaType);
	//ReadFeatures(imageAreaFile.c_str(), CheckimageFeatureList, ossimDFeature::ossimDPointType);
	//ReadFeatures(groundAreaFile.c_str(), CheckgroundFeatureList, ossimGFeature::ossimGPointType);
	vector<ossimTieFeature> ChecktieFeatureList;
	AppendTieFeatures(ChecktieFeatureList, CheckimageFeatureList, CheckgroundFeatureList);

	prj.SaveFeaturetoShape(point_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPointType);
	prj.SaveFeaturetoShape(line_shape_input, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(area_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.UpdateFeatures(prj.m_sensorModel, ChecktieFeatureList);
	prj.SaveFeaturetoShape(area_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.SaveFeaturetoShape(line_shape_output, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(point_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPointType);

	prj.Orthograph(outfile);
}

bool Tiff2RspfDem()
{
	ossimFilename tifFileName = "C:\\Users\\Administrator\\Desktop\\dem301.img";
	ossimFilename ossimDemFileName = "D:\\workspace\\dem_bj\\test.ras";
	ossimInit::instance()->initialize();
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim_2.0\\project\\Release\\ossimgdal_plugin.dll");
	ossimImageHandler* handler = ossimImageHandlerRegistry::instance()->open(tifFileName);
	if (!handler)
	{
		return false;
	}
	ossimKeywordlist kwl;
	handler->getImageGeometry()->saveState(kwl);
	ProjectionParameters proParam;
	proParam.DatumName = "北京54";
	proParam.TrueOriginLongitude = 117.0;
	proParam.TrueOriginLatitude = 0.0;
	proParam.EastingFalse = 500000;
	proParam.NorthingFalse = 0;
	proParam.PixelSizeX = 30.0;
	proParam.PixelSizeY = 30.0;
	ossimMapProjection* projection = CreateProjection(proParam, kwl);
	ossimImageGeometry imageGeom;
	imageGeom.setProjection(projection);
	handler->setImageGeometry(&imageGeom);


	ossimImageFileWriter* writer =
		ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("image/ras"));
	ossimImageRenderer* renderer = new ossimImageRenderer;
	renderer->connectMyInputTo(handler);

	ossimFilename Tmp = ossimDemFileName;
	Tmp.setExtension(".ras");
	writer->setFilename(Tmp);
	//emit addListener(writer);
	ossimStdOutProgress progress(0,true);
	writer->addListener(&progress);
	writer->connectMyInputTo(0,renderer);
	if(!writer->execute())
	{
		writer->removeListener(&progress);
		//delete writer;
		handler->close();
		delete handler;
		return false;
	}
	handler->close();
	writer->disableListener();
	writer->removeListener(&progress);
	writer->disconnectAllInputs();
	renderer->disconnectAllInputs();

	return true;
}


void areas_alos()
{

	//ossimFilename workfold = "E:\\testdata\\AlosData\\d1001111-001";
	////ossimFilename workfold = "C:\\LD2010003816";
	//ossimFilename sourcefile = workfold + "\\IMG-ALPSMW220853075_O1B2R_UW.tif";
	//ossimFilename elevationpath = "D:\\workspace\\dem";
	//ossimFilename gcpfile = workfold + "\\features\\points.txt";
	//ossimFilename imagePointFile = workfold + "\\features\\imagePoints.txt";
	//ossimFilename groundPointFile = workfold + "\\features\\groundPoints.txt";
	//ossimFilename imageLineFile = workfold + "\\features\\imageLines.txt";
	//ossimFilename groundLineFile = workfold + "\\features\\groundLines.txt";
	//ossimFilename imageAreaFile = workfold + "\\features\\imageAreas.txt";
	//ossimFilename groundAreaFile = workfold + "\\features\\groundAreas.txt";
	//ossimFilename reportfile = workfold + "\\reports\\report.txt";
	//ossimFilename featureReport = workfold + "\\reports\\feature_report.txt";
	//ossimFilename outfile = workfold + "\\rect.tif";

	//ossimFilename point_shape_input = workfold + "\\shapes\\points_input.shp";
	//ossimFilename line_shape_input = workfold + "\\shapes\\lines_input.shp";
	//ossimFilename area_shape_input = workfold + "\\shapes\\areas_input.shp";
	//ossimFilename point_shape_output = workfold + "\\shapes\\points_output.shp";
	//ossimFilename line_shape_output = workfold + "\\shapes\\lines_output.shp";
	//ossimFilename area_shape_output = workfold + "\\shapes\\areas_output.shp";


	//MyProject prj;
	//AlosPRISM alosUti(workfold);
	//if(!alosUti.getInitState())
	//{
	//	cout<<"warning: Alos PRISM数据\""<<workfold<<"\"初始化失败！"<<endl;
	//	return;
	//}
	//prj.theMgr = ossimElevManager::instance();
	//if(!prj.theMgr->loadElevationPath(ossimFilename(elevationpath)))
	//{
	//	cout<<"warning: 加载DEM失败！"<<endl;
	//	return;
	//}
	//prj.ReadGcpAndProjection(gcpfile);
	//prj.theMgr = ossimElevManager::instance();
	//prj.theMgr->loadElevationPath(ossimFilename(elevationpath));
	//prj.m_DemPath = ossimFilename(elevationpath);
	//prj.GetElevations(prj.m_CtrlGptSet);
	//prj.GetElevations(prj.m_ChkGptSet);

	//prj.m_MapProjection = alosUti.getKeywordlist();
	//prj.m_MapPar = alosUti.getMapProjection();

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


	//vector<ossimDFeature> imageFeatureList;
	//vector<ossimGFeature> groundFeatureList;
	//// 分别读取点线面特征
	////prj.ReadImagePoints(imagePointFile, imageFeatureList);
	////prj.ReadGroundPoints(groundPointFile, groundFeatureList);
	////prj.ReadImageLines(imageLineFile, imageFeatureList);
	////prj.ReadGroundLines(groundLineFile, groundFeatureList);
	////prj.ReadImageAreas(imageAreaFile, imageFeatureList);
	////prj.ReadGroundAreas(groundAreaFile, groundFeatureList);
	////ReadFeatures(imagePointFile.c_str(), imageFeatureList, ossimDFeature::ossimDUnknown);
	////ReadFeatures(groundPointFile.c_str(), groundFeatureList, ossimGFeature::ossimGUnknown);
	////ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDPointType);
	////ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGPointType);
	//ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), imageFeatureList, ossimFeatureType::ossimPolygonType);
	//ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), groundFeatureList, ossimFeatureType::ossimPolygonType);

	//prj.m_ImgFileNameUnc = sourcefile;
	//prj.m_OutBandList.clear();
	//prj.m_OutBandList.push_back(0);

	//if( NULL == prj.m_sensorModel )
	//	return;
	//vector<ossimTieFeature> tieFeatureList;
	////根据序号连接对应的特征
	//AppendTieFeatures(tieFeatureList, imageFeatureList, groundFeatureList);
	//prj.GetElevations(tieFeatureList);



	//// 使用RPC模型
	////ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	////if(!handler) return;   //应该弹出警告对话框
	////prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	////prj.m_sensorModel->m_proj = prj.m_MapPar;
	////handler->close();
	////delete handler;


	//prj.UpdateSensorModel(tieFeatureList, prj.m_sensorModel, prj.geom);
	////prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	////prj.UpdateSensorModel(tieLineList, prj.m_sensorModel, prj.geom);


	////生成虚拟控制点
	////ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	////if(!handler) return;   //应该弹出警告对话框
	////ossimTieGptSet *gptSet = createTieGptSet(handler->getImageRectangle(), *prj.m_sensorModel);
	////ossimFilename virtualGcpfile = workfold + "\\features\\virtualGcps.txt";
	////prj.SavePointToFile(virtualGcpfile, gptSet, prj.m_ChkGptSet);
	////handler->close();
	////delete handler;

	//prj.m_sensorModel->loadState(prj.geom);

	////NEWMAT::ColumnVector residue0 = prj.getResidue(prj.m_sensorModel, tieFeatureList);
	////double ki0=residue0.SumSquare();
	////cout<<residue0<<endl;
	////cout<<endl<<ki0<<endl<<endl;

	//prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	//FeatureReport(featureReport, prj.m_sensorModel, tieFeatureList);

	////
	//vector<ossimDFeature> CheckimageFeatureList;
	//vector<ossimGFeature> CheckgroundFeatureList;
	//ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	//ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	//vector<ossimTieFeature> ChecktieFeatureList;

	//AppendTieFeatures(ChecktieFeatureList, CheckimageFeatureList, CheckgroundFeatureList);

	//prj.SaveFeaturetoShape(point_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPointType);
	//prj.SaveFeaturetoShape(line_shape_input, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	//prj.SaveFeaturetoShape(area_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	//prj.UpdateFeatures(prj.m_sensorModel, ChecktieFeatureList);
	//prj.SaveFeaturetoShape(area_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	//prj.SaveFeaturetoShape(line_shape_output, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	//prj.SaveFeaturetoShape(point_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPointType);

	//return;
	//prj.OutputReport1(featureReport, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	////prj.Orthograph(outfile);
}

void areas_quickbird()
{
	ossimFilename workfold = "E:\\testdata\\QB\\Src\\1\\052407492010_01_P002";
	ossimFilename sourcefile = workfold + "\\10OCT21042139-P2AS-052407492010_01_P002.TIL";
	ossimFilename elevationpath = "D:\\workspace\\dem_xn";

	ossimFilename gcpfile = workfold + "\\chekpt.txt";
	ossimFilename imageAreaFile = workfold + "\\features\\origin\\imageAreas.txt";
	ossimFilename groundAreaFile = workfold + "\\features\\origin\\groundAreas.txt";
	ossimFilename reportfile = workfold + "\\reports\\report.txt";
	ossimFilename featureReport = workfold + "\\reports\\feature_report.txt";
	ossimFilename outfile = workfold + "\\rect.tif";

	ossimFilename point_shape_input = workfold + "\\shapes\\points_input.shp";
	ossimFilename line_shape_input = workfold + "\\shapes\\lines_input.shp";
	ossimFilename area_shape_input = workfold + "\\shapes\\areas_input.shp";
	ossimFilename point_shape_output = workfold + "\\shapes\\points_output.shp";
	ossimFilename line_shape_output = workfold + "\\shapes\\lines_output.shp";
	ossimFilename area_shape_output = workfold + "\\shapes\\areas_output.shp";

	MyProject prj;
	QuickbirdRpcModel *qbRpcModel = new QuickbirdRpcModel(sourcefile);
	prj.ReadGcpAndProjection(gcpfile, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(elevationpath));//
	prj.m_DemPath=ossimFilename(elevationpath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.m_sensorModel = qbRpcModel->m_sensorModel;
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	// 分别读取点线面特征
	//prj.ReadImagePoints(imagePointFile, imageFeatureList);
	//prj.ReadGroundPoints(groundPointFile, groundFeatureList);
	//prj.ReadImageLines(imageLineFile, imageFeatureList);
	//prj.ReadGroundLines(groundLineFile, groundFeatureList);
	//prj.ReadImageAreas(imageAreaFile, imageFeatureList);
	//prj.ReadGroundAreas(groundAreaFile, groundFeatureList);
	//ReadFeatures(imagePointFile.c_str(), imageFeatureList, ossimDFeature::ossimDUnknown);
	//ReadFeatures(groundPointFile.c_str(), groundFeatureList, ossimGFeature::ossimGUnknown);
	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDPointType);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGPointType);
	ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), imageFeatureList, ossimFeatureType::ossimPolygonType);
	ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), groundFeatureList, ossimFeatureType::ossimPolygonType);

	prj.m_ImgFileNameUnc = sourcefile;
	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);

	if( NULL == prj.m_sensorModel )
		return;
	vector<ossimTieFeature> tieFeatureList;
	//根据序号连接对应的特征
	AppendTieFeatures(tieFeatureList, imageFeatureList, groundFeatureList);
	prj.GetElevations(tieFeatureList);



	// 使用RPC模型
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//prj.m_sensorModel = createRpcModelFromProjection(handler->getImageRectangle(), *prj.m_sensorModel);
	//prj.m_sensorModel->m_proj = prj.m_MapPar;
	//handler->close();
	//delete handler;


	prj.UpdateSensorModel(tieFeatureList, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	//prj.UpdateSensorModel(tieLineList, prj.m_sensorModel, prj.geom);


	//生成虚拟控制点
	//ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(sourcefile);
	//if(!handler) return;   //应该弹出警告对话框
	//ossimTieGptSet *gptSet = createTieGptSet(handler->getImageRectangle(), *prj.m_sensorModel);
	//ossimFilename virtualGcpfile = workfold + "\\features\\virtualGcps.txt";
	//prj.SavePointToFile(virtualGcpfile, gptSet, prj.m_ChkGptSet);
	//handler->close();
	//delete handler;

	prj.m_sensorModel->loadState(prj.geom);

	//NEWMAT::ColumnVector residue0 = prj.getResidue(prj.m_sensorModel, tieFeatureList);
	//double ki0=residue0.SumSquare();
	//cout<<residue0<<endl;
	//cout<<endl<<ki0<<endl<<endl;

	prj.OutputReport(reportfile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	FeatureReport(featureReport, prj.m_sensorModel, tieFeatureList);

	//
	vector<ossimDFeature> CheckimageFeatureList;
	vector<ossimGFeature> CheckgroundFeatureList;
	ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), CheckimageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), CheckgroundFeatureList, ossimFeatureType::ossimUnknown);
	vector<ossimTieFeature> ChecktieFeatureList;

	AppendTieFeatures(ChecktieFeatureList, CheckimageFeatureList, CheckgroundFeatureList);
	prj.GetElevations(ChecktieFeatureList);
	//FeatureReport(featureReport, prj.m_sensorModel, ChecktieFeatureList);

	prj.SaveFeaturetoShape(point_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPointType);
	prj.SaveFeaturetoShape(line_shape_input, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(area_shape_input, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.UpdateFeatures(prj.m_sensorModel, ChecktieFeatureList);
	prj.SaveFeaturetoShape(area_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPolygonType);
	prj.SaveFeaturetoShape(line_shape_output, ChecktieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(point_shape_output, ChecktieFeatureList, ossimFeatureType::ossimPointType);

	return;
	prj.OutputReport1(featureReport, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
}

void spotTest()
{
	ossimFilename workfold = "G:\\testdata\\spot5\\beijing";
	workfold = "D:\\workspace\\Spot5\\tianjin";
	//ossimFilename workfold = "E:\\beijing";
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\281268-20040602-2.5\\scene01\\imagery.tif";
	sourceFile = workfold + "\\282270-20040518-2.5\\scene01\\imagery.tif";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";

	//ossimFilename demPath = "I:\\spot5\\beijing";
	//ossimFilename demPath = "D:\\workspace\\dem";
	//ossimFilename demPath = "E:\\dem_bj";
	//ossimFilename chkFile = workfold + "\\281268_2m5_gcp09033_hgt.txt";
	ossimFilename chkFile = workfold + "\\wgs84_hgt.txt";
	chkFile = workfold + "\\gcp_wgs84_tm.txt";
	//ossimFilename chkFile = workfold + "\\gcp_all_hgt.txt";
	ossimFilename demPath = "D:\\workspace\\dem\\aster30";
	
	ossimFilename gcp_hgt_File = workfold + "\\wgs84_gcp_hgt.txt";
	//ossimFilename chkFile = workfold + "\\gcp.txt";
	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report_orbit.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(chkFile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.SavePointToFile(gcp_hgt_File, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;


	fstream fs;
	fs.open(reportFile.c_str(), ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(2);

	//prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);

	////prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet);
	////prj.m_sensorModel->updateModel();
	////prj.m_sensorModel->saveState(prj.geom);
	//prj.m_sensorModel->loadState(prj.geom);

	//prj.OutputReport(reportFile, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	prj.OutputReport(fs, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	fs.close();

	//prj.Orthograph(outFile);
}


void Hj1Test()
{
	ossimInit::instance()->initialize();
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401050432_201401050443";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_SYC_201401050300_201401050309";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401070348_201401070354";
	//ossimFilename workfold = "E:\\beijing";
	ossimFilename workfold = "E:\\HJ1\\test";
	MyProject prj;
	//ossimFilename sourceFile = workfold + "\\IMAGE.tif";
	ossimFilename sourceFile = workfold + "\\IMAGE_1.tif";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";

	ossimFilename gcpFile = workfold + "\\gcp.txt";
	ossimFilename demPath = "D:\\workspace\\dem\\aster30";

	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpFile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	fstream fs;
	fs.open(reportFile.c_str(), ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(2);

	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.m_sensorModel->loadState(prj.geom);

	prj.OutputReport(fs, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	fs.close();

	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);
	prj.Orthograph(outFile);
}

void Hj1CalibrationTest()
{
	ossimInit::instance()->initialize();
	//ossimFilename workfold = "E:\\beijing";
	//ossimFilename workfold = "E:\\HJ1\\test";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401050432_201401050443";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_SYC_201401050300_201401050309";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_KSC_201401070348_201401070354";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-1_MYC_201403240251_201403240259";
	//ossimFilename workfold = "E:\\HJ1\\HJ-1A_CCD-2_MYC_201403240251_201403240259";
	ossimFilename workfold = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201403230317_201403230320";
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\IMAGE.tif";
	//ossimFilename sourceFile = workfold + "\\IMAGE_1.tif";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";

	ossimFilename gcpFile = workfold + "\\gcp.txt";
	ossimFilename demPath = "D:\\workspace\\dem\\srtm90";

	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report.txt";
	ossimFilename residualFile = workfold + "\\residual.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpFile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);

	prj.SavePointToFile(workfold + "\\gcp-hgt.txt", prj.m_CtrlGptSet, NULL);

	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;


	ossimDpt imagepoint,cimagepoint;
	ossimGpt goundpoint,tGpt;

	int num = static_cast<int>(prj.m_CtrlGptSet->size());

	vector<ossimRefPtr<ossimTieGpt> >& theTPV = prj.m_CtrlGptSet->refTiePoints();

	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = prj.m_MapPar->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
		tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	//setup initail values
	int iter=0;
	int iter_max = 200;
	double minDelta = 1e-5;
	double optimizer_delta = 1e10;
	vector<int> exteriorList;
	vector<int> innerList;
	exteriorList.push_back(ossimHj1Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	exteriorList.push_back(ossimHj1Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	exteriorList.push_back(ossimHj1Model::AdjustParamIndex::CCD_YAW_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::ROLL_A1_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::PITCH_A1_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::YAW_A1_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::ROLL_A2_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::PITCH_A2_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::YAW_A2_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::LINE_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A0_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A1_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A2_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A3_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A4_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::CCDLEN_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::AX0_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::AY0_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::AY1_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::AY2_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::AY3_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::AY4_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::X_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::Y_OFFSET);
	//exteriorList.push_back(ossimHj1Model::AdjustParamIndex::Z_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::CCD_ROLL_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::CCD_PITCH_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::CCD_YAW_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A0_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A1_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A2_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A3_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::LINE_A4_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::AX0_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::AY0_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::AY1_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::AY2_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::AY3_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::AY4_OFFSET);
	//innerList.push_back(ossimHj1Model::AdjustParamIndex::LINE_OFFSET);
	ossimAdjustmentInfo cadj;
	prj.m_sensorModel->getAdjustment(cadj);
	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	int np = (int)parmlist.size();
	NEWMAT::ColumnVector old_parm(np), new_parm(np);
	for(int n=0;n<np;++n)
	{
		old_parm(n+1) = parmlist[n].getParameter();
	}
	while(optimizer_delta > minDelta && iter < iter_max)
	{
		// do exterior orientation
		prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, exteriorList);
		prj.m_sensorModel->updateModel();
		prj.m_sensorModel->saveState(prj.geom);
		
		// do inner orientation
		prj.m_sensorModel->optimizeFit(*prj.m_CtrlGptSet, innerList);
		prj.m_sensorModel->updateModel();
		prj.m_sensorModel->saveState(prj.geom);

		ossimAdjustmentInfo cadj;
		prj.m_sensorModel->getAdjustment(cadj);
		std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
		for(int n=0;n<np;++n)
		{
			new_parm(n+1) = parmlist[n].getParameter();
			cout<<new_parm(n+1)<<endl;
		}
		// then calculate the change value of the optimizers
		optimizer_delta = (new_parm - old_parm).NormInfinity();
		old_parm = new_parm;
		cout<<"iteration "<<1+iter++<<" :"<<optimizer_delta<<endl;
	}// end while

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		ossimDpt dpt = prj.m_sensorModel->m_proj->forward(*(*tit));
		ossimGpt gpt(dpt.x,dpt.y);
		(*tit)->setGroundPoint(ossimGpt(dpt.x,dpt.y,(*tit)->hgt));
	}

	fstream fs;
	fs.open(reportFile.c_str(), ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(2);
	prj.OutputReport(fs, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	fs.close();

	fstream ofs;
	ofs.open(residualFile.c_str(), ios_base::out);
	for (int i = 0;i < prj.m_CtrlGptSet->getTiePoints().size();++i)
	{
		ossimGpt ll;
		prj.m_sensorModel->lineSampleToWorld(prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint(), ll);
		ossimDpt gpt = prj.m_MapPar->forward(ll);
		double delta_lat = gpt.x - prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint().lat;
		double delta_lon = gpt.y - prj.m_CtrlGptSet->getTiePoints()[i]->getGroundPoint().lon;
		if (0 != i)
		{
			ofs<<endl;
		}
		ofs<<i+1<<"\t"
			<<prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().x<<"\t"
			<<prj.m_CtrlGptSet->getTiePoints()[i]->getImagePoint().y<<"\t"
			//<<gptSet->getTiePoints()[i]->getGroundPoint().lat<<"\t"
			//<<gptSet->getTiePoints()[i]->getGroundPoint().lon<<"\t"
			//<<gptSet->getTiePoints()[i]->getGroundPoint().hgt<<"\t"
			<<delta_lat<<"\t"
			<<delta_lon;
	}
	ofs.close();

	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);
	system("pause");
	prj.Orthograph(outFile);
}

void TheosTest()
{
	ossimFilename workfold = "D:\\workspace\\Theos\\112577-116719-TH2014002255(1A1116)";
	//ossimFilename workfold = "E:\\HJ1\\test";
	//ossimFilename workfold = "E:\\beijing";
	MyProject prj;
	ossimFilename sourceFile = workfold + "\\THEOS1_20140221065605507_14002255-0_01_1_1\\imagery.tif";
	//ossimFilename sourceFile = workfold + "\\IMAGE_1.tif";
	//ossimFilename demPath = "D:\\workspace\\beijingDem";

	ossimFilename gcpFile = workfold + "\\gcps.txt";
	ossimFilename demPath = "D:\\workspace\\dem\\aster30";

	ossimFilename outFile = workfold + "\\rect.tif";
	ossimFilename reportFile = workfold + "\\report.txt";
	ModelType modelType = OrbitType;

	prj.m_ImgFileNameUnc = sourceFile;
	prj.m_ImgFileNameout = outFile;
	prj.m_ModelType = modelType;
	prj.ReadGcpAndProjection(ossimFilename(gcpFile));

	prj.theMgr = ossimElevManager::instance();
	prj.theMgr->loadElevationPath(ossimFilename(demPath));
	prj.m_DemPath=ossimFilename(demPath);

	prj.GetElevations(prj.m_CtrlGptSet);
	prj.GetElevations(prj.m_ChkGptSet);
	
	prj.InitiateSensorModel();
	prj.m_sensorModel->m_proj = prj.m_MapPar;

	fstream fs;
	fs.open(reportFile.c_str(), ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(2);

	prj.UpdateSensorModel(*prj.m_CtrlGptSet, prj.m_sensorModel, prj.geom);
	prj.m_sensorModel->loadState(prj.geom);

	prj.OutputReport(fs, prj.m_sensorModel, prj.m_CtrlGptSet, prj.m_ChkGptSet);
	fs.close();

	prj.m_OutBandList.clear();
	prj.m_OutBandList.push_back(0);
	prj.Orthograph(outFile);
}

//void AlosBatchOrth(ossimFilename inputFold, ossimFilename outputFold)
//{
//	AlosBatch *m_alosBatch = new AlosBatch(inputFold, outputFold, 0);
//	m_alosBatch->Run();
//}


void GCP2Shape()
{
	ossimFilename gcpfile = "E:\\testdata\\QB\\Src\\gcp_pci.txt";
	ossimFilename filenametosave = "E:\\testdata\\QB\\Src\\gcp_pci.shp";
	//ossimFilename gcpfile = "E:\\testdata\\QB\\Src\\052407492010_01_P002\\gcp.txt";
	//ossimFilename filenametosave = "E:\\testdata\\QB\\Src\\052407492010_01_P002\\gcp.shp";
	MyProject prj;
	prj.ReadGcpAndProjection(gcpfile);
	prj.SavePoint2Shape(filenametosave, prj.m_CtrlGptSet);
}

void Feature2Shapefile()
{
	ossimFilename workfold = "D:\\test\\qb";

	ossimFilename gcpfile = workfold + "\\chkPoints.txt";
	ossimFilename imageAreaFile = workfold + "\\features\\imageAreas.txt";
	ossimFilename groundAreaFile = workfold + "\\features\\groundAreas.txt";

	ossimFilename point_shape = workfold + "\\shapes\\points.shp";
	ossimFilename line_shape = workfold + "\\shapes\\lines.shp";
	ossimFilename area_shape = workfold + "\\shapes\\areas.shp";
	ossimFilename chk_shape = workfold + "\\shapes\\chk.shp";

	MyProject prj;
	prj.ReadGcpAndProjection(gcpfile, prj.m_CtrlGptSet, prj.m_ChkGptSet);

	vector<ossimDFeature> imageFeatureList;
	vector<ossimGFeature> groundFeatureList;
	// 分别读取点线面特征
	ReadFeatures<ossimDFeature, ossimDpt>(imageAreaFile.c_str(), imageFeatureList, ossimFeatureType::ossimUnknown);
	ReadFeatures<ossimGFeature, ossimGpt>(groundAreaFile.c_str(), groundFeatureList, ossimFeatureType::ossimUnknown);
	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDPointType);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGPointType);
	//ReadFeatures(imageAreaFile.c_str(), imageFeatureList, ossimDFeature::ossimDAreaType);
	//ReadFeatures(groundAreaFile.c_str(), groundFeatureList, ossimGFeature::ossimGAreaType);

	vector<ossimTieFeature> tieFeatureList;
	//根据序号连接对应的特征
	AppendTieFeatures(tieFeatureList, imageFeatureList, groundFeatureList);
	
	prj.SaveFeaturetoShape(point_shape, tieFeatureList, ossimFeatureType::ossimPointType);
	prj.SaveFeaturetoShape(line_shape, tieFeatureList, ossimFeatureType::ossimStraightLineType);
	prj.SaveFeaturetoShape(area_shape, tieFeatureList, ossimFeatureType::ossimPolygonType);

	prj.SavePoint2Shape(chk_shape, prj.m_CtrlGptSet);

	return;
}
