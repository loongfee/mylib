#include <iostream>
#include <iterator>
#include <fstream>
using namespace std;

#include <func.h>
#include <mprojectdefine.h>
#include "..\QuickbirdRpcModel.h"

#include <strUtil.h>
#include <fileUtil.h>

#include <ossim/base/ossimProcessInterface.h>
#include <ossim/base/ossimObjectFactoryRegistry.h>
using namespace mylib;

//void matchPoint(ossimFilename sourcePath, ossimFilename stdPath, ossimFilename outPath)
//{
//	//ossimInit::instance()->initialize();	
//	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimsurf_plugin.dll");
//	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimreg_plugin.dll");
//	bool result;
//	ossimKeywordlist geomp;
//	vector<ossimFilename> inputfile;
//	vector<ossimFilename> allmatchfile;
//	ossimFilename master,slave,file_tmp,outimgfile;
//
//
//	// 标准影像
//	ossimDirectory dirsrc(stdPath);
//	//dirsrc.findAllFilesThatMatch(inputfile,ossimString("TM_[0-9]*.tif"));//取得标准影像
//	dirsrc.findAllFilesThatMatch(inputfile,ossimString(".tif"));//取得标准影像
//
//	// 待校正影像
//	// landsat5
//	allmatchfile.clear();
//	char* sx="header.dat";
//	_chdir(sourcePath);
//	Search_Directory(sx, allmatchfile);
//
//	// landsat 7
//	char* sLx="L*_HRF.FST";
//	//_chdir("E:\\data_share\\error");
//	_chdir(sourcePath);
//	Search_Directory(sLx, allmatchfile);
//
//	//for(int h=0;h<fileList.size();h++)
//	//{
//	//	allmatchfile.push_back(fileList[h]);
//	//}
//
//
//	ifstream is;
//	char thePathRowNumber[10],theAcquisitionDate[9],orderNUM[13];;
//	ossimString str_tmp1,str_tmp2,str_tmp3,str_tmp4;
//	ossimString strpath,strrow,strcenterll,strcenterllTif,orderstr;
//	ossimFilename newdir;
//	bool flag=false;
//	for(int i = 0;i < allmatchfile.size();i++)
//	{
//		flag = false;
//
//		is.open(allmatchfile[i].c_str());
//
//		if( allmatchfile[i].ext().contains("FST") || allmatchfile[i].ext().contains("fst")) 
//		{
//			is.seekg(34, ios::beg);
//			is.get(thePathRowNumber, 10,' ');
//			is.seekg(70, ios::beg);
//			is.get(theAcquisitionDate, 9);
//			ossimFfL7 *L7Header = new ossimFfL7(allmatchfile[i].c_str());
//			double ll = L7Header->getParam(4);
//			strcenterll = ossimString::toString(ll);
//			L7Header->unref_nodelete();
//		}
//		else {
//			is.seekg(26, ios::beg);
//			is.get(thePathRowNumber, 10,' ');
//			is.seekg(54, ios::beg);
//			is.get(theAcquisitionDate, 9);
//
//			for (int ii=0; ii < 15; ii++)
//			{
//				is.seekg(PROJ_PARAM_OFFSET[ii], ios::beg);
//				is.get(theUsgsProjParam[ii], 25);
//			}
//			strcenterll=theUsgsProjParam[4];
//		}
//		is.close();
//
//
//		//if ((strcenterll.toDouble()-strcenterllTif.toDouble()) >0.5) continue;
//
//		str_tmp1=thePathRowNumber;
//		strpath = str_tmp1.substr(0, 3);
//		strrow = str_tmp1.substr(4, 3);
//		//strpath=inputfile[i].file().substr(1,3);
//		//strrow=inputfile[i].file().substr(5,3);
//		//str_tmp2=strpath+"/"+strrow;
//
//		ossimString strFilter = "p" + strpath + "r" + strrow;
//
//		//if(str_tmp1.substr(0,7)==(str_tmp2))
//		for(int j = 0;j < static_cast<int>(inputfile.size());++j)
//		{
//			if(inputfile[j].file().contains(strFilter))
//			{
//				slave = allmatchfile[i];
//				ifstream isordernum;
//				ossimFilename fill = slave.path()+"\\"+"ProductDescription.self";
//				if (!fill.exists()) continue;
//				isordernum.open(fill.c_str());
//				isordernum.seekg(37, ios::beg);
//				isordernum.get(orderNUM, 13);
//				orderstr=orderNUM;
//				isordernum.close();
//
//				master = inputfile[j];
//
//				ossimString ACQUISITION_DATE = theAcquisitionDate;
//
//				cout<<slave.c_str()<<endl;
//				if( slave.ext().contains("FST") || slave.ext().contains("fst")) 
//					file_tmp= outPath+"\\L7-ETM+-"+strpath+"-"+strrow+"-"+ACQUISITION_DATE+"-"+orderstr+"-L4-GCP-ALL.TXT";
//				else
//					file_tmp= outPath+"\\L5-TM-"+strpath+"-"+strrow+"-"+ACQUISITION_DATE+"-"+orderstr+"-L4-GCP-ALL.TXT";
//				if (file_tmp.exists()) continue;
//
//				ossimRefPtr<ossimObject> icObject = ossimObjectFactoryRegistry::instance()->createObject(ossimString("ossimImageCorrelator"));//自定义自动选取同名点的类，为dll加载形式。
//				ossimOutputSource* icSource = PTR_CAST(ossimOutputSource, icObject.get());
//				ossimStdOutProgress progress(0,true);
//				icSource->addListener(&progress);
//				ossimProcessInterface* icProcessInterface = PTR_CAST(ossimProcessInterface, icObject.get());
//				ossimPropertyInterface* icPropertyInterface = PTR_CAST(ossimPropertyInterface, icObject.get());
//				ossimRefPtr<ossimProperty> masterBand        = icSource->getProperty("master_band");//参考影像选用的波段
//				ossimRefPtr<ossimProperty> slaveBand         = icSource->getProperty("slave_band");//待校正影像选用的波段
//				ossimRefPtr<ossimProperty> scaleRatio        = icSource->getProperty("scale_ratio");//尺度因子
//				ossimRefPtr<ossimProperty> cornerDensity     = icSource->getProperty("corner_density");//harriscorner算子选取点的密度
//				ossimRefPtr<ossimProperty> minCorrel         = icSource->getProperty("min_correl");//最小相关度
//				ossimRefPtr<ossimProperty> templateRadius    = icSource->getProperty("template_radius");//模板半径
//				ossimRefPtr<ossimProperty> slaveAccuracy     = icSource->getProperty("slave_accuracy");//待校正影像的误差
//				ossimRefPtr<ossimProperty> projectionType    = icSource->getProperty("projection_type");//采用的投影方式
//				ossimRefPtr<ossimProperty> outputFilename    = icSource->getProperty("output_filename");//输出同名点txt文件
//				str_tmp1= "3";
//				icPropertyInterface->setProperty("master_band", str_tmp1);
//
//				str_tmp1= "3";
//				icPropertyInterface->setProperty("slave_band", str_tmp1);
//
//				//str_tmp=mStaticTextEparameter->GetValue().c_str();
//				str_tmp1="150";
//
//				icPropertyInterface->setProperty("slave_accuracy", str_tmp1);//1102
//
//				str_tmp1="0.7";
//
//				icPropertyInterface->setProperty("min_correl", str_tmp1);//1102
//
//				str_tmp1="50";
//				icPropertyInterface->setProperty("template_radius", str_tmp1);//1102
//
//				icPropertyInterface->setProperty("master_filename", master);//参考影像
//
//				icPropertyInterface->setProperty("slave_filename", slave);//待校正影像
//
//
//
//
//				icPropertyInterface->setProperty("output_filename", file_tmp); 
//				// if (!file_tmp.exists()) 
//
//				result = icProcessInterface->execute();
//
//				icObject=NULL;
//
//				flag = true;
//				break;
//			}
//		}
//	}
//}

void matchPoint(ossimFilename sourcePath, ossimFilename stdPath, ossimFilename outPath)
{
	//ossimInit::instance()->initialize();	
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossim-registration.dll");
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimreg_plugin.dll");
	ossimInit::instance()->loadPlugins("D:\\opensource\\ossim\\2.0\\project\\Release\\ossimgdal_plugin.dll");
	bool result;
	ossimKeywordlist geomp;
	vector<ossimFilename> inputfile;
	ossimFilename master,slave,file_tmp,outimgfile;

	
	ossimRefPtr<ossimObject> icObject = ossimObjectFactoryRegistry::instance()->createObject(ossimString("ossimImageCorrelator"));//自定义自动选取同名点的类，为dll加载形式。
	ossimOutputSource* icSource = PTR_CAST(ossimOutputSource, icObject.get());
	ossimStdOutProgress progress(0,true);
	icSource->addListener(&progress);
	ossimProcessInterface* icProcessInterface = PTR_CAST(ossimProcessInterface, icObject.get());
	ossimPropertyInterface* icPropertyInterface = PTR_CAST(ossimPropertyInterface, icObject.get());
	ossimRefPtr<ossimProperty> masterBand        = icSource->getProperty("master_band");//参考影像选用的波段
	ossimRefPtr<ossimProperty> slaveBand         = icSource->getProperty("slave_band");//待校正影像选用的波段
	ossimRefPtr<ossimProperty> scaleRatio        = icSource->getProperty("scale_ratio");//尺度因子
	ossimRefPtr<ossimProperty> cornerDensity     = icSource->getProperty("corner_density");//harriscorner算子选取点的密度
	ossimRefPtr<ossimProperty> minCorrel         = icSource->getProperty("min_correl");//最小相关度
	ossimRefPtr<ossimProperty> templateRadius    = icSource->getProperty("template_radius");//模板半径
	ossimRefPtr<ossimProperty> slaveAccuracy     = icSource->getProperty("slave_accuracy");//待校正影像的误差
	ossimRefPtr<ossimProperty> projectionType    = icSource->getProperty("projection_type");//采用的投影方式
	ossimRefPtr<ossimProperty> outputFilename    = icSource->getProperty("output_filename");//输出同名点txt文件

	ossimString master_band = "3";
	ossimString slave_band = "3";
	ossimString slave_accuracy = "100";
	ossimString point_number = "200";
	ossimString tile_size = "64";

	ossimString sift_nfeatures = "0";
	ossimString sift_noctavelayers = "3";
	ossimString sift_contrastthreshold = "0.01";
	ossimString sift_edgethreshold = "10";
	ossimString sift_sigma = "1.6";
	
	icPropertyInterface->setProperty("master_band", master_band);
	icPropertyInterface->setProperty("slave_band", slave_band);
	icPropertyInterface->setProperty("slave_accuracy", slave_accuracy);
	icPropertyInterface->setProperty("point_number", point_number);
	icPropertyInterface->setProperty("slave_filename", sourcePath);//待校正影像
	icPropertyInterface->setProperty("master_filename", stdPath);//参考影像
	icPropertyInterface->setProperty("tile_size", tile_size);

	icPropertyInterface->setProperty("sift_nfeatures", sift_nfeatures);
	icPropertyInterface->setProperty("sift_noctavelayers", sift_noctavelayers);
	icPropertyInterface->setProperty("sift_contrastthreshold", sift_contrastthreshold);
	icPropertyInterface->setProperty("sift_edgethreshold", sift_edgethreshold);
	icPropertyInterface->setProperty("sift_sigma", sift_sigma);


	icPropertyInterface->setProperty("output_filename", "tmp.xml");

	result = icProcessInterface->execute();

	icObject=NULL;

	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	readGcpFile("tmp.xml", gcpSet, chkSet, &prjKwl);

	ossimImageHandler* handlerM = ossimImageHandlerRegistry::instance()->open(stdPath);
	handlerM->getImageGeometry()->getProjection()->saveState(prjKwl);
	saveGcpFile(outPath, gcpSet, NULL, &prjKwl, false);
}