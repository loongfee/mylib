#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <gcpUtil.h>
#include <time.h>
#include "gdal_priv.h"
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>

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
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "mlpack.lib")

#include <mprojectdefine.h>
#include <ossim/base/ossimPreferences.h>
#include <ossim_plugin/radi/radiRpcSolver.h>
#include <ossim_plugin/radi/radiRpcModel.h>
#include <ossim_plugin/radi/radiBlockAdjustment.h>
#include <ossim_plugin/radi/radiBlockTieGptSet.h>
#include <ossim/parallel/ossimMultiThreadSequencer.h>
#include <func.h>

using namespace std;
using namespace mylib;
using namespace ossimplugins;

const char *pszBlockTieGptSetFile = "";
const char *pszProjectFile = "";
const char *pszPreferenceFile = "preference.txt";
const char *pszReportFile = "";
const char *pszOutputFolder = "";
int nThreads = 0;
bool bOverwrite = false;
bool bOverview = false;
bool bRpc = true;
bool bUseL1 = false;
radiBlockAdjustment::RobustMode robustModel = radiBlockAdjustment::RobustMode::NONE;

void outputResidue(const radiBlockAdjustment& ba, radiBlockTieGptSet* blockTieGptSet, const char* outputFilename,
				   ossimRefPtr<ossimMapProjection> mapProjection)
{
	FILE *pf = fopen(outputFilename, "w+");
	int counter = 0;
	for (vector<ossimRefPtr<radiBlockTieGpt> >::const_iterator iter =  blockTieGptSet->getTiePoints().begin();
		iter !=  blockTieGptSet->getTiePoints().end();++iter,++counter)
	{
		if (radiBlockTieGpt::point_type::unknown_tie_image_points == (*iter)->m_nType)
		{
			ossimGpt gpt;
			ba.lineSampleToWorld((*iter)->m_DptList[0], gpt);
			ossimDpt dpt1 = mapProjection->forward(gpt);
			ba.lineSampleToWorld((*iter)->m_DptList[1], gpt);
			ossimDpt dpt2 = mapProjection->forward(gpt);
			ossimDpt res = dpt1 - dpt2;
			double rmse = sqrt(res.x*res.x + res.y*res.y);
			fprintf(pf, "%4d%10.3lf%10.3lf%10.3lf%8s%8.3lf%5.03d%10.3lf%10.3lf%5.03d%10.3lf%10.3lf\n", counter, rmse, res.x, res.y, "TIE",
				(*iter)->score,
				(*iter)->m_DptList[0].first,
				(*iter)->m_DptList[0].second.x, (*iter)->m_DptList[0].second.y,
				(*iter)->m_DptList[1].first,
				(*iter)->m_DptList[1].second.x, (*iter)->m_DptList[1].second.y);
		}
		else
		{
			pair<int, ossimDpt> dptPair((*iter)->getSlaveId(), (*iter)->getImagePoint());
			ossimGpt gpt;
			ba.lineSampleToWorld(dptPair, gpt);
			ossimDpt dpt1 = mapProjection->forward(gpt);
			gpt = (*iter)->getGroundPoint();
			if (ossimElevManager::instance())
			{
				gpt.hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(gpt);
			}
			ossimDpt dpt2 = mapProjection->forward( (*iter)->getGroundPoint());
			ossimDpt res = dpt1 - dpt2;
			double rmse = sqrt(res.x*res.x + res.y*res.y);
			cout<<rmse<<endl;
			//fprintf(pf, "%4d%10.3lf%10.3lf%10.3lf%8s%10.3lf%10.3lf%10.3lf%10.3lf\n", counter, rmse, res.x, res.y, "CONTROL",
			//	(*iter)->m_DptList[0].second.x, (*iter)->m_DptList[0].second.y,
			//	(*iter)->m_DptList[1].second.x, (*iter)->m_DptList[1].second.y);
			const char* tp = NULL;
			if (radiBlockTieGpt::point_type::known_ground_control_points == (*iter)->m_nType)
			{
				tp = "CONTROL";
			}
			else
			{
				tp = "CHECK";
			}
			fprintf(pf, "%4d%10.3lf%10.3lf%10.3lf%8s%8.3lf%5.03d%10.3lf%10.3lf%5.03d%10.3lf%10.3lf\n", counter, rmse, res.x, res.y, tp,
				(*iter)->score,
				(*iter)->getSlaveId(),
				dpt1.x, dpt1.y,
				(*iter)->getMasterId(),
				dpt2.x, dpt2.y);
		}
	}

	fclose(pf);
}
//
//void block_test_l5()
//{
//	radiBlockAdjustment ba;
//	radiBlockTieGptSet* blockTieGptSet = new radiBlockTieGptSet;
//	ossimFilename blockTieGptSetFile = "I:\\testdata\\GF1-zhongya\\blockadjustment\\tie_points.xml";
//	ossimXmlDocument gmlDoc(blockTieGptSetFile);
//	ossimRefPtr<ossimXmlNode> aGmlNode = gmlDoc.getRoot();
//	blockTieGptSet->importFromGmlNode(aGmlNode);
//	int nImage = (int)blockTieGptSet->getImageList().size();
//	std::vector< ossimRefPtr<ossimSensorModel> > sensorModelList(nImage, NULL);
//	bool bRpc = true;
//	bool bUseL1 = false;
//	for (int i = 0;i < nImage;++i)
//	{
//		int imageId = blockTieGptSet->getImageList()[i].first;
//		ossimFilename imageFile = blockTieGptSet->getImageList()[i].second;
//		ossimRefPtr<ossimImageHandler> handler   = ossimImageHandlerRegistry::instance()->open(imageFile);
//		if(!handler)	
//		{
//			cerr<<"Failed to open image: "<<imageFile<<endl;
//			return;
//		}
//		if (bRpc)
//		{
//			// 强制使用RPC模型
//			ossimplugins::radiRpcModel *rpcModel = new ossimplugins::radiRpcModel;
//			if(!rpcModel->parseRpcFile(imageFile))
//			{
//				//
//				cout<<"Failed to find rpc file for rpc model."<<endl;
//				return;
//			}
//			rpcModel->setUseL1(bUseL1);
//			sensorModelList[i] = PTR_CAST(ossimSensorModel, rpcModel);
//		}
//		else
//		{
//			sensorModelList[i] = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());
//		}
//		//sensorModelList[i] = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());
//		ba.addSensorModel(imageId, sensorModelList[i].get());
//	}
//	ba.adjustment(*blockTieGptSet, robustModel);
//
//	ossimKeywordlist prjKwl;
//	ossimFilename projectionFile(pszProjectFile);
//	if (projectionFile.exists())
//	{
//		prjKwl.clear();
//		prjKwl.addFile(projectionFile);
//	}
//	ossimRefPtr<ossimMapProjection> mapProjection = PTR_CAST(ossimMapProjection,
//		ossimMapProjectionFactory::instance()->createProjection(prjKwl));
//	const char* outputFilename = "I:\\testdata\\GF1-zhongya\\blockadjustment\\report.txt";
//	outputResidue(ba, blockTieGptSet, outputFilename,  mapProjection);
//
//	for (int i = 0;i < nImage;++i)
//	{
//		int imageId = blockTieGptSet->getImageList()[i].first;
//		ossimFilename modelFile = blockTieGptSet->getImageList()[i].second;
//		modelFile.setExtension("geom");
//		ossimKeywordlist model_geom;
//		sensorModelList[i]->saveState(model_geom);
//		fstream fs;
//		fs.open(modelFile.c_str(), ios_base::out);
//		fs<<model_geom;
//		fs.close();
//	}
//	return;
//
//	vector<ossim_uint32> outBandList;
//	for (int i = 0;i < nImage;++i)
//	{
//		int imageId = blockTieGptSet->getImageList()[i].first;
//		ossimFilename imageFile = blockTieGptSet->getImageList()[i].second;
//		ossimFilename outputImageFile = ossimFilename(pszOutputFolder) + "\\o" + imageFile.fileNoExtension() + ".tif";
//
//
//		ossimRefPtr<ossimImageHandler> handler   = ossimImageHandlerRegistry::instance()->open(imageFile);
//		if(!handler)	
//		{
//			cerr<<"Failed to open image: "<<imageFile<<endl;
//			return;
//		}
//
//		sensorModelList[i]->m_proj = PTR_CAST(ossimMapProjection,
//			ossimMapProjectionFactory::instance()->createProjection(prjKwl));
//
//		ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
//		imageGeom->setProjection(sensorModelList[i].get());
//		handler->setImageGeometry(imageGeom.get());//1128
//		//cout<<geom<<endl;
//		ossimKeywordlist tt_geom;
//
//		ossimDpt imagesize(handler->getImageRectangle().width(), handler->getImageRectangle().height());
//		ossimImageRenderer* renderer = new ossimImageRenderer;
//		ossimPolyCutter* theCutter;
//		ossimBandSelector* theBandSelector;
//		vector<ossimDpt> polygon;
//		ossimIrect bound,boundw;
//		theCutter = new ossimPolyCutter;
//		///////////////////以下四行应该从界面取得//
//		int startpixel = 0;
//		int startline = 0;
//		int endpixel = 0;
//		int endline = 0;
//		if (endpixel == 0)
//		{
//			endpixel = imagesize.x - 1;
//		}
//		if (endline == 0)
//		{
//			endline = imagesize.y - 1;
//		}
//		theBandSelector = new ossimBandSelector;
//		theBandSelector->connectMyInputTo(0, handler.get());
//		if (outBandList.empty())
//		{
//			int nBands = handler->getNumberOfInputBands();
//			for (int i=0;i<nBands;++i)
//			{
//				outBandList.push_back(i);
//			}
//		}
//
//		theBandSelector->setOutputBandList(outBandList);
//
//
//		ossimDpt ps(startpixel,startline),p2(endpixel,startline),p3(endpixel,endline),p4(startpixel,endline),p5(startpixel,startline);
//
//		polygon.push_back(ps);
//		polygon.push_back(p2);
//		polygon.push_back(p3);
//		polygon.push_back(p4);
//		polygon.push_back(p5);
//		theCutter->connectMyInputTo(theBandSelector);
//		theCutter->setPolygon(polygon);
//		theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
//		theCutter->setNumberOfPolygons(1);
//
//		renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
//		ossimRefPtr<ossimImageFileWriter> writer;
//		{
//			writer = ossimImageWriterFactoryRegistry::instance()->
//				createWriterFromExtension( outputImageFile.ext() );
//			//writer->setOverviewCompressType(COMPRESSION_JPEG);
//			//writer->setOverviewCompressType(COMPRESSION_LZW);
//			//ossimRefPtr<ossimProperty> prop =  new ossimStringProperty("gdal_overview_type", "nearest", false);
//			//writer->setProperty(prop);
//		}
//		if (bOverview)
//		{
//			writer->setWriteOverviewFlag(true);
//		}
//		if(writer==NULL) return;
//		tt_geom.clear();
//		writer->saveState(tt_geom);
//		tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA", true);
//		writer->loadState(tt_geom);
//
//		renderer->connectMyInputTo(theCutter);
//		renderer->setView(sensorModelList[i]->m_proj);
//
//		myOutProgress *m_progress = new myOutProgress(0,true);
//		writer->addListener(m_progress);
//
//		writer->setFilename(outputImageFile);
//
//		if (nThreads < 1)
//		{
//			nThreads = OpenThreads::GetNumberOfProcessors() * 2;
//		}
//		cout<<"using "<<nThreads<<" threads:"<<endl;
//		ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
//		writer->changeSequencer(sequencer.get());
//		writer->connectMyInputTo(0,renderer);
//		//ossimIrect bounding = handler->getBoundingRect();
//		//writer->setAreaOfInterest(bounding);
//
//		//bool bResult = writer->execute();
//
//		if(!writer->execute())
//		{
//			writer->removeListener(m_progress);
//			writer->disconnectAllInputs();
//			renderer->disconnectAllInputs();
//			theCutter->disconnectAllInputs();
//			//delete handler;
//			//handler->close();
//			//delete handler;
//			return;
//		}
//
//		writer->disableListener();
//		writer->removeListener(m_progress);
//		writer->disconnectAllInputs();
//		renderer->disconnectAllInputs();
//		theCutter->disconnectAllInputs();
//	}
//
//	cout<<"OK"<<endl;
//	//int i;
//	//for (i = 0;i < num;++i)
//	//{
//	//	InitializePrj(prj[i], sourceFile[i], demPath[i], chkFile[i]); 
//	//	ba.addSensorModel(prj[i].m_sensorModel);
//	//}
//
//	//vector< radiBlockTieGpt > gptList;
//	//// 加载（已知）控制点
//	//AddGpts("D:\\workspace\\testdata\\BlockAdjustment\\gpt.txt", gptList);
//	//// 加载同名点
//	//AddCpts("D:\\workspace\\testdata\\BlockAdjustment\\cpt.txt", gptList);
//	//// 进行区域网平差
//	//ba.adjustment(gptList);
//	//for(i = 0;i < num;i++)
//	//{
//	//	prj[i].m_sensorModel->saveState(prj[i].geom);
//	//}
//
//	//// 用优化后的模型输出检查点精度
//	//fstream fs;
//	//fs.open("D:\\workspace\\testdata\\BlockAdjustment\\residue.txt",ios_base::out);
//	//fs.setf(ios::fixed, ios::floatfield);
//	//fs.precision(6);
//	//for (i = 0;i < num;++i)
//	//{
//	//	AppendResidue2File(prj[i].m_sensorModel, prj[i].m_CtrlGptSet, fs);
//	//	fs<<endl;
//	//}
//	//fs.close();
//
//	//for(i = 0;i < num;i++)
//	//{
//	//	prj[i].m_OutBandList.clear();
//	//	prj[i].m_OutBandList.push_back(6);
//	//	prj[i].m_OutBandList.push_back(3);
//	//	prj[i].m_OutBandList.push_back(1);
//	//	prj[i].Orthograph(outFile[i]);
//	//}
//}

void ba()
{
	radiBlockAdjustment ba;
	radiBlockTieGptSet* blockTieGptSet = new radiBlockTieGptSet;
	ossimFilename blockTieGptSetFile(pszBlockTieGptSetFile);
	ossimXmlDocument gmlDoc(blockTieGptSetFile);
	ossimRefPtr<ossimXmlNode> aGmlNode = gmlDoc.getRoot();
	blockTieGptSet->importFromGmlNode(aGmlNode);
	int nImage = (int)blockTieGptSet->getImageList().size();
	std::vector< ossimRefPtr<ossimSensorModel> > sensorModelList(nImage, NULL);
	for (int i = 0;i < nImage;++i)
	{
		int imageId = blockTieGptSet->getImageList()[i].first;
		ossimFilename imageFile = blockTieGptSet->getImageList()[i].second;
		ossimRefPtr<ossimImageHandler> handler   = ossimImageHandlerRegistry::instance()->open(imageFile);
		if(!handler)	
		{
			cerr<<"Failed to open image: "<<imageFile<<endl;
			return;
		}
		if (bRpc)
		{
			// 强制使用RPC模型
			ossimplugins::radiRpcModel *rpcModel = new ossimplugins::radiRpcModel;
			if(!rpcModel->parseRpcFile(imageFile))
			{
				//
				cout<<"Failed to find rpc file for rpc model."<<endl;
				return;
			}
			rpcModel->setUseL1(bUseL1);
			sensorModelList[i] = PTR_CAST(ossimSensorModel, rpcModel);
		}
		else
		{
			sensorModelList[i] = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());
		}
		//sensorModelList[i] = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());
		ba.addSensorModel(imageId, sensorModelList[i].get());
	}
	//ba.adjustment(*blockTieGptSet);
	ba.adjustment(*blockTieGptSet, robustModel);
	//ba.adjustment(*blockTieGptSet, radiBlockAdjustment::RobustMode::Bisquare);

	if (0 != strcmp(pszReportFile, ""))
	{
		ossimKeywordlist prjKwl;
		ossimFilename projectionFile(pszProjectFile);
		if (projectionFile.exists())
		{
			prjKwl.clear();
			prjKwl.addFile(projectionFile);
		}
		ossimRefPtr<ossimMapProjection> mapProjection = PTR_CAST(ossimMapProjection,
			ossimMapProjectionFactory::instance()->createProjection(prjKwl));
		outputResidue(ba, blockTieGptSet, pszReportFile,  mapProjection);
	}

	for (int i = 0;i < nImage;++i)
	{
		int imageId = blockTieGptSet->getImageList()[i].first;
		ossimFilename modelFile = blockTieGptSet->getImageList()[i].second;
		modelFile.setExtension("model");
		ossimKeywordlist model_geom;
		sensorModelList[i]->saveState(model_geom);
		fstream fs;
		fs.open(modelFile.c_str(), ios_base::out);
		fs<<model_geom;
		fs.close();
	}
}

//int main( int argc, char** argv )
//{
//	ossimFilename preferences_file = ossimFilename("pszPreferenceFile");
//	if (preferences_file.exists())
//	{
//		ossimPreferences::instance()->loadPreferences(preferences_file);
//	}
//	else
//	{
//		preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//		if (preferences_file.exists())
//		{
//			ossimPreferences::instance()->loadPreferences(preferences_file);
//		}
//	}
//	std::string	tempString;
//	ossimArgumentParser::ossimParameter	stringParam(tempString);
//	ossimArgumentParser argumentParser(&argc, argv);
//	ossimInit::instance()->addOptions(argumentParser);
//	argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//	ossimInit::instance()->initialize(argumentParser);
//	ossimElevManager::instance()->initialize();
//	pszProjectFile = "I:\\testdata\\GF1-zhongya\\MSS\\projection_mss.txt";
//	pszOutputFolder = "I:\\testdata\\GF1-zhongya\\MSS\\orth";
//	block_test_l5();
//	return 0;
//}

/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: ba -i blockTiePointSetFile [-report reportFile] [-overwrite] \n"
		"\t[-prj projectionfile] [-pref preferencefile] [-rm robustmodel]\n"
		"\t[-rpc] [-l1 {0|1}]\n");

	printf( 
		"-rm robustmodel: NONE(default), HUBER, BISQUARE\n");

	if( pszErrorMsg != NULL )
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

int main( int argc, char** argv )
{
	GDALAllRegister();
	//OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持
	
	if (argc > 1)
	{
		/* -------------------------------------------------------------------- */
		/*      Parse arguments.                                                */
		/* -------------------------------------------------------------------- */
		for( int i = 1; i < argc; i++ )
		{
			if( 0 == _stricmp(argv[i],"-i") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszBlockTieGptSetFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-report") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszReportFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-rpc") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
				bRpc = true ;
			}
			else if( 0 == _stricmp(argv[i],"-l1") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				int flag = atoi( argv[++i] );	
				if (0 == flag)
				{
					bUseL1 = false;
				}
				else
				{
					bUseL1 = true;
				}
			}
			else if( 0 == _stricmp(argv[i],"-prj") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszProjectFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-pref") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszPreferenceFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-overwrite") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
				bOverwrite = true ;
			}
			else if( 0 == _stricmp(argv[i],"-rm") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				string rm = argv[++i];
				std::transform(rm.begin(), rm.end(),rm.begin(), ::toupper);
				if (0 == strcmp(rm.c_str(), "NONE"))
				{
					robustModel = radiBlockAdjustment::RobustMode::NONE;
				}
				else if(0 == strcmp(rm.c_str(), "HUBER"))
				{
					robustModel = radiBlockAdjustment::RobustMode::HUBER;
				}
				else if(0 == strcmp(rm.c_str(), "BISQUARE"))
				{
					robustModel = radiBlockAdjustment::RobustMode::BISQUARE;
				}
				else
				{
					cerr<<"unknown robust model:"<<rm<<endl;
					robustModel = radiBlockAdjustment::RobustMode::NONE;
				}
			}
			else
			{
				Usage();
			}
		}

		if (0 == strcmp(pszBlockTieGptSetFile, ""))
		{
			printf("blockTiePointSetFile can not be empty!\n");
			Usage();
		}
		else
		{
			ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
			if (preferences_file.exists())
			{
				ossimPreferences::instance()->loadPreferences(preferences_file);
			}
			else
			{
				preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
				if (preferences_file.exists())
				{
					ossimPreferences::instance()->loadPreferences(preferences_file);
				}
			}

			std::string	tempString;
			ossimArgumentParser::ossimParameter	stringParam(tempString);
			ossimArgumentParser argumentParser(&argc, argv);
			ossimInit::instance()->addOptions(argumentParser);
			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
			ossimInit::instance()->initialize(argumentParser);
			//pszInputFile = "I:\\GF1-Testdata\\1\\PAN\\GF1_PMS1_E76.2_N42.5_20130918_L1A0000085664-PAN1.tiff";
			//pszGcpFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\gcp.txt";
			//pszOutFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\rpb_orth.tif";

			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
			ba();
			clockEnd = clock();
			printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
	}
	else
	{
		ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
		if (preferences_file.exists())
		{
			ossimPreferences::instance()->loadPreferences(preferences_file);
		}
		else
		{
			preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
			if (preferences_file.exists())
			{
				ossimPreferences::instance()->loadPreferences(preferences_file);
			}
		}

		std::string	tempString;
		ossimArgumentParser::ossimParameter	stringParam(tempString);
		ossimArgumentParser argumentParser(&argc, argv);
		ossimInit::instance()->addOptions(argumentParser);
		argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
		ossimInit::instance()->initialize(argumentParser);

		// 这里是调试时使用
		pszBlockTieGptSetFile = "I:\\testdata\\GF1\\xinjiang-jierjisisitan\\blockadjustment\\tie_points.xml";
		pszReportFile = "I:\\testdata\\GF1\\xinjiang-jierjisisitan\\blockadjustment\\report.txt";
		pszProjectFile = "I:\\testdata\\GF1\\xinjiang-jierjisisitan\\data\\projection_mss.txt";
		bOverwrite = true;

		ba();
		Usage(0);
	}
	return 0;
}