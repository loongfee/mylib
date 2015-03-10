#include "radiImageRegistration.h"
#include <ossim/base/ossimProcessInterface.h>
#include <ossim/base/ossimStdOutProgress.h>
#include <ossim/init/ossimInit.h>
#include <ossim/base/ossimPreferences.h>
#include <gcpUtil.h>
#include <time.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>
#include <ossim/imaging/ossimImageGeometry.h>

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
//#pragma comment(lib, "ossimopencv_plugin.lib")
#pragma comment(lib, "OpenThreads.lib")
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")
#pragma comment(lib, "opencv_ts300.lib")
#pragma comment(lib, "opencv_world300.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "mlpack.lib")
#pragma comment(lib, "libfftw-3.3.lib")

using namespace mylib;

const char* pszSourceFile = "";
const char* pszReferenceFile = "";
const char* pszGcpFile = "";
const char *pszProjectFile = "";
const char *pszPreferenceFile = "preference.txt";

// master band
int mBand = 0;
// slave band
int sBand = 0;
// accuracy: pixels
double sAccuracy = 50;
// tile size
int tileSize = 64;
// required point number
int nPoint = 64;
// area of interesting: upper left point and lower right point (pixels)
int aoi_x0 = 0;
int aoi_x1 = 0;
int aoi_y0 = 0;
int aoi_y1 = 0;
// number of threads
int nThreads = 0;
// type of points: control points or tie points
radiImageRegistration::point_type pointType = radiImageRegistration::point_type::control;
// whether to use geographic coordinates if the output projection is not set
bool bGeographic = false;
int slaveId = 1;
int masterId = 2;
bool bDebug = true;
bool bAppend = false;

void matchPoint(ossimFilename sourcePath, ossimFilename stdPath, ossimFilename outPath)
{
	//ossimInit::instance()->initialize();	
	bool result;
	ossimKeywordlist geomp;
	vector<ossimFilename> inputfile;
	ossimFilename master,slave,file_tmp,outimgfile;

	radiImageRegistration *img_reg = new radiImageRegistration;
	ossimStdOutProgress progress(0,true);
	img_reg->addListener(&progress);

	//ossimString master_band = "0";
	//ossimString slave_band = "2";
	//ossimString slave_accuracy = "50";
	//ossimString point_number = "100";
	//ossimString tile_size = "64";

	//ossimString sift_nfeatures = "0";
	//ossimString sift_noctavelayers = "3";
	//ossimString sift_contrastthreshold = "0.01";
	//ossimString sift_edgethreshold = "10";
	//ossimString sift_sigma = "1.6";

	img_reg->setMasterBand(mBand);
	img_reg->setSlaveBand(sBand);
	img_reg->setSlaveAccuracy(sAccuracy);
	img_reg->setPointNumber(nPoint);
	img_reg->setSlave(sourcePath);
	img_reg->setMaster(stdPath);
	img_reg->setTileSize(tileSize);
	img_reg->setThreadNum(nThreads);
	img_reg->setPointType(pointType);
	img_reg->setUseGeographic(bGeographic);
	img_reg->setSlaveId(slaveId);
	img_reg->setMasterId(masterId);
	img_reg->setDebug(bDebug);
	//img_reg->setPointType(radiImageRegistration::point_type::tie);

	//img_reg->setSiftNfeatures(0);
	//img_reg->setSiftNOctaveLayers(3);
	//img_reg->setSiftContrastThreshold(0.01);
	//img_reg->setSiftEdgeThreshold(10);
	//img_reg->setSiftSigma(1.6);

	img_reg->setSiftNfeatures(0);
	img_reg->setSiftNOctaveLayers(2);
	img_reg->setSiftContrastThreshold(0.001);
	img_reg->setSiftEdgeThreshold(10);
	img_reg->setSiftSigma(1.6);

	img_reg->setOutputName(outPath);
	img_reg->setAreaOfInterest(ossimIrect(aoi_x0, aoi_y0, aoi_x1, aoi_y1));


	result = img_reg->execute();

	//img_reg->appendControlPoints(outPath);
	//img_reg->appendTiePoints(outPath);
	//ossimRefPtr<ossimImageHandler> handler = ossimImageHandlerRegistry::instance()->open(stdPath);
	//ossimKeywordlist outPrj;
	//img_reg->writePoints(outPath, PTR_CAST(ossimMapProjection, handler->getImageGeometry()->getProjection()));
	ossimImageHandler* handlerM = ossimImageHandlerRegistry::instance()->open(stdPath);
	ossimMapProjection* prj = PTR_CAST(ossimMapProjection, handlerM->getImageGeometry()->getProjection());
	prj = NULL;
	if (bAppend)
	{
		img_reg->appendPoints(outPath, prj);
	}
	else
	{
		img_reg->writePoints(outPath, prj);
	}

	delete img_reg;
	img_reg = NULL;

	//ossimTieGptSet* gcpSet = new ossimTieGptSet;
	//ossimTieGptSet* chkSet = new ossimTieGptSet;
	//ossimKeywordlist prjKwl;
	//readGcpFile("tmp.xml", gcpSet, chkSet, &prjKwl);

	//ossimImageHandler* handlerM = ossimImageHandlerRegistry::instance()->open(stdPath);
	//handlerM->getImageGeometry()->getProjection()->saveState(prjKwl);
	//saveGcpFile(outPath, gcpSet, NULL, &prjKwl, false);
}

//#include <ossim/base/ossimProcessInterface.h>
//#include <ossim/base/ossimObjectFactoryRegistry.h>
//void matchPoint1(ossimFilename sourcePath, ossimFilename stdPath, ossimFilename outPath)
//{
//	bool result;
//	ossimKeywordlist geomp;
//	vector<ossimFilename> inputfile;
//	ossimFilename master,slave,file_tmp,outimgfile;
//
//
//	ossimRefPtr<ossimObject> icObject = ossimObjectFactoryRegistry::instance()->createObject(ossimString("ossimImageCorrelator"));//自定义自动选取同名点的类，为dll加载形式。
//	ossimOutputSource* icSource = PTR_CAST(ossimOutputSource, icObject.get());
//	ossimStdOutProgress progress(0,true);
//	icSource->addListener(&progress);
//	ossimProcessInterface* icProcessInterface = PTR_CAST(ossimProcessInterface, icObject.get());
//	ossimPropertyInterface* icPropertyInterface = PTR_CAST(ossimPropertyInterface, icObject.get());
//	ossimRefPtr<ossimProperty> masterBand        = icSource->getProperty("master_band");//参考影像选用的波段
//	ossimRefPtr<ossimProperty> slaveBand         = icSource->getProperty("slave_band");//待校正影像选用的波段
//	ossimRefPtr<ossimProperty> scaleRatio        = icSource->getProperty("scale_ratio");//尺度因子
//	ossimRefPtr<ossimProperty> cornerDensity     = icSource->getProperty("corner_density");//harriscorner算子选取点的密度
//	ossimRefPtr<ossimProperty> minCorrel         = icSource->getProperty("min_correl");//最小相关度
//	ossimRefPtr<ossimProperty> templateRadius    = icSource->getProperty("template_radius");//模板半径
//	ossimRefPtr<ossimProperty> slaveAccuracy     = icSource->getProperty("slave_accuracy");//待校正影像的误差
//	ossimRefPtr<ossimProperty> projectionType    = icSource->getProperty("projection_type");//采用的投影方式
//	ossimRefPtr<ossimProperty> outputFilename    = icSource->getProperty("output_filename");//输出同名点txt文件
//
//	ossimString master_band = "2";
//	ossimString slave_band = "2";
//	ossimString slave_accuracy = "50";
//	ossimString point_number = "100";
//	ossimString tile_size = "64";
//
//	ossimString sift_nfeatures = "0";
//	ossimString sift_noctavelayers = "3";
//	ossimString sift_contrastthreshold = "0.01";
//	ossimString sift_edgethreshold = "10";
//	ossimString sift_sigma = "1.6";
//
//	icPropertyInterface->setProperty("master_band", master_band);
//	icPropertyInterface->setProperty("slave_band", slave_band);
//	icPropertyInterface->setProperty("slave_accuracy", slave_accuracy);
//	icPropertyInterface->setProperty("point_number", point_number);
//	icPropertyInterface->setProperty("slave_filename", sourcePath);
//	icPropertyInterface->setProperty("master_filename", stdPath);
//	icPropertyInterface->setProperty("template_radius", "64");
//	icPropertyInterface->setProperty("min_correl", "0.70");
//
//	//icPropertyInterface->setProperty("sift_nfeatures", sift_nfeatures);
//	//icPropertyInterface->setProperty("sift_noctavelayers", sift_noctavelayers);
//	//icPropertyInterface->setProperty("sift_contrastthreshold", sift_contrastthreshold);
//	//icPropertyInterface->setProperty("sift_edgethreshold", sift_edgethreshold);
//	//icPropertyInterface->setProperty("sift_sigma", sift_sigma);
//
//
//	icPropertyInterface->setProperty("output_filename", "tmp.xml");
//	//MPI_Comm_rank(MPI_COMM_WORLD,&myid);
//	//MPI_Get_processor_name(processor_name, &namelen);
//	result = icProcessInterface->execute();
//
//	icObject=NULL;
//
//	ossimTieGptSet* gcpSet = new ossimTieGptSet;
//	ossimTieGptSet* chkSet = new ossimTieGptSet;
//	ossimKeywordlist prjKwl;
//	readGcpFile("tmp.xml", gcpSet, chkSet, &prjKwl);
//
//	ossimImageHandler* handlerM = ossimImageHandlerRegistry::instance()->open(stdPath);
//	handlerM->getImageGeometry()->getProjection()->saveState(prjKwl);
//	saveGcpFile(outPath, gcpSet, NULL, &prjKwl, false);
//}

//int main()
//{
//	string dropbox = getenv("DROPBOX");
//	string preferences_file = string(getenv("OSSIM2_DIR")) + "\\preference.txt";
//
//	ossimPreferences::instance()->loadPreferences(ossimFilename(preferences_file));
//	ossimInit::instance()->initialize();
//
//	clock_t  clockBegin, clockEnd;
//	clockBegin = clock();
//
//	ossimElevManager::instance()->loadElevationPath("D:\\workspace\\dem\\srtm90");
//	//matchPoint("E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\IMAGE.TIF", 
//	//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\reference\\referenceL8.tif", 
//	//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\ossim_matching.txt");
//	//matchPoint("E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01\\IMAGE.tif", 
//	//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\111\\1.pix", 
//	//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01\\ossim_matching.txt");
//	matchPoint("E:\\HJ1\\HJ-1B_CCD-1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\IMAGE.tif", 
//		"E:\\HJ1\\HJ-1B_CCD-1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\reference\\referenceL8.tif", 
//		"E:\\HJ1\\HJ-1B_CCD-1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\ossim_matching.txt");
//	//matchPoint("D:\\workspace\\Spot5\\beijing\\281268-20040602-2.5\\scene01\\imagery.tif", 
//	//"D:\\workspace\\Spot5\\beijing\\281268_20040602.tif", "D:\\workspace\\Spot5\\beijing\\ossim-matching.txt");
//
//	//matchPoint("I:\\20140523\\MSS\\GF1_PMS1_E100.8_N24.4_20140202_L1A0000160810-MSS1.tiff", 
//	//	"I:\\20140523\\fusion\\GF1_PMS1_E100.8_N24.4_20140202_L1A0000160810-PAN1.tif", 
//	//	"I:\\20140523\\MSS\\ossim_matching.txt");
//
//
//	//matchPoint("D:\\workspace\\Landsat\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat", 
//	//	"D:\\workspace\\Landsat\\146030_UTM.tif", "D:\\workspace\\Landsat\\ossim_matching.txt");
//
//	//matchPoint1("E:\\HJ1\\HJ-1B_CCD-1_MYC_201405190109_201405190117\\Scene01\\IMAGE_B1.tif", 
//	//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201405190109_201405190117\\Scene01\\IMAGE_B3.tif", 
//	//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201405190109_201405190117\\Scene01\\1-3.txt");
//	//matchPoint1("E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\IMAGE.TIF", 
//	//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\reference\\referenceL8.tif", 
//	//	"E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\ossim_matching.txt");
//	//matchPoint1("D:\\workspace\\Landsat\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat", 
//	//	"D:\\workspace\\Landsat\\146030_UTM.tif", "D:\\workspace\\Landsat\\ossim_matching.txt");
//	clockEnd = clock();
//	printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
//	return 1;
//}


/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: img-reg -i inputFile -r refFile [-gcp gcpFile]\n"
		"\t[-prj projectionfile] [-pref preferencefile] [-nt nThreads]\n"
		"\t[-mb mband] [-sb sband] [-ts tilesize] [-np npoint] [-accuracy tolerance]"
		"\t[-pt {\"control\"|\"tie\"}] [-ll {0|1}] [-debug {0|1}] [-append]\n");

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
	
	int null_value = 0;

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
				pszSourceFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-r") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszReferenceFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-gcp") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszGcpFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-prj") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszProjectFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-mb") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				mBand = atoi(argv[++i]) - 1 ;
			}
			else if( 0 == _stricmp(argv[i],"-sb") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				sBand = atoi(argv[++i]) - 1 ;
			}
			else if( 0 == _stricmp(argv[i],"-ts") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				tileSize = atoi(argv[++i]) ;
			}
			else if( 0 == _stricmp(argv[i],"-np") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				nPoint = atof(argv[++i]) ;
			}
			else if( 0 == _stricmp(argv[i],"-accuracy") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				sAccuracy = atof(argv[++i]) ;
			}
			else if( 0 == _stricmp(argv[i],"-aoi") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(4);
				aoi_x0 = atoi(argv[++i]) ;
				aoi_y0 = atoi(argv[++i]) ;
				aoi_x1 = atoi(argv[++i]) ;
				aoi_y1 = atoi(argv[++i]) ;
			}
			else if(0 == _stricmp(argv[i],"-nt") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				nThreads = atoi( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-pt") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				const char* strPointType = argv[++i];
				if (0 == _stricmp(strPointType,"control"))
				{
					pointType = radiImageRegistration::point_type::control;
				}
				else if (0 == _stricmp(strPointType,"tie"))
				{
					pointType = radiImageRegistration::point_type::tie;
				}
				else
				{
					cerr<<"the value of \"-pt\" option can only be \"control\" or \"tie\""<<endl;
					Usage();
				}
			}
			else if(0 == _stricmp(argv[i],"-ll") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				int ll = atoi( argv[++i] );
				if (0 == ll)
				{
					bGeographic = false;
				}
				else
				{
					bGeographic = true;
				}
			}
			else if(0 == _stricmp(argv[i],"-sid") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				slaveId = atoi(argv[++i]);	
			}
			else if(0 == _stricmp(argv[i],"-mid") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				masterId = atoi(argv[++i]);	
			}
			else if(0 == _stricmp(argv[i],"-append") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
				bAppend = true;	
			}
			else if(0 == _stricmp(argv[i],"-debug") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				if (atoi(argv[++i]) == 0)
				{
					bDebug = false;
				}
				else
				{
					bDebug = true;
				}
			}
			else if( 0 == _stricmp(argv[i],"-pref") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszPreferenceFile = argv[++i] ;
			}
			else
			{
				Usage();
			}
		}

		if (0 == strcmp(pszSourceFile, ""))
		{
			printf("input file can not be empty!\n");
			Usage();
		}
		else if (0 == strcmp(pszReferenceFile, ""))
		{
			printf("reference file can not be empty!\n");
			Usage();
		}
		else
		{
			ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
			if (preferences_file.exists())
			{
				//cout<<"found1"<<endl;
				ossimPreferences::instance()->loadPreferences(preferences_file);
			}
			else
			{
				preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
				if (preferences_file.exists())
				{
					//cout<<"found2"<<endl;
					ossimPreferences::instance()->loadPreferences(preferences_file);
				}
			}
			std::string	tempString;
			ossimArgumentParser::ossimParameter	stringParam(tempString);
			ossimArgumentParser argumentParser(&argc, argv);
			ossimInit::instance()->addOptions(argumentParser);
			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
			ossimInit::instance()->initialize(argumentParser);
			ossimElevManager::instance()->initialize();

			ossimFilename sourceFile(pszSourceFile);
			ossimFilename referenceFile(pszReferenceFile);
			ossimFilename gcpFile(pszGcpFile);			
			if (0 == strcmp(pszGcpFile, ""))
			{
				gcpFile = sourceFile.expand().path() + "\\gcp_auto.txt";
			}

			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
			matchPoint(sourceFile, referenceFile, gcpFile);
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
		ossimElevManager::instance()->initialize();

		//pszSourceFile = "E:\\HJ1\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201405250250_201405250300\\Scene07\\IMAGE_B3.tif";
		pszSourceFile = "I:\\testdata\\GF1\\xinjiang-jierjisisitan\\data\\GF1_PMS2_E76.9_N40.7_20140626_L1A0000260221-MSS2.tiff";
		pszSourceFile = "D:\\workspace\\LD2010003816\\header.dat";
		//pszReferenceFile = "D:\\workspace\\fuhai\\FCGC600205282\\IMG_PHR1A_P_001\\rect.tif";
		//pszReferenceFile = "I:\\testdata\\std\\index.shp";
		//pszReferenceFile = "I:\\testdata\\tmchina\\china.vrt";
		pszReferenceFile = "G:\\testdata\\GF1\\xinjiang-jierjisisitan\\data\\GF1_PMS2_E77.0_N41.0_20140806_L1A0000296787-MSS2.tiff";
		pszReferenceFile = "D:\\workspace\\LD2010003816\\TM_123036.TIF";
		pszGcpFile = "G:\\testdata\\GF1\\xinjiang-jierjisisitan\\GCP_vrt.xml";
		pszGcpFile = "D:\\workspace\\LD2010003816\\test.txt";

		//pszSourceFile = "D:\\workspace\\Landsat\\feature_test\\LS5_TM_20100717_050010_050035_145033_FASTB_L2\\header.dat";
		//pszReferenceFile = "D:\\workspace\\Landsat\\feature_test\\LS5_TM_20080625_000000_000000_145033_GEOTIFF_L4\\L5-TM-145-033-20080625-L4.TIF";
		//pszGcpFile = "D:\\workspace\\Landsat\\feature_test\\100.txt";
		//pszSourceFile = "E:\\HJ1\\HJ-L2\\HJ\\2010.tif";
		////pszReferenceFile = "D:\\workspace\\fuhai\\FCGC600205282\\IMG_PHR1A_P_001\\rect.tif";
		//pszReferenceFile = "E:\\HJ1\\HJ-L2\\TM\\2010.tif";

		//pszSourceFile = "D:\\workspace\\Spot5\\tianjin\\282270-20040518-2.5\\scene01\\imagery.tif";
		//pszReferenceFile = "D:\\workspace\\Spot5\\tianjin\\282270_20040518.tif";
		//pszGcpFile = "D:\\workspace\\Spot5\\tianjin\\gcp_auto.txt";

		//pszSourceFile = "D:\\workspace\\HJ\\西北\\b3.TIF";
		//pszReferenceFile = "G:\\testdata\\tmchina\\china.vrt";
		////pszReferenceFile = "D:\\workspace\\HJ\\西北\\HJ1B-CCD1-38-68-20100831-L20000384864-3.TIF";
		//pszGcpFile = "D:\\workspace\\HJ\\西北\\gcp_auto.txt";

		//pszSourceFile = "E:\\Kuaipan\\Programs\\发给龙师兄\\L.tif";
		//pszReferenceFile = "E:\\Kuaipan\\Programs\\发给龙师兄\\R.tif";
		//pszGcpFile = "E:\\Kuaipan\\Programs\\发给龙师兄\\gcp_auto.txt";

		pszSourceFile = "D:\\workspace\\LD2010003816\\header.dat";
		//pszReferenceFile = "D:\\workspace\\LD2010003816\\TM_123036.TIF";
		pszReferenceFile = "G:\\testdata\\china\\china2005-20141115.tif";
		pszGcpFile = "D:\\workspace\\LD2010003816\\auto.txt";

		sBand = 0;
		mBand = 0;
		pointType = radiImageRegistration::point_type::control;
		//pointType = radiImageRegistration::point_type::tie;
		nPoint = 25;
		nThreads = 1;
		bAppend = false;
		bDebug = true;
		sAccuracy = 30;
		sBand = 3;
		mBand = 0;
		tileSize = 128;
		bGeographic = true;
		////tileSize = 256;
		//nPoint = 100;
		//aoi_x0 = 4500;
		//aoi_y0 = 0;
		//aoi_x1 = 4600;
		//aoi_y1 = 11999;
		
		ossimFilename sourceFile(pszSourceFile);
		ossimFilename referenceFile(pszReferenceFile);
		ossimFilename gcpFile(pszGcpFile);
		if (0 == strcmp(pszGcpFile, ""))
		{
			gcpFile = sourceFile.path() + "\\gcp_auto.xml";
		}

		clock_t  clockBegin, clockEnd;
		clockBegin = clock();
		matchPoint(sourceFile, referenceFile, gcpFile);
		clockEnd = clock();
		printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);

		Usage(0);
	}
	return 0;
}