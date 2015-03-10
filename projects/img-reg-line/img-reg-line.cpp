#include "radiLineRegistration.h"
#include <ossim/base/ossimProcessInterface.h>
#include <ossim/base/ossimStdOutProgress.h>
#include <ossim/init/ossimInit.h>
#include <ossim/base/ossimPreferences.h>
#include <gcpUtil.h>
#include <time.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>
#include <ossim/imaging/ossimImageGeometry.h>

#include "lineConstant.h"


#pragma comment(lib, "ossim20.lib")
#pragma comment(lib, "blas_win32_MT.lib")
#pragma comment(lib, "lapack_win32_MT.lib")

#pragma comment(lib, "ossim_plugin.lib")
#pragma comment(lib, "OpenThreads.lib")
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")
#pragma comment(lib, "opencv_ts300.lib")
#pragma comment(lib, "opencv_world300.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "mlpack.lib")
#pragma comment(lib, "levmar.lib")

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
int tileSize = 256;
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
radiLineRegistration::point_type pointType = radiLineRegistration::point_type::control;
// whether to use geographic coordinates if the output projection is not set
bool bGeographic = false;
int slaveId = 1;
int masterId = 2;
bool bDebug = true;
bool bAppend = false;

void saveTieLinesSlave(ossimFilename filenametosave,
	vector<LinePolarPair> linePairs,
	ossimKeywordlist* prjKwl = NULL,
	bool extern_file = true)
{
	FILE *pf = fopen(filenametosave.c_str(), "w+");
	fclose(pf);

	// 输出txt
	pf = fopen(filenametosave.c_str(), "a+");

	for (int i = 0; i<(int)linePairs.size(); i++)
	{
		if (i != 0)
		{
			fprintf(pf, "\n");
		}
		fprintf(pf, "4%03d\t%s\t%d", i+1, "StraightLine", 2);
		fprintf(pf, "\n%20.8lf\t%20.8lf", linePairs[i].line.pt1.x, linePairs[i].line.pt1.y);
		fprintf(pf, "\n%20.8lf\t%20.8lf", linePairs[i].line.pt2.x, linePairs[i].line.pt2.y);
	}
	fclose(pf);
}


void saveTieLinesMaster(ossimFilename filenametosave,
	vector<LinePolarPair> linePairs,
	ossimKeywordlist* prjKwl = NULL,
	bool extern_file = true)
{
	FILE *pf = fopen(filenametosave.c_str(), "w+");
	fclose(pf);

	// 输出txt
	pf = fopen(filenametosave.c_str(), "a+");

	for (int i = 0; i<(int)linePairs.size(); i++)
	{
		if (i != 0)
		{
			fprintf(pf, "\n");
		}
		fprintf(pf, "4%03d\t%s\t%d", i + 1, "StraightLine", 2);
		fprintf(pf, "\n%20.8lf\t%20.8lf", linePairs[i].line_prime.pt1.x, linePairs[i].line_prime.pt1.y);
		fprintf(pf, "\n%20.8lf\t%20.8lf", linePairs[i].line_prime.pt2.x, linePairs[i].line_prime.pt2.y);
	}
	fclose(pf);
}

void matchLine(ossimFilename sourcePath, ossimFilename stdPath, ossimFilename outPath)
{
	//ossimInit::instance()->initialize();	
	bool result;
	ossimKeywordlist geomp;
	vector<ossimFilename> inputfile;
	ossimFilename master,slave,file_tmp,outimgfile;

	radiLineRegistration *img_reg = new radiLineRegistration;
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
	//img_reg->setPointType(radiLineRegistration::point_type::tie);

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

	//////img_reg->appendControlPoints(outPath);
	//////img_reg->appendTiePoints(outPath);
	ossimRefPtr<ossimImageHandler> handler = ossimImageHandlerRegistry::instance()->open(stdPath);
	ossimKeywordlist outPrj;
	//////img_reg->writePoints(outPath, PTR_CAST(ossimMapProjection, handler->getImageGeometry()->getProjection()));
	//if (bAppend)
	//{
	//	img_reg->appendPoints(outPath);
	//}
	//else
	//{
	//	img_reg->writePoints(outPath);
	//}

	//ossimFilename strProjectionFile(pszProjectFile);
	//ossimKeywordlist outPrj;
	//ossimRefPtr < ossimMapProjection > MapProjection;
	//if (strProjectionFile.exists())
	//{
	//	outPrj.addFile(strProjectionFile);
	//	// 测试投影有效性
	//	if (!(MapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(outPrj))))
	//	{
	//		// 如果无效，则清空
	//		cout << "Warning: 无效的投影文件" << endl;
	//	}
	//}
	img_reg->writeTieLines(outPath.path(), PTR_CAST(ossimMapProjection, handler->getImageGeometry()->getProjection()));
	//img_reg->writeTieLines(outPath.path());

	//saveTieLinesSlave("source.txt", img_reg->m_TieLines, NULL);
	//saveTieLinesMaster("reference.txt", img_reg->m_TieLines, NULL);

	delete img_reg;
	img_reg = NULL;
}

#include <ossim/base/ossimProcessInterface.h>
#include <ossim/base/ossimObjectFactoryRegistry.h>
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
					pointType = radiLineRegistration::point_type::control;
				}
				else if (0 == _stricmp(strPointType,"tie"))
				{
					pointType = radiLineRegistration::point_type::tie;
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
					//ossimPreferences::instance()->loadPreferences(preferences_file);
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
			matchLine(sourceFile, referenceFile, gcpFile);
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
		pszGcpFile = "I:\\testdata\\GF1\\xinjiang-jierjisisitan\\GCP_vrt.xml";

		//pszSourceFile = "E:\\HJ1\\HJ-L2\\HJ\\2010.tif";
		////pszReferenceFile = "D:\\workspace\\fuhai\\FCGC600205282\\IMG_PHR1A_P_001\\rect.tif";
		//pszReferenceFile = "E:\\HJ1\\HJ-L2\\TM\\2010.tif";

		pszSourceFile = "D:\\workspace\\Landsat\\xinjiang\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\header.dat";
		pszReferenceFile = "D:\\workspace\\Landsat\\xinjiang\\features\\reference.tif";
		pszGcpFile = "D:\\workspace\\Landsat\\xinjiang\\200lines.txt";
		//pszProjectFile = "D:\\workspace\\Landsat\\xinjiang\\gcps\\gcps_hgt_26.geom";

		//pszSourceFile = "D:\\workspace\\Landsat\\feature_test\\LS5_TM_20100717_050010_050035_145033_FASTB_L2\\header.dat";
		//pszReferenceFile = "D:\\workspace\\Landsat\\feature_test\\LS5_TM_20080625_000000_000000_145033_GEOTIFF_L4\\L5-TM-145-033-20080625-L4.TIF";
		//pszGcpFile = "D:\\workspace\\Landsat\\feature_test\\200lines.txt";

		//pszSourceFile = "D:\\workspace\\Radarsat2\\imagery_LeeH.pix";
		//pszReferenceFile = "D:\\workspace\\jiading\\image.pix";
		//pszGcpFile = "D:\\workspace\\Radarsat2\\200lines.txt";
		sBand = 4;
		mBand = 0;
		//pointType = radiLineRegistration::point_type::tie;
		pointType = radiLineRegistration::point_type::control;
		nPoint = 9;
		nThreads = 1;
		sAccuracy = 25;
		bAppend = false;
		bDebug = true;
		sBand = 0;
		mBand = 0;
		tileSize = 256*1;
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
		matchLine(sourceFile, referenceFile, gcpFile);
		clockEnd = clock();
		printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);

		Usage(0);
	}
	return 0;
}