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
#include <ossim/parallel/ossimMultiThreadSequencer.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>

#include <QDir>
#include <fstream>
#include <func.h>



#include <strUtil.h>
#include <fileUtil.h>
#include <mprojectdefine.h>
#include <time.h>

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

const char *pszInputFile = "";
const char *pszGcpFile = "";
const char *pszOutFile = "";
const char *pszInputPath = "";
const char *pszFilter = "*.rpb";
const char *pszOutPath = "orth";
const char *pszDemPath = "";
const char *pszPluginPath = "";
bool bOverwrite = false;

const char *pszProjectFile = "";
const char *pszPreferenceFile = "preference.txt";

bool searchRpcFiles(QString imageFileName, ossimRpcModel::rpcModelStruct &rpcStruct)
{
	// rpb file
	QString rpcFile = QBeforeLast(imageFileName, '.') + ".rpb";
	if (QFileInfo(rpcFile).exists())
	{
		mylib::readRPBFile(ossimFilename(rpcFile.toLatin1()), rpcStruct);
		return true;
	}

	rpcFile = QBeforeLast(imageFileName, '.') + "_rpc.txt";
	if (QFileInfo(rpcFile).exists())
	{
		mylib::readRPCFile(ossimFilename(rpcFile.toLatin1()), rpcStruct);
		return true;
	}

	rpcFile = QBeforeLast(imageFileName, '.') + ".rpc";
	if (QFileInfo(rpcFile).exists())
	{
		mylib::readRPCFile(ossimFilename(rpcFile.toLatin1()), rpcStruct);
		return true;
	}

	return false;
}

void rpc_orth(ossimFilename inputFile, ossimFilename outputFile,
			  ossimFilename gcpFile, ossimFilename projectionFile)
{
	ossimplugins::radiRpcModel *rpcModel = new ossimplugins::radiRpcModel;
	if(!rpcModel->parseRpcFile(inputFile))
	{
		std::cerr<<"cannot find the rpc file."<<endl;
	}
	//ossimRpcModel::rpcModelStruct rpcStruct;
	//if(!searchRpcFiles(QString(pszImageFile), rpcStruct))
	//{
	//	std::cerr<<"cannot find the rpc file."<<endl;
	//}

	//ossimRpcModel *rpcModel = new ossimRpcModel;
	//rpcModel->setAttributes(rpcStruct);

	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	if (gcpFile.exists())
	{
		mylib::readGcpFile(gcpFile, gcpSet, chkSet, &prjKwl);
		mylib::get_elevation(gcpSet, prjKwl, 0.0);
		mylib::get_elevation(chkSet, prjKwl, 0.0);
		mylib::projection2ll(gcpSet, prjKwl);
		mylib::projection2ll(chkSet, prjKwl);
	}

	if (projectionFile.exists())
	{
		prjKwl.clear();
		prjKwl.addFile(projectionFile);
	}

	rpcModel->m_proj = PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(prjKwl));
	rpcModel->setUseL1(false);
	rpcModel->optimizeFit(*gcpSet);
	rpcModel->updateModel();
	rpcModel->print(cout);
	mylib::OutputReport("report.txt", rpcModel, gcpSet, chkSet, true, true);

	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(ossimFilename(inputFile));
	if(!handler)
	{
		std::cerr<<"cannot open file : "<< inputFile<<endl;
		return;   //应该弹出警告对话框
	}

	ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
	imageGeom->setProjection(rpcModel);

	handler->setImageGeometry(imageGeom.get());

	ossimKeywordlist tt_geom;
	ossimDpt imagesize(handler->getImageRectangle().width(), handler->getImageRectangle().height());
	ossimImageRenderer* renderer = new ossimImageRenderer;
	ossimPolyCutter* theCutter;
	ossimBandSelector* theBandSelector;
	vector<ossimDpt> polygon;
	ossimIrect bound,boundw;
	theCutter = new ossimPolyCutter;
	int starline = 0;
	int starpixel =	0;
	int endline = 0;
	int endpixel = 0;

	//starline = 8500;
	//starpixel =	11000;
	//endline = 20650;
	//endpixel = 21500;

	//starline = 2010;
	//starpixel =	3768;
	//endline = 10368;
	//endpixel = 12557;
	if (endpixel == 0)
	{
		endpixel = imagesize.x - 1;
	}
	if (endline == 0)
	{
		endline = imagesize.y - 1;
	}

	theBandSelector = new ossimBandSelector;
	theBandSelector->connectMyInputTo(0, handler);
	int nBands = handler->getNumberOfInputBands();
	std::vector<ossim_uint32> bandList;
	for (int i=0;i<nBands;++i)
	{
		bandList.push_back(i);
	}
	theBandSelector->setOutputBandList(bandList);


	ossimDpt ps(starpixel,starline),p2(endpixel,starline),p3(endpixel,endline),p4(starpixel,endline),p5(starpixel,starline);
	polygon.push_back(ps);
	polygon.push_back(p2);
	polygon.push_back(p3);
	polygon.push_back(p4);
	polygon.push_back(p5);
	theCutter->connectMyInputTo(theBandSelector);
	theCutter->setPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);


	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
	////renderer->getResampler()->setFilterType(filter_type.c_str());////////////////应该从界面取得//
	//if(m_SampleType.contains("NEAREST_NEIGHBOR"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_NEAREST_NEIGHBOR);
	//if(m_SampleType.contains("BILINEAR"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	//if(m_SampleType.contains("BICUBIC"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);

	ossimRefPtr<ossimImageFileWriter> writer = ossimImageWriterFactoryRegistry::instance()->
		createWriterFromExtension( outputFile.ext() );
	if(writer==NULL) return;
	tt_geom.clear();
	writer->saveState(tt_geom);
	tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA",true);
	writer->loadState(tt_geom);

	int nThreads = OpenThreads::GetNumberOfProcessors() * 2;
	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
	writer->changeSequencer(sequencer.get());
	renderer->connectMyInputTo(theCutter);
	renderer->setView(rpcModel->m_proj);

	myOutProgress *progress = new myOutProgress(0,true);
	writer->addListener(progress);

	writer->setFilename(outputFile);
	writer->connectMyInputTo(0,renderer);

	if(!writer->execute())
	{
		writer->removeListener(progress);
		writer->disconnectAllInputs();
		renderer->disconnectAllInputs();
		theCutter->disconnectAllInputs();
		std::cerr<<"orth rectification failed."<<endl;
		return;
	}

	writer->disableListener();
	writer->removeListener(progress);
	writer->disconnectAllInputs();
	renderer->disconnectAllInputs();
	theCutter->disconnectAllInputs();
	handler->close();
}

void test()
{

	//pszImageFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\13AUG03025915-P2AS-053553869020_01_P001.TIF";
	//pszGcpFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\gcp.txt";
	//pszOutFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\orth.tif";

	//pszImageFile = "D:\\workspace\\WV01\\53553869030_01-beijing\\053553869030_01\\053553869030_01_P001_MUL\\13SEP21033152-M2AS-053553869030_01_P001.TIF";
	//pszGcpFile = "D:\\workspace\\WV01\\53553869030_01-beijing\\053553869030_01\\053553869030_01_P001_MUL\\gcp.txt";
	//pszOutFile = "D:\\workspace\\WV01\\53553869030_01-beijing\\053553869030_01\\053553869030_01_P001_MUL\\orth.tif";


	//pszImageFile = "D:\\workspace\\WV01\\K3_20130914052922_07080_08591299_L1R\\K3_20130914052922_07080_08591299_L1R_P.tif";
	//pszGcpFile = "D:\\workspace\\WV01\\K3_20130914052922_07080_08591299_L1R\\gcp.txt";
	//pszOutFile = "D:\\workspace\\WV01\\K3_20130914052922_07080_08591299_L1R\\orth.tif";

	//pszImageFile = "D:\\workspace\\WV01\\K3_20130914052922_07080_08591299_L1R\\mul\\K3_20130914052922_07080_08591299_L1R_N.tif";
	//pszGcpFile = "D:\\workspace\\WV01\\K3_20130914052922_07080_08591299_L1R\\mul\\gcp.txt";
	//pszOutFile = "D:\\workspace\\WV01\\K3_20130914052922_07080_08591299_L1R\\mul\\orth_K3_20130914052922_07080_08591299_L1R_N.tif";

	////pszImageFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\MSC_110526013121_25751_10031299PN12_1R.tif";
	////pszGcpFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\gcp.txt";
	////pszOutFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\orth.tif";

	//pszImageFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\MSC_110526013121_25751_10031299M1N12G_1R.tif";
	//pszGcpFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\mul_gcp.txt";
	//pszOutFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\orth_MSC_110526013121_25751_10031299M1N12G_1R.tif";
	//rpc_orth();

	//pszImageFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\MSC_110526013121_25751_10031299M2N12B_1R.tif";
	//pszGcpFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\mul_gcp.txt";
	//pszOutFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\orth_MSC_110526013121_25751_10031299M2N12B_1R.tif";
	//rpc_orth();

	//pszImageFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\MSC_110526013121_25751_10031299M3N12N_1R.tif";
	//pszGcpFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\mul_gcp.txt";
	//pszOutFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\orth_MSC_110526013121_25751_10031299M3N12N_1R.tif";
	//rpc_orth();

	//pszImageFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\MSC_110526013121_25751_10031299M4N12R_1R.tif";
	//pszGcpFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\mul_gcp.txt";
	//pszOutFile = "D:\\workspace\\WV01\\K2_MSC_110526013121_25751_10031299\\orth_MSC_110526013121_25751_10031299M4N12R_1R.tif";

	//pszImageFile = "I:\\20140523\\MSS\\GF1_PMS1_E100.8_N24.4_20140202_L1A0000160810-MSS1.tiff";
	//pszGcpFile = "I:\\20140523\\fusion\\gcp1.txt";
	//pszOutFile = "I:\\20140523\\fusion\\GF1_PMS1_E100.8_N24.4_20140202_L1A0000160810-MSS1.tif";

	//pszImageFile = "I:\\20140523\\MSS\\GF1_PMS1_E100.8_N24.7_20140202_L1A0000160809-MSS1.tiff";
	//pszGcpFile = "I:\\20140523\\fusion\\gcp2.txt";
	//pszOutFile = "I:\\20140523\\fusion\\GF1_PMS1_E100.8_N24.7_20140202_L1A0000160809-MSS1.tif";

	//pszImageFile = "I:\\20140523\\MSS\\GF1_PMS2_E101.1_N24.3_20140202_L1A0000160875-MSS2.tiff";
	//pszGcpFile = "I:\\20140523\\fusion\\gcp3.txt";
	//pszOutFile = "I:\\20140523\\fusion\\GF1_PMS2_E101.1_N24.3_20140202_L1A0000160875-MSS2.tif";
	//rpc_orth();

	//pszImageFile = "I:\\20140523\\MSS\\GF1_PMS2_E101.2_N24.6_20140202_L1A0000160874-MSS2.tiff";
	//pszGcpFile = "I:\\20140523\\fusion\\gcp4.txt";
	//pszOutFile = "I:\\20140523\\fusion\\GF1_PMS2_E101.2_N24.6_20140202_L1A0000160874-MSS2.tif";


	//pszImageFile = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\IMAGE.TIF";
	//pszGcpFile = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\gcp1111.txt";
	//pszOutFile = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404110131_201404110140\\Scene02\\rpb_orth.tif";

	//pszImageFile = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01\\SatAngle.TIF";
	//pszGcpFile = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01\\gcp_angle.txt";
	//pszOutFile = "E:\\HJ1\\HJ-1B_CCD-1_MYC_201404100243_201404100250\\Scene01\\rpb_orth_SatAngle.tif";

	//pszImageFile = "I:\\GF1_test\\GF1_WFV1\\SatAngle.tif";
	//pszGcpFile = "I:\\GF1_test\\GF1_WFV1\\gcp1111.txt";
	//pszOutFile = "I:\\GF1_test\\GF1_WFV1\\rpb_orth_SatAngle.tif";

	pszInputFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\13AUG03025915-P2AS-053553869020_01_P001.TIF";
	pszGcpFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\gcp.txt";
	pszOutFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\rpb_orth.tif";
	//rpc_orth();
}

/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: RpcOrth -i inputFile [-o outputFile] [-gcp gcpFile] [-overwrite] \n"
		"\t[-prj projectionfile] [-pref preferencefile]\n");

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

	const char *pszInputPath = "";
	const char *pszInputFile = "";
	int null_value = 0;
	//ossimRefPtr<ossimEquDistCylProjection> llProjection = new ossimEquDistCylProjection;
	//ossimKeywordlist kwl;
	//llProjection->saveState(kwl);
	//cout<<kwl<<endl;
	//fstream fs;
	//fs.open("ll.txt", ios_base::out);
	//fs<<kwl;
	//fs.close();


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
				pszInputFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-o") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutFile = argv[++i] ;
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
			else
			{
				Usage();
			}
		}

		if (0 == strcmp(pszInputFile, ""))
		{
			printf("inputFile can not be empty!\n");
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
			
			
			ossimFilename inputFile(QDir::toNativeSeparators(QDir(pszInputFile).absolutePath()).toStdString());
			//cout<<inputFile<<endl;
			ossimFilename gcpFile(pszGcpFile);
			ossimFilename projectionFile(pszProjectFile);
			ossimFilename outputFile(pszOutFile);
			//ossimFilename gcpFile(QDir::toNativeSeparators(QDir(pszGcpFile).absolutePath()).toStdString());
			//ossimFilename projectionFile(QDir::toNativeSeparators(QDir(pszProjectFile).absolutePath()).toStdString());
			//ossimFilename outputFile(QDir::toNativeSeparators(QDir(pszOutFile).absolutePath()).toStdString());
			if (0 == strcmp(pszOutFile, ""))
			{
				outputFile = inputFile.path() + "//o" + inputFile.fileNoExtension() + ".tif";
			}
			if (outputFile.exists() && !bOverwrite)
			{
				printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
					" \"-overwrite\" option.\n", outputFile);
				Usage(0);
			}
			else
			{
				gcpFile.setExtension(gcpFile.ext());
				//cout<<gcpFile<<endl;
				clock_t  clockBegin, clockEnd;
				clockBegin = clock();

				rpc_orth(inputFile, outputFile,
					gcpFile, projectionFile);
				clockEnd = clock();
				printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
			}

		}
	}
	else
	{
		pszInputFile = "I:\\GF1-Testdata\\wuzhishan\\GF1_PMS2_E113.9_N37.7_20140526_L1A0000234945-MSS2.tiff";
		pszGcpFile = "I:\\GF1-Testdata\\wuzhishan\\gcp_auto.txt";
		pszInputFile = "D:\\workspace\\fuhai\\FCGC600205282\\IMG_PHR1A_P_001\\IMG_PHR1A_P_201309040246085_SEN_922557101-001_R1C1.TIF";
		pszGcpFile = "D:\\workspace\\fuhai\\FCGC600205282\\IMG_PHR1A_P_001\\g.txt";
		bOverwrite = true;
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
		//ossimElevManager::instance()->loadElevationPath("D:\\workspace\\dem\\srtm90");

		//pszInputFile = "I:\\GF1-Testdata\\1\\PAN\\GF1_PMS1_E76.2_N42.5_20130918_L1A0000085664-PAN1.tiff";
		//pszGcpFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\gcp.txt";
		//pszOutFile = "D:\\workspace\\WV01\\53553869020_01-beijing\\053553869020_01\\053553869020_01_P001_PAN\\rpb_orth.tif";


		ossimFilename inputFile(QDir::toNativeSeparators(QDir(pszInputFile).absolutePath()).toStdString());
		//cout<<inputFile<<endl;
		ossimFilename gcpFile(pszGcpFile);
		ossimFilename projectionFile(pszProjectFile);
		ossimFilename outputFile(pszOutFile);
		//ossimFilename gcpFile(QDir::toNativeSeparators(QDir(pszGcpFile).absolutePath()).toStdString());
		//ossimFilename projectionFile(QDir::toNativeSeparators(QDir(pszProjectFile).absolutePath()).toStdString());
		//ossimFilename outputFile(QDir::toNativeSeparators(QDir(pszOutFile).absolutePath()).toStdString());
		if (0 == strcmp(pszOutFile, ""))
		{
			outputFile = inputFile.path() + "//o" + inputFile.fileNoExtension() + ".tif";
		}
		if (outputFile.exists() && !bOverwrite)
		{
			printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
				" \"-overwrite\" option.\n", outputFile);
			Usage(0);
		}
		else
		{
			gcpFile.setExtension(gcpFile.ext());
			//cout<<gcpFile<<endl;
			clock_t  clockBegin, clockEnd;
			clockBegin = clock();

			rpc_orth(inputFile, outputFile,
				gcpFile, projectionFile);
			clockEnd = clock();
			printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
		Usage(0);
	}
	return 0;
}