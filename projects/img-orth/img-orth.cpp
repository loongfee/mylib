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
#include <ossim/base/ossimStringProperty.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>

#include <QDir>
#include <QDebug>
#include <fstream>
#include <func.h>

#include <QApplication>  
#include <QTextCodec>  
#include <QLabel>  

#include <strUtil.h>
#include <fileUtil.h>
#include <mprojectdefine.h>
#include <time.h>
#include <tiff.h>


#include <argtable2/argtable2.h>
#include <assert.h>
#include <boost/filesystem.hpp>

/* for memory leak debugging */
#ifdef DMALLOC
#include "dmalloc.h"
#endif

using namespace std;
using namespace mylib;
namespace fs = boost::filesystem;

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
#pragma comment(lib, "argtable2.lib")
//#pragma comment(lib, "mpi.lib")

//// 代码一定要是: UTF-8(BOM)  
////qt版本高于等于qt5,VS版本高于或等于VS2010
//#if _MSC_VER >= 1600  
//#pragma execution_character_set("utf-8")  
//#endif


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
const char *pszModelFile = "";
const char *pszFilter = "*.rpb";
const char *pszOutPath = "orth";
const char *pszDemPath = "";
const char *pszPluginPath = "";
const char *pszReportFile = "";
vector<ossim_uint32> outBandList;
int nThreads = 0;
bool bOverwrite = false;
bool bOverview = false;
bool bRpc = false;
bool bUseL1 = false;

const char *pszProjectFile = "";
const char *pszPreferenceFile = "preference.txt";

//bool searchRpcFiles(QString imageFileName, ossimRpcModel::rpcModelStruct &rpcStruct)
//{
//	// rpb file
//	QString rpcFile = QBeforeLast(imageFileName, '.') + ".rpb";
//	if (QFileInfo(rpcFile).exists())
//	{
//		mylib::readRPBFile(ossimFilename(rpcFile.toLatin1()), rpcStruct);
//		return true;
//	}
//
//	rpcFile = QBeforeLast(imageFileName, '.') + "_rpc.txt";
//	if (QFileInfo(rpcFile).exists())
//	{
//		mylib::readRPCFile(ossimFilename(rpcFile.toLatin1()), rpcStruct);
//		return true;
//	}
//
//	rpcFile = QBeforeLast(imageFileName, '.') + ".rpc";
//	if (QFileInfo(rpcFile).exists())
//	{
//		mylib::readRPCFile(ossimFilename(rpcFile.toLatin1()), rpcStruct);
//		return true;
//	}
//
//	return false;
//}


bool orth(ossimFilename inputFile, ossimFilename outputFile,
		  ossimFilename gcpFile, ossimFilename projectionFile)
{	
	ossimRefPtr<ossimImageHandler> handler   = ossimImageHandlerRegistry::instance()->open(inputFile);
	if(!handler) return false;   //应该弹出警告对话框

	ossimRefPtr<ossimSensorModel> sensorModel;
	if (bRpc || bUseL1)
	{
		// 强制使用RPC模型
		ossimplugins::radiRpcModel *rpcModel = new ossimplugins::radiRpcModel;
		if(!rpcModel->parseRpcFile(inputFile))
		{
			//
			cout<<"Failed to find rpc file for rpc model."<<endl;
			return false;
		}
		rpcModel->setUseL1(bUseL1);
		sensorModel = PTR_CAST(ossimSensorModel, rpcModel);
	}
	else
	{
		sensorModel = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());
	}

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
	ossimFilename modelFile(pszModelFile);
	if (modelFile.exists())
	{
		ossimKeywordlist model_kwl;
		model_kwl.addFile(modelFile);
		sensorModel->loadState(model_kwl);
	}

	sensorModel->m_proj = PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(prjKwl));
	
	sensorModel->optimizeFit(*gcpSet);
	sensorModel->updateModel();
	//sensorModel->print(cout);
	//ossimFilename rpcFile = inFile;
	//rpcFile.setExtension("rpd");
	//((ossimplugins::radiRpcModel*)sensorModel.get())->writeRpcFile(rpcFile);

	if (0 != strcmp(pszReportFile, ""))
	{
		mylib::OutputReport(pszReportFile, sensorModel.get(), gcpSet, chkSet, true, true);
	}

	//ossimKeywordlist prjKwl;
	//prjKwl.addFile(prjFile);
	//ossimMapProjection* MapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(prjKwl));

	//ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
	//imageGeom->setProjection(m_sensorModel);
	//imageGeom->loadState(geom);

	ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
	imageGeom->setProjection(sensorModel.get());
	handler->setImageGeometry(imageGeom.get());//1128
	//cout<<geom<<endl;
	ossimKeywordlist tt_geom;

	ossimDpt imagesize(handler->getImageRectangle().width(), handler->getImageRectangle().height());
	ossimImageRenderer* renderer = new ossimImageRenderer;
	ossimPolyCutter* theCutter;
	ossimBandSelector* theBandSelector;
	vector<ossimDpt> polygon;
	ossimIrect bound,boundw;
	theCutter = new ossimPolyCutter;
	///////////////////以下四行应该从界面取得//
	int startpixel = 0;
	int startline = 0;
	int endpixel = 0;
	int endline = 0;
	if (endpixel == 0)
	{
		endpixel = imagesize.x - 1;
	}
	if (endline == 0)
	{
		endline = imagesize.y - 1;
	}
	theBandSelector = new ossimBandSelector;
	theBandSelector->connectMyInputTo(0, handler.get());
	if (outBandList.empty())
	{
		int nBands = handler->getNumberOfInputBands();
		for (int i=0;i<nBands;++i)
		{
			outBandList.push_back(i);
		}
	}

	theBandSelector->setOutputBandList(outBandList);


	ossimDpt ps(startpixel,startline),p2(endpixel,startline),p3(endpixel,endline),p4(startpixel,endline),p5(startpixel,startline);

	polygon.push_back(ps);
	polygon.push_back(p2);
	polygon.push_back(p3);
	polygon.push_back(p4);
	polygon.push_back(p5);
	theCutter->connectMyInputTo(theBandSelector);
	theCutter->setPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);



	//renderer->getResampler()->setFilterType(filter_type.c_str());////////////////应该从界面取得//
	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
	ossimRefPtr<ossimImageFileWriter> writer;
	//if (0 == strcmp(outputFile.ext().upcase(), "TIF"))
	//{
	//	writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_gtiff"));

	//}
	//else if (0 == strcmp(outputFile.ext().upcase(), "PIX"))
	//{
	//	writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_PCIDSK"));
	//}
	//else if (0 == strcmp(outputFile.ext().upcase(), "PIX"))
	//{
	//	writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_PCIDSK"));
	//}
	//else
	{
		writer = ossimImageWriterFactoryRegistry::instance()->
			createWriterFromExtension( outputFile.ext() );
		//writer->setOverviewCompressType(COMPRESSION_JPEG);
		//writer->setOverviewCompressType(COMPRESSION_LZW);
		//ossimRefPtr<ossimProperty> prop =  new ossimStringProperty("gdal_overview_type", "nearest", false);
		//writer->setProperty(prop);
	}
	if (bOverview)
	{
		writer->setWriteOverviewFlag(true);
	}
	if(writer==NULL) return false;
	tt_geom.clear();
	writer->saveState(tt_geom);
	tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA", true);
	writer->loadState(tt_geom);

	renderer->connectMyInputTo(theCutter);
	renderer->setView(sensorModel->m_proj);

	myOutProgress *m_progress = new myOutProgress(0,true);
	writer->addListener(m_progress);

	writer->setFilename(outputFile);

	if (nThreads < 1)
	{
		nThreads = OpenThreads::GetNumberOfProcessors() * 2;
	}
	cout<<"using "<<nThreads<<" threads:"<<endl;
	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
	writer->changeSequencer(sequencer.get());
	writer->connectMyInputTo(0,renderer);
	//ossimIrect bounding = handler->getBoundingRect();
	//writer->setAreaOfInterest(bounding);

	//bool bResult = writer->execute();

	if(!writer->execute())
	{
		writer->removeListener(m_progress);
		writer->disconnectAllInputs();
		renderer->disconnectAllInputs();
		theCutter->disconnectAllInputs();
		//delete handler;
		//handler->close();
		//delete handler;
		return false;
	}

	writer->disableListener();
	writer->removeListener(m_progress);
	writer->disconnectAllInputs();
	renderer->disconnectAllInputs();
	theCutter->disconnectAllInputs();
	return true;
}


///*************************************************************************/
///*                               Usage()                                */
///************************************************************************/
//
//static void Usage(const char* pszErrorMsg = NULL)
//
//{
//	printf( 
//		"Usage: img-orth -i inputFile [-o outputFile] [-gcp gcpFile] [-overwrite] \n"
//		"\t[-prj projectionfile] [-pref preferencefile] [-ov {0|1}] [-l1 {0|1}]\n"
//		"\t[-b nb b1 b2 ...] [-nt nthreads] [-rpc] [-model modelfile] [-report reportFile]\n");
//
//	printf( "   -ov {0|1}: whether create overview, default value is 0.\n");
//
//	if( pszErrorMsg != NULL )
//		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
//	exit(1);
//}
//
//#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
//	do { if (i + nExtraArg >= argc) \
//	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)



int mymain(struct arg_file *gcp_file,
	struct arg_file *report_file,
	struct arg_file *projection_file,
	struct arg_file *preference_file,
	struct arg_file *model_file,
	struct arg_int *out_bands,
	struct arg_int *use_threads,
	struct arg_lit *overview,
	struct arg_lit *rpc,
	struct arg_lit *l1,
	struct arg_lit *overwrite,
	struct arg_lit  *help,
	struct arg_lit  *version,
	struct arg_file *infile,
	struct arg_file *outfile)
{

	if (0 == strcmp(infile->filename[0], ""))
	{
		printf("inputFile can not be empty!\n");
		return 0;
	}
	else
	{
		ossimFilename preferences_file = ossimFilename(preference_file->filename[0]);
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

		ossimInit::instance()->initialize();


		if (overview->count > 0)
		{
			bOverview = true;
		}
		if (overwrite->count > 0)
		{
			bOverwrite = true;
		}
		if (rpc->count > 0)
		{
			bRpc = true;
		}
		if (l1->count > 0)
		{
			bUseL1 = true;
		}
		nThreads = use_threads->ival[0];


		outBandList.clear();
		for (int ib = 0; ib < out_bands->count; ++ib)
		{
			outBandList.push_back(out_bands->ival[ib] - 1);
		}

		if (!fs::exists(infile->filename[0]))
		{
			printf("inputFile does not exist!\n");
			return 0;
		}
		//fs::path inputPathname = fs::system_complete(fs::path(infile->filename[0], fs::native));	//将相对路径转换为绝对路径
		fs::path inputPathname = fs::absolute(infile->filename[0]);

		ossimFilename inputFile(inputPathname.string().c_str());
		ossimFilename gcpFile(gcp_file->filename[0]);
		ossimFilename projectionFile(projection_file->filename[0]);
		ossimFilename outputFile(outfile->filename[0]);
		pszReportFile = report_file->filename[0];
		if (0 == strcmp(outfile->filename[0], ""))
		{
			outputFile = inputFile.path() + "\\o" + inputFile.fileNoExtension() + ".tif";
		}
		if (outputFile.exists() && !bOverwrite)
		{
			printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
				" \"-overwrite\" option.\n", outputFile.c_str());
			return 0;
		}
		else
		{
			gcpFile.setExtension(gcpFile.ext());
			//cout<<gcpFile<<endl;
			cout << inputFile << endl;
			clock_t  clockBegin, clockEnd;
			clockBegin = clock();

			orth(inputFile, outputFile,
				gcpFile, projectionFile);
			clockEnd = clock();
			printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
		}

	}
	return 0;
}

//	if (0 != _strcmpi(sensorType[0], strHJ1)
//		&& 0 != _strcmpi(sensorType[0], strZY3))
//	{
//		printf("Unknown Sensor Type.\n");
//		return 0;
//	}
//
//	strInputFile = infile[0];
//	if (0 == _strcmpi(strInputFile.c_str(), ""))
//	{
//		printf("inputFile can not be empty!\n");
//		return 0;
//	}
//
//	if (!fs::exists(strInputFile))
//	{
//		printf("inputFile does not exist!\n");
//	}
//	else
//	{
//		fs::path inputPathname = fs::absolute(strInputFile);
//		//char fullPath[_MAX_PATH];
//		//_fullpath(fullPath, pszInputFile, _MAX_PATH);
//
//		if (0 == _strcmpi(outpath[0], ""))
//		{
//			strOutputDir = inputPathname.parent_path().string();
//			//pszOutputDir = SBeforeLast(string(fullPath), '\\').c_str();
//		}
//		else
//		{
//			strOutputDir = fs::absolute(outpath[0]).string();
//		}
//
//		if (!fs::exists(strOutputDir))
//		{
//			_mkdir(strOutputDir.c_str());
//		}
//
//		clock_t  clockBegin, clockEnd;
//		clockBegin = clock();
//		if (0 == _strcmpi(sensorType[0], strHJ1))
//		{
//			// HJ1
//			QVHJ::QVProcReader qv;
//			qv.read_by_scene(inputPathname.string().c_str(), strOutputDir.c_str());
//		}
//		else if (0 == _strcmpi(sensorType[0], strZY3))
//		{
//			// ZY3
//			QVProc::ZY3QVProcReader zy3_qv;
//			zy3_qv.read_by_scene(inputPathname.string().c_str(), strOutputDir.c_str());
//		}
//		else
//		{
//			// Unknown Type
//			printf("Unknown Sensor Type.\n");
//			return 0;
//		}
//
//		clockEnd = clock();
//		printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
//	}
//	return 0;
//}

int main( int argc, char** argv )
{
	GDALAllRegister();
	//OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持


	struct arg_file *gcp_file = arg_file0(NULL, "gcp", NULL, "gcp file");
	struct arg_file *report_file = arg_file0(NULL, "report", NULL, "report file");
	struct arg_file *projection_file = arg_file0(NULL, "prj", NULL, "preference file");
	struct arg_file *preference_file = arg_file0(NULL, "pref,preference", NULL, "projection file");
	struct arg_file *model_file = arg_file0(NULL, "model", NULL, "model file");
	struct arg_int *out_bands = arg_intn("b", NULL, NULL, 0, argc + 3, "output bands (begin with 1)");
	struct arg_int *use_threads = arg_int0(NULL, "nt,nthreads", NULL, "the number of threads (if 0, automatically choose upon the processors)");
	struct arg_lit *overview = arg_lit0(NULL, "ov,overview", "create overview");
	struct arg_lit *rpc = arg_lit0(NULL, "rpc", "mandatorily use rpc ");
	struct arg_lit *l1 = arg_lit0(NULL, "l1", "use l1 norm regularizer");
	struct arg_lit *overwrite = arg_lit0("ow", "overwrite", "overwrite the existing file");
	struct arg_lit  *help = arg_lit0("h", "help", "print this help and exit");
	struct arg_lit  *version = arg_lit0("v", "version", "print version information and exit");
	struct arg_file *infile = arg_file1(NULL, NULL, "<input file>", "input file");
	struct arg_file *outfile = arg_file0(NULL, NULL, "<output file>", "output file");
	struct arg_end  *end = arg_end(20);


	//struct arg_lit  *list = arg_lit0("lL", NULL, "list files");
	//struct arg_lit  *recurse = arg_lit0("R", NULL, "recurse through subdirectories");
	//struct arg_int  *repeat = arg_int0("k", "scalar", NULL, "define scalar value k (default is 3)");
	//struct arg_str  *defines = arg_strn("D", "define", "MACRO", 0, argc + 2, "macro definitions");
	////struct arg_file *outfile = arg_file0("o", NULL, "<output>", "output file (default is \"-\")");
	//struct arg_lit  *verbose = arg_lit0("v", "verbose,debug", "verbose messages");
	//struct arg_lit  *help = arg_lit0(NULL, "help", "print this help and exit");
	//struct arg_lit  *version = arg_lit0(NULL, "version", "print version information and exit");
	//struct arg_file *infiles = arg_filen(NULL, NULL, "<input>", 1, argc + 2, "input file(s)");
	//struct arg_file *outfile = arg_file0(NULL, NULL, "<output>", "output file");
	//struct arg_end  *end = arg_end(20);
	////void* argtable[] = { list, recurse, repeat, defines, outfile, verbose, help, version, infiles, end };
	//void* argtable[] = { list, recurse, repeat, defines, verbose, help, version, infiles, outfile, end };
	void* argtable[] = { gcp_file, report_file, projection_file, preference_file,
		model_file, out_bands, use_threads, overview, rpc, l1, overwrite, help, version, infile, outfile, end };
	const char* progname = "img-orth";
	int nerrors;
	int exitcode = 0;

	/* verify the argtable[] entries were allocated sucessfully */
	if (arg_nullcheck(argtable) != 0)
	{
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n", progname);
		exitcode = 1;
		goto exit;
	}

	// default value
	use_threads->ival[0] = 0;

	///* set any command line default values prior to parsing */
	//repeat->ival[0] = 3;
	//outfile->filename[0] = "-";

	/* Parse the command line as defined by argtable[] */
	nerrors = arg_parse(argc, argv, argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0)
	{
		printf("Usage: %s", progname);
		arg_print_syntax(stdout, argtable, "\n");
		printf("This program ortho-rectifies remotely sensed image.\n");
		//printf("for parsing command line arguments. Argtable accepts integers\n");
		//printf("in decimal (123), hexadecimal (0xff), octal (0o123) and binary\n");
		//printf("(0b101101) formats. Suffixes KB, MB and GB are also accepted.\n");
		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
		exitcode = 0;
		goto exit;
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0)
	{
		printf("'%s' Ortho-rectification.\n", progname);
		printf("December 2014, Long Tengfei.\n");
		exitcode = 0;
		goto exit;
	}

	/* If the parser returned any errors then display them and exit */
	if (nerrors > 0)
	{
		/* Display the error details contained in the arg_end struct.*/
		arg_print_errors(stdout, end, progname);
		printf("Try '%s -h or --help' for more information.\n", progname);
		exitcode = 1;
		goto exit;
	}

	/* special case: uname with no command line options induces brief help */
	if (argc == 1)
	{
		printf("Try '%s -h or --help' for more information.\n", progname);
		exitcode = 0;
		goto exit;
	}

	/* normal case: take the command line options at face value */
	exitcode = mymain(gcp_file, report_file, projection_file, preference_file,
		model_file, out_bands, use_threads, overview, rpc, l1, overwrite, help, version, infile, outfile);

exit:
	/* deallocate each non-null entry in argtable[] */
	arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

	return exitcode;

////	//qt版本低于qt5
////#if QT_VERSION < QT_VERSION_CHECK(5,0,0)
////	//VS版本低于VS2010
////#if defined(_MSC_VER) && (_MSC_VER < 1600)
////	QTextCodec::setCodecForTr(QTextCodec::codecForName("GBK"));
////#else    
////	QTextCodec::setCodecForTr(QTextCodec::codecForName("UTF-8"));
////#endif    
////
////#endif
//
//	const char *pszInputPath = "";
//	const char *pszInputFile = "";
//	int null_value = 0;
//	//ossimRefPtr<ossimEquDistCylProjection> llProjection = new ossimEquDistCylProjection;
//	//ossimKeywordlist kwl;
//	//llProjection->saveState(kwl);
//	//cout<<kwl<<endl;
//	//fstream fs;
//	//fs.open("ll.txt", ios_base::out);
//	//fs<<kwl;
//	//fs.close();
//
//
//	if (argc > 1)
//	{
//		/* -------------------------------------------------------------------- */
//		/*      Parse arguments.                                                */
//		/* -------------------------------------------------------------------- */
//		for( int i = 1; i < argc; i++ )
//		{
//			if( 0 == _stricmp(argv[i],"-i") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszInputFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-o") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszOutFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-gcp") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszGcpFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-report") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszReportFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-prj") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszProjectFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-pref") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszPreferenceFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-model") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				pszModelFile = argv[++i] ;
//			}
//			else if( 0 == _stricmp(argv[i],"-b") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				int nOutBands = atoi(argv[++i]);
//				if (nOutBands < 1)
//				{
//					printf("output band number cannot less than 1\n");
//					exit(0);
//				}
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nOutBands);
//				for (int ib = 0;ib < nOutBands;++ib)
//				{
//					outBandList.push_back(atoi(argv[++i])-1);
//				}
//			}
//			else if(0 == _stricmp(argv[i],"-nt") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				nThreads = atoi( argv[++i] );	
//			}
//			else if( 0 == _stricmp(argv[i],"-ov") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				if(atoi( argv[++i] ) > 0)
//				{
//					bOverview = true;
//				}
//				else
//				{
//					bOverview = false;
//				}
//			}
//			else if( 0 == _stricmp(argv[i],"-rpc") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
//				bRpc = true ;
//			}
//			else if( 0 == _stricmp(argv[i],"-l1") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
//				int flag = atoi( argv[++i] );	
//				if (0 == flag)
//				{
//					bUseL1 = false;
//				}
//				else
//				{
//					bUseL1 = true;
//				}
//			}
//			else if( 0 == _stricmp(argv[i],"-overwrite") )
//			{
//				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
//				bOverwrite = true ;
//			}
//			else
//			{
//				Usage();
//			}
//		}
//
//		if (0 == strcmp(pszInputFile, ""))
//		{
//			printf("inputFile can not be empty!\n");
//			Usage();
//		}
//		else
//		{
//			ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//			if (preferences_file.exists())
//			{
//				ossimPreferences::instance()->loadPreferences(preferences_file);
//			}
//			else
//			{
//				preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//				if (preferences_file.exists())
//				{
//					ossimPreferences::instance()->loadPreferences(preferences_file);
//				}
//			}
//
//
//			std::string	tempString;
//			ossimArgumentParser::ossimParameter	stringParam(tempString);
//			ossimArgumentParser argumentParser(&argc, argv);
//			ossimInit::instance()->addOptions(argumentParser);
//			argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//			ossimInit::instance()->initialize(argumentParser);
//			
//			ossimString strApplicationPath = ossimFilename(argumentParser.getApplicationName()).path();
//			ossimFilename inputFile(pszInputFile);
//			ossimFilename gcpFile(pszGcpFile);
//			ossimFilename projectionFile(pszProjectFile);
//			ossimFilename outputFile(pszOutFile);
//			if (0 == strcmp(pszOutFile, ""))
//			{
//				outputFile = inputFile.path() + "\\o" + inputFile.fileNoExtension() + ".tif";
//			}
//			if (outputFile.exists() && !bOverwrite)
//			{
//				printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
//					" \"-overwrite\" option.\n", outputFile.c_str());
//				Usage(0);
//			}
//			else
//			{
//				gcpFile.setExtension(gcpFile.ext());
//				//cout<<gcpFile<<endl;
//				cout << inputFile << endl;
//				clock_t  clockBegin, clockEnd;
//				clockBegin = clock();
//
//				orth(inputFile, outputFile,
//					gcpFile, projectionFile);
//				clockEnd = clock();
//				printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
//			}
//
//		}
//	}
//	else
//	{
//		//// 这里是调试时使用
//		//pszInputFile = "I:\\GF1-Testdata\\wuzhishan\\GF1_PMS2_E113.9_N37.7_20140526_L1A0000234945-MSS2.tiff";
//		////pszGcpFile = "I:\\GF1-Testdata\\wuzhishan\\gcp_auto.txt";
//		//pszInputFile = "I:\\testdata\\qb\\cangshan\\data\\052407462010_02_P001_PAN\\10OCT04030208-P2AS-052407462010_02_P001.tif";
//		//pszGcpFile = "I:\\testdata\\qb\\cangshan\\gcp_all.txt";
//		//pszProjectFile = "I:\\testdata\\qb\\cangshan\\gcp_all.geom";
//		//pszOutFile = "I:\\testdata\\qb\\cangshan\\out.tif";
//		////pszOutFile = "蛤.tif";
//		//bOverwrite = true;
//		//bRpc = true;
//		//ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
//		//if (preferences_file.exists())
//		//{
//		//	ossimPreferences::instance()->loadPreferences(preferences_file);
//		//}
//		//else
//		//{
//		//	preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
//		//	if (preferences_file.exists())
//		//	{
//		//		ossimPreferences::instance()->loadPreferences(preferences_file);
//		//	}
//		//}
//		//std::string	tempString;
//		//ossimArgumentParser::ossimParameter	stringParam(tempString);
//		//ossimArgumentParser argumentParser(&argc, argv);
//		//ossimInit::instance()->addOptions(argumentParser);
//		//argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
//		//ossimInit::instance()->initialize(argumentParser);
//
//		//ossimString strApplicationPath = ossimFilename(argumentParser.getApplicationName()).path();
//		//ossimFilename inputFile(pszInputFile);
//		//ossimFilename gcpFile(pszGcpFile);
//		//ossimFilename projectionFile(pszProjectFile);
//		//ossimFilename outputFile(pszOutFile);
//		//nThreads = 4;
//		//if (0 == strcmp(pszOutFile, ""))
//		//{
//		//	outputFile = inputFile.path() + "\\o" + inputFile.fileNoExtension() + ".tif";
//		//}
//		//if (outputFile.exists() && !bOverwrite)
//		//{
//		//	printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
//		//		" \"-overwrite\" option.\n", outputFile);
//		//	Usage(0);
//		//}
//		//else
//		//{
//		//	gcpFile.setExtension(gcpFile.ext());
//		//	clock_t  clockBegin, clockEnd;
//		//	clockBegin = clock();
//
//		//	orth(inputFile, outputFile,
//		//		gcpFile, projectionFile);
//		//	clockEnd = clock();
//		//	printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
//		//}
//		//Usage(0);
//	}
	return 0;
}