
/////// ossim
#include <ossim/base/ossimDirectory.h>
#include <ossim/base/ossimDirectoryTree.h>
#include <ossim/base/ossimString.h>
#include "ossim/base/ossimFilename.h"
#include "ossim/base/ossimString.h"
#include "ossim/base/ossimerrorcodes.h"
#include "ossim/imaging/ossimImageHandlerRegistry.h"
#include "ossim/imaging/ossimImageHandler.h"
#include "ossim/imaging/ossimImageFileWriter.h"
#include "ossim/imaging/ossimImageWriterFactoryRegistry.h"
#include <ossim/base/ossimDirectory.h>
#include <ossim/base/ossimThreeParamDatum.h>
#include <ossim/imaging/ossimFilterResampler.h>
#include "ossim/projection/ossimProjection.h"
#include <ossim/projection/ossimRpcProjection.h> 
#include "ossim/projection/ossimUtmProjection.h"
#include "ossim/projection/ossimTransMercatorProjection.h"
#include "ossim/projection/ossimProjectionFactoryRegistry.h"
#include "ossim/imaging/ossimImageRenderer.h"
#include "ossim/init/ossimInit.h"
#include "ossim/projection/ossimIkonosRpcModel.h"
#include "ossim/projection/ossimquickbirdrpcmodel.h"
#include "ossim/projection/ossimLandSatModel.h"
#include "ossim/support_data/ossimFfL5.h"
#include <ossim/support_data/ossimSpotDimapSupportData.h>
#include <ossim/projection/ossimSpot5Model.h> 
#include <ossim/projection/ossimProjection.h>
#include <ossim/projection/ossimMapProjectionFactory.h>
#include <ossim/projection/ossimProjectionFactoryRegistry.h>
#include <ossim/base/ossimStdOutProgress.h>
#include "ossim/base/ossimGpt.h"
#include "ossim/base/ossimDpt.h"
#include <ossim/imaging/ossimBandSelector.h>
#include <ossim/imaging/ossimPolyCutter.h>
#include <ossim/plugin/ossimSharedPluginRegistry.h>
#include <ossim\parallel\ossimMultiThreadSequencer.h>
#include <ossim/base/ossimPreferences.h>
#include <ossim/projection/ossimEquDistCylProjection.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>

/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"

#include <func.h>

#include <ossim/projection/ossimMapProjection.h>

#include <QDir>

#include <fstream>

using namespace std;
using namespace mylib;

#include <strUtil.h>

#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")
#pragma comment(lib, "ossim20.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "OpenThreads.lib")
#pragma comment(lib, "mlpack.lib")


const char *pszFilter = "*.tif";
const char *pszOutPath = "..\reprojection";
const char *pszDemPath = "";
const char *pszProjectFile = "projection.txt";
const char *pszPreferenceFile = "preference.txt";

int ImageReproject(const char* m_In, const char* m_Out, const char* geomFile,
				   vector<ossim_uint32> outBandList = vector<ossim_uint32>(0))
{
	ossimFilename outfile(m_Out);
	ossimKeywordlist in_geom_kwl,out_geom_kwl,tt_geom_kwl;

	
	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(ossimFilename(m_In));
	if(!handler)
	{
		//	cout << "Unable to open input image: "<< endl;
		return 0;
	}

	handler->getImageGeometry()->saveState(in_geom_kwl);
	//in_geom_kwl.add("pixel_type","area",true);
	//handler->loadState(in_geom_kwl);

	// 指定输出投影
	out_geom_kwl.addFile(geomFile);
	ossimRefPtr<ossimProjection> proj;
	proj = ossimProjectionFactoryRegistry::instance()->createProjection(out_geom_kwl);
	ossimMapProjection* pNewProjection = PTR_CAST(ossimMapProjection, proj.get());

	
	// 选择输出波段	
	ossimBandSelector* theBandSelector;
	theBandSelector = new ossimBandSelector;
	theBandSelector->connectMyInputTo(0, handler);	
	int nbands = handler->getNumberOfInputBands();
	//vector<ossim_uint32> outBandList(nbands);
	if (outBandList.size() < 1)
	{
		outBandList.clear();
		for (int i = 0; i < nbands; i++)
		{
			outBandList.push_back(i);
		}
	}
	theBandSelector->setOutputBandList(outBandList);
	
	// 输出指定范围图像
	vector<ossimDpt> polygon;
	ossimIrect bound,boundw;
	ossimPolyCutter* theCutter = new ossimPolyCutter;
	ossimDpt imagesize(handler->getImageRectangle().width(), handler->getImageRectangle().height());
	int starline,starpixel,endpixel,endline;
	starline	=	0;
	starpixel	=	0;
	endpixel	=	0;
	endline		=	0;
	if (0 == endpixel)
	{
		endpixel = imagesize.x - 1;
	}
	if (0 == endline)
	{
		endline = imagesize.y - 1;
	}/////////////////////////////////////////////////////////////////////
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
	

	ossimImageRenderer* renderer = new ossimImageRenderer;	
	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);
	renderer->setView(pNewProjection);
	renderer->connectMyInputTo(theCutter);
	
	ossimImageFileWriter* writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
	if(!writer)
	{
		//	cout << "Unable to writer input image: " << endl;
		return 0;
	}
	writer->saveState(tt_geom_kwl);
	tt_geom_kwl.add("pixel_type","area",true);
	writer->loadState(tt_geom_kwl);
	// cout<<"writer  "<<endl;
	//cout<<tt_geom;
	writer->setFilename(outfile);
	//ossimStdOutProgress progress(0,true);
	//writer->addListener(&progress);
	// 开始监听
	ossimStdOutProgress progress(0,true);
	writer->addListener(&progress);

	// 写图像

	int nThreads = OpenThreads::GetNumberOfProcessors() * 2;
	ossimRefPtr<ossimMultiThreadSequencer> sequencer = new ossimMultiThreadSequencer(0, nThreads);
	writer->changeSequencer(sequencer.get());
	writer->connectMyInputTo(0,renderer);
	writer->execute();

	// 结束监听

	// 清理
	//delete progress;
	//delete writer;
	writer->close();
	handler->close();
	//delete theBandSelector;
	in_geom_kwl.clear();
	out_geom_kwl.clear();
	tt_geom_kwl.clear();

	return 1;
}

void batchReprojection(QString inputPath, QStringList filters = QStringList("*.tif"),
					   QString outputPath = "", QString projectionFile = "projection.txt")
{	
	ossimInit::instance()->initialize();
	ossimFilename gdal_plugin = "ossimgdal_plugin.dll";
	if(!ossimSharedPluginRegistry::instance()->getPlugin(gdal_plugin))
	{
		ossimInit::instance()->loadPlugins(gdal_plugin);
	}

	QStringList TifFiles;
	mylib::QFindFile(inputPath, filters, TifFiles);

	int n = (int)TifFiles.size();
	if (n > 0)
	{
		if (!QDir(pszOutPath).exists())
		{
			_mkdir(pszOutPath);
		}
	}
	for (int i = 0; i < n; i++)
	{
		QString inFile = TifFiles[i];
		QString outFile = outputPath + "\\" + QFileInfo(inFile).fileName();
		cout<<i+1<<" / "<<n<<"..."<<endl;
		vector<ossim_uint32> outBandList;
		outBandList.push_back(4);
		outBandList.push_back(3);
		outBandList.push_back(2);
		if (!QFileInfo(outFile).exists())
		{
			ImageReproject(inFile.toLatin1(), outFile.toLatin1(), projectionFile.toLatin1());
		}
	}
}


/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: image_reprojection  [-o outpath] [-f filter] \n"
		"\t[-prj projectionfile] [-pref preferencefile]\n"
		"  -o outpath\t: output directory\n"
		"  -f filter\t: output directory\n"
		"  -prj projectionfile: output directory\n"
		"  -pref preferencefile: preference file\n");

	if( pszErrorMsg != NULL )
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

int main( int argc, char** argv )
{
	//ossimInit::instance()->loadPlugins("ossimgdal_plugin.dll");
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
	if (argc == 1)
	{
		batchReprojection("", QStringList(pszFilter), pszOutPath, pszProjectFile);
	}

	if (argc > 1)
	{
		/* -------------------------------------------------------------------- */
		/*      Parse arguments.                                                */
		/* -------------------------------------------------------------------- */
		for( int i = 1; i < argc; i++ )
		{
			if( 0 == _stricmp(argv[i],"-o") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutPath = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-i") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszInputFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-f") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszFilter = argv[++i] ;
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
			else
			{
				Usage();
			}
		}

		//batchReprojection("", QStringList(pszFilter), pszOutPath, pszProjectFile);
		if (!QDir(pszOutPath).exists())
		{
			_mkdir(pszOutPath);
		}

		QString outFile = QString(pszOutPath) + "\\" + QFileInfo(pszInputFile).fileName();

		ossimPreferences::instance()->loadPreferences(ossimFilename(pszPreferenceFile));
		std::string	tempString;
		ossimArgumentParser::ossimParameter	stringParam(tempString);
		ossimArgumentParser argumentParser(&argc, argv);
		ossimInit::instance()->addOptions(argumentParser);
		argumentParser.getApplicationUsage()->setApplicationName(argumentParser.getApplicationName());
		ossimInit::instance()->initialize(argumentParser);
		ImageReproject(pszInputFile, outFile.toLatin1(), pszProjectFile);
	}
	return 0;
}