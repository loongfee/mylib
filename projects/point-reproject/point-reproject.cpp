
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
#include  <io.h>
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


const char *pszInFile = "";
const char *pszOutFile = "";
const char *pszProjectFile = "";
const char *pszOutProjectFile = "";
const char *pszPreferenceFile = "";
bool bOverwrite = true;
bool bGetElevation = false;

void reproject()
{
	//ossimTieGptSet* gcpSet = new ossimTieGptSet;
	//ossimTieGptSet* chkSet = new ossimTieGptSet;
	//ossimKeywordlist* prjKwl;
	//readGcpFile(pszInFile, gcpSet, chkSet, prjKwl);
	//ossimTieGptSet* gptSet = mylib::ReadGptFromFile(pszInFile);
	//projection2ll(pszInFile, pszOutFile, bGetElevation);

	ossimKeywordlist outPrjKwl;
	ossimFilename strOutProjectFile(pszOutProjectFile);
	if (strOutProjectFile.exists())
	{
		outPrjKwl.addFile(strOutProjectFile);
		//const char* out_ellipse_code = outPrjKwl.find(ossimKeywordNames::DATUM_KW);
		ossimRefPtr<ossimMapProjection> MapProjection;
		// 测试投影有效性
		if (!(MapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(outPrjKwl))))
		{
			// 如果无效，则清空
			cout << "Warning: 无效的投影文件" << endl;
			outPrjKwl.clear();
		}
	}
	else
	{
		projection2ll(pszInFile, pszOutFile, bGetElevation);
		return;
	}

	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist inPrjKwl;
	readGcpFile(pszInFile, gcpSet, chkSet, &inPrjKwl);

	get_elevation(gcpSet, inPrjKwl, 0.0);
	get_elevation(chkSet, inPrjKwl, 0.0);

	gcpSet = reprojectionPoints(gcpSet, inPrjKwl, outPrjKwl);
	chkSet = reprojectionPoints(chkSet, inPrjKwl, outPrjKwl);

	if (bGetElevation)
	{
		get_elevation(gcpSet, inPrjKwl, 0.0);
		get_elevation(chkSet, inPrjKwl, 0.0);
	}

	gcpSet = projection2ll(gcpSet, inPrjKwl);
	chkSet = projection2ll(chkSet, inPrjKwl);

	saveGcpFile(pszOutFile, gcpSet, chkSet, &outPrjKwl);	
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

	if (pszErrorMsg != NULL)
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

int main(int argc, char** argv)
{
	//ossimInit::instance()->loadPlugins("ossimgdal_plugin.dll");
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
		for (int i = 1; i < argc; i++)
		{
			if (0 == _stricmp(argv[i], "-o"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-i"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszInFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-oprj"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutProjectFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-hgt"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
				bGetElevation = true;
			}
			else if (0 == _stricmp(argv[i], "-pref"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszPreferenceFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-overwrite"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(0);
				bOverwrite = true;
			}
			else
			{
				Usage();
			}
		}


		if (0 == strcmp(pszInFile, ""))
		{
			printf("input file can not be empty!\n");
			Usage();
		}
		else if (0 == strcmp(pszOutFile, ""))
		{
			printf("output file can not be empty!\n");
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

			if (strcmp(pszInFile, pszOutFile) == 0)
			{
				Usage("Source and destination datasets must be different.");
			}

			if (_access(pszOutFile, 0) != -1 && !bOverwrite)
			{
				printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
					" \"-overwrite\" option.\n", pszOutFile);
				Usage(0);
			}


			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
			reproject();
			clockEnd = clock();
			printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
	}
	else
	{
		pszInFile = "E:\\Dropbox\\Document\\Notes\\Research\\L1 norm\\TGRS\\revise1\\test\\spot\\terrain-independent\\bj54_tm.txt";
		pszOutFile = "E:\\Dropbox\\Document\\Notes\\Research\\L1 norm\\TGRS\\revise1\\test\\spot\\terrain-independent\\wgs84.txt";
		pszOutProjectFile = "E:\\Dropbox\\Document\\Notes\\Research\\L1 norm\\TGRS\\revise1\\test\\spot\\terrain-independent\\wgs84_tm_prj.geom";
		pszPreferenceFile = "D:\\opensource\\ossim\\preference.txt";
		bGetElevation = true;


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

		if (strcmp(pszInFile, pszOutFile) == 0)
		{
			Usage("Source and destination datasets must be different.");
		}

		if (_access(pszOutFile, 0) != -1 && !bOverwrite)
		{
			printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
				" \"-overwrite\" option.\n", pszOutFile);
			Usage(0);
		}


		reproject();
		Usage(0);
	}
	return 0;
}