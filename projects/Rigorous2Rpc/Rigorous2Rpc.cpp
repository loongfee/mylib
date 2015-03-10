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
#include <ogr_spatialref.h>

#include <ossim/projection/ossimMapProjection.h>
#include <ossim/base/ossimPreferences.h>
#include <ossim/parallel/ossimMultiThreadSequencer.h>
#include <ossim/base/ossimStringProperty.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>
#include <ossim/base/ossimXmlDocument.h>
#include <ossim/base/ossimXmlAttribute.h>
#include <ossim/base/ossimXmlNode.h>
#include <ossim_plugin/radi/radiRpcSolver.h>

#include <fstream>
#include <func.h>

#include <fileUtil.h>
#include <mprojectdefine.h>
#include <time.h>
#include <tiff.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>


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
#pragma comment(lib, "mlpack.lib")


const char *pszPreferenceFile = "preference.txt";

string strInputImageFile = "";
string strOutputRpcFile = "";


void relative2absolute(ossimFilename &strPath, boost::filesystem::path current_path)
{
	boost::filesystem::path p(strPath.c_str());
	if (p.is_relative())
	{
		boost::filesystem::path absolute = boost::filesystem::absolute(p, current_path).normalize();
		strPath = absolute.string().c_str();
	}
}

void relative2absolute(string &strPath, boost::filesystem::path current_path)
{
	boost::filesystem::path p(strPath.c_str());
	if (p.is_relative())
	{
		boost::filesystem::path absolute = boost::filesystem::absolute(p, current_path).normalize();
		strPath = absolute.string().c_str();
	}
}

bool Rigorous2Rpc()
{
	ossimRefPtr<ossimImageHandler> handler = ossimImageHandlerRegistry::instance()->open(ossimFilename(strInputImageFile.c_str()));
	if (!handler)
	{
		cout << "Failed to open image file: " << strInputImageFile << endl;
		return false;
	}

	ossimRefPtr<ossimSensorModel> sensorModel = PTR_CAST(ossimSensorModel, handler->getImageGeometry()->getProjection());

	if (!sensorModel)
	{
		cout << "Failed to find the sensor model of image file: " << strInputImageFile << endl;
		return false;
	}


	ossimTieGptSet* gptSet = new ossimTieGptSet;
	double max_height = 1100.0;
	double min_height = 0.0;

	ossimIrect imageRect(handler->getImageRectangle().ul(),
		handler->getImageRectangle().lr());
	int nLevels = 5;
	for (int i = 0; i < nLevels; ++i)
	{
		double hgt = min_height + i*(max_height - min_height) / max(1, (nLevels - 1));
		mylib::create3DGridPoints(imageRect, *sensorModel, hgt, gptSet, 15, 15, true, false);
	}

	ossimTieGptSet* chkSet = new ossimTieGptSet;
	nLevels = 10;
	for (int i = 0; i < nLevels; ++i)
	{
		double hgt = min_height + i*(max_height - min_height) / max(1, (nLevels - 1));
		mylib::create3DGridPoints(imageRect, *sensorModel, hgt, chkSet, 30, 30, true, false);
	}
	handler->close();

	ossimRefPtr<ossimplugins::radiRpcSolver> rpcSolver = new ossimplugins::radiRpcSolver(true);
	rpcSolver->solveCoefficients(gptSet, true, ossimplugins::radiRpcSolver::EstimationMethodIndex::LASSO, 1e-5);
	ossimRefPtr<ossimplugins::radiRpcModel> rpcModel = new ossimplugins::radiRpcModel();
	rpcSolver->setRpcModel(rpcModel.get());

	rpcModel->writeRpcFile(ossimFilename(strOutputRpcFile.c_str()));

	MyProject prj;
	boost::filesystem::path p(strOutputRpcFile);
	string reportfile = p.remove_filename().string() + "\\rpcReport.txt";
	prj.OutputReport(ossimFilename(reportfile.c_str()), rpcModel.get(), gptSet, chkSet, true);

	return true;
}


/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf(
		"Usage: Rigorous2Rpc -i inputFile [-o refFile] [-pref preferencefile]\n");

	if (pszErrorMsg != NULL)
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

int main(int argc, char** argv)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持
	
	if (argc > 1)
	{
		/* -------------------------------------------------------------------- */
		/*      Parse arguments.                                                */
		/* -------------------------------------------------------------------- */
		for (int i = 1; i < argc; i++)
		{
			if (0 == _stricmp(argv[i], "-i"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				strInputImageFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-o"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				strOutputRpcFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-pref"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszPreferenceFile = argv[++i];
			}
			else
			{
				Usage();
			}
		}

		if (0 == strcmp(strInputImageFile.c_str(), ""))
		{
			printf("input file can not be empty!\n");
			Usage();
		}
		else
		{
			relative2absolute(strInputImageFile, boost::filesystem::current_path());
			if (0 == strcmp(strOutputRpcFile.c_str(), ""))
			{
				boost::filesystem::path p(strInputImageFile);
				strOutputRpcFile = p.replace_extension("rpb").string();
			}
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
			ossimInit::instance()->initialize();
			Rigorous2Rpc();
		}
	}
	else
	{
		Usage();
	}
	return 0;
}