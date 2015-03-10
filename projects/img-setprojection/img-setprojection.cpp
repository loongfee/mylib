
/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"

#include <fileUtil.h>
using namespace mylib;

#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "mylib.lib")

int fix(const char* destImage, const char* refImage)
{	
	GDALDataset *poDataset;
	GDALDriver *poDriver;

	GDALAllRegister();	

	if (!QFileInfo(destImage).exists() || !QFileInfo(refImage).exists())
	{
		return 0;
	}


	GDALDataset *poRefDataset = (GDALDataset *) GDALOpen( refImage, GA_ReadOnly );
	double adfThisGeoTransform[6];
	poRefDataset->GetGeoTransform( adfThisGeoTransform );
	poDataset = (GDALDataset *) GDALOpen( destImage, GA_Update );
	poDataset->SetGeoTransform(adfThisGeoTransform);
	poDataset->SetProjection(poRefDataset->GetProjectionRef());

	GDALClose(poDataset);
	GDALClose(poRefDataset);

	return 1;
}

/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	//printf( 
	//	"Usage: img-reg -i inputFile -o outputFile [-trim percentage] [-g gamma]\n"
	//	" \n"
	//	"\n");

	printf( "Usage: img-setprojection -i inputFile -r refFile\n");

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

	const char* pszInputFile = "";
	const char* pszReferenceFile = "";
	if (argc > 1)
	{
		/* -------------------------------------------------------------------- */
		/*      Parse arguments.                                                */
		/* -------------------------------------------------------------------- */
		for( int i = 1; i < argc; i++ )
		{
			if( 0 == _stricmp(argv[i],"-i"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszInputFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-r") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszReferenceFile = argv[++i] ;
			}
			else
			{
				Usage();
			}
		}
		fix(pszInputFile, pszReferenceFile);
	}
	else
	{
		Usage();
	}
	return 0;
}