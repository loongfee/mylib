#include  <io.h>
#include  <stdio.h>
#include  <stdlib.h>
/////////////// gdal
#include "gdal_priv.h"
#include "ogr_srs_api.h"
#include <ogrsf_frmts.h>
#include <gdal.h>
#include <ogr_api.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"
#include <time.h>

#include <GdalRasterApp.h>

#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "ossim20.lib")
#pragma comment(lib, "opencv_ts300.lib")
#pragma comment(lib, "opencv_world300.lib")
#pragma comment(lib, "blas_win32_MT.lib")
#pragma comment(lib, "lapack_win32_MT.lib")

const char* pszInputShpFile;
const char* pszOutputShpFile;
const char* pszInputRasterFile;
double raster_min = 0;
double raster_max = 255;
bool bOverwrite = true;

struct fPoint{
	double x;
	double y;
};

void raster_filter()
{
#define GETVAL(papoSource, eSrcType) \
      (eSrcType == GDT_Byte ? \
          ((GByte *)papoSource)[0] : \
      (eSrcType == GDT_Float32 ? \
          ((float *)papoSource)[0] : \
      (eSrcType == GDT_Float64 ? \
          ((double *)papoSource)[0] : \
      (eSrcType == GDT_Int32 ? \
          ((GInt32 *)papoSource)[0] : \
      (eSrcType == GDT_UInt16 ? \
          ((GUInt16 *)papoSource)[0] : \
      (eSrcType == GDT_Int16 ? \
          ((GInt16 *)papoSource)[0] : \
      (eSrcType == GDT_UInt32 ? \
          ((GUInt32 *)papoSource)[0] : \
      (eSrcType == GDT_CInt16 ? \
          ((GInt16 *)papoSource)[0 * 2] : \
      (eSrcType == GDT_CInt32 ? \
          ((GInt32 *)papoSource)[0 * 2] : \
      (eSrcType == GDT_CFloat32 ? \
          ((float *)papoSource)[0 * 2] : \
      (eSrcType == GDT_CFloat64 ? \
          ((double *)papoSource)[0 * 2] : 0)))))))))))

	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");//得到shp文件的处理器
	OGRDataSource* poDS = poDriver->Open(pszInputShpFile, NULL);//打开文件
	if (!poDS)
	{
		printf("Failed to open the input shapefile:%s.\n", pszInputShpFile);
		return;
	}

	OGRLayer* poLayer = poDS->GetLayer(0);//获取shp图层
	OGRSpatialReference *inOSRS = poLayer->GetSpatialRef();
	OGRSpatialReference poLatLong;
	poLatLong.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *poTransform = OGRCreateCoordinateTransformation(inOSRS, &poLatLong);
	if (poTransform == NULL)
	{
		cout << "Projection Error!" << endl;
		system("Pause");
	}

	// create the output shapefile
	//GDAL写shp文件
	const char* nm = poLayer->GetName();
	if (_access(pszOutputShpFile, 0) != -1)
	{
		poDriver->DeleteDataSource(pszOutputShpFile);
	}
	OGRDataSource* poDS_out = poDriver->CreateDataSource(pszOutputShpFile, NULL);//创建shp文件
	if (poDS_out == NULL)
	{
		printf("Failed to create the output shapefile:%s.\n", pszOutputShpFile);
		return;
	}
	OGRLayer* poLayer_out = poDS_out->CreateLayer(poLayer->GetName(), inOSRS, wkbPoint, NULL);//创建图层

	// open the raster file
	mylib::GdalRasterApp gdalApp;
	gdalApp.open(pszInputRasterFile);
	gdalApp.setTileWidth(512);
	gdalApp.setTileHeight(512);

	//读取几何和属性值
	OGRFeature *poFeature;
	poLayer->ResetReading();//确保是从该层的开头开始
	int counter = 0;
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry *poGeometry;
		poGeometry = poFeature->GetGeometryRef();
		if (0 == counter)
		{
			OGRwkbGeometryType geoType = poGeometry->getGeometryType();
			if (poGeometry == NULL
				|| wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
			{
				printf("Failed: the shapefile must be point feature layer.\n");
				gdalApp.close();
				return;
			}

			OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
			int iField;
			for (iField = 0; iField < poFDefn->GetFieldCount(); iField++)
			{
				OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
				printf("%s=", poFieldDefn->GetNameRef());
				if (poFieldDefn->GetType() == OFTInteger)
				{
					printf("%d\n", poFeature->GetFieldAsInteger(iField));
				}
				else if (poFieldDefn->GetType() == OFTReal)
				{
					printf("%.3f\n", poFeature->GetFieldAsDouble(iField));
				}
				else if (poFieldDefn->GetType() == OFTString)
				{
					printf("%s\n", poFeature->GetFieldAsString(iField));
				}
				else
				{
					printf("%s\n", poFeature->GetFieldAsString(iField));
				}
				if (poLayer_out->CreateField(poFieldDefn) != OGRERR_NONE)
				{
					printf( "Creating Name field failed.\n" );
					return;
				}
			}

		}
		if (poGeometry != NULL
			&& wkbFlatten(poGeometry->getGeometryType()) == wkbPoint)//点文件
		{
			OGRPoint *poPoint = (OGRPoint *)poGeometry;
			//printf("%.3f,%.3f\n", poPoint->getX(), poPoint->getY());
			fPoint pt;
			pt.x = poPoint->getX();
			pt.y = poPoint->getY();
			if (!poTransform->Transform(1, &pt.x, &pt.y))
			{
				cout << "x=" << pt.x << "\t" << "y=" << pt.y << endl;
			}
			gdalApp.lonlat2world(pt);
			fPoint linesample;
			gdalApp.world2linesample(pt, linesample);
			GDALDataType dataType = gdalApp.getDataType();
			void *buf = NULL;
			buf = new GByte[64];
			gdalApp.getPixel(linesample.x, linesample.y, &buf);
			double data_value = GETVAL(buf, dataType);
			delete[]buf;
			buf = NULL;

			if (data_value >= raster_min && data_value <= raster_max)
			{
				cout << "lon=" << pt.x << "\t" << "lat=" << pt.y << endl;
				cout << "x=" << linesample.x << "\t" << "y=" << linesample.y << endl;
				cout << "height=" << data_value << endl;
				printf("\n");
				//创建几何和Feature
				OGRFeature *poFeature_out = new OGRFeature( poLayer_out->GetLayerDefn() );
				poFeature_out = poFeature->Clone();
				if (poLayer_out->CreateFeature(poFeature_out) != OGRERR_NONE)
				{
					printf("Failed to create feature in shapefile.\n");
					return;
				}
				OGRFeature::DestroyFeature(poFeature);
				OGRFeature::DestroyFeature(poFeature_out);
			}
			counter++;
		}
		else
		{
			printf("no geometry\n");
		}
		//printf("\n");
	}
	OGRDataSource::DestroyDataSource(poDS);
	OGRDataSource::DestroyDataSource(poDS_out);
	OGRCleanupAll();//资源清理
	return;
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


	printf("Usage: RasterFilter -i inputFile -o outputFile\n"
		"       -raster rasterFile -range min max \n");

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
	//OGRRegisterAll();//注册所有的文件格式驱动
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
				pszInputShpFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-o"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutputShpFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-raster"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszInputRasterFile = argv[++i];
			}
			else if (0 == _stricmp(argv[i], "-range"))
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(2);
				raster_min = atof(argv[++i]);
				raster_max = atof(argv[++i]);
			}
			else
			{
				Usage();
			}
		}

		if (0 == strcmp(pszInputShpFile, ""))
		{
			printf("input file can not be empty!\n");
			Usage();
		}
		else if (0 == strcmp(pszOutputShpFile, ""))
		{
			printf("output file can not be empty!\n");
			Usage();
		}
		else if (0 == strcmp(pszInputRasterFile, ""))
		{
			printf("output file can not be empty!\n");
			Usage();
		}
		else
		{
			if (strcmp(pszInputShpFile, pszOutputShpFile) == 0)
			{
				Usage("Source and destination datasets must be different.");
			}

			if (_access(pszOutputShpFile, 0) != -1 && !bOverwrite)
			{
				printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
					" \"-overwrite\" option.\n", pszOutputShpFile);
				Usage(0);
			}

			//if (!bFormatExplicitelySet)
			//	CheckExtensionConsistency(pszOutputFile, pszFormat);

			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
			raster_filter();
			clockEnd = clock();
			printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
	}
	else
	{
		pszInputShpFile = "D:\\workspace\\dem\\b4\\point_all.shp";
		pszInputRasterFile = "D:\\workspace\\dem\\b4\\clip1.tif";
		pszOutputShpFile = "D:\\workspace\\dem\\b4\\1000-2000.shp";
		raster_min = 1000;
		raster_max = 2000;
		raster_filter();
		Usage(0);
	}
	return 0;
}
