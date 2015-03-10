#include  <io.h>
#include  <stdio.h>
#include  <stdlib.h>
/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"
#include "commonutils.h"
#include "vrtdataset.h"
#include "jpeglib.h"

#include <gcpUtil.h>
#include <time.h>

#pragma comment(lib, "OpenThreads.lib")
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "jpeg.lib")

using namespace mylib;

const char* pszInputFile = "";
const char* pszOutputFile = "";
double trim_percentage = 1.0;
int blockSize = 1024;
bool bOverwrite = false;
int nBins = 256;
double null_value = 0.0;
GByte out_null_value = 0;
bool bIgnorNull = true;
bool bStretching = true;
//const int nBins = 65536;
double gamma = 1.0;
double statistics_scale = 0;
int nThreads = 0;
int maxThreads = 10;
const char* pszOXSize = NULL;
const char* pszOYSize = NULL;
int outSizeX = 0;
int outSizeY = 0;
int nOutBands = 0;
vector<int> outBandList;
double out_scale = 1.0;
const char* pszFormat = "GTiff";
bool bFormatExplicitelySet = false;
char** papszCreateOptions = NULL;
char** papszMetadataOptions = NULL;

vector<double> nBinStepList;
vector<double> gMinValueList;
vector<double> gMaxValueList;
vector<vector<long long> > gHistogramList;
vector<double> vLowList;
vector<double> vHighList;
vector<double> vLenList;

static GDALDataset  *pSrcDataset = NULL;
static GDALDataset  *pDstDataset = NULL;
static jpeg_compress_struct sCInfo;
static jpeg_error_mgr sJErr;

static FILE *fpJpgImage = NULL;
static bool bJpg = false;
//static VRTDataset *poVDS = NULL;
static int totalBlocks;
static int finishedBlocks;


#include <OpenThreads/Thread>
#include <OpenThreads/Mutex>
#include <OpenThreads/Barrier>

struct _RECT
{
	int startSample;
	int startLine;
	int endSample;
	int endLine;
	_RECT(int sSample, int sLine, int eSample, int eLine)
		:startSample(sSample),
		startLine(sLine),
		endSample(eSample),
		endLine(eLine)
	{
		
	};
};
static OpenThreads::Barrier bar;
static int GLOBAL_NUM_THREADS;

enum processType{
	ExtremeStatistics = 0,
	HistogramStatistcs,
	Stretching,
	Translate,
};

char* getTempImageFilename()
{
	char temp_img_file_name[MAX_PATH];
	char temp_path[MAX_PATH];

	GetTempPath(sizeof(temp_path), (wchar_t*)temp_path);
	GetTempFileName((wchar_t*)temp_path, (wchar_t*)"", 0, (wchar_t*)temp_img_file_name);
	//GetTempPath(sizeof(temp_path), temp_path);
	//GetTempFileName(temp_path, "", 0, temp_img_file_name);

	//QString absolutePath = QFileInfo(temp_img_file_name).absolutePath();

	return temp_img_file_name;
}

bool initializeJpg(int nXSize, int nYSize)
{
	int  nBands = (int)outBandList.size();

	int  nQuality = 75;
	int  bProgressive = FALSE;

	// --------------------------------------------------
	//      Some somerudimentary checks                  
	// --------------------------------------------------
	if( nBands != 1&& nBands != 3 )
	{
		CPLError( CE_Failure,CPLE_NotSupported,
			"JPEG driver doesn't support %d bands.  Must be 1 (grey) "
			"or 3 (RGB) bands.\n", nBands );

		return false;
	}
	// --------------------------------------------------
	//      What optionshas the user selected?                            
	// --------------------------------------------------
	if( CSLFetchNameValue(papszCreateOptions,"QUALITY")!= NULL )
	{
		nQuality = atoi(CSLFetchNameValue(papszCreateOptions,"QUALITY"));
		if( nQuality <10 || nQuality > 100 )
		{
			CPLError( CE_Failure,CPLE_IllegalArg,
				"QUALITY=%s is not a legal value in the range10-100.",
				CSLFetchNameValue(papszCreateOptions,"QUALITY") );
			return false;
		}
	}

	if( CSLFetchNameValue(papszCreateOptions,"PROGRESSIVE")!= NULL )
	{
		bProgressive = TRUE;
	}

	// --------------------------------------------------
	//      Create thedataset.                                            
	// --------------------------------------------------

	fpJpgImage = VSIFOpen(pszOutputFile, "wb");
	if( fpJpgImage == NULL )
	{
		CPLError( CE_Failure,CPLE_OpenFailed,
			"Unable to create jpeg file %s.\n",
			pszOutputFile );
		return false;
	}

	// --------------------------------------------------
	//      InitializeJPG access to the file.                             
	// --------------------------------------------------
	struct jpeg_compress_structsCInfo;
	struct jpeg_error_mgrsJErr;

	sCInfo.err = jpeg_std_error( &sJErr);
	jpeg_create_compress( &sCInfo );

	jpeg_stdio_dest( &sCInfo,fpJpgImage );

	sCInfo.image_width= nXSize;
	sCInfo.image_height= nYSize;
	sCInfo.input_components= nBands;

	if( nBands == 1 )
	{
		sCInfo.in_color_space= JCS_GRAYSCALE;
	}
	else
	{
		sCInfo.in_color_space= JCS_RGB;
	}

	jpeg_set_defaults( &sCInfo);

	jpeg_set_quality( &sCInfo,nQuality, TRUE);

	if( bProgressive )
		jpeg_simple_progression( &sCInfo );

	jpeg_start_compress( &sCInfo,TRUE );

	return true;
}

//static GDALDataset *
//	JPEGCreateCopy( const char * pszFilename, GDALDataset*poSrcDS,
//	int bStrict, char ** papszOptions, void* pProgressData )
//{
//	jpeg_compress_struct sCInfo;
//	jpeg_error_mgr sJErr;
//	int  nBands = poSrcDS->GetRasterCount();
//	int  nXSize = poSrcDS->GetRasterXSize();
//	int  nYSize = poSrcDS->GetRasterYSize();
//	int  nQuality = 75;
//	int  bProgressive = FALSE;
//
//	// --------------------------------------------------
//	//      Some somerudimentary checks                  
//	// --------------------------------------------------
//	if( nBands != 1&& nBands != 3 )
//	{
//		CPLError( CE_Failure,CPLE_NotSupported,
//			"JPEG driver doesn't support %d bands.  Must be 1 (grey) "
//			"or 3 (RGB) bands.\n", nBands );
//
//		return NULL;
//	}
//
//	if( poSrcDS->GetRasterBand(1)->GetRasterDataType()!= GDT_Byte && bStrict )
//	{
//		CPLError( CE_Failure,CPLE_NotSupported,
//			"JPEG driver doesn't support data type %s. "
//			"Only eight bit byte bands supported.\n",
//			GDALGetDataTypeName(
//			poSrcDS->GetRasterBand(1)->GetRasterDataType()) );
//
//		return NULL;
//	}
//
//	// --------------------------------------------------
//	//      What optionshas the user selected?                            
//	// --------------------------------------------------
//	if( CSLFetchNameValue(papszOptions,"QUALITY")!= NULL )
//	{
//		nQuality = atoi(CSLFetchNameValue(papszOptions,"QUALITY"));
//		if( nQuality <10 || nQuality > 100 )
//		{
//			CPLError( CE_Failure,CPLE_IllegalArg,
//				"QUALITY=%s is not a legal value in the range10-100.",
//				CSLFetchNameValue(papszOptions,"QUALITY") );
//			return NULL;
//		}
//	}
//
//	if( CSLFetchNameValue(papszOptions,"PROGRESSIVE")!= NULL )
//	{
//		bProgressive = TRUE;
//	}
//
//	// --------------------------------------------------
//	//      Create thedataset.                                            
//	// --------------------------------------------------
//	FILE    *fpImage;
//
//	fpImage = VSIFOpen(pszFilename, "wb");
//	if( fpImage == NULL )
//	{
//		CPLError( CE_Failure,CPLE_OpenFailed,
//			"Unable to create jpeg file %s.\n",
//			pszFilename );
//		return NULL;
//	}
//
//	// --------------------------------------------------
//	//      InitializeJPG access to the file.                             
//	// --------------------------------------------------
//	struct jpeg_compress_structsCInfo;
//	struct jpeg_error_mgrsJErr;
//
//	sCInfo.err = jpeg_std_error( &sJErr);
//	jpeg_create_compress( &sCInfo );
//
//	jpeg_stdio_dest( &sCInfo,fpImage );
//
//	sCInfo.image_width= nXSize;
//	sCInfo.image_height= nYSize;
//	sCInfo.input_components= nBands;
//
//	if( nBands == 1 )
//	{
//		sCInfo.in_color_space= JCS_GRAYSCALE;
//	}
//	else
//	{
//		sCInfo.in_color_space= JCS_RGB;
//	}
//
//	jpeg_set_defaults( &sCInfo);
//
//	jpeg_set_quality( &sCInfo,nQuality, TRUE);
//
//	if( bProgressive )
//		jpeg_simple_progression( &sCInfo );
//
//	jpeg_start_compress( &sCInfo,TRUE );
//
//	// --------------------------------------------------
//	//      Loop overimage, copying image data.                           
//	// --------------------------------------------------
//	GByte   *pabyScanline;
//	CPLErr      eErr;
//
//	pabyScanline = (GByte*) CPLMalloc( nBands* nXSize );
//
//	for( int iLine = 0; iLine< nYSize; iLine++)
//	{
//		JSAMPLE     *ppSamples;
//
//		for( int iBand = 0; iBand< nBands; iBand++)
//		{
//			GDALRasterBand * poBand= poSrcDS->GetRasterBand(iBand+1 );
//			eErr = poBand->RasterIO( GF_Read,0, iLine, nXSize,1,
//				pabyScanline + iBand,nXSize, 1, GDT_Byte,
//				nBands, nBands* nXSize );
//		}
//
//		ppSamples = pabyScanline;
//		jpeg_write_scanlines( &sCInfo, &ppSamples, 1 );
//	}
//
//	CPLFree( pabyScanline);
//
//	jpeg_finish_compress( &sCInfo );
//	jpeg_destroy_compress( &sCInfo );
//
//	VSIFClose( fpImage);
//
//	return (GDALDataset*) GDALOpen( pszFilename,GA_ReadOnly );
//}

/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	//printf( 
	//	"Usage: img-reg -i inputFile -o outputFile [-trim percentage] [-g gamma]\n"
	//	" \n"
	//	"\n");


	printf( "Usage: img-translate -i inputFile -o outputFile\n"
		"       [-trim percentage] [-g gamma] [-overwrite]\n"
		"       [-nv null_value] [-onv out_null_value] [-in {0|1}]\n"
		"       [-b nb b1 b2 ...] [-stretching {0|1}]\n"
		"       [-nt nThreads] [-bs blocksize] [-ss statistics_scale]\n"
		"       [-osx xsize[%%]] [-osy ysize[%%]]\n");
	printf( "   -trim percentage: defualt value is %lf\n", trim_percentage);
	printf( "   -g gamma: defualt value is %lf.\n"
		"         if gamma is greater than 1 increase contrast, otherwise decrease contrast.\n", gamma);
	printf( "   -stretching: defualt value is 1.\n");
	printf( "   -in {0|1}: whether ignor null value, default value is 1.\n");

	if( pszErrorMsg != NULL )
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)


/************************************************************************/
/*                           AttachMetadata()                           */
/************************************************************************/

static void AttachMetadata( GDALDatasetH hDS, char **papszMetadataOptions )

{
	int nCount = CSLCount(papszMetadataOptions);
	int i;

	for( i = 0; i < nCount; i++ )
	{
		char    *pszKey = NULL;
		const char *pszValue;

		pszValue = CPLParseNameValue( papszMetadataOptions[i], &pszKey );
		GDALSetMetadataItem(hDS,pszKey,pszValue,NULL);
		CPLFree( pszKey );
	}

	CSLDestroy( papszMetadataOptions );
}


class processThread: public OpenThreads::Thread
{
public:
	processThread()
		: OpenThreads::Thread(){};
	virtual  ~processThread()
	{
		//m_pRegistration = NULL;
	};
	void setData(vector<_RECT> rectList, processType tp)
	{
		m_rectList = rectList;
		m_processType = tp;
	}
	virtual void run()
	{
		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.
		if (processType::ExtremeStatistics == m_processType)
		{
			for (int i = 0;i < (int)m_rectList.size();++i)
			{
				extremeStatistics(m_rectList[i]);
			}
		}
		else if (processType::HistogramStatistcs == m_processType)
		{
			for (int i = 0;i < (int)m_rectList.size();++i)
			{
				histogramStatistics(m_rectList[i]);
			}
		}
		else if (processType::Stretching == m_processType)
		{
			for (int i = 0;i < (int)m_rectList.size();++i)
			{
				stretching(m_rectList[i]);
			}
		}
		else if (processType::Translate == m_processType)
		{
			for (int i = 0;i < (int)m_rectList.size();++i)
			{
				translate(m_rectList[i]);
			}
		}
		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.

		//---------------------------------------------------------------------
		// Now that we've done our work, wait for a sign that we should quit.
		//
		while (true) {

			_quitmutex.lock();
			if(_quitflag == true) break;
			_quitmutex.unlock();

			OpenThreads::Thread::YieldCurrentThread();
		}
	};

	void extremeStatistics(_RECT rect)
	{
		GDALDatasetH tmpHSrcDs = GDALOpenShared(pszInputFile, GA_ReadOnly );

		//int nBands = pSrcDataset->GetRasterCount();
		int nBands = (int)outBandList.size();
		GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
		int dataTypeSize = GDALGetDataTypeSize(eDT)/8;

		int nRasterXSize = rect.endSample - rect.startSample + 1;
		int nRasterYSize = rect.endLine - rect.startLine + 1;
		int nStatisticsRasterXSize = nRasterXSize * statistics_scale;
		int nStatisticsRasterYSize = nRasterYSize * statistics_scale;
		GByte* pBuf = new GByte[nStatisticsRasterXSize * nStatisticsRasterYSize * dataTypeSize];

		vector<double> minValue;
		vector<double> maxValue;
		for (int ib = 0;ib < nBands;++ib)
		{
			minValue.push_back(1e10);
			maxValue.push_back(-1e10);
		}

		// find min and max values
		for (int ib = 0;ib < nBands;++ib)
		{
			//((GDALDataset*)tmpHSrcDs)->RasterIO(GF_Read, 0, rect.startLine, nRasterXSize, nRasterYSize, 
			//	pBuf, nStatisticsRasterXSize, nStatisticsRasterYSize, eDT, nBands,
			//	0, 0, 0, 0);
			((GDALDataset*)tmpHSrcDs)->GetRasterBand(outBandList[ib])->RasterIO(GF_Read, 0, rect.startLine, nRasterXSize, nRasterYSize, 
				pBuf, nStatisticsRasterXSize, nStatisticsRasterYSize, eDT, 0, 0);
			for (int i = 0; i < nStatisticsRasterYSize; i++)
			{
				for (int j = 0;j < nStatisticsRasterXSize;++j)
				{
					int pos = i * nStatisticsRasterXSize + j;
					double data = SRCVAL(pBuf, eDT, pos);
					if (bIgnorNull && data == null_value)
					{
						continue;
					}
					if (data < minValue[ib])
					{
						minValue[ib] = data;
					}
					if (data > maxValue[ib])
					{
						maxValue[ib] = data;
					}
				}
			}
		}

		CPLFree( pBuf );
		GDALClose(tmpHSrcDs);

		// reduce
		_quitmutex.lock();
		for (int ib = 0;ib < nBands;++ib)
		{
			gMinValueList[ib] = min(minValue[ib], gMinValueList[ib]);
			gMaxValueList[ib] = max(maxValue[ib], gMaxValueList[ib]);
		}
		finishedBlocks++;
		_quitmutex.unlock();

		int percent = (int)(finishedBlocks / (double)(totalBlocks) * 25 + 0.5);
		printf("\r%3d%%", percent);
		fflush( stdout );
	};


	void histogramStatistics(_RECT rect)
	{
		GDALDatasetH tmpHSrcDs = GDALOpenShared(pszInputFile, GA_ReadOnly );
		//int nBands = pSrcDataset->GetRasterCount();
		int nBands = (int)outBandList.size();
		GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
		int dataTypeSize = GDALGetDataTypeSize(eDT)/8;

		int nRasterXSize = rect.endSample - rect.startSample + 1;
		int nRasterYSize = rect.endLine - rect.startLine + 1;
		int nStatisticsRasterXSize = nRasterXSize * statistics_scale;
		int nStatisticsRasterYSize = nRasterYSize * statistics_scale;
		GByte* pBuf = new GByte[nStatisticsRasterXSize * nStatisticsRasterYSize * dataTypeSize];


		//((GDALDataset*)tmpHSrcDs)->RasterIO(GF_Read, 0, rect.startLine, nRasterXSize, nRasterYSize, pBuf, 
		//	nStatisticsRasterXSize, nStatisticsRasterYSize, eDT, nBands,
		//	0, 0, 0, 0);

		// histogram statistics
		vector<vector<long long> > histogramList;
		for (int ib = 0;ib < nBands;++ib)
		{
			((GDALDataset*)tmpHSrcDs)->GetRasterBand(outBandList[ib])->RasterIO(GF_Read, 0, rect.startLine, nRasterXSize, nRasterYSize, 
				pBuf, nStatisticsRasterXSize, nStatisticsRasterYSize, eDT, 0, 0);
			vector<long long> histogram(nBins, 0.0);
			for (int i = 0; i < nStatisticsRasterYSize; i++)
			{
				for (int j = 0;j < nStatisticsRasterXSize;++j)
				{
					int pos = i * nStatisticsRasterXSize + j;
					//int pos = ib*nStatisticsRasterXSize*nStatisticsRasterYSize + i * nStatisticsRasterXSize + j;
					double data = SRCVAL(pBuf, eDT, pos);
					if (bIgnorNull && data == null_value)
					{
						continue;
					}
					// scale to bins
					int v = int((data - gMinValueList[ib])/nBinStepList[ib]+0.5);
					v = v > (nBins-1) ? (nBins-1) : v;
					v = v < 0 ? 0 : v;
					histogram[v]++;
				}
			}
			histogramList.push_back(histogram);
		}

		CPLFree( pBuf );
		GDALClose(tmpHSrcDs);

		// reduce
		_quitmutex.lock();
		for (int ib = 0;ib < nBands;++ib)
		{
			for (int i = 0;i < nBins;++i)
			{
				gHistogramList[ib][i] += histogramList[ib][i];
			}
		}
		finishedBlocks++;
		_quitmutex.unlock();

		int percent = 25+(int)(finishedBlocks / (double)(totalBlocks) * 25 + 0.5);
		printf("\r%3d%%", percent);
		fflush( stdout );
	};


	void stretching(_RECT rect)
	{
		GDALDatasetH tmpHSrcDs = GDALOpenShared(pszInputFile, GA_ReadOnly );
		//int nBands = pSrcDataset->GetRasterCount();
		int nBands = (int)outBandList.size();
		GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
		int dataTypeSize = GDALGetDataTypeSize(eDT)/8;

		int nRasterXSize = rect.endSample - rect.startSample + 1;
		int nRasterYSize = rect.endLine - rect.startLine + 1;
		int nOutRasterXSize = int(nRasterXSize * out_scale + 0.5);
		int nOutRasterYSize = int(nRasterYSize * out_scale + 0.5);
		int nOutBlockSize = int(blockSize * out_scale + 0.5);
		int nOutStartLine = rect.startLine / blockSize * nOutBlockSize;
		
		//((GDALDataset*)tmpHSrcDs)->RasterIO(GF_Read, 0, rect.startLine, nRasterXSize, nRasterYSize, pBuf,
		//	nOutRasterXSize, nOutRasterYSize, eDT, nBands,
		//	0, 0, 0, 0);


		// stretch to byte

		GByte* pBuf = new GByte[nOutRasterXSize * nOutRasterYSize * dataTypeSize];
		GByte* pNewBuf = new GByte[nOutRasterXSize * nOutRasterYSize];
		for (int ib = 0;ib < nBands;++ib)
		{
			((GDALDataset*)tmpHSrcDs)->GetRasterBand(outBandList[ib])->RasterIO(GF_Read, 0, rect.startLine, nRasterXSize, nRasterYSize, 
				pBuf, nOutRasterXSize, nOutRasterYSize, eDT, 0, 0);
			for (int i = 0; i < nOutRasterYSize; i++)
			{
				for (int j = 0;j < nOutRasterXSize;++j)
				{
					int pos = i * nOutRasterXSize + j;
					//int pos = ib*nOutRasterXSize*nOutRasterYSize + i * nOutRasterXSize + j;
					double data = SRCVAL(pBuf, eDT, pos);

					if (bIgnorNull && data == null_value)
					{
						pNewBuf[pos] = out_null_value;
						continue;
					}
					int oData = (int)(pow((data - vLowList[ib]) / vLenList[ib], gamma) * 255.0 + 0.5);
					oData = oData > 255 ? 255 : oData;
					oData = oData < 0 ? 0 : oData;
					pNewBuf[pos] = (GByte)oData;
				}
			}
			_quitmutex.lock();
			pDstDataset->GetRasterBand( ib+1 )->RasterIO(GF_Write, 0, nOutStartLine, nOutRasterXSize, nOutRasterYSize, pNewBuf, 
				nOutRasterXSize, nOutRasterYSize, GDT_Byte, 0, 0);
			_quitmutex.unlock();
			//poVRTBand = NULL;
		}
		CPLFree( pBuf );
		CPLFree( pNewBuf );
		// reduce
		_quitmutex.lock();
		finishedBlocks++;
		_quitmutex.unlock();

		GDALClose(tmpHSrcDs);

		int percent = 50+(int)(finishedBlocks / (double)(totalBlocks) * 50 + 0.5);
		printf("\r%3d%%", percent);
		fflush( stdout );
	};


	void translate(_RECT rect)
	{
		GDALDatasetH tmpHSrcDs = GDALOpenShared(pszInputFile, GA_ReadOnly );
		//int nBands = pSrcDataset->GetRasterCount();
		int nBands = (int)outBandList.size();
		GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
		int dataTypeSize = GDALGetDataTypeSize(eDT)/8;

		int nRasterXSize = rect.endSample - rect.startSample + 1;
		int nRasterYSize = rect.endLine - rect.startLine + 1;
		int nOutRasterXSize = int(nRasterXSize * out_scale + 0.5);
		int nOutRasterYSize = int(nRasterYSize * out_scale + 0.5);
		int nOutBlockSize = int(blockSize * out_scale + 0.5);
		int nOutStartLine = rect.startLine / blockSize * nOutBlockSize;
		GByte* pBuf = new GByte[nOutRasterXSize * nOutRasterYSize * dataTypeSize];


		// stretch to byte
		for (int ib = 0;ib < nBands;++ib)
		{
			((GDALDataset*)tmpHSrcDs)->GetRasterBand(outBandList[ib])->RasterIO(GF_Read, 0, rect.startLine, nRasterXSize, nRasterYSize, 
				pBuf, nOutRasterXSize, nOutRasterYSize, eDT, 0, 0);

			_quitmutex.lock();
			pDstDataset->GetRasterBand( ib+1 )->RasterIO(GF_Write, 0, nOutStartLine, nOutRasterXSize, nOutRasterYSize, pBuf, 
				nOutRasterXSize, nOutRasterYSize, eDT, 0, 0);
			_quitmutex.unlock();
		}

		// reduce
		_quitmutex.lock();
		finishedBlocks++;
		_quitmutex.unlock();

		CPLFree( pBuf );
		GDALClose(tmpHSrcDs);

		int percent = 50+(int)(finishedBlocks / (double)(totalBlocks) * 50 + 0.5);
		printf("\r%3d%%", percent);
		fflush( stdout );
	};

	void quit() {
		_quitmutex.lock();
		_quitflag = true;
		_quitmutex.unlock();
	};
private:
	vector<_RECT> m_rectList;
	//int m_startLine;
	//int m_endLine;
	//vector<double> m_minValue;
	//vector<double> m_maxValue;
	processType m_processType;
	volatile bool _quitflag;
	OpenThreads::Mutex _quitmutex;
};

void img_stretching()
{
	GDALDatasetH hSrcDs = GDALOpenShared(pszInputFile, GA_ReadOnly );
	pSrcDataset = (GDALDataset  *)hSrcDs;
	int nRasterXSizeRead = pSrcDataset->GetRasterXSize();
	int nRasterYSizeRead = pSrcDataset->GetRasterYSize();
	int nBandTotal = pSrcDataset->GetRasterCount();
	GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();

	// get output size
	if( pszOXSize != NULL)
	{
		outSizeX = (int) ((pszOXSize[strlen(pszOXSize)-1]=='%' 
			? CPLAtofM(pszOXSize)/100*nRasterXSizeRead : atoi(pszOXSize)));
	}
	if(pszOYSize != NULL)
	{
		outSizeY = (int) ((pszOYSize[strlen(pszOYSize)-1]=='%' 
			? CPLAtofM(pszOYSize)/100*nRasterYSizeRead : atoi(pszOYSize)));
	}
	if (0 < outSizeX)
	{
		out_scale = outSizeX / (double)nRasterXSizeRead;
	}
	else if (0 < outSizeY)
	{
		out_scale = outSizeY / (double)nRasterYSizeRead;
	}

	// get statistics scale
	if (0.0 == statistics_scale)
	{
		statistics_scale = min(1024/(double)nRasterXSizeRead, 1024/(double)nRasterYSizeRead);
	}
	if (statistics_scale > 1.0)
	{
		statistics_scale = 1.0;
	}

	// check band list
	nOutBands = (int)outBandList.size();
	for (int i = 0;i < nOutBands;i++)
	{
		if (outBandList[i] < 1 || outBandList[i] > nBandTotal)
		{
			printf("warning: band %d is out of range (1~%d).\n", outBandList[i], nBandTotal);
			outBandList.erase(outBandList.begin()+i);
			i--;
		}
	}
	nOutBands = (int)outBandList.size();
	if (nOutBands < 1)
	{
		//printf("warning: output band number cannot less than 1.\n"
		//	"Therefore, all the bands of the input file will be outputted.\n");
		
		for (int i = 0;i < nBandTotal;++i)
		{
			outBandList.push_back(i+1);
		}
		nOutBands = (int)outBandList.size();
	}

	int num_threads = nThreads;
	if (num_threads < 1)
	{
		num_threads = OpenThreads::GetNumberOfProcessors() * 2;
	}
	if (num_threads > maxThreads)
	{
		num_threads = maxThreads;
	}
	//num_threads = 3;
	cout<<"using "<<num_threads<<" threads..."<<endl;
	GLOBAL_NUM_THREADS = num_threads + 1;
	std::vector<processThread *> threads(num_threads);
	OpenThreads::Thread::SetConcurrency(num_threads);
	OpenThreads::Thread::Init();

	//设置光标位置
	printf("\r%3d%%", 0);

	// 0: initialize
	for (int ib = 0;ib < nOutBands;++ib)
	{
		gMinValueList.push_back(1e10);
		gMaxValueList.push_back(-1e10);
	}

	// counter start
	finishedBlocks = 0;

	totalBlocks = ceil(nRasterYSizeRead / (double)blockSize );
	vector<int> blockList;
	for (int i = 0;i < totalBlocks;++i)
	{
		blockList.push_back(i);
	}

	// 1: extreme statistics
	for(int i = 0;i < num_threads;++i)
	{

		threads[i] = new processThread();

		vector<_RECT> rectList;
		for (int j = i;j < totalBlocks; j += num_threads)
		{
			int iblock = (int)blockList[j];

			int nLines = blockSize;
			if (iblock == totalBlocks - 1)
			{
				// last block
				nLines = min(blockSize, nRasterYSizeRead-iblock*blockSize);
			}
			rectList.push_back(_RECT(0, iblock*blockSize, nRasterXSizeRead-1, iblock*blockSize+nLines-1));
		}
		threads[i]->setData(rectList, processType::ExtremeStatistics);
		threads[i]->start();
	}
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished

	nBinStepList.clear();
	for (int ib = 0;ib < nOutBands;++ib)
	{
		nBinStepList.push_back((gMaxValueList[ib]-gMinValueList[ib]) / (double)nBins);
	}

	// 2: histogram statistics
	// counter refresh
	finishedBlocks = 0;
	gHistogramList.clear();
	for (int ib = 0;ib < nOutBands;++ib)
	{
		vector<long long> histogram(nBins, 0.0);
		gHistogramList.push_back(histogram);
	}
	for(int i = 0;i < num_threads;++i)
	{
		vector<_RECT> rectList;
		for (int j = i;j < totalBlocks; j += num_threads)
		{
			int iblock = (int)blockList[j];

			int nLines = blockSize;
			if (iblock == totalBlocks - 1)
			{
				// last block
				nLines = min(blockSize, nRasterYSizeRead-iblock*blockSize);
			}
			rectList.push_back(_RECT(0, iblock*blockSize, nRasterXSizeRead-1, iblock*blockSize+nLines-1));
		}
		threads[i]->setData(rectList, processType::HistogramStatistcs);
		threads[i]->start();
	}
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished

	// 3: find the stretch percent position
	//vector<double> vLowList;
	//vector<double> vHighList;
	//vector<double> vLenList;
	vLowList.clear();
	vHighList.clear();
	vLenList.clear();
	
	for (int ib = 0;ib < nOutBands;++ib)
	{
		long long counter = 0;
		for (int i = 0;i < nBins;++i)
		{
			counter += gHistogramList[ib][i];
		}
		int stretch_count = int(trim_percentage*0.01 * counter + 0.5);

		int low_pos = 0;
		int low_count = 0;
		for (int i = 0;i < nBins;++i)
		{
			low_count += gHistogramList[ib][i];
			if (low_count > stretch_count)
			{
				low_pos = i;
				break;
			}
		}
		int high_pos = nBins - 1;
		int high_count = 0;
		for (int i = nBins - 1;i > 0;--i)
		{
			high_count += gHistogramList[ib][i];
			if (high_count > stretch_count)
			{
				high_pos = i;
				break;
			}
		}

		// high and low to data
		if (low_pos == high_pos)
		{
			if (low_pos > 0)
			{
				low_pos--;
			}
			else
			{
				high_pos++;
			}
		}
		double vlow = low_pos*nBinStepList[ib] + gMinValueList[ib];
		double vhigh = high_pos*nBinStepList[ib] + gMinValueList[ib];
		double vlen = vhigh - vlow;
		if (vlen < 1e-6)
		{
			vlen += 1e-6;
		}

		vLowList.push_back(vlow);
		vHighList.push_back(vhigh);
		vLenList.push_back(vlen);
	}


	// 4: stretching
	// counter refresh
	finishedBlocks = 0;
	int nOutRasterXSize = int(nRasterXSizeRead * out_scale + 0.5);
	int nOutRasterYSize = int(nRasterYSizeRead * out_scale + 0.5);

	//const char* outTiffFile = pszOutputFile;
	char* pszDestExtension = CPLStrdup(CPLGetExtension(pszOutputFile));
	bJpg = (0 == strcmp(strupr(pszDestExtension), "JPG"));
	if (bJpg)
	{
		// a single thread is used for jpeg writer
		if (initializeJpg(nOutRasterXSize, nOutRasterYSize))
		{
			//int nBands = pSrcDataset->GetRasterCount();
			int nBands = (int)outBandList.size();
			GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
			int dataTypeSize = GDALGetDataTypeSize(eDT)/8;
			int outBlockSize = int(blockSize * out_scale+0.5);

			for (int j = 0;j < totalBlocks; j += 1)
			{
				int iblock = (int)blockList[j];

				int outLines = outBlockSize;
				int inLines = blockSize;
				if (iblock == totalBlocks - 1)
				{
					// last block
					outLines = min(outBlockSize, nOutRasterYSize-iblock*outBlockSize);
					inLines = min(blockSize, nRasterYSizeRead-iblock*blockSize);
				}

				int nBlockOutRasterXSize = nOutRasterXSize;
				int nBlockOutRasterYSize = outLines;

				int nBlockRasterXSize = nRasterXSizeRead;
				int nBlockRasterYSize = inLines;

				int inStartLine = iblock * blockSize;
				int outStartLine = iblock * outBlockSize;
				

				// stretch to byte
				// --------------------------------------------------
				//      Loop overimage, copying image data.                           
				// --------------------------------------------------
				CPLErr      eErr;
				GByte* pBuf = new GByte[nOutRasterXSize * nOutRasterYSize * dataTypeSize];
				GByte* pNewBuf = new GByte[nOutRasterXSize * nOutRasterYSize * nBands];
				for (int ib = 0;ib < nBands;++ib)
				{
					pSrcDataset->GetRasterBand(outBandList[ib])->RasterIO(GF_Read, 0, inStartLine, nBlockRasterXSize, nBlockRasterYSize, 
						pBuf, nBlockOutRasterXSize, nBlockOutRasterYSize, eDT, 0, 0);
					for (int i = 0; i < nBlockOutRasterYSize; i++)
					{
						for (int j = 0;j < nBlockOutRasterXSize;++j)
						{
							int pos = i * nBlockOutRasterXSize + j;
							int newpos = i * nBlockOutRasterXSize * nBands + j * nBands + ib;
							double data = SRCVAL(pBuf, eDT, pos);

							if (bIgnorNull && data == null_value)
							{
								pNewBuf[newpos] = out_null_value;
								continue;
							}
							int oData = (int)(pow((data - vLowList[ib]) / vLenList[ib], gamma) * 255.0 + 0.5);
							oData = oData > 255 ? 255 : oData;
							oData = oData < 0 ? 0 : oData;
							pNewBuf[newpos] = (GByte)oData;
						}
					}
				}

				for( int iLine = 0; iLine< nBlockOutRasterYSize; iLine++)
				{
					JSAMPLE     *ppSamples;
					ppSamples = pNewBuf + iLine * nBlockOutRasterXSize * nBands;

					jpeg_write_scanlines( &sCInfo, &ppSamples, 1 );
				}
				CPLFree( pBuf );
				CPLFree( pNewBuf );

				finishedBlocks++;
				int percent = 50+(int)(finishedBlocks / (double)(totalBlocks) * 50 + 0.5);
				printf("\r%3d%%", percent);
				fflush( stdout );
			}

			jpeg_finish_compress( &sCInfo );
			jpeg_destroy_compress( &sCInfo );

			VSIFClose( fpJpgImage);
		}
		else
		{
			printf("create jpeg file %s failed.\n", pszOutputFile);
		}
	}
	else
	{
		GDALDriverH hDriver = GDALGetDriverByName( "GTiff" );
		pDstDataset = (GDALDataset*)GDALCreate( hDriver, pszOutputFile, nOutRasterXSize, nOutRasterYSize, nOutBands, GDALDataType::GDT_Byte, NULL );
		double adfThisGeoTransform[6];
		GDALGetGeoTransform( hSrcDs, adfThisGeoTransform );
		adfThisGeoTransform[1] /= out_scale;
		adfThisGeoTransform[5] /= out_scale;
		pDstDataset->SetGeoTransform(adfThisGeoTransform);
		pDstDataset->SetProjection(GDALGetProjectionRef( hSrcDs ));

		for(int i = 0;i < num_threads;++i)
		{
			vector<_RECT> rectList;
			for (int j = i;j < totalBlocks; j += num_threads)
			{
				int iblock = (int)blockList[j];

				int nLines = blockSize;
				if (iblock == totalBlocks - 1)
				{
					// last block
					nLines = min(blockSize, nRasterYSizeRead-iblock*blockSize);
				}
				rectList.push_back(_RECT(0, iblock*blockSize, nRasterXSizeRead-1, iblock*blockSize+nLines-1));
			}
			threads[i]->setData(rectList, processType::Stretching);
			threads[i]->start();
		}
		bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
		bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished
		GDALClose(pDstDataset);
	}

	///* -------------------------------------------------------------------- */
	///*      Find the output driver.                                         */
	///* -------------------------------------------------------------------- */
	//GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	//if( hDriver == NULL )
	//{
	//	int	iDr;

	//	printf( "Output driver `%s' not recognised.\n", pszFormat );
	//	printf( "The following format drivers are configured and support output:\n" );
	//	for( iDr = 0; iDr < GDALGetDriverCount(); iDr++ )
	//	{
	//		GDALDriverH hDriver = GDALGetDriver(iDr);

	//		if( GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATE, NULL ) != NULL
	//			|| GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATECOPY,
	//			NULL ) != NULL )
	//		{
	//			printf( "  %s: %s\n",
	//				GDALGetDriverShortName( hDriver  ),
	//				GDALGetDriverLongName( hDriver ) );
	//		}
	//	}
	//	printf( "\n" );
	//	Usage();

	//	GDALClose( hSrcDs );
	//	GDALDestroyDriverManager();
	//	CSLDestroy( papszCreateOptions );
	//	exit( 1 );
	//}

	///* ==================================================================== */
	///*      Create a virtual dataset.                                       */
	///* ==================================================================== */

	///* -------------------------------------------------------------------- */
	///*      Make a virtual clone.                                           */
	///* -------------------------------------------------------------------- */
	//poVDS = (VRTDataset *) VRTCreate( nOutRasterXSize, nOutRasterYSize );


	//double adfThisGeoTransform[6];
	//GDALGetGeoTransform( hSrcDs, adfThisGeoTransform );
	//adfThisGeoTransform[1] /= out_scale;
	//adfThisGeoTransform[5] /= out_scale;
	//poVDS->SetGeoTransform(adfThisGeoTransform);
	//poVDS->SetProjection(GDALGetProjectionRef( hSrcDs ));

	///* -------------------------------------------------------------------- */
	///*      Transfer generally applicable metadata.                         */
	///* -------------------------------------------------------------------- */
	//char** papszMetadata = CSLDuplicate(((GDALDataset*)hSrcDs)->GetMetadata());
	//poVDS->SetMetadata( papszMetadata );
	//CSLDestroy( papszMetadata );
	//AttachMetadata( (GDALDatasetH) poVDS, papszMetadataOptions );

	//const char* pszInterleave = GDALGetMetadataItem(hSrcDs, "INTERLEAVE", "IMAGE_STRUCTURE");
	//if (pszInterleave)
	//	poVDS->SetMetadataItem("INTERLEAVE", pszInterleave, "IMAGE_STRUCTURE");

	//for (int i = 0;i < (int)outBandList.size();++i)
	//{
	//	poVDS->AddBand( GDT_Byte, NULL );
	//}
	
	


	///* -------------------------------------------------------------------- */
	///*      Write to the output file using CopyCreate().                    */
	///* -------------------------------------------------------------------- */
	//GDALDatasetH hOutDS = GDALCreateCopy( hDriver, pszOutputFile, (GDALDatasetH) poVDS,
	//	FALSE, papszCreateOptions, 
	//	NULL, NULL );
	//if( hOutDS != NULL )
	//{
	//	int bHasGotErr = FALSE;
	//	CPLErrorReset();
	//	GDALFlushCache( hOutDS );
	//	if (CPLGetLastErrorType() != CE_None)
	//		bHasGotErr = TRUE;
	//	GDALClose( hOutDS );
	//	if (bHasGotErr)
	//		hOutDS = NULL;
	//}

	//GDALClose( (GDALDatasetH) poVDS );


	// 5: clean
	//if (bJpg)
	//{
	//	//GDALDriverH hJpgDriver = GDALGetDriverByName( "JPEG" );
	//	//GDALCreateCopy(hJpgDriver, pszOutputFile, pDstDataset, FALSE, papszCreateOptions, 0, 0);
	//	GDALClose( JPEGCreateCopy( pszOutputFile, pDstDataset,
	//		FALSE, papszCreateOptions, 0));
	//}
	
	GDALClose( pSrcDataset );

	printf("\r%3d%%\n", 100);
}


void img_translate()
{
	GDALDatasetH hSrcDs = GDALOpenShared(pszInputFile, GA_ReadOnly );
	pSrcDataset = (GDALDataset  *)hSrcDs;
	int nRasterXSizeRead = pSrcDataset->GetRasterXSize();
	int nRasterYSizeRead = pSrcDataset->GetRasterYSize();
	int nBandTotal = pSrcDataset->GetRasterCount();
	GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
	GDALDataType eOutDT = eDT;
	char* pszDestExtension = CPLStrdup(CPLGetExtension(pszOutputFile));
	bJpg = (0 == strcmp(strupr(pszDestExtension), "JPG"));

	if (eDT != GDALDataType::GDT_Byte && bJpg)
	{
		printf("Only eight bit byte bands supported for jpeg.\n"
			"the input bands will be stretched to 0~255\n");
		trim_percentage = 0.0;
		gamma = 1.0;
		img_stretching();
		return;
	}

	// get output size
	if( pszOXSize != NULL)
	{
		outSizeX = (int) ((pszOXSize[strlen(pszOXSize)-1]=='%' 
			? CPLAtofM(pszOXSize)/100*nRasterXSizeRead : atoi(pszOXSize)));
	}
	if(pszOYSize != NULL)
	{
		outSizeY = (int) ((pszOYSize[strlen(pszOYSize)-1]=='%' 
			? CPLAtofM(pszOYSize)/100*nRasterYSizeRead : atoi(pszOYSize)));
	}
	if (0 < outSizeX)
	{
		out_scale = outSizeX / (double)nRasterXSizeRead;
	}
	else if (0 < outSizeY)
	{
		out_scale = outSizeY / (double)nRasterYSizeRead;
	}

	// check band list
	nOutBands = (int)outBandList.size();
	for (int i = 0;i < nOutBands;i++)
	{
		if (outBandList[i] < 1 || outBandList[i] > nBandTotal)
		{
			printf("warning: band %d is out of range (1~%d).\n", outBandList[i], nBandTotal);
			outBandList.erase(outBandList.begin()+i);
			i--;
		}
	}
	nOutBands = (int)outBandList.size();
	if (nOutBands < 1)
	{
		//printf("warning: output band number cannot less than 1.\n"
		//	"Therefore, all the bands of the input file will be outputted.\n");

		for (int i = 0;i < nBandTotal;++i)
		{
			outBandList.push_back(i+1);
		}
		nOutBands = (int)outBandList.size();
	}

	int num_threads = nThreads;
	if (num_threads < 1)
	{
		num_threads = OpenThreads::GetNumberOfProcessors() * 2;
	}
	if (num_threads > maxThreads)
	{
		num_threads = maxThreads;
	}
	//num_threads = 3;
	cout<<"using "<<num_threads<<" threads..."<<endl;
	GLOBAL_NUM_THREADS = num_threads + 1;
	std::vector<processThread *> threads(num_threads);
	OpenThreads::Thread::SetConcurrency(num_threads);
	OpenThreads::Thread::Init();

	//设置光标位置
	printf("\r%3d%%", 0);

	// counter start
	finishedBlocks = 0;

	totalBlocks = ceil(nRasterYSizeRead / (double)blockSize );
	vector<int> blockList;
	for (int i = 0;i < totalBlocks;++i)
	{
		blockList.push_back(i);
	}

	// translate
	// counter refresh
	finishedBlocks = 0;
	int nOutRasterXSize = int(nRasterXSizeRead * out_scale + 0.5);
	int nOutRasterYSize = int(nRasterYSizeRead * out_scale + 0.5);

	if (bJpg)
	{
		// a single thread is used for jpeg writer
		if (initializeJpg(nOutRasterXSize, nOutRasterYSize))
		{
			//int nBands = pSrcDataset->GetRasterCount();
			int nBands = (int)outBandList.size();
			GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
			int dataTypeSize = GDALGetDataTypeSize(eDT)/8;
			int outBlockSize = int(blockSize * out_scale+0.5);

			for (int j = 0;j < totalBlocks; j += 1)
			{
				int iblock = (int)blockList[j];

				int outLines = outBlockSize;
				int inLines = blockSize;
				if (iblock == totalBlocks - 1)
				{
					// last block
					outLines = min(outBlockSize, nOutRasterYSize-iblock*outBlockSize);
					inLines = min(blockSize, nRasterYSizeRead-iblock*blockSize);
				}

				int nBlockOutRasterXSize = nOutRasterXSize;
				int nBlockOutRasterYSize = outLines;

				int nBlockRasterXSize = nRasterXSizeRead;
				int nBlockRasterYSize = inLines;

				int inStartLine = iblock * blockSize;
				int outStartLine = iblock * outBlockSize;


				// stretch to byte
				// --------------------------------------------------
				//      Loop overimage, copying image data.                           
				// --------------------------------------------------
				CPLErr      eErr;
				GByte* pBuf = new GByte[nOutRasterXSize * nOutRasterYSize * dataTypeSize * nBands];

				for (int ib = 0;ib < nBands;++ib)
				{
					pSrcDataset->GetRasterBand(outBandList[ib])->RasterIO(GF_Read, 0, inStartLine, nBlockRasterXSize, nBlockRasterYSize, 
						pBuf+ib*dataTypeSize, nBlockOutRasterXSize, nBlockOutRasterYSize, eDT, nBands*dataTypeSize, nBlockOutRasterXSize*nBands*dataTypeSize);
				}

				for( int iLine = 0; iLine< nBlockOutRasterYSize; iLine++)
				{
					JSAMPLE     *ppSamples;
					ppSamples = pBuf + iLine * nBlockOutRasterXSize * nBands * dataTypeSize;

					jpeg_write_scanlines( &sCInfo, &ppSamples, 1 );
				}
				CPLFree( pBuf );

				finishedBlocks++;
				int percent = 50+(int)(finishedBlocks / (double)(totalBlocks) * 50 + 0.5);
				printf("\r%3d%%", percent);
				fflush( stdout );
			}

			jpeg_finish_compress( &sCInfo );
			jpeg_destroy_compress( &sCInfo );

			VSIFClose( fpJpgImage);
		}
		else
		{
			printf("create jpeg file %s failed.\n", pszOutputFile);
		}
	}
	else
	{
		GDALDriverH hDriver = GDALGetDriverByName( "GTiff" );
		pDstDataset = (GDALDataset*)GDALCreate( hDriver, pszOutputFile, nOutRasterXSize, nOutRasterYSize, nOutBands, eOutDT, NULL );
		double adfThisGeoTransform[6];
		GDALGetGeoTransform( hSrcDs, adfThisGeoTransform );
		adfThisGeoTransform[1] /= out_scale;
		adfThisGeoTransform[5] /= out_scale;
		pDstDataset->SetGeoTransform(adfThisGeoTransform);
		pDstDataset->SetProjection(GDALGetProjectionRef( hSrcDs ));

		for(int i = 0;i < num_threads;++i)
		{

			threads[i] = new processThread();
			vector<_RECT> rectList;
			for (int j = i;j < totalBlocks; j += num_threads)
			{
				int iblock = (int)blockList[j];

				int nLines = blockSize;
				if (iblock == totalBlocks - 1)
				{
					// last block
					nLines = min(blockSize, nRasterYSizeRead-iblock*blockSize);
				}
				rectList.push_back(_RECT(0, iblock*blockSize, nRasterXSizeRead-1, iblock*blockSize+nLines-1));
			}
			threads[i]->setData(rectList, processType::Translate);
			threads[i]->start();
		}
		bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
		bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished
		GDALClose(pDstDataset);
	}


	///* -------------------------------------------------------------------- */
	///*      Write to the output file using CopyCreate().                    */
	///* -------------------------------------------------------------------- */
	//GDALDatasetH hOutDS = GDALCreateCopy( hDriver, pszOutputFile, (GDALDatasetH) poVDS,
	//	FALSE, papszCreateOptions, 
	//	NULL, NULL );
	//if( hOutDS != NULL )
	//{
	//	int bHasGotErr = FALSE;
	//	CPLErrorReset();
	//	GDALFlushCache( hOutDS );
	//	if (CPLGetLastErrorType() != CE_None)
	//		bHasGotErr = TRUE;
	//	GDALClose( hOutDS );
	//	if (bHasGotErr)
	//		hOutDS = NULL;
	//}

	//GDALClose( (GDALDatasetH) poVDS );


	// 5: clean
	GDALClose( pSrcDataset );

	printf("\r%3d%%\n", 100);
}

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
				pszInputFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-o") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutputFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-trim") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				trim_percentage = atof(argv[++i]) ;
			}
			else if( 0 == _stricmp(argv[i],"-g") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				gamma = atof(argv[++i]) ;
			}
			else if( 0 == _stricmp(argv[i],"-stretching") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				if(atoi( argv[++i] ) > 0)
				{
					bStretching = true;
				}
				else
				{
					bStretching = false;
				}
			}
			else if(0 == _stricmp(argv[i],"-nv") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				null_value = CPLAtofM( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-in") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				if(atoi( argv[++i] ) > 0)
				{
					bIgnorNull = true;
				}
				else
				{
					bIgnorNull = false;
				}
			}
			else if(0 == _stricmp(argv[i],"-onv") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				out_null_value = atoi( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-nt") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				nThreads = atoi( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-bs") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				blockSize = atoi( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-nb") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				nBins = atoi( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-ss") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				statistics_scale = atof( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-osx") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOXSize = argv[++i];
				//outSizeX = atoi( argv[++i] );	
			}
			else if(0 == _stricmp(argv[i],"-osy") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				//outSizeY = atoi( argv[++i] );
				pszOYSize = argv[++i];
			}
			else if(0 == _stricmp(argv[i],"-of") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszFormat = argv[++i];
				bFormatExplicitelySet = true;
			}
			else if( EQUAL(argv[i],"-co") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				papszCreateOptions = CSLAddString( papszCreateOptions, argv[++i] );
			} 
			else if( EQUAL(argv[i],"-mo") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				papszMetadataOptions = CSLAddString( papszMetadataOptions, argv[++i] );
			} 
			else if( 0 == _stricmp(argv[i],"-b") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				nOutBands = atoi(argv[++i]);
				if (nOutBands < 1)
				{
					printf("output band number cannot less than 1\n");
					exit(0);
				}
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nOutBands);
				for (int ib = 0;ib < nOutBands;++ib)
				{
					outBandList.push_back(atoi(argv[++i]));
				}
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
			printf("input file can not be empty!\n");
			Usage();
		}
		else if (0 == strcmp(pszOutputFile, ""))
		{
			printf("output file can not be empty!\n");
			Usage();
		}
		else
		{
			if ( strcmp(pszInputFile, pszOutputFile) == 0)
			{
				Usage("Source and destination datasets must be different.");
			}

			if (_access( pszOutputFile, 0 ) != -1 && !bOverwrite)
			{
				printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
					" \"-overwrite\" option.\n", pszOutputFile);
				Usage(0);
			}

			//if (!bFormatExplicitelySet)
			//	CheckExtensionConsistency(pszOutputFile, pszFormat);

			clock_t  clockBegin, clockEnd;
			clockBegin = clock();
			if (bStretching)
			{
				img_stretching();
			}
			else
			{
				img_translate();
			}
			clockEnd = clock();
			printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
	}
	else
	{
		//pszInputFile = "d:\\workspace\\fuhai\\FCGC600205282\\IMG_PHR1A_MS_002\\rect.tif";
		//pszOutputFile = "d:\\workspace\\fuhai\\FCGC600205282\\IMG_PHR1A_MS_002\\2.jpg";
		//bIgnorNull = true;
		//null_value = 0.0;
		//pszOXSize = "2048";
		//trim_percentage = 0.0;
		//outBandList.push_back(3);
		//outBandList.push_back(4);
		//outBandList.push_back(2);
		//clock_t  clockBegin, clockEnd;
		//clockBegin = clock();
		//img_stretching();
		//clockEnd = clock();
		//printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
		Usage(0);
	}
	return 0;
}