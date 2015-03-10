#include <stdlib.h>
#include <math.h>
#include <direct.h>
#include  <io.h>
#include  <stdio.h>

/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"

#include <fstream>

#include <strUtil.h>
#include <fileUtil.h>
#include <time.h>

using namespace std;
using namespace mylib;

#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "OpenThreads.lib")

static GDALDataset  *pDstDataset = NULL;
static int totalBlocks;
static int finishedBlocks;
const char* pszInputFile = "";
const char* pszOutputFile = "";
int nBands;
vector<int> bandList;
bool bOverwrite = false;
int nThreads = 0;
int maxThreads = 10;
int blockSize = 1024;

#include <OpenThreads/Thread>
#include <OpenThreads/Mutex>
#include <OpenThreads/Barrier>


static OpenThreads::Barrier bar;
static int GLOBAL_NUM_THREADS;

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

class processThread: public OpenThreads::Thread
{
public:
	processThread()
		: OpenThreads::Thread(){};
	virtual  ~processThread()
	{
		//m_pRegistration = NULL;
	};
	void setData(vector<_RECT> rectList)
	{
		m_rectList = rectList;
	}
	virtual void run()
	{
		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.
		for (int i = 0;i < (int)m_rectList.size();++i)
		{
			rebands(m_rectList[i]);
			// reduce
			_quitmutex.lock();
			finishedBlocks++;
			_quitmutex.unlock();

			int percent = (int)(finishedBlocks / (double)(totalBlocks) * 100 + 0.5);
			printf("\r%3d%%", percent);
			fflush( stdout );
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

	void rebands(_RECT rect)
	{
		GDALDataset  * pSrcDataset = (GDALDataset  *)GDALOpen(pszInputFile, GA_ReadOnly );
		int nRasterXSizeRead = pSrcDataset->GetRasterXSize();
		int nRasterYSizeRead = pSrcDataset->GetRasterYSize();
		int nBandTotal = pSrcDataset->GetRasterCount();
		GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
		int dataTypeSize = GDALGetDataTypeSize(eDT)/8;

		int nRasterXSize = rect.endSample - rect.startSample + 1;
		int nRasterYSize = rect.endLine - rect.startLine + 1;

		int outBands = (int)bandList.size();
		for (int i = 0;i < outBands;i++)
		{
			if (bandList[i] < 1 || bandList[i] > nBandTotal)
			{
				printf("warning: band %d is out of range (1~%d).\n", bandList[i], nBandTotal);
				bandList.erase(bandList.begin()+i);
				i--;
			}
		}
		outBands = (int)bandList.size();
		if (outBands < 1)
		{
			printf("output band number cannot less than 1\n");
			return;
		}


		GByte *oldBuf = new GByte[nRasterXSizeRead*nRasterYSizeRead*dataTypeSize];
		for (int i = 0;i < outBands;++i)
		{
			pSrcDataset->GetRasterBand(bandList[i])->RasterIO(GF_Read, 0, rect.startLine, nRasterXSize, nRasterYSize,
				oldBuf, nRasterXSize, nRasterYSize, eDT, NULL, NULL);
			_quitmutex.lock();
			pDstDataset->GetRasterBand(i+1)->RasterIO(GF_Write, 0, rect.startLine, nRasterXSize, nRasterYSize,
				oldBuf, nRasterXSize, nRasterYSize, eDT, NULL, NULL);
			_quitmutex.unlock();
		}

		delete []oldBuf;
		GDALClose(pSrcDataset);
	}
	
	void quit() {
		_quitmutex.lock();
		_quitflag = true;
		_quitmutex.unlock();
	};
private:
	vector<_RECT> m_rectList;
	volatile bool _quitflag;
	OpenThreads::Mutex _quitmutex;
};

void img_rebands()
{
	double adfThisGeoTransform[6];
	GDALDataset *pOldDataset = (GDALDataset *) GDALOpen( pszInputFile, GA_ReadOnly );
	int iWidth = pOldDataset->GetRasterXSize();
	int iHeight = pOldDataset->GetRasterYSize();
	pOldDataset->GetGeoTransform(adfThisGeoTransform);
	const char* projectString = pOldDataset->GetProjectionRef();
	int nBandTotal = pOldDataset->GetRasterCount();
	GDALDataType eDT = pOldDataset->GetRasterBand(1)->GetRasterDataType();
	int dataTypeSize = GDALGetDataTypeSize(eDT) / 8;

	int outBands = (int)bandList.size();
	for (int i = 0;i < outBands;i++)
	{
		if (bandList[i] < 1 || bandList[i] > nBandTotal)
		{
			printf("warning: band %d is out of range (1~%d).\n", bandList[i], nBandTotal);
			bandList.erase(bandList.begin()+i);
			i--;
		}
	}
	outBands = (int)bandList.size();
	if (outBands < 1)
	{
		printf("output band number cannot less than 1\n");
		return;
	}


	const char *pszFormat = "GTiff";
	char *pszTargetSRS = NULL;
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	int nBands = 1;
	pDstDataset = (GDALDataset*)GDALCreate( hDriver, pszOutputFile, iWidth, iHeight, outBands, eDT, NULL );
	pDstDataset->SetProjection(projectString);
	pDstDataset->SetGeoTransform(adfThisGeoTransform);


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

	// 0%
	printf("\r%3d%%", 0);

	// counter start
	finishedBlocks = 0;

	totalBlocks = ceil(iHeight / (double)blockSize );
	vector<int> blockList;
	for (int i = 0;i < totalBlocks;++i)
	{
		blockList.push_back(i);
	}

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
				nLines = min(blockSize, iHeight-iblock*blockSize);
			}
			rectList.push_back(_RECT(0, iblock*blockSize, iWidth-1, iblock*blockSize+nLines-1));
		}
		threads[i]->setData(rectList);
		threads[i]->start();
	}
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished

	GDALClose(pDstDataset);
	// 0%
	printf("\r%3d%%\n", 100);
}

/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: img-bands -i inputFile -o outputFile [-b nb b1 b2 ...] \n"
		" [-overwrite] [-nt nThreads] [-bs blocksize] \n");

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
			else if( 0 == _stricmp(argv[i],"-b") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				nBands = atoi(argv[++i]);
				if (nBands < 1)
				{
					printf("output band number cannot less than 1\n");
					exit(0);
				}
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nBands);
				bandList.clear();
				for (int ib = 0;ib < nBands;++ib)
				{
					bandList.push_back(atoi(argv[++i]));
				}
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
		else if(0 == strcmp(pszOutputFile, ""))
		{
			printf("outputFile can not be empty!\n");
			Usage();
		}
		else
		{

			if (_access( pszOutputFile, 0 ) != -1 && !bOverwrite)
			{
				printf("process skipped: %s is already exist, to overwrite the existing result file, please use"
					" \"-overwrite\" option.\n", pszOutputFile);
				Usage(0);
			}
			clock_t  clockBegin, clockEnd;
			clockBegin = clock();

			img_rebands();
			clockEnd = clock();
			printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
		}
	}
	else
	{
		Usage();
	}
	return 0;
}