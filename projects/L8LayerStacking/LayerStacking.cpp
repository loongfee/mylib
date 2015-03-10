#include <iostream>
#include <string>
#include <vector>
using namespace std;

//////gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"

#include "strUtil.h"
#include "FileAndDirFinder.h"

bool parseInputFiles(char* papszSrcFiles, std::vector<string>& inFileList, bool bFilter);

COORD GetConsoleCursorPosition(HANDLE hHandle)
{  
	CONSOLE_SCREEN_BUFFER_INFO info={0};  
	GetConsoleScreenBufferInfo( hHandle , &info );  
	return info.dwCursorPosition;  
} 

/************************************************************************/
/*                               Usage()                                */
/************************************************************************/
static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: LayerStacking infiles dstfile\n");

	if( pszErrorMsg != NULL )
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);

	exit( 0 );
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

int main(int argc, char* argv[])
{
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持


	GDALDatasetH hDstDS;
	//char *pszDstFilename = NULL;
	char **papszSrcFiles = NULL;
	bool bOverwrite = FALSE;
	bool bFilter = FALSE;
	double null_value = 0;
	char *pszDstFile = NULL;
	int i;
	/* -------------------------------------------------------------------- */
	/*      Parse arguments.                                                */
	/* -------------------------------------------------------------------- */
	for( i = 1; i < argc; i++ )
	{
		if ((0 == _stricmp(argv[i],"-h")) || (0 == _stricmp(argv[i],"--help"))) {
			Usage(argv[0]);
			return 0;
		}
		else if( 0 == _stricmp(argv[i],"-overwrite") )
			bOverwrite = TRUE;
		else if( (0 == _stricmp(argv[i],"-f")) || (0 == _stricmp(argv[i],"--filter")) )
			bFilter = TRUE;
		else if( argv[i][0] == '-' )
			Usage(CPLSPrintf("Unkown option name '%s'", argv[i]));

		else 
			papszSrcFiles = CSLAddString( papszSrcFiles, argv[i] );
			//printf(argv[i]);
	}
	/* -------------------------------------------------------------------- */
	/*      The last filename in the file list is really our destination    */
	/*      file.                                                           */
	/* -------------------------------------------------------------------- */
	if( CSLCount(papszSrcFiles) > 1 )
	{
		pszDstFile = papszSrcFiles[CSLCount(papszSrcFiles)-1];
		papszSrcFiles[CSLCount(papszSrcFiles)-1] = NULL;
	}

	if( pszDstFile == NULL )
	{
		Usage("No target filename specified.");
		return 0;
	}

	std::vector<string> inFileList;
	parseInputFiles(papszSrcFiles[0], inFileList, bFilter);


	/* FIXME ? source filename=target filename and -overwrite is definitely */
	/* an error. But I can't imagine of a valid case (without -overwrite), */
	/* where it would make sense. In doubt, let's keep that dubious possibility... */
	for (int i=0;i<(int)inFileList.size();++i)
	{
		if ( CSLCount(papszSrcFiles) == 1 &&
			strcmp(inFileList[i].c_str(), pszDstFile) == 0 && bOverwrite)
		{
			fprintf(stderr, "Source and destination datasets must be different.\n");
			return 0;
		}

	}

	CPLPushErrorHandler( CPLQuietErrorHandler );
	hDstDS = GDALOpen( pszDstFile, GA_Update );
	//hDstDS = GDALOpen( pszDstFile, GA_ReadOnly );
	CPLPopErrorHandler();

	if( hDstDS != NULL && bOverwrite )
	{
		GDALClose(hDstDS);
		hDstDS = NULL;
	}

	// 已存在，不覆盖
	if( hDstDS != NULL && !bOverwrite )
	{
		return 0;
	}

	/* Avoid overwriting an existing destination file that cannot be opened in */
	/* update mode with a new GTiff file */
	if ( hDstDS == NULL && !bOverwrite )
	{
		CPLPushErrorHandler( CPLQuietErrorHandler );
		hDstDS = GDALOpen( pszDstFile, GA_ReadOnly );
		CPLPopErrorHandler();

		if (hDstDS)
		{
			fprintf( stderr, 
				"Output dataset %s exists, but cannot be opened in update mode\n",
				pszDstFile );
			GDALClose(hDstDS);
			return 0;
		}
	}
	//GDALClose(hDstDS);

	GDALDataset  *hSrcDataset, *outtiff;
	GDALDriverH		hDriver;
	GDALDatasetH	hSrcDs;
	CPLErr      eErr = CE_None,eErr1=CE_None;
	int nRasterXSizeRead,nRasterYSizeRead;
	GDALDataType eDT;
	GDALRasterBand *rasterband,*rasterband_IMG;
	GByte   *pabyLine,*bufout;

	int nBands = inFileList.size();
	if(nBands < 1) return false;

	hSrcDs=GDALOpen(inFileList[0].c_str(), GA_ReadOnly );
	hSrcDataset = (GDALDataset  *)hSrcDs;

	rasterband_IMG=hSrcDataset->GetRasterBand(1);
	eDT=rasterband_IMG->GetRasterDataType();
	nRasterXSizeRead=hSrcDataset->GetRasterXSize();
	nRasterYSizeRead=hSrcDataset->GetRasterYSize();

	const char         *pszFormat = "GTiff";
	char               *pszTargetSRS = NULL;
	double adfThisGeoTransform[6];

	hDriver = GDALGetDriverByName( pszFormat );
	hDstDS = GDALCreate( hDriver, pszDstFile, nRasterXSizeRead, nRasterYSizeRead, 
		nBands, eDT, NULL );

	if( hDstDS == NULL ) {	
		GDALClose( hSrcDataset );
		cout<<pszDstFile<<"创建失败！"; 
		return false;
	}


	GDALGetGeoTransform( hSrcDs, adfThisGeoTransform );
	pszTargetSRS = CPLStrdup(GDALGetProjectionRef(hSrcDs));

	GDALSetProjection( hDstDS, pszTargetSRS );
	GDALSetGeoTransform( hDstDS, adfThisGeoTransform );
	outtiff=(GDALDataset  *)hDstDS;
	GDALClose( hSrcDataset );

	int percent = 0;
	//设置光标位置
	HANDLE hOut; 
	hOut = GetStdHandle(STD_OUTPUT_HANDLE); 
	//得到当前光标位置
	COORD pos= GetConsoleCursorPosition(hOut);
	SetConsoleCursorPosition(hOut, pos);
	cout<<percent<<"%";

	////////////////////////
	// 一次处理多行数据，提高速度
	int linenumber = 2000;
	int nyNum=(nRasterYSizeRead-1)/linenumber+1;//计算列方向块数
	for (int nYI=0;nYI<nyNum;nYI++)//块循环
	{
		//cout<<"波段合并： "<<nYI+1<<" of "<<nyNum<<"..."<<endl;
		int Bufsizex = nRasterXSizeRead;
		int Bufsizey = linenumber;


		//列末尾小块处理
		if (nYI == nyNum-1)
		{
			Bufsizey = nRasterYSizeRead - (nyNum-1) * linenumber;//得出当前块的宽度Bufsizex，高度Bufsizey
			Bufsizey = min(Bufsizey, linenumber);
		}

		pabyLine=(GByte *) CPLMalloc (Bufsizex * Bufsizey * GDALGetDataTypeSize( eDT ) / 8);
		bufout=(GByte *) CPLMalloc (Bufsizex * Bufsizey * GDALGetDataTypeSize( eDT ) / 8);
		for (i=0;i<nBands;i++) {
			hSrcDs=GDALOpen(inFileList[i].c_str(), GA_ReadOnly );
			hSrcDataset = (GDALDataset  *)hSrcDs;
			rasterband_IMG=hSrcDataset->GetRasterBand(1);
			rasterband=outtiff->GetRasterBand(i+1);
			eErr = rasterband_IMG->RasterIO( GF_Read, 0, nYI*linenumber, Bufsizex , Bufsizey, pabyLine, Bufsizex, Bufsizey, eDT,  0, 0 );
			eErr1 = rasterband->RasterIO( GF_Write, 0, nYI*linenumber, Bufsizex , Bufsizey,pabyLine, Bufsizex, Bufsizey, eDT,  0, 0 ); 
			GDALClose( hSrcDataset );

			// 计算当前进度
			int tmpPercent = (int)((nYI+i/(double)(nBands)) / (double)(nyNum) * 100 + 0.5);
			if(tmpPercent > percent)
			{
				percent = tmpPercent;
				// 在保留光标位置输出，覆盖原输出内容
				SetConsoleCursorPosition(hOut, pos);
				cout<<percent<<"%";
			}
		}
		CPLFree( pabyLine );
		CPLFree( bufout );
	}

	GDALClose( outtiff );

	SetConsoleCursorPosition(hOut, pos);
	cout<<100<<"%"<<endl;

	return 1;
}

bool parseInputFiles(char* papszSrcFiles, std::vector<string>& inFileList, bool bFilter)
{
	//string srcFiles( papszSrcFiles[0], strnlen(papszSrcFiles[0], 10240) );
	std::string srcFiles((const char *)(papszSrcFiles));
	//srcFiles.assign( papszSrcFiles[0], 10240 );
	//srcFiles.erase( std::find( srcFiles.begin(), srcFiles.end(), '\0' ), srcFiles.end() );
	if (!bFilter)
	{
		// 直接是以";"分割的文件列表
		SplitString(srcFiles, ";", inFileList, false);
		for (int i = 0; i < inFileList.size(); i++)
		{
			printf("%s\n", inFileList[i]);
		}
	}
	else{
		// 按通配符查找
		string strFold = SBeforeLast(srcFiles, '\\');
		if (strFold.empty())
		{
			//当前目录
			string filter = srcFiles;
			CFileAndDirFinder cFileAndDirFinder;
			cFileAndDirFinder.FindAllFileHere(filter.c_str(), inFileList);
		}
		else{
			string filter = SAfterLast(srcFiles, '\\');
			CFileAndDirFinder cFileAndDirFinder;
			cFileAndDirFinder.FindAllFile(strFold.c_str(), filter.c_str(), inFileList);
		}
	}
	return true;
}