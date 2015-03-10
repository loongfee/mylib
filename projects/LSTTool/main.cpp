#pragma once

//////gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"
#include "gdalwarper.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "ogr_api.h"
#include "commonutils.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "gdalwarp.h"
#include "raster_clip.h"
#include "util.h"

#include <fileUtil.h>
#include <strUtil.h>
#include "FileAndDirFinder.h"

using namespace std;
using namespace mylib;

#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")

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

/************************************************************************/
/*                               GDALExit()                             */
/*  This function exits and cleans up GDAL and OGR resources            */
/*  Perhaps it should be added to C api and used in all apps?           */
/************************************************************************/

static int GDALExit( int nCode )
{
  const char  *pszDebug = CPLGetConfigOption("CPL_DEBUG",NULL);
  if( pszDebug && (_stricmp(pszDebug,"ON") || _stricmp(pszDebug,"") ) )
  {  
    GDALDumpOpenDatasets( stderr );
    CPLDumpSharedList( NULL );
  }

  GDALDestroyDriverManager();

#ifdef OGR_ENABLED
  OGRCleanupAll();
#endif

  exit( nCode );
}

/************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
    printf( 
        "Usage: image_average  [-f filter] [-nv nullvalue]\n"
        "    [-cutline datasource] [-cl layer] [-cwhere expression]\n"
		" dstfile\n"
		"\n");

    if( pszErrorMsg != NULL )
        fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);

    GDALExit( 1 );
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
    do { if (i + nExtraArg >= argc) \
        Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

/************************************************************************/
/*                               main()                                */
/************************************************************************/
int main( int argc, char ** argv )
{	
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持
	
	char **papszSrcFiles = NULL;
	const char *pszFilter = "*.tif";
	const char *pszInputPath = "";
	double null_value = 0;
	char *pszDstFile = NULL;
    char *pszCutlineDSName = NULL;
	char *pszCLayer = NULL, *pszCWHERE = NULL;
	bool bFilter = FALSE;
	int i;

	const char* tmpName = "average_tmp.tif";
/* -------------------------------------------------------------------- */
/*      Parse arguments.                                                */
/* -------------------------------------------------------------------- */
    for( i = 1; i < argc; i++ )
    {
        if( 0 == _stricmp(argv[i],"-i") )
        {
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			papszSrcFiles = CSLAddString( papszSrcFiles, argv[++i] );
            //pszFilter = argv[++i] ;
		}
		if( 0 == _stricmp(argv[i],"-f") )
		{
			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
			bFilter = true;
			//pszFilter = argv[++i] ;
		}
		else if( 0 == _stricmp(argv[i],"-nv") )
		{
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
            null_value = CPLAtofM( argv[++i] );	
		}
		else if( EQUAL(argv[i],"-cutline") )
        {
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
            pszCutlineDSName = argv[++i];
        }
        else if( EQUAL(argv[i],"-cwhere") )
        {
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
            pszCWHERE = argv[++i];
        }
        else if( EQUAL(argv[i],"-cl") )
        {
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
            pszCLayer = argv[++i];
        }
        else if( argv[i][0] == '-' )
            Usage(CPLSPrintf("Unkown option name '%s'", argv[i]));
        else 
            pszDstFile = argv[i];
    }
	if (!pszDstFile)
	{
		Usage();
		return 0;
	}

	std::vector<string> files;
	parseInputFiles(papszSrcFiles[0], files, bFilter);

	//TCHAR *tchar_Filter = char2TCHAR(pszFilter);
	//FilesVec files;
	//FindFiles( GetAbsolutePathName(pszFilter),files );

	//FindFiles( _T("E:\\Dropbox\\Programs\\image_average\\Release\\*.tif") ,files );
	//cout<<"Filter:"<<pszFilter<<endl;

	vector<char*> fileList;
	int nWidth, nHeight;
	GDALDataType eDT;
	for (int i = 0; i < (int)files.size(); i++)
	{
		char * filename = TCHAR2char(GetAbsolutePathName(files[i].c_str()));
		if(0 == strcmp(filename, TCHAR2char(GetAbsolutePathName(tmpName))))
		{
			//printf("%s结果图像不作为输入图像.\n", filename);
			continue;
		}
		//char *buf = NULL;
		//Tchar2Char (GetAbsolutePathName(files[i]), filename);
		//filename = new char[1024];
		//sprintf_s(filename, 1024, "E:\\Dropbox\\Programs\\image_average\\Release\\%s", buf);
		GDALDatasetH hDataset = GDALOpen(filename, GA_ReadOnly );
		int w = ((GDALDataset  *)hDataset)->GetRasterXSize();
		int h = ((GDALDataset  *)hDataset)->GetRasterYSize();
		GDALDataType _eDT = ((GDALDataset  *)hDataset)->GetRasterBand(1)->GetRasterDataType();
		if (0 == i)
		{
			nWidth = w;
			nHeight = h;
			eDT = _eDT;
		}

		if (_eDT != eDT)
		{
			printf("%s图像数据位数不一致.\n", filename);
		}
		else if(w != nWidth || h != nHeight)
		{
			printf("%s图像大小不一致.\n", filename);
		}else{
			fileList.push_back(filename);
		}
		GDALClose( hDataset );
	}
	//cout<<fileList.size()<<" images found."<<endl;
	
	// 创建一个TIFF
	const char *pszFormat = "GTiff";
	char *pszTargetSRS = NULL;
	double adfThisGeoTransform[6];
	int nBands = 1;	// 只处理单波段

	GDALDataType outeDT = GDT_Float32;
	char* tmpPath = TCHAR2char(GetAbsolutePathName(tmpName));
	remove(tmpPath);

	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALDatasetH hDstDS = GDALCreate( hDriver, tmpPath, nWidth, nHeight, 
                        nBands, outeDT, NULL );
    
	if( hDstDS == NULL ){
		cout<<tmpPath<<"创建失败！"; 
		return 0;
	}
       
	
	GDALDatasetH hDataset = GDALOpen(fileList[0], GA_ReadOnly );
	GDALGetGeoTransform( hDataset, adfThisGeoTransform );
	pszTargetSRS = CPLStrdup(GDALGetProjectionRef(hDataset));
	GDALClose( hDataset );

    GDALSetProjection( hDstDS, pszTargetSRS );
    GDALSetGeoTransform( hDstDS, adfThisGeoTransform );
    GDALDataset  *pDstDataset = (GDALDataset  *)hDstDS;

	////////////////////////
	// 一次处理多行数据，提高速度
	int linenumber = 2000;
	int nyNum=(nHeight-1)/linenumber+1;//计算列方向块数
	for (int nYI=0;nYI<nyNum;nYI++)//块循环
	{
		cout<<"计算均值： "<<nYI+1<<" of "<<nyNum<<"..."<<endl;
		int Bufsizex = nWidth;
		int Bufsizey = linenumber;


		//列末尾小块处理
		if (nYI == nyNum-1)
		{
			Bufsizey = nHeight - (nyNum-1) * linenumber;//得出当前块的宽度Bufsizex，高度Bufsizey
			Bufsizey = min(Bufsizey, linenumber);
		}

		int *counter = new int[Bufsizex * Bufsizey];
		memset(counter, 0, Bufsizex * Bufsizey * sizeof(int));
		
		float *outBuf = (float *) CPLMalloc (Bufsizex * Bufsizey * GDALGetDataTypeSize( outeDT ) / 8);
		//memset(outBuf, 0, Bufsizex * Bufsizey * GDALGetDataTypeSize( eDT ) / 8);

		for (int iFile = 0; iFile < (int)fileList.size(); iFile++)
		{
	
			GDALDatasetH hDataset = GDALOpen(fileList[iFile], GA_ReadOnly );
			if( hDataset == NULL ){
				cout<<hDataset<<"创建失败！"; 
				continue;
			}
			GByte *inBuf = (GByte *) CPLMalloc (Bufsizex * Bufsizey * GDALGetDataTypeSize( eDT ) / 8);
			CPLErr eErr = ((GDALDataset  *)hDataset)->GetRasterBand(1)->RasterIO( GF_Read, 0, nYI*linenumber, Bufsizex , Bufsizey, inBuf, Bufsizex, Bufsizey, eDT,  0, 0 );
			for (int iHeight = 0; iHeight < Bufsizey; iHeight++)
			{
				for (int iWidth = 0; iWidth < Bufsizex; iWidth++)
				{					
					switch(eDT) {
					case GDT_Byte: 
						if (null_value == *(((GByte*)inBuf) + iHeight*Bufsizex+iWidth))
						{
							break;
						}
						if (0 == iFile)
						{
							*(outBuf + iHeight*Bufsizex+iWidth) = (float)*(((GByte*)inBuf) + iHeight*Bufsizex+iWidth);
						}else
						{
							*(outBuf + iHeight*Bufsizex+iWidth) += (float)*(((GByte*)inBuf) + iHeight*Bufsizex+iWidth);
						}
						counter[iHeight*Bufsizex+iWidth]++;
						break;
					case GDT_UInt16:
						if (null_value == *(((GUInt16*)inBuf) + iHeight*Bufsizex+iWidth))
						{
							break;
						}
						if (0 == iFile)
						{
							*(outBuf + iHeight*Bufsizex+iWidth) = (float)*(((GInt16*)inBuf) + iHeight*Bufsizex+iWidth);
						}else
						{
							*(outBuf + iHeight*Bufsizex+iWidth) += (float)*(((GUInt16*)inBuf) + iHeight*Bufsizex+iWidth);
						}
						counter[iHeight*Bufsizex+iWidth]++;
						break;
					case GDT_Int16: 
						if (null_value == *(((GInt16*)inBuf) + iHeight*Bufsizex+iWidth))
						{
							break;
						}
						if (0 == iFile)
						{
							*(outBuf + iHeight*Bufsizex+iWidth) = (float)*(((GInt16*)inBuf) + iHeight*Bufsizex+iWidth);
						}else
						{
							*(outBuf + iHeight*Bufsizex+iWidth) += (float)*(((GInt16*)inBuf) + iHeight*Bufsizex+iWidth);
						}
						counter[iHeight*Bufsizex+iWidth]++;
						break;
					case GDT_UInt32: 
						if (null_value == *(((GUInt32*)inBuf) + iHeight*Bufsizex+iWidth))
						{
							break;
						}
						if (0 == iFile)
						{
							*(outBuf + iHeight*Bufsizex+iWidth) = (float)*(((GUInt32*)inBuf) + iHeight*Bufsizex+iWidth);
						}else
						{
							*(outBuf + iHeight*Bufsizex+iWidth) += (float)*(((GUInt32*)inBuf) + iHeight*Bufsizex+iWidth);
						}
						counter[iHeight*Bufsizex+iWidth]++;
						break;
					case GDT_Int32: 
						if (null_value == *(((GInt32*)inBuf) + iHeight*Bufsizex+iWidth))
						{
							break;
						}
						if (0 == iFile)
						{
							*(outBuf + iHeight*Bufsizex+iWidth) = (float)*(((GInt32*)inBuf) + iHeight*Bufsizex+iWidth);
						}else
						{
							*(outBuf + iHeight*Bufsizex+iWidth) += (float)*(((GInt32*)inBuf) + iHeight*Bufsizex+iWidth);
						}
						counter[iHeight*Bufsizex+iWidth]++;
						break;
					case GDT_Float32: 
						if (null_value == *(((float*)inBuf) + iHeight*Bufsizex+iWidth))
						{
							break;
						}
						if (0 == iFile)
						{
							*(outBuf + iHeight*Bufsizex+iWidth) = (float)*(((float*)inBuf) + iHeight*Bufsizex+iWidth);
						}else
						{
							*(outBuf + iHeight*Bufsizex+iWidth) += (float)*(((float*)inBuf) + iHeight*Bufsizex+iWidth);
						}
						counter[iHeight*Bufsizex+iWidth]++;
						break;
					case GDT_Float64: 
						if (null_value == *(((double*)inBuf) + iHeight*Bufsizex+iWidth))
						{
							break;
						}
						if (0 == iFile)
						{
							*(outBuf + iHeight*Bufsizex+iWidth) = (float)*(((double*)inBuf) + iHeight*Bufsizex+iWidth);
						}else
						{
							*(outBuf + iHeight*Bufsizex+iWidth) += (float)*(((double*)inBuf) + iHeight*Bufsizex+iWidth);
						}
						counter[iHeight*Bufsizex+iWidth]++;
						break;
					}
				}
			}
			
			GDALClose( hDataset );
			CPLFree( inBuf );
		}

		for (int iHeight = 0; iHeight < Bufsizey; iHeight++)
		{
			for (int iWidth = 0; iWidth < Bufsizex; iWidth++)
			{
				int n = counter[iHeight*Bufsizex+iWidth];
				if (0==n)
				{
					*(outBuf + iHeight*Bufsizex+iWidth) = 0;
				}else
				{
					*(outBuf + iHeight*Bufsizex+iWidth) = *(outBuf + iHeight*Bufsizex+iWidth)/(double)n*0.02-273.15;
				}
			}
		}

		
		for (i=0;i<nBands;i++) {
			CPLErr eErr = pDstDataset->GetRasterBand(1)->RasterIO( GF_Write, 0, nYI*linenumber, Bufsizex , Bufsizey, outBuf, Bufsizex, Bufsizey, outeDT,  0, 0 );
		}
		CPLFree( outBuf );
	}
	// 清理
	GDALClose( hDstDS );

	
	cut_by_shp(tmpName, TCHAR2char(GetAbsolutePathName(pszCutlineDSName)), TCHAR2char(GetAbsolutePathName(pszDstFile)));

	// 删除临时文件
	remove(tmpName);
	return 1;
}