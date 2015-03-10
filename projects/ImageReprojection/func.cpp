#include "func.h"
#include <QDateTime>


bool leap_year(int iyear)
{
	if ((iyear%4==0 && iyear%100!=0) || iyear%400==0)
	{
		return true;
	}
	return false;
}

void parseDate(QString strDay, QDate& dt)
{
	int nYear = atoi(strDay.left(4).toAscii());
	//int months_start_table[12];
	//memcpy(months_start_table, months_start, sizeof(months_start));
	int months_start_table[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 303, 304, 334};
	if (leap_year(nYear))
	{
		for (int i = 2; i < 12; i++)
		{
			months_start_table[i]++;
		}
	}
	int days = atoi(strDay.right(3).toAscii());
	int nMonth = 0;
	for (int i = 0; i < 12; i++)
	{
		if (months_start_table[i] < days && months_start_table[i+1] >= days)
		{
			nMonth = i+1;
			break;
		}
	}
	int nDay = days - months_start_table[nMonth-1];
	dt.setDate(nYear, nMonth, nDay);
}


bool FindFile(const QString& path, const QStringList& filter, QStringList &vFiles)
{
	QDir dir(path);
	if(!dir.exists())
		return false;
	dir.setFilter(QDir::Dirs | QDir::Files);
	dir.setSorting(QDir::DirsFirst);

	QFileInfoList list = dir.entryInfoList(filter, QDir::Files);
	int i = 0;
	for(i = 0;i < list.size();i++)
	{
		vFiles.push_back(list.at(i).absoluteFilePath());
	}

	i=0;
	list = dir.entryInfoList(QDir::Dirs);
	do{
		QFileInfo fileInfo=list.at(i);
		if(fileInfo.fileName()=="."||fileInfo.fileName()=="..")
		{
			i++;
			continue;
		}
		bool bisDir=fileInfo.isDir();
		if(bisDir)
		{
			FindFile(fileInfo.filePath(), filter, vFiles);
		}
		i++;
	}while(i<list.size());
	return true;
}

bool FindFile(const QString& path, const QStringList& filter, list<QString> &vFiles)
{
	QDir dir(path);
	if(!dir.exists())
		return false;
	dir.setFilter(QDir::Dirs | QDir::Files);
	dir.setSorting(QDir::DirsFirst);

	QFileInfoList list = dir.entryInfoList(filter, QDir::Files);
	int i = 0;
	for(i = 0;i < list.size();i++)
	{
		vFiles.push_back(list.at(i).absoluteFilePath());
	}

	i=0;
	list = dir.entryInfoList(QDir::Dirs);
	do{
		QFileInfo fileInfo=list.at(i);
		if(fileInfo.fileName()=="."||fileInfo.fileName()=="..")
		{
			i++;
			continue;
		}
		bool bisDir=fileInfo.isDir();
		if(bisDir)
		{
			FindFile(fileInfo.filePath(), filter, vFiles);
		}
		i++;
	}while(i<list.size());
	return true;
}

QString QBeforeLast(const QString& str, QChar ch)
{
	QStringList fields = str.split(ch);
	fields.removeLast();
	return fields.join(ch);
}

QString QAfterFirst(const QString& str, QChar ch)
{
	QStringList fields = str.split(ch);
	fields.removeFirst();
	return fields.join(ch);
}

QString QAfterLast(const QString& str, QChar ch)
{
	QStringList fields = str.split(ch);
	return fields.takeLast();
}

QString QBeforeFirst(const QString& str, QChar ch)
{
	QStringList fields = str.split(ch);
	return fields.takeFirst();
}

bool CopyFolder(QString pathFrom, QString pathTo, QStringList filters/* = QStringList("*.*")*/, bool bOverWrite/* = true*/)
{
	QStringList FileList;
	FindFile(pathFrom, filters, FileList);	//找出对应的全色文件	
	
	if (!QDir(pathTo).exists())
	{
		_mkdir(pathTo.toUtf8());
	}
	for(int i = 0;i < FileList.count();i++)
	{
		QFileInfo fileinfo(FileList.at(i));
		QString newpath = QDir::toNativeSeparators(pathTo + "\\" + fileinfo.fileName());
		if(QFile(newpath).exists() && bOverWrite)
		{
			QFile(newpath).remove();
		}
		QFile::copy(FileList.at(i), newpath);
	}
	return true;
}


bool bandMerge(std::vector<string>& fileList,string createtiffile)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDataset  *hSrcDataset, *outtiff;
	GDALDriverH		hDriver;
	GDALDatasetH	hSrcDs, hDstDS;
	CPLErr      eErr = CE_None,eErr1=CE_None;
	int i;
	int nRasterXSizeRead,nRasterYSizeRead;
	GDALDataType eDT;
	GDALRasterBand *rasterband,*rasterband_IMG;
	GByte   *pabyLine,*bufout;

	int nBands = fileList.size();
	if(nBands < 1) return false;
	
	hSrcDs=GDALOpen(fileList[0].c_str(), GA_ReadOnly );
	hSrcDataset = (GDALDataset  *)hSrcDs;

	rasterband_IMG=hSrcDataset->GetRasterBand(1);
	eDT=rasterband_IMG->GetRasterDataType();
	nRasterXSizeRead=hSrcDataset->GetRasterXSize();
	nRasterYSizeRead=hSrcDataset->GetRasterYSize();

	const char         *pszFormat = "GTiff";
	char               *pszTargetSRS = NULL;
	double adfThisGeoTransform[6];

	hDriver = GDALGetDriverByName( pszFormat );
	hDstDS = GDALCreate( hDriver, createtiffile.c_str(), nRasterXSizeRead, nRasterYSizeRead, 
                        nBands, eDT, NULL );
    
	if( hDstDS == NULL ) {	
		GDALClose( hSrcDataset );
		cout<<createtiffile<<"创建失败！"; 
		return false;
	}
       
	
	GDALGetGeoTransform( hSrcDs, adfThisGeoTransform );
	pszTargetSRS = CPLStrdup(GDALGetProjectionRef(hSrcDs));

    GDALSetProjection( hDstDS, pszTargetSRS );
    GDALSetGeoTransform( hDstDS, adfThisGeoTransform );
    outtiff=(GDALDataset  *)hDstDS;
	GDALClose( hSrcDataset );

	////////////////////////
	// 一次处理多行数据，提高速度
	int linenumber = 2000;
	int nyNum=(nRasterYSizeRead-1)/linenumber+1;//计算列方向块数
	for (int nYI=0;nYI<nyNum;nYI++)//块循环
	{
		cout<<"波段合并： "<<nYI+1<<" of "<<nyNum<<"..."<<endl;
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
			hSrcDs=GDALOpen(fileList[i].c_str(), GA_ReadOnly );
			hSrcDataset = (GDALDataset  *)hSrcDs;
			rasterband_IMG=hSrcDataset->GetRasterBand(1);
			rasterband=outtiff->GetRasterBand(i+1);
				 eErr = rasterband_IMG->RasterIO( GF_Read, 0, nYI*linenumber, Bufsizex , Bufsizey, pabyLine, Bufsizex, Bufsizey, eDT,  0, 0 );
				 eErr1 = rasterband->RasterIO( GF_Write, 0, nYI*linenumber, Bufsizex , Bufsizey,pabyLine, Bufsizex, Bufsizey, eDT,  0, 0 ); 
				 GDALClose( hSrcDataset );
		}
		CPLFree( pabyLine );
		CPLFree( bufout );
	}

	//int iCount = 0;

	//for (i=0;i<nBands;i++) {
	//	hSrcDs=GDALOpen(fileList[i].c_str(), GA_ReadOnly );
	//	hSrcDataset = (GDALDataset  *)hSrcDs;
	//	rasterband_IMG=hSrcDataset->GetRasterBand(1);
	//	rasterband=outtiff->GetRasterBand(i+1);
	//		 eErr = rasterband_IMG->RasterIO( GF_Read, 0, 0,nRasterXSizeRead , linenumber*nRasterYSizeRead, pabyLine, nRasterXSizeRead, linenumber*nRasterYSizeRead, eDT,  0, 0 );
	//		 eErr1 = rasterband->RasterIO( GF_Write, 0, 0,nRasterXSizeRead , linenumber*nRasterYSizeRead,pabyLine, nRasterXSizeRead, linenumber*nRasterYSizeRead, eDT,  0, 0 ); 
	//		 GDALClose( hSrcDataset );
	//}

	//CPLFree( pabyLine );
	//CPLFree( bufout );
	GDALClose( outtiff );
	return true;
}