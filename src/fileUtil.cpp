#include "fileUtil.h"
#include "strUtil.h"

namespace mylib
{

#ifdef USE_QT
bool QFindFile(const QString& path, const QStringList& filter, QStringList &vFiles)
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
			QFindFile(fileInfo.filePath(), filter, vFiles);
		}
		i++;
	}while(i<list.size());
	return true;
}

bool QFindFile(const QString& path, const QStringList& filter, list<QString> &vFiles)
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
			QFindFile(fileInfo.filePath(), filter, vFiles);
		}
		i++;
	}while(i<list.size());
	return true;
}


bool QCopyFolder(QString pathFrom, QString pathTo, QStringList filters/* = QStringList("*.*")*/, bool bOverWrite/* = true*/)
{
	QStringList FileList;
	QFindFile(pathFrom, filters, FileList);

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
#endif

char* getExePath()
{
	HMODULE hModule = GetModuleHandleW(NULL);
	TCHAR buf[MAX_PATH];
	GetModuleFileNameW(hModule, buf, MAX_PATH);
	return TCHAR2CHAR(buf);
}

} // end of namespace mylib