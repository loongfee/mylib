#include "FileAndDirFinder.h"
#include "util.h"
#include <stdio.h>


CFileAndDirFinder::CFileAndDirFinder(void)
{
}


CFileAndDirFinder::~CFileAndDirFinder(void)
{
}


void CFileAndDirFinder::FindAllDir( const char* pCurDir, vector<string>& vtDirs )
{
	// 当前目录
	char szDir[MAX_PATH] = {0};
	
	sprintf_s(szDir, MAX_PATH, "%s\\*.*", pCurDir);


	WIN32_FIND_DATA findFileData = {0};
	HANDLE hFind = FindFirstFile(char2TCHAR(szDir), &findFileData);


	if (INVALID_HANDLE_VALUE == hFind)
	{
		return ;
	}


	do
	{
		/* 返回的文件名中会包含"."和".."。“.'代表本目录，".."代表上一层目录。
		一般情况下需要把这两个名称过滤掉。比如要进行文件删除操作
		*/
		if (findFileData.cFileName[0] != '.')//不是当前路径或者父目录的快捷方式
		{
			if(findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
			{
				// 这是一个普通目录
				char tmpDir[MAX_PATH] = {0};
				sprintf_s(tmpDir, MAX_PATH, "%s\\%s", pCurDir, TCHAR2char(findFileData.cFileName));
				vtDirs.push_back(tmpDir);


				// 递归调用查找子目录
				FindAllDir(tmpDir, vtDirs);
			}
		}
	}while (FindNextFile(hFind, &findFileData));


	FindClose(hFind);
	hFind = INVALID_HANDLE_VALUE;
}


void CFileAndDirFinder::FindAllFile( const char* pCurDir, const char* pFilter, vector<string>& vtFiles )
{
	// 当前目录
	char szDir[MAX_PATH] = {0};
	sprintf_s(szDir, MAX_PATH, "%s\\%s", pCurDir, pFilter);


	WIN32_FIND_DATA findFileData = {0};
	HANDLE hFind = FindFirstFile(char2TCHAR(szDir), &findFileData);


	if (INVALID_HANDLE_VALUE == hFind)
	{
		return ;
	}


	do
	{
		/* 返回的文件名中会包含"."和".."。“.'代表本目录，".."代表上一层目录。
		一般情况下需要把这两个名称过滤掉。比如要进行文件删除操作
		*/
		if (findFileData.cFileName[0] != '.')//不是当前路径或者父目录的快捷方式
		{
			// 这是一个文件
			char tmpFile[MAX_PATH] = {0};
			sprintf_s(tmpFile, MAX_PATH, "%s\\%s", pCurDir, TCHAR2char(findFileData.cFileName));


			vtFiles.push_back(tmpFile);
		}
	}while (FindNextFile(hFind, &findFileData));


	FindClose(hFind);
	hFind = INVALID_HANDLE_VALUE;
}

void CFileAndDirFinder::FindAllFileHere(const char* pFilter, vector<string>& vtFiles )
{
	// 当前目录
	WIN32_FIND_DATA findFileData = {0};
	HANDLE hFind = FindFirstFile(char2TCHAR(pFilter), &findFileData);


	if (INVALID_HANDLE_VALUE == hFind)
	{
		return ;
	}


	do
	{
		/* 返回的文件名中会包含"."和".."。“.'代表本目录，".."代表上一层目录。
		一般情况下需要把这两个名称过滤掉。比如要进行文件删除操作
		*/
		if (findFileData.cFileName[0] != '.')//不是当前路径或者父目录的快捷方式
		{
			// 这是一个文件
			char tmpFile[MAX_PATH] = {0};
			sprintf_s(tmpFile, MAX_PATH, "%s", TCHAR2char(findFileData.cFileName));


			vtFiles.push_back(tmpFile);
		}
	}while (FindNextFile(hFind, &findFileData));


	FindClose(hFind);
	hFind = INVALID_HANDLE_VALUE;
}


void CFileAndDirFinder::FindAllFileE( const char* pCurDir, const char* pFilter, vector<string>& vtFiles )
{
	// 获取子目录
	vector<string> vtDirs;
	FindAllDir(pCurDir, vtDirs);


	// 加入当前目录
	vtDirs.push_back(pCurDir);


	vector<string>::iterator itDir = vtDirs.begin();
	for (; itDir!=vtDirs.end(); ++itDir)
	{
		FindAllFile(itDir->c_str(), pFilter, vtFiles);
	}
}