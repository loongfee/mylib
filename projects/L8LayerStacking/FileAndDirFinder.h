#pragma once
#include <string>
#include <vector>
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

class CFileAndDirFinder
{
public:
	CFileAndDirFinder(void);
	~CFileAndDirFinder(void);


	// 查找当前目录下的所有目录(不包括当前目录)
	void FindAllDir(const char* pCurDir, vector<string>& vtDirs);
	// 查找当前目录下的所有文件(不包括子目录)，制定查找文件类型，如：*.txt,*.lua,*.*
	void FindAllFile(const char* pCurDir, const char* pFileType, vector<string>& vtFiles);
	// 查找当前目录下的所有文件(不包括子目录)
	void FindAllFileHere(const char* pFilter, vector<string>& vtFiles );
	// 查找当前目录下的所有文件(包括子目录)，制定查找文件类型，如：*.txt,*.lua,*.*
	void FindAllFileE(const char* pCurDir, const char* pFileType, vector<string>& vtFiles);


};