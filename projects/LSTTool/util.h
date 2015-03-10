#pragma once
// 使用UNICODE字符集
#ifndef UNICODE
#define UNICODE
#endif
// 如果不使用Unicode可以解锁下面的：
//#undef UNICODE
//#undef _UNICODE
#include <windows.h>
#include <string>
#include <vector>
// 版本控制
// 因为我不知道VC7.1以前的版本是否有tchar.h
#if _MSC_VER < 1310
#error "注意，请您使用VC7.1或者以上版本编译此程序"
#endif
#include <tchar.h>
#pragma message("请注意：如果您的文件夹下面有中文文件名的文件的话，最好使用UNICODE字符集")
// 提供字符串的定义有中文的时候最好使用wstring
#ifdef UNICODE
typedef std::wstring String;
#else
typedef std::string String;
#endif

/* 
*********************************************************************** 
* 函数： TCHAR2Char 
* 描述：将TCHAR* 转换为 char* 
* 日期：
*********************************************************************** 
*/ 
static char* TCHAR2char(TCHAR* tchStr) 
{ 
	int iLen = 2*wcslen(tchStr);//CString,TCHAR汉字算一个字符，因此不用普通计算长度 
	char* chRtn = new char[iLen+1];
	size_t converted = 0;
	wcstombs_s(&converted, chRtn, iLen+1, tchStr, _TRUNCATE);//转换成功返回为非负值 
	return chRtn; 
} 

static char* TCHAR2char(const TCHAR* tchStr) 
{ 
	int iLen = 2*wcslen(tchStr);//CString,TCHAR汉字算一个字符，因此不用普通计算长度 
	char* chRtn = new char[iLen+1];
	size_t converted = 0;
	wcstombs_s(&converted, chRtn, iLen+1, tchStr, _TRUNCATE);//转换成功返回为非负值 
	return chRtn; 
} 

/*
*********************************************************************** 
* 函数： char2tchar
* 描述：将 char* 转换为 TCHAR*
* 日期：
*********************************************************************** 
*/ 
static TCHAR *char2TCHAR(char *str)
{
	int iLen = strlen(str);
	TCHAR *chRtn = new TCHAR[iLen+1];
	//mbstowcs(chRtn, str, iLen+1); return chRtn;
	size_t converted = 0;
	mbstowcs_s(&converted, chRtn, iLen+1, str, _TRUNCATE);//转换成功返回为非负值
	return chRtn;
}

static TCHAR *char2TCHAR(const char *str)
{
	int iLen = strlen(str);
	TCHAR *chRtn = new TCHAR[iLen+1];
	//mbstowcs(chRtn, str, iLen+1); return chRtn;
	size_t converted = 0;
	mbstowcs_s(&converted, chRtn, iLen+1, str, _TRUNCATE);//转换成功返回为非负值
	return chRtn;
}


static LPWSTR GetAbsolutePathName(const char *relativePathName)
{	
	TCHAR *fullPath = new TCHAR[MAX_PATH];
	TCHAR *tchar_ = char2TCHAR(relativePathName);
	GetFullPathName(tchar_, MAX_PATH, fullPath, NULL);
	return fullPath;
}

static LPWSTR GetAbsolutePathName(char *relativePathName)
{	
	TCHAR *fullPath = new TCHAR[MAX_PATH];
	TCHAR *tchar_ = char2TCHAR(relativePathName);
	GetFullPathName(tchar_, MAX_PATH, fullPath, NULL);
	return fullPath;
}

static LPWSTR GetAbsolutePathName(TCHAR *relativePathName)
{	
	TCHAR *fullPath = new TCHAR[MAX_PATH];
	TCHAR *tchar_ = NULL;
	GetFullPathName(relativePathName, MAX_PATH, fullPath, NULL);
	return fullPath;
}

static LPWSTR GetAbsolutePathName(const TCHAR *relativePathName)
{	
	TCHAR *fullPath = new TCHAR[MAX_PATH];
	TCHAR *tchar_ = NULL;
	GetFullPathName(relativePathName, MAX_PATH, fullPath, NULL);
	return fullPath;
}

// 定义文件名的定义
typedef std::vector<LPWSTR> FilesVec;
// 查找当前目录下的所有文件返回所有文件的文件名请以FilesVec的引用的方式传递
// pathName 为路径名
static void FindFiles( const TCHAR * pszFilter ,FilesVec &files )
{
   WIN32_FIND_DATA FindFileData;
   HANDLE hFind = INVALID_HANDLE_VALUE;
   //TCHAR PathBuffer[ _MAX_PATH ];
  
   //_tcscpy_s( PathBuffer,_MAX_PATH, pathName );
   //_tcscat_s( PathBuffer,_MAX_PATH, _T("\\*") ); // 我也不知道为什么要添加，反正MSDN是这么做的，我只好照着做了
   hFind = ::FindFirstFile( pszFilter ,&FindFileData );
   if( INVALID_HANDLE_VALUE == hFind )      // 如果出现了某种异常就直接抛出便可
   {
      //char buffer[56];
      //sprintf_s( buffer,"Invalid File Handle.Error is %u\n",GetLastError() );
      //throw std::exception( buffer );
	   // 没有找到
	   return;
   }
   else // 然后再接着查找下去
   {
      files.push_back( GetAbsolutePathName(FindFileData.cFileName) );  
     
      while( 0 != FindNextFile( hFind,&FindFileData ) )
      {
		  LPWSTR fullPath = GetAbsolutePathName(FindFileData.cFileName);
         files.push_back( fullPath );
      }
      FindClose( hFind );
   }
}