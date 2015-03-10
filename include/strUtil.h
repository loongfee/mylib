#pragma once
#ifndef STR_UTIL_HPP
#define STR_UTIL_HPP

#include <ossim/base/ossimString.h>

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <wtypes.h>
#include <atltypes.h>

#ifdef USE_QT
#include <QStringList>
#include <QString>
#endif

#include <string>
#include <vector>
#include <iostream>

#include <sstream>
#include <algorithm>
#include <iterator>

#include <locale>
#include <codecvt>
#include <string>

using namespace std;

namespace mylib{
//static wxString ossimString2wxString(const ossimString& rs)
//{
//	return wxString::FromUTF8(rs.c_str());
//}
//
//static ossimString wxString2ossimString(const wxString& ws)
//{
//	return ossimString(ws.mb_str());
//}

/************************************************************************/
/*           string                                                     */
/************************************************************************/
int SplitString(const string& input, 
				const string& delimiter, vector<string>& results, 
				bool includeEmpties);


/************************************************************************/
/* 
 说明：查找字符串中第一次出现分隔符的位置
 参数：
const string& str 待查找的字符串
const vector<string>& delimiterList为分隔符列表
int* indexOfDelimiter 为第一次出现的分隔符序号
 返回值：
 第一次出现分隔符的位置，如果没有找到任何一个分隔符，则返回-1
*/
/************************************************************************/
int findDelimiter(const string& str, const vector<string>& delimiterList, int* indexOfDelimiter = NULL);

/************************************************************************/
/* 
 说明：用多个分割符拆分字符串
 参数：
const string& input 待拆分字符串
const vector<string>& delimiterList 为分隔符列表
 输出：vector<string>& results，存放拆分结果
*/
/************************************************************************/
void splitString(const string& input, 
	const vector<string>& delimiterList, vector<string>& results);


/************************************************************************/
/*           wstring                                                     */
/************************************************************************/
int SplitString(const wstring& input,
	const wstring& delimiter, vector<wstring>& results,
	bool includeEmpties);


/************************************************************************/
/*
说明：查找字符串中第一次出现分隔符的位置
参数：
const wstring& str 待查找的字符串
const vector<wstring>& delimiterList为分隔符列表
int* indexOfDelimiter 为第一次出现的分隔符序号
返回值：
第一次出现分隔符的位置，如果没有找到任何一个分隔符，则返回-1
*/
/************************************************************************/
int findDelimiter(const wstring& str, const vector<wstring>& delimiterList, int* indexOfDelimiter = NULL);

/************************************************************************/
/*
说明：用多个分割符拆分字符串
参数：
const wstring& input 待拆分字符串
const vector<wstring>& delimiterList 为分隔符列表
输出：vector<wstring>& results，存放拆分结果
*/
/************************************************************************/
void splitString(const wstring& input,
	const vector<wstring>& delimiterList, vector<wstring>& results);

//std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
//	std::stringstream ss(s);
//	std::string item;
//	while (std::getline(ss, item, delim)) {
//		elems.push_back(item);
//	}
//	return elems;
//}
//
//
//std::vector<std::string> split(const std::string &s, char delim) {
//	std::vector<std::string> elems;
//	split(s, delim, elems);
//	return elems;
//}

string SBeforeLast(const string& str, char ch);

string SAfterFirst(const string& str, char ch);

string SAfterLast(const string& str, char ch);

string SBeforeFirst(const string& str, char ch);


/************************************************************************/
/* ossimString                                                           */
/************************************************************************/
/************************************************************************/
/* 
 说明：查找字符串中第一次出现分隔符的位置
 参数：
 CString str 待查找的字符串
 vector<char> chList为分隔符列表
 返回值：
 第一次出现分隔符的位置，如果没有找到任何一个分隔符，则返回-1
*/
/************************************************************************/
int findChar(ossimString str, vector<char> chList);


/************************************************************************/
/* 
 说明：用多个分割符拆分字符串
 参数：
 ossimString str 待拆分字符串
 vector<char> chList为分隔符列表
 输出：vector<ossimString>& strArray，存放拆分结果
*/
/************************************************************************/
void splitString(ossimString str, vector<char> chList, vector<ossimString>& strArray);

ossimString ossimSBeforeLast(const ossimString& str, char ch);

ossimString ossimSAfterFirst(const ossimString& str, char ch);

ossimString ossimSAfterLast(const ossimString& str, char ch);

ossimString ossimSBeforeFirst(const ossimString& str, char ch);

#ifdef USE_QT
/************************************************************************/
/* QString                                                              */
/************************************************************************/
QString QBeforeLast(const QString& str, QChar ch);

QString QAfterFirst(const QString& str, QChar ch);

QString QAfterLast(const QString& str, QChar ch);

QString QBeforeFirst(const QString& str, QChar ch);
#endif

std::string ConvertFromUtf16ToUtf8(const std::wstring& wstr);

std::wstring ConvertFromUtf8ToUtf16(const std::string& str);

std::string TCHAR2STRING(TCHAR *STR);

char* TCHAR2CHAR(TCHAR *STR);

TCHAR* CHAR2TCHAR(CHAR *STR);


char* WCHAR2CHAR(wchar_t *STR);
const wchar_t* CHAR2WCHAR(const char *STR);

}; // end of namespace mylib
#endif