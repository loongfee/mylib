#pragma once
#ifndef strUtil_H
#define strUtil_H

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

#include <stdio.h>
#include <iostream>
#include <windows.h>
using namespace std;

static int SplitString(const string& input, 
				const string& delimiter, vector<string>& results, 
				bool includeEmpties)
{
	int iPos = 0;
	int newPos = -1;
	int sizeS2 = (int)delimiter.size();
	int isize = (int)input.size();

	if( 
		( isize == 0 )
		||
		( sizeS2 == 0 )
		)
	{
		return 0;
	}

	vector<int> positions;

	newPos = input.find (delimiter, 0);

	if( newPos < 0 )
	{ 
		return 0; 
	}

	int numFound = 0;

	while( newPos >= iPos )
	{
		numFound++;
		positions.push_back(newPos);
		iPos = newPos;
		newPos = input.find (delimiter, iPos+sizeS2);
	}

	if( numFound == 0 )
	{
		return 0;
	}

	for( int i=0; i <= (int)positions.size(); ++i )
	{
		string s("");
		if( i == 0 ) 
		{ 
			s = input.substr( i, positions[i] ); 
			if( includeEmpties || ( s.size() > 0 ) )
			{
				results.push_back(s);
			}
			continue;
		}
		int offset = positions[i-1] + sizeS2;
		if( offset < isize )
		{
			if( i == positions.size() )
			{
				s = input.substr(offset);
			}
			else if( i > 0 )
			{
				s = input.substr( positions[i-1] + sizeS2, 
					positions[i] - positions[i-1] - sizeS2 );
			}
		}
		if( includeEmpties || ( s.size() > 0 ) )
		{
			results.push_back(s);
		}
	}
	return numFound;
}


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
static int findDelimiter(const string& str, const vector<string>& delimiterList, int* indexOfDelimiter = NULL)
{
	// 获取分隔符的个数
	int num = static_cast<int>(delimiterList.size());
	int iPos = -1;	//定义一个游标
	int index_ = 0;
	for(int i = 0;i < num;i++)
	{
		//依次次查找各分隔符
		int tmp;
		if((tmp = str.find(delimiterList[i])) != -1)
		{
			//如果找到某分隔符
			if(-1 == iPos || tmp < iPos)
			{
				index_ = i;
				iPos = tmp;
			}
		}
	}

	if (indexOfDelimiter)
	{
		*indexOfDelimiter = index_;
	}
	//返回第一个分隔符的位置，，如果没有找到任何一个分隔符，则返回-1
	return iPos;
}

/************************************************************************/
/* 
 说明：用多个分割符拆分字符串
 参数：
const string& input 待拆分字符串
const vector<string>& delimiterList 为分隔符列表
 输出：vector<string>& results，存放拆分结果
*/
/************************************************************************/
static void splitString(const string& input, 
	const vector<string>& delimiterList, vector<string>& results)
{
	results.clear();

	int iPos=-1;   //定义一个游标
	int tmpPos = -1;

	string str = input;
	//删除字符串首的分隔符
	for (int i = 0;i < (int)delimiterList.size();++i)
	{
		string delimiter = delimiterList[i];
		int pos;
		while((pos = str.find (delimiter, 0)) == 0)
		{
			str = str.substr( pos+delimiter.size(), str.size() ); 
		}
	}

	while((iPos=findDelimiter(str, delimiterList))!=-1) //找到一个delimiter，索引时从0开始的
	{
		results.push_back(str.substr(0, iPos));//获取一个元素，并插入数组

		//删除该元素
		//str.erase(0,iPos+1);
		str = str.substr(iPos+1, str.size());

		//删除多余的分隔符
		int indexOfDelimiter;
		while(findDelimiter(str, delimiterList, &indexOfDelimiter) == 0)
		{
			// 如果字符串首存在分隔符
			// 则删除该分隔符
			str = str.substr(delimiterList[indexOfDelimiter].size(),  str.size());
		}
	}
	if(str != "")
	{
		// 如果最后的一个元素不为空
		// 则将最后的一个元素加入数组
		results.push_back(str);
	}
}

static string SBeforeLast(const string& str, char ch)
{
	int pos = str.find_last_of(ch);
	if (-1 == pos)
	{
		return "";
	}
	return str.substr(0, pos);
}

static string SAfterFirst(const string& str, char ch)
{
	int pos = str.find_first_of(ch);
	if (-1 == pos)
	{
		return "";
	}
	return str.substr(pos + 1, str.size() - 1);
}

static string SAfterLast(const string& str, char ch)
{
	int pos = str.find_last_of(ch);
	if (-1 == pos)
	{
		return "";
	}
	return str.substr(pos + 1, str.size() - 1);
}

static string SBeforeFirst(const string& str, char ch)
{
	int pos = str.find_first_of(ch);
	if (-1 == pos)
	{
		return "";
	}
	return str.substr(0, pos);
}

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

#endif