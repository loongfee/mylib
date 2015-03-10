#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

#include "QVProcReader.h"
#include "ZY3QVProcReader.h"

#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")
#pragma comment(lib, "gdal_i.lib")

int main()
{
	//string strPath = "E:\\HJ1\\test";
	//string qvFile = strPath + "\\HJ1A_C1_027802_20131105_MY020_S1_1.dat";
	string strPath = "E:\\HJ1";
	//strPath = "I:\\testdata\\HJ\\QV";
	//string qvFile = strPath + "\\HJ-1A_CCD-1_KSC_201401070348_201401070354\\HJ-1A_CCD-1_KSC_201401070348_201401070354.dat";
	//string qvFile = strPath + "\\201403240251\\HJ-1A_CCD-1_MYC_201403240251_201403240259.dat";
	//string qvFile = strPath + "\\201403240114\\HJ-1A_CCD-1_MYC_201403240114_201403240125.dat";
	//string qvFile = strPath + "\\HJ-1B_CCD-1_MYC_201403230317_201403230320\\data.dat";
	//string qvFile = strPath + "\\HJ-1A_CCD-2\\HJ-1A_CCD-2_MYC_201406030321_201406030329\\data.dat";
	//string qvFile = strPath + "\\HJ-1B_CCD-1_SYC_201404020108_201404020114\\data.dat";
	//string qvFile = strPath + "\\HJ-1A_CCD-2_KSC_201403250456_201403250500\\data.dat";
	//string qvFile = strPath + "\\ZY03_BWD_MYC_201401111529_201401111535\\data.dat";
	string qvFile = strPath + "\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201407290219_201407290226.dat";
	QVHJ::QVProcReader qv;
	//qv.read(qvFile.c_str());
	//qv.read_band_seperate(qvFile.c_str());
	//qv.read_by_scene_band_seperate(qvFile.c_str());
	//qv.read_by_scene(qvFile.c_str());

	//qvFile = strPath + "\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201403260030_201403260036\\data.dat";
	//qv.read_band_seperate(qvFile.c_str());
	//qvFile = strPath + "\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201404110217_201404110229\\data.dat";
	//qv.read_band_seperate(qvFile.c_str());
	//qvFile = strPath + "\\HJ-1A_CCD-1\\HJ-1A_CCD-1_MYC_201405250250_201405250300\\data.dat";
	//qv.read_band_seperate(qvFile.c_str());

	qvFile = strPath + "\\ZY-3\\ZY-3_BWD_MYC_201409221253_201409221255.dat";
	string qvHeader = strPath + "\\ZY-3\\qvheader.dat";
	QVProc::ZY3QVProcReader zy3_qv;
	//zy3_qv.extractHeader(qvFile.c_str());
	zy3_qv.readHeader(qvHeader.c_str());
	//zy3_qv.read_by_scene(qvFile.c_str());
	return 1;
}		