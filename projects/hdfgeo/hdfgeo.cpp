#include "gdalwarper.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "ogr_api.h"
#include "gdal_priv.h"
#include "gdal.h"
#include <vector>
#include <iostream>
#include <string>
#include <func.h>
#include <stdlib.h>
#include <boost/filesystem.hpp>

#include <ossim/imaging/ossimImageHandlerFactory.h>
#include <ossim/projection/ossimImageViewProjectionTransform.h>
#include <ossim/base/ossimPreferences.h>
#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
namespace fs = boost::filesystem;
namespace po = boost::program_options;

#undef BOOST_WINDOWS_API


using namespace std;
using namespace mylib;


const char * pszHdfFile = "";
const char * pszRefFile = "";
const char * pszOutFile = "";
const char * pszHegPath = "C:\\HEGtools\\HEG_Win";

#ifndef _WIN64
#pragma comment(lib, "ossim20.lib")
#else
#pragma comment(lib, "ossim20x64.lib")
#endif

#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "mlpack.lib")

//const char* TAG_LAT = "Geodedic_Latitude";
//const char* TAG_LON = "Geodedic_Longitude";
const char* TAG_LAT = "Latitude";
const char* TAG_LON = "Longitude";
//const char* TAG_VAPOR = "Water_Vapor_Near_Infrared";
const char* pszPreferenceFile = "preference.txt";
//const char* TAG = "Water_Vapor_Near_Infrared";
string dimensionName = "Water_Vapor_Near_Infrared";

int NULL_VALUE = -9999;
double add_offset = 0;
//const float value_scale = 1e-3f;
double scale_factor = 0.00100000004749745;

const int percentList[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

string getTempImageFilename(string ext = "tif")
{
	char temp_img_file_name[MAX_PATH];
	char temp_path[MAX_PATH];

	//GetTempPath(sizeof(temp_path), (wchar_t*)temp_path);
	//GetTempFileName((wchar_t*)temp_path, (wchar_t*)"", 0, (wchar_t*)temp_img_file_name);
	GetTempPath(sizeof(temp_path), temp_path);
	GetTempFileName(temp_path, "", 0, temp_img_file_name);

	//QString absolutePath = QFileInfo(temp_img_file_name).absolutePath();
	//QString absolutePath = QBeforeLast(QString(temp_img_file_name), '.') + "." + QString(ext.c_str());
	string absolutePath = SBeforeLast(temp_img_file_name, '.') + "." + ext;

	//return absolutePath.toStdString();
	return absolutePath;
}

string getTempName(string hdfFile)
{
	string tmp_vapor_file = SBeforeLast(string(hdfFile), '.') + "_tmp.tif";
	return tmp_vapor_file;
}

string getOutName(string hdfFile)
{
	string out_file = SBeforeLast(hdfFile, '\\') + "\\MOD05_" + SAfterLast(hdfFile, '\\');
	out_file = SBeforeLast(out_file, '.') + ".tif";
	//string out_file = SBeforeLast(string(hdfFile), '.') + "_MOD05.tif";
	return out_file;
}

string getPrmName(string hdfFile)
{
	string out_file = SBeforeLast(string(hdfFile), '.') + ".prm";
	return out_file;
}

string getBatchName(string hdfFile)
{
	string out_file = SBeforeLast(string(hdfFile), '.') + ".bat";
	return out_file;
}

void create_prm(string hdfFile, string refFile, string outFile)
{
	GDALDataset *poDataset = (GDALDataset *)GDALOpen(hdfFile.c_str(), GA_ReadOnly);
	if (!poDataset) {
		std::cerr << "cannot open " << hdfFile << std::endl;
		exit(1);
	}

	//char **metadata = poDataset->GetMetadata("SUBDATASETS");
	//for(int  i = 0; metadata[i] != NULL; i++ )
	//{
	//	cout<<metadata[i]<<endl;
	//}


	// read dimension sizes
	int xsize;
	int ysize;

	GDALDriver * driver = poDataset->GetDriver ();
	string papszMetadata = GDALGetDriverShortName((GDALDriverH)poDataset);
	int index=papszMetadata.find_first_of("hdf");
	if(index<0)
		index=papszMetadata.find_first_of("HDF");
	if(!index)
	{
		std::cerr << "invalid MOD05_L2 file. " << std::endl;
		exit(1);
	}

	string swathName = "mod04";

	char ** SUBDATASETS = GDALGetMetadata( (GDALDatasetH)poDataset, "SUBDATASETS" );
	//if (CSLCount(SUBDATASETS) > 0)
	//{
	//	for (int i = 0; SUBDATASETS[i] != NULL; i++)
	//	{
	//		cout << SUBDATASETS[i] << endl;
	//	}
	//}
	GDALDataset * lat_dt;
	GDALDataset * lon_dt;
	GDALDataset * data_dt;
	if( CSLCount(SUBDATASETS) > 0 )
	{
		for(int  i = 0; SUBDATASETS[i] != NULL; i++ )
		{
			if(i%2==0)
			{
				string tmpfilename = string(SUBDATASETS[i]);
				tmpfilename = tmpfilename.substr(tmpfilename.find_first_of("=") + 1);
				string DESC = string(SUBDATASETS[i+1]);
				DESC = DESC.substr(DESC.find_last_of("=") + 1);
				//string tag_str = tmpfilename.substr(tmpfilename.find_last_of(":") + 1);
				//if (0 == tag_str.compare(TAG_LAT))
				//{
				//	lat_dt = (GDALDataset *) GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
				//}
				//if (0 == tag_str.compare(TAG_LON))
				//{
				//	lon_dt = (GDALDataset *) GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
				//}
				vector<string> strListTmp;
				mylib::splitString(DESC, { " " }, strListTmp);
				if (strListTmp.size() < 3)
				{
					cerr << "description of " << tmpfilename << " is not valid." << endl;
					return;
				}
				if (0 == strListTmp[1].compare(TAG_LAT))
				{
					lat_dt = (GDALDataset *)GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
				}
				if (0 == strListTmp[1].compare(TAG_LON))
				{
					lon_dt = (GDALDataset *)GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
				}
				if (0 == strListTmp[1].compare(dimensionName)
					&& tmpfilename.find(dimensionName) != std::string::npos)
				{
					data_dt = (GDALDataset *)GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
					swathName = strListTmp[2];
				}
			}            
		}
	}

	int xsize_lat = lat_dt->GetRasterXSize();
	int ysize_lat = lat_dt->GetRasterYSize();
	int nb_lat = lat_dt->GetRasterCount();
	GDALDataType eDataType_lat = lat_dt->GetRasterBand(1)->GetRasterDataType();
	int datatypeSize_lat = GDALGetDataTypeSize(eDataType_lat) / 8;
	GByte *buffer_lat = new GByte [xsize_lat * ysize_lat * nb_lat * datatypeSize_lat];
	GByte *buffer_lon = new GByte [xsize_lat * ysize_lat * nb_lat * datatypeSize_lat];

	vector<int> bandList;
	for (int iband = 0; iband < nb_lat;++iband)
	{
		bandList.push_back(iband+1);
	}
	lat_dt->RasterIO(GF_Read, 0, 0, xsize_lat, ysize_lat, buffer_lat, xsize_lat, ysize_lat, eDataType_lat, nb_lat, &bandList[0], 0, 0, 0);
	lon_dt->RasterIO(GF_Read, 0, 0, xsize_lat, ysize_lat, buffer_lon, xsize_lat, ysize_lat, eDataType_lat, nb_lat, &bandList[0], 0, 0, 0);

	vector<float> latList;
	vector<float> lonList;
	latList.push_back(((float*)buffer_lat)[0 * xsize_lat + 0]);
	latList.push_back(((float*)buffer_lat)[0 * xsize_lat + (xsize_lat-1)]);
	latList.push_back(((float*)buffer_lat)[(ysize_lat-1) * xsize_lat + (xsize_lat-1)]);
	latList.push_back(((float*)buffer_lat)[(ysize_lat-1) * xsize_lat + 0]);

	lonList.push_back(((float*)buffer_lon)[0 * xsize_lat + 0]);
	lonList.push_back(((float*)buffer_lon)[0 * xsize_lat + (xsize_lat-1)]);
	lonList.push_back(((float*)buffer_lon)[(ysize_lat-1) * xsize_lat + (xsize_lat-1)]);
	lonList.push_back(((float*)buffer_lon)[(ysize_lat-1) * xsize_lat + 0]);

	float minLat = *std::min_element(latList.begin(), latList.end());
	float maxLat = *std::max_element(latList.begin(), latList.end());
	float minLon = *std::min_element(lonList.begin(), lonList.end());
	float maxLon = *std::max_element(lonList.begin(), lonList.end());

	//float maxLat = ((float*)buffer_lat)[0 * xsize_lat + 0];
	//float minLon = ((float*)buffer_lon)[0 * xsize_lat + 0];
	//float minLat = ((float*)buffer_lat)[(ysize_lat-1) * xsize_lat + (xsize_lat-1)];
	//float maxLon = ((float*)buffer_lon)[(ysize_lat-1) * xsize_lat + (xsize_lat-1)];

	//float OUTPUT_PIXEL_SIZE_X = (-lonList[0]-lonList[3]+lonList[1]+lonList[2]) * 0.1 / xsize_lat;
	//float OUTPUT_PIXEL_SIZE_Y = (latList[0]+latList[1]-latList[2]-latList[3]) * 0.1 / ysize_lat;
	//float OUTPUT_PIXEL_SIZE_X = (latList[0]+latList[1]-latList[2]-latList[3]) * 0.1 / xsize_lat;
	//float OUTPUT_PIXEL_SIZE_Y = (-lonList[0]-lonList[3]+lonList[1]+lonList[2]) * 0.1 / ysize_lat;
	//float OUTPUT_PIXEL_SIZE_X = (min(lonList[1], lonList[2]) - max(lonList[0], lonList[3])) * 0.2 / xsize_lat;
	//float OUTPUT_PIXEL_SIZE_Y = (min(latList[0], latList[1]) - max(latList[2], latList[3])) * 0.2 / ysize_lat;
	float mean_lat = (minLat + maxLat) * 0.5;
	float OUTPUT_PIXEL_SIZE_X = fabs((-lonList[0]+lonList[2]) * 0.2 / xsize_lat * cos(mean_lat*M_PI/180.0));
	float OUTPUT_PIXEL_SIZE_Y = fabs((latList[0] - latList[2]) * 0.2 / ysize_lat);
	//float OUTPUT_PIXEL_SIZE_X = (min(lonList[1], lonList[2]) - max(lonList[0], lonList[3])) * 0.2 / xsize_lat;
	//float OUTPUT_PIXEL_SIZE_Y = (min(latList[0], latList[1]) - max(latList[2], latList[3])) * 0.2 / ysize_lat;

	if (NULL == data_dt)
	{
		cerr << "Dimension \"" << dimensionName << "\" is found." << endl;
		exit(1);
	}
	//get data dataset info
	char ** metaData = data_dt->GetMetadata();
	if (CSLCount(metaData) > 0)
	{
		for (int i = 0; metaData[i] != NULL; i++)
		{
			vector<string> strListTmp;
			mylib::splitString(metaData[i], { "=" }, strListTmp);
			string strTagName = strListTmp[0];
			boost::algorithm::trim_if(strTagName, boost::algorithm::is_any_of(" "));
			string strTagValue = strListTmp[1];
			boost::algorithm::trim_if(strTagValue, boost::algorithm::is_any_of(" "));
			if (0 == strTagName.compare("_FillValue"))
			{
				NULL_VALUE = stoi(strTagValue);
			}
			else if (0 == strTagName.compare("add_offset"))
			{
				add_offset = stod(strTagValue);
			}
			else if (0 == strTagName.compare("scale_factor"))
			{
				scale_factor = stod(strTagValue);
			}
		}
	}

	//string tmp_vapor_file = getTempName(hdfFile);
	string tmp_data_file =  getTempImageFilename("tif");
	int task_num = 1;
	//string prm_file = getPrmName(hdfFile);
	string prm_file =  getTempImageFilename("prm");
	fstream fs;
	fs.open(prm_file.c_str(), ios_base::out);
	fs<<"\n";
	fs<<"NUM_RUNS = "<<task_num<<"\n";
	fs<<"\n";
	fs<<"BEGIN"<<"\n";
	fs<<"INPUT_FILENAME = "<<hdfFile<<"\n";
	fs << "OBJECT_NAME = " << swathName << "\n";
	fs << "FIELD_NAME = " << dimensionName << "|" << "\n";
	fs<<"BAND_NUMBER = "<<"1"<<"\n";
	fs<<"OUTPUT_PIXEL_SIZE_X = "<<OUTPUT_PIXEL_SIZE_X<<"\n";
	fs<<"OUTPUT_PIXEL_SIZE_Y = "<<OUTPUT_PIXEL_SIZE_Y<<"\n";
	fs<<"SPATIAL_SUBSET_UL_CORNER = ( "<<maxLat<<" "<<minLon<<" )"<<"\n";
	fs<<"SPATIAL_SUBSET_LR_CORNER = ( "<<minLat<<" "<<maxLon<<" )"<<"\n";
	fs<<"RESAMPLING_TYPE = "<<"NN"<<"\n";
	fs << "OUTPUT_FILENAME = " << tmp_data_file.c_str() << "\n";
	fs << "OUTPUT_PROJECTION_TYPE = GEO"<< "\n";
	fs << "ELLIPSOID_CODE = WGS84" << "\n";		
	fs<<"OUTPUT_TYPE = "<<"GEO"<<"\n";
	fs<<"END"<<"\n";
	fs<<"\n";
	fs.close();
	//FILE *pf = fopen(prm_file.c_str(), "w+");
	//fprintf(pf, "\n");
	//fprintf(pf, "NUM_RUNS = %d\n", task_num);
	//fclose(pf);
	
	if (fs::exists(tmp_data_file))
	{
		DeleteFileA(tmp_data_file.c_str());
	}

	//string batFile = getBatchName(hdfFile);
	string batFile = getTempImageFilename("bat");
	fstream batfs;
	batfs.open(batFile.c_str(), ios_base::out);
	batfs<<"@echo off"<<endl;
	batfs<<"set CYGWIN=nodosfilewarning"<<endl;
	batfs<<"set LD_LIBRARY_PATH="<<pszHegPath<<"\\bin"<<endl;
	batfs<<"set MRTDATADIR="<<pszHegPath<<"\\data"<<endl;
	batfs<<"set MRTBINDIR="<<pszHegPath<<"\\bin"<<endl;
	batfs<<"set PGSHOME="<<pszHegPath<<"\\TOOLKIT_MTD"<<endl;


	batfs<<"dos2unix "<<prm_file<<endl;
	batfs<<"swtif.exe -p "<<prm_file<<endl;
	batfs.close();

	system(batFile.c_str());



	//////////////////////////////////////////////////////////////////////////
	// clip
	printf( "cliping and resampling... \n" );
	ossimRefPtr<ossimImageHandler> refHandler = ossimImageHandlerFactory::instance()->open(ossimFilename(refFile));
	if (!refHandler)
	{
		std::cerr << "cannot open " << refFile << std::endl;
		exit(1);
	}
	ossimIrect refBoundary = refHandler->getBoundingRect();
	ossimRefPtr<ossimMapProjection> refProjection = PTR_CAST(ossimMapProjection,
		refHandler->getImageGeometry()->getProjection());

	ossimRefPtr<ossimImageHandler> tmpHandler = ossimImageHandlerFactory::instance()->open(ossimFilename(tmp_data_file.c_str()));
	if (!tmpHandler)
	{
		std::cerr << "cannot open " << tmp_data_file.c_str() << std::endl;
		exit(1);
	}
	ossimRefPtr<ossimMapProjection> srcProjection = PTR_CAST(ossimMapProjection,
		tmpHandler->getImageGeometry()->getProjection());
	
	ossimRefPtr<ossimImageRenderer> renderer = new ossimImageRenderer;
	//renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_NEAREST_NEIGHBOR);
	renderer->connectMyInputTo(tmpHandler.get());
	renderer->setView(refProjection.get());

	vector<ossimDpt> polygon;
	ossimDpt dpt;
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ul()), dpt);
	polygon.push_back(dpt);
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ur()), dpt);
	polygon.push_back(dpt);
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.lr()), dpt);
	polygon.push_back(dpt);
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ll()), dpt);
	polygon.push_back(dpt);
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ul()), dpt);
	polygon.push_back(dpt);
	ossimRefPtr<ossimPolyCutter> theCutter = new ossimPolyCutter;
	theCutter->connectMyInputTo(tmpHandler.get());
	theCutter->setPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, renderer->getImageViewTransform());

	//统计外接矩形
	dpt = renderer->getImageViewTransform()->imageToView(polygon[0]);
	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
	for (int i = 1; i < (int) polygon.size(); i++)
	{
		dpt = renderer->getImageViewTransform()->imageToView(polygon[i]);
		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
	}
	int offset = 0;
	dpt = renderer->getImageViewTransform()->viewToImage(UL1);
	ossimDpt UL(dpt.x - offset, dpt.y - offset);
	dpt = renderer->getImageViewTransform()->viewToImage(LR1);
	ossimDpt LR(dpt.x + offset, dpt.y + offset);
	UL1 = renderer->getImageViewTransform()->imageToView(UL);
	LR1 = renderer->getImageViewTransform()->imageToView(LR);
	//ossimDrect viewRegion = calcRegion(clipRect, inmapinfo, renderer->getImageViewTransform());
	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
	ossimDpt UR = renderer->getImageViewTransform()->viewToImage(ossimDpt(LR1.lon, UL1.lat));
	ossimDpt LL = renderer->getImageViewTransform()->viewToImage(ossimDpt(UL1.lon, LR1.lat));

	vector<ossimDpt> clipRect;
	clipRect.clear();
	clipRect.push_back(UL);
	clipRect.push_back(UR);
	clipRect.push_back(LR);
	clipRect.push_back(LL);
	clipRect.push_back(UL);

	theCutter->setNumberOfPolygons(0);
	theCutter->addPolygon(clipRect);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	theCutter->connectMyInputTo(tmpHandler.get());

	renderer->connectMyInputTo(theCutter.get());

	ossimFilename clipTempFile = getTempImageFilename("tif");

	ossimRefPtr<ossimImageFileWriter> writer = ossimImageWriterFactoryRegistry::instance()->
		createWriterFromExtension( clipTempFile.ext() );
	if(writer==NULL) return;

	ossimKeywordlist tt_geom;
	tt_geom.clear();
	writer->saveState(tt_geom);
	tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA",true);
	writer->loadState(tt_geom);


	ossimStdOutProgress* progress = new ossimStdOutProgress(0,true);
	writer->addListener(progress);

	writer->setFilename(clipTempFile);
	writer->connectMyInputTo(0, renderer.get());
	writer->setAreaOfInterest(viewRegion);
	//ossimIrect bounding = handler->getBoundingRect();
	//writer->setAreaOfInterest(bounding);

	//bool bResult = writer->execute();

	if(!writer->execute())
	{
		writer->removeListener(progress);
		std::cerr << "geometric correction failed." << std::endl;
		exit(1);
	}
	refHandler->close();
	tmpHandler->close();

	//printf( "it's ok! please enter any key! \n" );

	//////////////////////////////////////////////////////////////////////////
	printf( "fill null pixels... \n" );

	string outFilename;
	if ("" == outFile)
	{
		outFilename = getOutName(refFile);
	}
	else
	{
		outFilename = outFile;
	}
	GDALDataset *pSrcDataset;
	pSrcDataset = (GDALDataset *) GDALOpen( clipTempFile, GA_ReadOnly );
	double adfThisGeoTransform[6];
	pSrcDataset->GetGeoTransform(adfThisGeoTransform);
	const char* pszTargetSRS = pSrcDataset->GetProjectionRef();
	int nX = pSrcDataset->GetRasterXSize();
	int nY = pSrcDataset->GetRasterYSize();
	int nBand = pSrcDataset->GetRasterCount();
	GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
	int dataTypeSize = GDALGetDataTypeSize(eDT) / 8;
	if (1 != nBand)
	{
		std::cerr << "water vapor file should have 1 band." << std::endl;
		exit(1);
	}
	bandList.clear();
	for (int iband = 0; iband < nBand;++iband)
	{
		bandList.push_back(iband+1);
	}
	GByte* pSrcData = new GByte[nX * nY * nBand * dataTypeSize];
	float* pDstData = new float[nX * nY * nBand];
	int nCount = 0;
	double mean_value = 0.0;
	pSrcDataset->RasterIO(GF_Read, 0, 0, nX, nY, pSrcData, nX, nY, eDT, nBand, &bandList[0], 0, 0, 0);

	int fill_pos = 0;
	cout<<percentList[fill_pos++]<<"%";
	int percent = 0;
	float total_size = nX * nY * 2;
	for (int i = 0; i < nX; i++)
	{
		for (int j = 0;j < nY;j++)
		{
			GInt16 d = ((GInt16*)pSrcData)[j * nX + i];
			if (NULL_VALUE != d && d > 0.0)
			{
				mean_value += (d - add_offset) * scale_factor;
				nCount++;
			}
			percent += 1;
			if (int(percent*100.0f/(float)total_size+0.5) == percentList[fill_pos])
			{
				cout<<" "<<percentList[fill_pos]<<"%";
				fill_pos++;
			}			
		}
	}
	if (0 == nCount)
	{
		// all pixels are null 
		cout << " 100%";
		cout << endl;

		GDALClose(pSrcDataset);
		CPLFree(pSrcData);
		CPLFree(pDstData);

		//////////////////////////////////////////////////////////////////////////
		// clean temp files
		printf("clean temp fiels... \n");
		if (fs::exists(tmp_data_file))
		{
			DeleteFileA(tmp_data_file.c_str());
		}
		if (fs::exists(tmp_data_file + ".met"))
		{
			DeleteFileA((tmp_data_file + ".met").c_str());
		}

		if (fs::exists(clipTempFile.c_str()))
		{
			DeleteFileA(clipTempFile.c_str());
		}
		std::cerr << "all pixels are null" << std::endl;
		exit(1);
		mean_value = 0.0;
	}
	else
	{
		mean_value /= nCount;
	}
	// new values
	for (int i = 0; i < nX; i++)
	{
		for (int j = 0;j < nY;j++)
		{
			GInt16 d = ((GInt16*)pSrcData)[j * nX + i];
			if (NULL_VALUE != d && d > 0.0)
			{
				pDstData[j * nX + i] = (d - add_offset) * scale_factor;
			}
			else
			{
				pDstData[j * nX + i] = mean_value;
			}
			percent += 1;
			if (int(percent*100.0f/(float)total_size+0.5) == percentList[fill_pos])
			{
				cout<<" "<<percentList[fill_pos]<<"%";
				fill_pos++;
			}	
		}
	}
	if (fill_pos < sizeof(percentList)/sizeof(percentList[0]) - 1)
	{
		cout<<" 100%";
	}
	cout<<endl;

	GDALDatasetH hDstDS = new GDALDatasetH;
	GDALDriverH hDriver = GDALGetDriverByName( "GTiff" );
	hDstDS = GDALCreate( hDriver, outFilename.c_str(), nX, nY, nBand, GDT_Float32, NULL );
	((GDALDataset*)hDstDS)->SetGeoTransform(adfThisGeoTransform);
	((GDALDataset*)hDstDS)->SetProjection(pszTargetSRS);
	((GDALDataset*)hDstDS)->RasterIO(GF_Write, 0, 0, nX, nY, pDstData, nX, nY, GDT_Float32, nBand, &bandList[0], 0, 0, 0);

	GDALClose( pSrcDataset );
	GDALClose( hDstDS );
	CPLFree( pSrcData );
	CPLFree( pDstData );
	
	//////////////////////////////////////////////////////////////////////////
	// clean temp files
	printf( "clean temp fiels... \n" );
	if (fs::exists(tmp_data_file))
	{
		DeleteFileA(tmp_data_file.c_str());
	}
	if (fs::exists(tmp_data_file + ".met"))
	{
		DeleteFileA((tmp_data_file + ".met").c_str());
	}

	if (fs::exists(prm_file))
	{
		DeleteFileA(prm_file.c_str());
	}

	if (fs::exists(batFile))
	{
		DeleteFileA(batFile.c_str());
	}

	if (fs::exists(clipTempFile.c_str()))
	{
		DeleteFileA(clipTempFile.c_str());
	}
}

void hdfgeo(string hdfFile, string refFile, string outFile)
{

	GDALDataset *poDataset = (GDALDataset *)GDALOpen(hdfFile.c_str(), GA_ReadOnly);
	if (!poDataset) {
		std::cerr << "cannot open " << hdfFile << std::endl;
		return;
	}

	//char **metadata = poDataset->GetMetadata("SUBDATASETS");
	//for(int  i = 0; metadata[i] != NULL; i++ )
	//{
	//	cout<<metadata[i]<<endl;
	//}


	// read dimension sizes
	int xsize;
	int ysize;

	GDALDriver * driver = poDataset->GetDriver();
	string papszMetadata = GDALGetDriverShortName((GDALDriverH)poDataset);
	int index = papszMetadata.find_first_of("hdf");
	if (index<0)
		index = papszMetadata.find_first_of("HDF");
	if (!index)
	{
		std::cerr << "invalid MOD05_L2 file. " << std::endl;
		return;
	}

	string swathName = "mod04";

	char ** SUBDATASETS = GDALGetMetadata((GDALDatasetH)poDataset, "SUBDATASETS");
	//if (CSLCount(SUBDATASETS) > 0)
	//{
	//	for (int i = 0; SUBDATASETS[i] != NULL; i++)
	//	{
	//		cout << SUBDATASETS[i] << endl;
	//	}
	//}
	GDALDataset * lat_dt;
	GDALDataset * lon_dt;
	GDALDataset * data_dt;
	string dataset_filename;
	if (CSLCount(SUBDATASETS) > 0)
	{
		for (int i = 0; SUBDATASETS[i] != NULL; i++)
		{
			if (i % 2 == 0)
			{
				string tmpfilename = string(SUBDATASETS[i]);
				tmpfilename = tmpfilename.substr(tmpfilename.find_first_of("=") + 1);
				string DESC = string(SUBDATASETS[i + 1]);
				DESC = DESC.substr(DESC.find_last_of("=") + 1);
				//string tag_str = tmpfilename.substr(tmpfilename.find_last_of(":") + 1);
				//if (0 == tag_str.compare(TAG_LAT))
				//{
				//	lat_dt = (GDALDataset *) GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
				//}
				//if (0 == tag_str.compare(TAG_LON))
				//{
				//	lon_dt = (GDALDataset *) GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
				//}
				vector<string> strListTmp;
				mylib::splitString(DESC, { " " }, strListTmp);
				if (strListTmp.size() < 3)
				{
					cerr << "description of " << tmpfilename << " is not valid." << endl;
					return;
				}
				if (0 == strListTmp[1].compare(TAG_LAT))
				{
					lat_dt = (GDALDataset *)GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
				}
				if (0 == strListTmp[1].compare(TAG_LON))
				{
					lon_dt = (GDALDataset *)GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
				}
				if (0 == strListTmp[1].compare(dimensionName)
					&& tmpfilename.find(dimensionName) != std::string::npos)
				{
					dataset_filename = tmpfilename;
					data_dt = (GDALDataset *)GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
					swathName = strListTmp[2];
				}
			}
		}
	}

	int xsize_lat = lat_dt->GetRasterXSize();
	int ysize_lat = lat_dt->GetRasterYSize();
	int nb_lat = lat_dt->GetRasterCount();
	GDALDataType eDataType_lat = lat_dt->GetRasterBand(1)->GetRasterDataType();
	int datatypeSize_lat = GDALGetDataTypeSize(eDataType_lat) / 8;
	GByte *buffer_lat = new GByte[xsize_lat * ysize_lat * nb_lat * datatypeSize_lat];
	GByte *buffer_lon = new GByte[xsize_lat * ysize_lat * nb_lat * datatypeSize_lat];

	vector<int> bandList;
	for (int iband = 0; iband < nb_lat; ++iband)
	{
		bandList.push_back(iband + 1);
	}

	if (NULL == data_dt)
	{
		cerr << "Dimension \"" << dimensionName << "\" is found." << endl;
		return;
	}
	//get data dataset info
	char ** metaData = data_dt->GetMetadata();
	if (CSLCount(metaData) > 0)
	{
		for (int i = 0; metaData[i] != NULL; i++)
		{
			vector<string> strListTmp;
			mylib::splitString(metaData[i], { "=" }, strListTmp);
			string strTagName = strListTmp[0];
			boost::algorithm::trim_if(strTagName, boost::algorithm::is_any_of(" "));
			string strTagValue = strListTmp[1];
			boost::algorithm::trim_if(strTagValue, boost::algorithm::is_any_of(" "));
			if (0 == strTagName.compare("_FillValue"))
			{
				NULL_VALUE = stoi(strTagValue);
			}
			else if (0 == strTagName.compare("add_offset"))
			{
				add_offset = stod(strTagValue);
			}
			else if (0 == strTagName.compare("scale_factor"))
			{
				scale_factor = stod(strTagValue);
			}
		}
	}

	string tmp_data_file = getTempImageFilename("tif");
	if (fs::exists(tmp_data_file))
	{
		DeleteFileA(tmp_data_file.c_str());
	}

	char cmdBuff[1024];
	sprintf_s(cmdBuff, "gdalwarp -srcnodata \"%d\" -t_srs \"+proj=longlat +datum=WGS84\" -r cubic %s %s\0", NULL_VALUE, dataset_filename.c_str(), tmp_data_file.c_str());
	system(cmdBuff);

	//////////////////////////////////////////////////////////////////////////
	// clip
	printf("cliping and resampling... \n");
	ossimRefPtr<ossimImageHandler> refHandler = ossimImageHandlerFactory::instance()->open(ossimFilename(refFile));
	if (!refHandler)
	{
		std::cerr << "cannot open " << refFile << std::endl;
		return;
	}
	ossimIrect refBoundary = refHandler->getBoundingRect();
	ossimRefPtr<ossimMapProjection> refProjection = PTR_CAST(ossimMapProjection,
		refHandler->getImageGeometry()->getProjection());

	ossimRefPtr<ossimImageHandler> tmpHandler = ossimImageHandlerFactory::instance()->open(ossimFilename(tmp_data_file.c_str()));
	if (!tmpHandler)
	{
		std::cerr << "cannot open " << tmp_data_file.c_str() << std::endl;
		return;
	}
	ossimRefPtr<ossimMapProjection> srcProjection = PTR_CAST(ossimMapProjection,
		tmpHandler->getImageGeometry()->getProjection());

	ossimRefPtr<ossimImageRenderer> renderer = new ossimImageRenderer;
	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	//renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_NEAREST_NEIGHBOR);
	renderer->connectMyInputTo(tmpHandler.get());
	renderer->setView(refProjection.get());

	vector<ossimDpt> polygon;
	ossimDpt dpt;
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ul()), dpt);
	polygon.push_back(dpt);
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ur()), dpt);
	polygon.push_back(dpt);
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.lr()), dpt);
	polygon.push_back(dpt);
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ll()), dpt);
	polygon.push_back(dpt);
	srcProjection->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ul()), dpt);
	polygon.push_back(dpt);
	ossimRefPtr<ossimPolyCutter> theCutter = new ossimPolyCutter;
	theCutter->connectMyInputTo(tmpHandler.get());
	theCutter->setPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, renderer->getImageViewTransform());

	//统计外接矩形
	dpt = renderer->getImageViewTransform()->imageToView(polygon[0]);
	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
	for (int i = 1; i < (int)polygon.size(); i++)
	{
		dpt = renderer->getImageViewTransform()->imageToView(polygon[i]);
		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
	}
	int offset = 0;
	dpt = renderer->getImageViewTransform()->viewToImage(UL1);
	ossimDpt UL(dpt.x - offset, dpt.y - offset);
	dpt = renderer->getImageViewTransform()->viewToImage(LR1);
	ossimDpt LR(dpt.x + offset, dpt.y + offset);
	UL1 = renderer->getImageViewTransform()->imageToView(UL);
	LR1 = renderer->getImageViewTransform()->imageToView(LR);
	//ossimDrect viewRegion = calcRegion(clipRect, inmapinfo, renderer->getImageViewTransform());
	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
	ossimDpt UR = renderer->getImageViewTransform()->viewToImage(ossimDpt(LR1.lon, UL1.lat));
	ossimDpt LL = renderer->getImageViewTransform()->viewToImage(ossimDpt(UL1.lon, LR1.lat));

	vector<ossimDpt> clipRect;
	clipRect.clear();
	clipRect.push_back(UL);
	clipRect.push_back(UR);
	clipRect.push_back(LR);
	clipRect.push_back(LL);
	clipRect.push_back(UL);

	theCutter->setNumberOfPolygons(0);
	theCutter->addPolygon(clipRect);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	theCutter->connectMyInputTo(tmpHandler.get());

	renderer->connectMyInputTo(theCutter.get());

	ossimFilename clipTempFile = getTempImageFilename("tif");

	ossimRefPtr<ossimImageFileWriter> writer = ossimImageWriterFactoryRegistry::instance()->
		createWriterFromExtension(clipTempFile.ext());
	if (writer == NULL) return;

	ossimKeywordlist tt_geom;
	tt_geom.clear();
	writer->saveState(tt_geom);
	tt_geom.add("", ossimKeywordNames::PIXEL_TYPE_KW, "PIXEL_IS_AREA", true);
	writer->loadState(tt_geom);


	ossimStdOutProgress* progress = new ossimStdOutProgress(0, true);
	writer->addListener(progress);

	writer->setFilename(clipTempFile);
	writer->connectMyInputTo(0, renderer.get());
	writer->setAreaOfInterest(viewRegion);
	//ossimIrect bounding = handler->getBoundingRect();
	//writer->setAreaOfInterest(bounding);

	//bool bResult = writer->execute();

	if (!writer->execute())
	{
		writer->removeListener(progress);
		std::cerr << "geometric correction failed." << std::endl;
		return;
	}
	refHandler->close();
	tmpHandler->close();

	//printf( "it's ok! please enter any key! \n" );

	//////////////////////////////////////////////////////////////////////////
	printf("fill null pixels... \n");

	string outFilename;
	if ("" == outFile)
	{
		outFilename = getOutName(refFile);
	}
	else
	{
		outFilename = outFile;
	}
	GDALDataset *pSrcDataset;
	pSrcDataset = (GDALDataset *)GDALOpen(clipTempFile, GA_ReadOnly);
	double adfThisGeoTransform[6];
	pSrcDataset->GetGeoTransform(adfThisGeoTransform);
	const char* pszTargetSRS = pSrcDataset->GetProjectionRef();
	int nX = pSrcDataset->GetRasterXSize();
	int nY = pSrcDataset->GetRasterYSize();
	int nBand = pSrcDataset->GetRasterCount();
	GDALDataType eDT = pSrcDataset->GetRasterBand(1)->GetRasterDataType();
	int dataTypeSize = GDALGetDataTypeSize(eDT) / 8;
	if (1 != nBand)
	{
		std::cerr << "water vapor file should have 1 band." << std::endl;
		return;
	}
	bandList.clear();
	for (int iband = 0; iband < nBand; ++iband)
	{
		bandList.push_back(iband + 1);
	}
	GByte* pSrcData = new GByte[nX * nY * nBand * dataTypeSize];
	float* pDstData = new float[nX * nY * nBand];
	int nCount = 0;
	double mean_value = 0.0;
	pSrcDataset->RasterIO(GF_Read, 0, 0, nX, nY, pSrcData, nX, nY, eDT, nBand, &bandList[0], 0, 0, 0);

	int fill_pos = 0;
	cout << percentList[fill_pos++] << "%";
	int percent = 0;
	float total_size = nX * nY * 2;
	for (int i = 0; i < nX; i++)
	{
		for (int j = 0; j < nY; j++)
		{
			GInt16 d = ((GInt16*)pSrcData)[j * nX + i];
			if (NULL_VALUE != d && d > 0.0)
			{
				mean_value += (d - add_offset) * scale_factor;
				nCount++;
			}
			percent += 1;
			if (int(percent*100.0f / (float)total_size + 0.5) == percentList[fill_pos])
			{
				cout << " " << percentList[fill_pos] << "%";
				fill_pos++;
			}
		}
	}
	if (0 == nCount)
	{
		// all pixels are null 
		cout << " 100%";
		cout << endl;

		GDALClose(pSrcDataset);
		CPLFree(pSrcData);
		CPLFree(pDstData);

		//////////////////////////////////////////////////////////////////////////
		// clean temp files
		printf("clean temp fiels... \n");
		if (fs::exists(tmp_data_file))
		{
			DeleteFileA(tmp_data_file.c_str());
		}
		if (fs::exists(tmp_data_file + ".met"))
		{
			DeleteFileA((tmp_data_file + ".met").c_str());
		}

		if (fs::exists(clipTempFile.c_str()))
		{
			DeleteFileA(clipTempFile.c_str());
		}
		throw "all pixels are null ";
		return;
		//goto exit_error;
		//mean_value = 0.0;
	}
	else
	{
		mean_value /= nCount;
	}
	// new values
	for (int i = 0; i < nX; i++)
	{
		for (int j = 0; j < nY; j++)
		{
			GInt16 d = ((GInt16*)pSrcData)[j * nX + i];
			if (NULL_VALUE != d && d > 0.0)
			{
				pDstData[j * nX + i] = (d - add_offset) * scale_factor;
			}
			else
			{
				pDstData[j * nX + i] = mean_value;
			}
			percent += 1;
			if (int(percent*100.0f / (float)total_size + 0.5) == percentList[fill_pos])
			{
				cout << " " << percentList[fill_pos] << "%";
				fill_pos++;
			}
		}
	}
	if (fill_pos < sizeof(percentList) / sizeof(percentList[0]) - 1)
	{
		cout << " 100%";
	}
	cout << endl;

	GDALDatasetH hDstDS = new GDALDatasetH;
	GDALDriverH hDriver = GDALGetDriverByName("GTiff");
	hDstDS = GDALCreate(hDriver, outFilename.c_str(), nX, nY, nBand, GDT_Float32, NULL);
	((GDALDataset*)hDstDS)->SetGeoTransform(adfThisGeoTransform);
	((GDALDataset*)hDstDS)->SetProjection(pszTargetSRS);
	((GDALDataset*)hDstDS)->RasterIO(GF_Write, 0, 0, nX, nY, pDstData, nX, nY, GDT_Float32, nBand, &bandList[0], 0, 0, 0);

	GDALClose(hDstDS);
exit_error:
	GDALClose(pSrcDataset);
	CPLFree(pSrcData);
	CPLFree(pDstData);

	//////////////////////////////////////////////////////////////////////////
	// clean temp files
	printf("clean temp fiels... \n");
	if (fs::exists(tmp_data_file))
	{
		DeleteFileA(tmp_data_file.c_str());
	}
	if (fs::exists(tmp_data_file + ".met"))
	{
		DeleteFileA((tmp_data_file + ".met").c_str());
	}

	if (fs::exists(clipTempFile.c_str()))
	{
		DeleteFileA(clipTempFile.c_str());
	}
}

//void hdfgeo(const char* hdfFile, const char* refFile, const char* outFile)
//{
//	GDALDataset *poDataset = (GDALDataset *)GDALOpen(hdfFile, GA_ReadOnly);
//	if (!poDataset) {
//		std::cerr << "cannot open " << hdfFile << std::endl;
//		return;
//	}
//
//	//char **metadata = poDataset->GetMetadata("SUBDATASETS");
//	//for(int  i = 0; metadata[i] != NULL; i++ )
//	//{
//	//	cout<<metadata[i]<<endl;
//	//}
//
//	string tmp_vapor_file = getTempName(hdfFile);
//
//	// read dimension sizes
//	int xsize;
//	int ysize;
//
//	GDALDriver * driver = poDataset->GetDriver ();
//	string papszMetadata = GDALGetDriverShortName((GDALDriverH)poDataset);
//	int index=papszMetadata.find_first_of("hdf");
//	if(index<0)
//		index=papszMetadata.find_first_of("HDF");
//	if(!index)
//	{
//		std::cerr << "invalid MOD05_L2 file. " << std::endl;
//		return;
//	}
//
//
//	char ** SUBDATASETS = GDALGetMetadata( (GDALDatasetH)poDataset, "SUBDATASETS" );
//	if( CSLCount(SUBDATASETS) > 0 )
//	{
//		//printf( "Subdatasets:\n" );
//		for(int  i = 0; SUBDATASETS[i] != NULL; i++ )
//		{
//			if(i%2==0)
//			{
//				string tmpfilename = string(SUBDATASETS[i]);
//				tmpfilename = tmpfilename.substr(tmpfilename.find_first_of("=") + 1);
//				string tag_str = tmpfilename.substr(tmpfilename.find_last_of(":") + 1);
//				if (0 == tag_str.compare(dimensionName))
//				{
//					GDALDataset * tmpdt = (GDALDataset *) GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
//					xsize = tmpdt->GetRasterXSize();
//					ysize = tmpdt->GetRasterYSize();
//					int nb = tmpdt->GetRasterCount();
//					GDALDataType eDataType = tmpdt->GetRasterBand(1)->GetRasterDataType();
//					int datatypeSize = GDALGetDataTypeSize(eDataType) / 8;
//					GByte *buffer = new GByte [xsize * ysize * nb * datatypeSize];
//
//					vector<int> bandList;
//					for (int iband = 0; iband < nb;++iband)
//					{
//						bandList.push_back(iband+1);
//					}
//					tmpdt->RasterIO(GF_Read, 0, 0, xsize, ysize, buffer, xsize, ysize, eDataType, nb, &bandList[0], 0, 0, 0);
//
//					GDALDatasetH hDstDS = new GDALDatasetH;
//					const char         *pszFormat = "GTiff";
//					//char               *pszTargetSRS = NULL;
//					//double adfThisGeoTransform[6];
//					GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
//					hDstDS = GDALCreate( hDriver, tmp_vapor_file.c_str(), xsize, ysize, nb, eDataType, NULL );
//					((GDALDataset*)hDstDS)->RasterIO(GF_Write, 0, 0, xsize, ysize, buffer, xsize, ysize, eDataType, nb, &bandList[0], 0, 0, 0);
//					GDALClose( hDstDS );
//					GDALClose( tmpdt );
//					CPLFree( buffer );
//				}
//			}            
//		}
//	}
//
//	//for(int y = 0; y < ySamples; ++y)
//	//{
//	//	for(int x = 0; x < xSamples; ++x)
//	//	{
//	//		ossimDpt imagePoint;
//	//		if(ySamples > 1)
//	//		{
//	//			ynorm = (double)y/(double)(ySamples - 1);
//	//		}
//	//		else
//	//		{
//	//			ynorm = 0.0;
//	//		}
//	//		if(xSamples > 1)
//	//		{
//	//			xnorm = (double)x/(double)(xSamples - 1);
//	//		}
//	//		else
//	//		{
//	//			xnorm = 0.0;
//	//		}
//
//	//		ossimDpt dpt(xsize*xnorm, ysize*ynorm);
//	//		char tmpStr[256];
//	//		sprintf(tmpStr, "%d", y * xSamples + x + 1);
//	//		ossimGpt gpt;
//	//		ossimTieGpt *aTiePt = new ossimTieGpt(gpt, dpt, 1.0, ossimString(tmpStr));
//	//		pTieGptSet->addTiePoint(aTiePt);
//	//	}
//	//}
//
//	GDALDataset * lat_dt;
//	GDALDataset * lon_dt;
//	if( CSLCount(SUBDATASETS) > 0 )
//	{
//		for(int  i = 0; SUBDATASETS[i] != NULL; i++ )
//		{
//			if(i%2==0)
//			{
//				string tmpfilename = string(SUBDATASETS[i]);
//				tmpfilename = tmpfilename.substr(tmpfilename.find_first_of("=") + 1);
//				string tag_str = tmpfilename.substr(tmpfilename.find_last_of(":") + 1);				
//				if (0 == tag_str.compare(TAG_LAT))
//				{
//					lat_dt = (GDALDataset *) GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
//				}
//				if (0 == tag_str.compare(TAG_LON))
//				{
//					lon_dt = (GDALDataset *) GDALOpen(tmpfilename.c_str(), GA_ReadOnly);
//				}
//			}            
//		}
//	}
//
//	int xsize_lat = lat_dt->GetRasterXSize();
//	int ysize_lat = lat_dt->GetRasterYSize();
//	int nb_lat = lat_dt->GetRasterCount();
//	GDALDataType eDataType_lat = lat_dt->GetRasterBand(1)->GetRasterDataType();
//	int datatypeSize_lat = GDALGetDataTypeSize(eDataType_lat) / 8;
//	GByte *buffer_lat = new GByte [xsize_lat * ysize_lat * nb_lat * datatypeSize_lat];
//	GByte *buffer_lon = new GByte [xsize_lat * ysize_lat * nb_lat * datatypeSize_lat];
//
//	vector<int> bandList;
//	for (int iband = 0; iband < nb_lat;++iband)
//	{
//		bandList.push_back(iband+1);
//	}
//	lat_dt->RasterIO(GF_Read, 0, 0, xsize_lat, ysize_lat, buffer_lat, xsize_lat, ysize_lat, eDataType_lat, nb_lat, &bandList[0], 0, 0, 0);
//	lon_dt->RasterIO(GF_Read, 0, 0, xsize_lat, ysize_lat, buffer_lon, xsize_lat, ysize_lat, eDataType_lat, nb_lat, &bandList[0], 0, 0, 0);
//
//	ossim_uint32 xSamples = 50;
//	ossim_uint32 ySamples = 50;
//
//	ossimTieGptSet *pTieGptSet = new ossimTieGptSet;
//	double ynorm;
//	double xnorm;
//	for(int y = 0; y < ySamples; ++y)
//	{
//		for(int x = 0; x < xSamples; ++x)
//		{
//			ossimDpt imagePoint;
//			if(ySamples > 1)
//			{
//				ynorm = (double)y/(double)(ySamples - 1);
//			}
//			else
//			{
//				ynorm = 0.0;
//			}
//			if(xSamples > 1)
//			{
//				xnorm = (double)x/(double)(xSamples - 1);
//			}
//			else
//			{
//				xnorm = 0.0;
//			}
//
//			ossimDpt dpt(xsize*xnorm, ysize*ynorm);
//			int ix = int(xsize*xnorm*0.2 + 0.5);
//			int iy = int(ysize*ynorm*0.2 + 0.5);
//			ix = ix < 0? 0 : ix;
//			iy = iy < 0? 0 : iy;
//			ix = ix > (xsize-1)? (xsize-1) : ix;
//			iy = iy > (ysize-1)? (ysize-1) : iy;
//			ossimGpt gpt(((float*)buffer_lat)[iy * xsize_lat + ix], ((float*)buffer_lon)[iy * xsize_lat + ix], 0.0);
//			char tmpStr[256];
//			sprintf(tmpStr, "%d", y * xSamples + x + 1);
//			ossimTieGpt *aTiePt = new ossimTieGpt(gpt, dpt, 1.0, ossimString(tmpStr));
//			pTieGptSet->addTiePoint(aTiePt);
//		}
//	}
//
//	GDALClose( lat_dt );
//	GDALClose( lon_dt );
//	CPLFree( buffer_lat );
//	CPLFree( buffer_lon );
//	GDALClose( poDataset );
//
//	ossimRefPtr<ossimImageHandler> refHandler = ossimImageHandlerFactory::instance()->open(ossimFilename(refFile));
//	if (!refHandler)
//	{
//		std::cerr << "cannot open " << refFile << std::endl;
//		return;
//	}
//	ossimIrect refBoundary = refHandler->getBoundingRect();
//	ossimRefPtr<ossimMapProjection> refProjection = PTR_CAST(ossimMapProjection,
//		refHandler->getImageGeometry()->getProjection());
//
//	ossimRefPtr<ossimPolynomProjection> sensorModel = new ossimPolynomProjection;
//	ossimString polyDegree="1 x y x2 xy y2";
//	//ossimString polyDegree="1 x y x2 xy y2 x3 y3 xy2 x2y";
//	sensorModel->setupOptimizer(polyDegree);
//	sensorModel->optimizeFit(*pTieGptSet);
//
//	mylib::saveGcpFile("1.txt", pTieGptSet);
//
//	ossimRefPtr<ossimImageHandler> tmpHandler = ossimImageHandlerFactory::instance()->open(ossimFilename(tmp_vapor_file.c_str()));
//	if (!tmpHandler)
//	{
//		std::cerr << "cannot open " << tmp_vapor_file.c_str() << std::endl;
//		return;
//	}
//	ossimRefPtr<ossimImageGeometry> imageGeom = new ossimImageGeometry;
//	imageGeom->setProjection(sensorModel.get());
//	tmpHandler->setImageGeometry(imageGeom.get());//1128
//
//
//	ossimRefPtr<ossimImageRenderer> renderer = new ossimImageRenderer;
//	renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_NEAREST_NEIGHBOR);
//	renderer->connectMyInputTo(tmpHandler.get());
//	renderer->setView(refProjection.get());
//
//	vector<ossimDpt> polygon;
//	ossimDpt dpt;
//	sensorModel->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ul()), dpt);
//	polygon.push_back(dpt);
//	sensorModel->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ur()), dpt);
//	polygon.push_back(dpt);
//	sensorModel->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.lr()), dpt);
//	polygon.push_back(dpt);
//	sensorModel->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ll()), dpt);
//	polygon.push_back(dpt);
//	sensorModel->worldToLineSample(refProjection->lineSampleToWorld(refBoundary.ul()), dpt);
//	polygon.push_back(dpt);
//	ossimRefPtr<ossimPolyCutter> theCutter = new ossimPolyCutter;
//	theCutter->connectMyInputTo(tmpHandler.get());
//	theCutter->setPolygon(polygon);
//	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
//	theCutter->setNumberOfPolygons(1);
//
//	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, renderer->getImageViewTransform());
//
//	//统计外接矩形
//	dpt = renderer->getImageViewTransform()->imageToView(polygon[0]);
//	ossimDpt UL1(dpt.x, dpt.y), LR1(dpt.x, dpt.y);
//	for (int i = 1; i < (int) polygon.size(); i++)
//	{
//		dpt = renderer->getImageViewTransform()->imageToView(polygon[i]);
//		UL1.x = (dpt.x < UL1.x) ? dpt.x : UL1.x;
//		UL1.y = (dpt.y < UL1.y) ? dpt.y : UL1.y;
//		LR1.x = (dpt.x > LR1.x) ? dpt.x : LR1.x;
//		LR1.y = (dpt.y > LR1.y) ? dpt.y : LR1.y;
//	}
//	int offset = 0;
//	dpt = renderer->getImageViewTransform()->viewToImage(UL1);
//	ossimDpt UL(dpt.x - offset, dpt.y - offset);
//	dpt = renderer->getImageViewTransform()->viewToImage(LR1);
//	ossimDpt LR(dpt.x + offset, dpt.y + offset);
//	UL1 = renderer->getImageViewTransform()->imageToView(UL);
//	LR1 = renderer->getImageViewTransform()->imageToView(LR);
//	//ossimDrect viewRegion = calcRegion(clipRect, inmapinfo, renderer->getImageViewTransform());
//	ossimDrect viewRegion(UL1.x, UL1.y, LR1.x, LR1.y);
//	ossimDpt UR = renderer->getImageViewTransform()->viewToImage(ossimDpt(LR1.lon, UL1.lat));
//	ossimDpt LL = renderer->getImageViewTransform()->viewToImage(ossimDpt(UL1.lon, LR1.lat));
//
//	vector<ossimDpt> clipRect;
//	clipRect.clear();
//	clipRect.push_back(UL);
//	clipRect.push_back(UR);
//	clipRect.push_back(LR);
//	clipRect.push_back(LL);
//	clipRect.push_back(UL);
//
//	theCutter->setNumberOfPolygons(0);
//	theCutter->addPolygon(clipRect);
//	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
//	theCutter->setNumberOfPolygons(1);
//
//	theCutter->connectMyInputTo(tmpHandler.get());
//
//	renderer->connectMyInputTo(theCutter.get());
//
//	ossimFilename outFilename;
//	if ("" == outFile)
//	{
//		outFilename = ossimFilename(getOutName(hdfFile).c_str());
//	}
//	else
//	{
//		outFilename = ossimFilename(outFile);
//	}
//
//	ossimRefPtr<ossimImageFileWriter> writer = ossimImageWriterFactoryRegistry::instance()->
//		createWriterFromExtension( outFilename.ext() );
//	if(writer==NULL) return;
//
//	ossimKeywordlist tt_geom;
//	tt_geom.clear();
//	writer->saveState(tt_geom);
//	tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA",true);
//	writer->loadState(tt_geom);
//
//
//	ossimStdOutProgress* progress = new ossimStdOutProgress(0,true);
//	writer->addListener(progress);
//
//	writer->setFilename(outFilename);
//	writer->connectMyInputTo(0, renderer.get());
//	writer->setAreaOfInterest(viewRegion);
//	//ossimIrect bounding = handler->getBoundingRect();
//	//writer->setAreaOfInterest(bounding);
//
//	//bool bResult = writer->execute();
//
//	if(!writer->execute())
//	{
//		writer->removeListener(progress);
//		std::cerr << "geometric correction failed." << std::endl;
//		return;
//	}
//
//	printf( "it's ok! please enter any key! \n" );
//
//	//char **metadata = poDataset->GetMetadata("");
//
//	//// read dimension sizes
//	//int xsize = poDataset->GetRasterXSize();
//	//int ysize = poDataset->GetRasterYSize();
//	//std::cout << "X size = " << xsize << ", " << "Y size = " << ysize << std::endl;
//
//	//// read elements
//	//float *buffer = new float[xsize * ysize];
//	//GDALRasterBand *rb = poDataset->GetRasterBand(1);
//	//GDALDataType eDataType = rb->GetRasterDataType();
//	//rb->RasterIO(GF_Read, 0, 0, xsize, ysize, buffer, xsize, ysize, GDT_Float32, 0, 0);
//	//for (int j = 0; j < ysize; ++j) {
//	//	for (int k = 0; k < xsize; ++k) {
//	//		std::cout << buffer[j * xsize + k] << " ";
//	//	}
//	//	std::cout << std::endl;
//	//}
//	//delete [] buffer;
//}


/*************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage(const char* pszErrorMsg = NULL)

{
	printf( 
		"Usage: hdfgeo  [-h hdffile] [-r reffile] \n"
		"\t[-o outfile] [-heg hegpath] \n"
		"  -h hdffile\t: hdf file\n"
		"  -r reffile\t: reference file\n"
		"  -o outfile\t: output file\n"
		"  -heg hegpath\t: heg path\n");

	if( pszErrorMsg != NULL )
		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
	exit(1);
}

#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
	do { if (i + nExtraArg >= argc) \
	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)

int main( int argc, char** argv )
{
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持

	//ossimInit::instance()->loadPlugins("ossimgdal_plugin.dll");

	if (argc == 1)
	{
	}

	if (argc > 1)
	{
		/* -------------------------------------------------------------------- */
		/*      Parse arguments.                                                */
		/* -------------------------------------------------------------------- */
		for( int i = 1; i < argc; i++ )
		{
			if( 0 == _stricmp(argv[i],"-h") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszHdfFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-r") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszRefFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-o") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszOutFile = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-tag") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				dimensionName = argv[++i];
			}
			else if( 0 == _stricmp(argv[i],"-heg") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszHegPath = argv[++i] ;
			}
			else if( 0 == _stricmp(argv[i],"-pref") )
			{
				CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
				pszPreferenceFile = argv[++i] ;
			}
			else
			{
				Usage();
			}
		}

		if ("" == pszHdfFile || "" == pszRefFile)
		{
			Usage();
		}
		//QDir(pszHdfFile).absoluteFilePath(pszHdfFile).toLatin1();
		//namespace fs = boost::filesystem;
		//fs::path full_p = fs::complete(fs::path(pszHdfFile)); // complete == absolute
		//const char* hdfFile = (const char *)full_p.c_str();
		//full_p = fs::complete(fs::path(pszRefFile)); // complete == absolute
		//const char* reffFile = (const char *)full_p.c_str();

		//	if (!fs::exists(strInputFile))
		//	{
		//		printf("inputFile does not exist!\n");
		//	}
		//	else
		//	{
		//		fs::path inputPathname = fs::absolute(strInputFile);

		//string hdfFile = QDir::toNativeSeparators(QDir(pszHdfFile).absolutePath()).toStdString();
		//string reffFile = QDir::toNativeSeparators(QDir(pszRefFile).absolutePath()).toStdString();
		string hdfFile = fs::absolute(string(pszHdfFile)).string();
		string reffFile = fs::absolute(string(pszRefFile)).string();

		ossimFilename preferences_file = ossimFilename(pszPreferenceFile);
		if (preferences_file.exists())
		{
			ossimPreferences::instance()->loadPreferences(preferences_file);
		}
		else
		{
			preferences_file = ossimFilename(getenv("OSSIM2_DIR")) + "\\preference.txt";
			if (preferences_file.exists())
			{
				ossimPreferences::instance()->loadPreferences(preferences_file);
			}
		}
		ossimInit::instance()->initialize();
		//hdfgeo(hdfFile, reffFile, pszOutFile);
		create_prm(hdfFile, reffFile, pszOutFile);
		//batchReprojection("", QStringList(pszFilter), pszOutPath, pszProjectFile);
	}
	else
	{
		pszHdfFile = "G:\\modis\\2014\\MOD04_L2\\LC81130262014257BJC00_B1_MOD04_L2.A2014257.0115.051.2014273103008.hdf";
		pszRefFile = "G:\\modis\\2014\\LC81130262014257BJC00_B4.TIF";
		dimensionName = "Image_Optical_Depth_Land_And_Ocean";
		pszOutFile = "G:\\modis\\2014\\tif\\LC81130262014257BJC00_B1_MOD04_L2.A2014257.0115.051.2014273103008.tif";
		//pszHdfFile = "E:\\modis\\1\\MOD05_L2.A2013244.0220.051.2013245223924.hdf";
		//pszRefFile = "E:\\modis\\1\\LC81230332013244LGN00_B4.TIF";
		////hdfgeo(pszHdfFile, pszRefFile);
		string hdfFile = fs::absolute(string(pszHdfFile)).string();
		string reffFile = fs::absolute(string(pszRefFile)).string();
		ossimInit::instance()->initialize();
		create_prm(hdfFile, reffFile, "");
		//hdfgeo(hdfFile, reffFile, pszOutFile);
		Usage();
	}

	return 0;
}