//#include <wx\iconbndl.h>
#include <ogrsf_frmts.h>
#include <gdal.h>
#include <ogr_api.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>
#include <string>
#include "ShpApp.h"
#include "GdalRasterApp.h"
using namespace std;

#pragma comment(lib, "opencv_ts300.lib")
#pragma comment(lib, "opencv_world300.lib")

bool Shp2Map(string srcfile, string dstfile)
{
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动

	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");//得到shp文件的处理器
	//OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("KML");//得到kml文件的处理器
	//OGRDataSource* poDS = poDriver->Open( srcfile, NULL );//打开文件

	OGRDataSource* poDS = poDriver->Open( srcfile.c_str(), NULL );//打开文件

	//int iLayerCount = poDS->GetLayerCount();

	OGRLayer* poLayer = poDS->GetLayer(0);//获取shp图层

	//读取几何和属性值
	OGRFeature *poFeature;
	poLayer->ResetReading();//确保是从该层的开头开始
	FILE *pf;
	pf = fopen(dstfile.c_str(),"w+");
	while ((poFeature=poLayer->GetNextFeature())!=NULL)
	{
		OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
		int iField;
		for( iField = 0; iField < poFDefn->GetFieldCount(); iField++ )
		{
			OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
			printf( "%s=", poFieldDefn->GetNameRef() );
			fprintf(pf, "%s=", poFieldDefn->GetNameRef() );
			if( poFieldDefn->GetType() == OFTInteger )
			{
				printf( "%d\n", poFeature->GetFieldAsInteger( iField ) );
				fprintf(pf, "%d\n", poFeature->GetFieldAsInteger( iField ) );
			}
			else if( poFieldDefn->GetType() == OFTReal )
			{
				printf( "%.3f\n", poFeature->GetFieldAsDouble(iField) );
				fprintf(pf, "%.3f\n", poFeature->GetFieldAsDouble(iField) );
			}
			else if( poFieldDefn->GetType() == OFTString )
			{
				printf( "%s\n", poFeature->GetFieldAsString(iField) );
				fprintf(pf, "%s\n", poFeature->GetFieldAsString(iField) );
			}
			else
			{
				printf( "%s\n", poFeature->GetFieldAsString(iField) );
				fprintf(pf, "%s\n", poFeature->GetFieldAsString(iField) );
			}
		}

		OGRGeometry *poGeometry;
		poGeometry = poFeature->GetGeometryRef();
		OGRwkbGeometryType geoType = poGeometry->getGeometryType();
		if( poGeometry != NULL 
			&& wkbFlatten(poGeometry->getGeometryType()) == wkbPoint )//点文件
		{
			OGRPoint *poPoint = (OGRPoint *) poGeometry;
			printf( "%.3f,%.3f\n", poPoint->getX(), poPoint->getY() );
			//fprintf(pf, "x=%lf,y=%lf\n",poPoint->getX(), poPoint->getY() );
			fprintf(pf, "%lf\t%lf\n",poPoint->getX(), poPoint->getY());
		}
		else if( poGeometry != NULL 
			&& wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon )//多边形文件
		{
			//OGRPolygon *poPolygon = (OGRPolygon *)poGeometry;//将该指针转化为OGRPolygon几何对象类型的指针
			//备份
			unsigned char *pData;//unsigned char*类型指针，用于指向数据缓冲区         
			pData = new unsigned char[poGeometry->WkbSize()];
			
			printf("WkbSize:%d\n",poGeometry->WkbSize());
			fprintf(pf, "WkbSize:%d\n",poGeometry->WkbSize());
			
			poGeometry->exportToWkb(wkbNDR ,pData);
			/*for(int i=0;i<poGeometry->WkbSize();i+=2)
			{
				printf("%d,%d\n",pData[i],pData[i+1]);
			}*/
			int pos = 0;
			char byteOrder;
			int wkbType;
			int num_wkbFeature;
			memcpy(&byteOrder, pData + pos, sizeof(byteOrder));
			pos += sizeof(byteOrder);
			memcpy(&wkbType, pData + pos, sizeof(wkbType));
			pos += sizeof(wkbType);
			memcpy(&num_wkbFeature, pData+pos, sizeof(num_wkbFeature));
			pos += sizeof(num_wkbFeature);
			printf("byteOrder=%d\nwkbType=%d\nnum_wkbFeature=%d\n",byteOrder,wkbType,num_wkbFeature);
			fprintf(pf, "byteOrder=%d\nwkbType=%d\nnum_wkbFeature=%d\n",byteOrder,wkbType,num_wkbFeature);
			for(int iFeature = 0;iFeature < num_wkbFeature;iFeature++)
			{
				int nPoint;
				memcpy(&nPoint, pData + pos, sizeof(nPoint));
				pos += sizeof(nPoint);
				printf("nPoint=%d\n",nPoint);
				fprintf(pf, "nPoint=%d\n",nPoint);
				for(int iPoint = 0;iPoint < nPoint;iPoint++)
				{
					double x;
					double y;
					memcpy(&x, pData + pos, sizeof(x));
					pos += sizeof(x);
					memcpy(&y, pData + pos, sizeof(y));
					pos += sizeof(y);
					printf("x=%lf,y=%lf\n",x,y);
					//fprintf(pf, "x=%lf,y=%lf\n",x,y);
					fprintf(pf, "%lf\t%lf\n",x,y);
				}
			}
			OGRFeature::DestroyFeature( poFeature );		
		}
		else if( poGeometry != NULL 
			&& wkbFlatten(poGeometry->getGeometryType()) == wkbLineString )//线文件
		{
			unsigned char *pData;//unsigned char*类型指针，用于指向数据缓冲区         
			pData = new unsigned char[poGeometry->WkbSize()];

			printf("WkbSize:%d\n",poGeometry->WkbSize());
			fprintf(pf, "WkbSize:%d\n",poGeometry->WkbSize());

			poGeometry->exportToWkb(wkbNDR ,pData);

			int pos = 0;
			char byteOrder;
			int wkbType;
			int num_wkbFeature;
			int nPoint;
			memcpy(&byteOrder, pData + pos, sizeof(byteOrder));
			pos += sizeof(byteOrder);
			memcpy(&wkbType, pData + pos, sizeof(wkbType));
			pos += sizeof(wkbType);
			memcpy(&nPoint, pData+pos, sizeof(nPoint));
			pos += sizeof(nPoint);
			printf("byteOrder=%d\nwkbType=%d\nPoint=%d\n",byteOrder,wkbType,nPoint);
			fprintf(pf, "byteOrder=%d\nwkbType=%d\nPoint=%d\n",byteOrder,wkbType,nPoint);

			for(int iPoint = 0;iPoint < nPoint;iPoint++)
			{
				double x;
				double y;
				memcpy(&x, pData + pos, sizeof(x));
				pos += sizeof(x);
				memcpy(&y, pData + pos, sizeof(y));
				pos += sizeof(y);
				printf("x=%lf,y=%lf\n",x,y);
				//fprintf(pf, "x=%lf,y=%lf\n",x,y);
				fprintf(pf, "%lf\t%lf\n",x,y);
			}
			OGRFeature::DestroyFeature( poFeature );
		}
		else
		{
			printf( "no geometry\n" );
		}
		printf("\n");
		fprintf(pf, "\n");
	}
	fclose(pf);
	OGRDataSource::DestroyDataSource( poDS );
	OGRCleanupAll();//资源清理
	return TRUE;
	
	//GDAL写shp文件
	//GDALAllRegister();
	//OGRRegisterAll();//注册所有的文件格式驱动	

	//OGRSFDriver* poDriver = Registrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");//得到shp文件的处理器
	//OGRDataSource* poDS = poDriver->CreateDataSource( "D:\\lakes.shp", NULL );//创建shp文件
	//OGRLayer* poLayer= poDS->CreateLayer( "tbLine", NULL, wkbLineString, NULL );//创建图层

	////创建字段
	//OGRFieldDefn oField1("GeoObjNum",OFTString);
	//oField1.SetWidth(8);// 字符串
	//if( poLayer->CreateField( &oField1 ) != OGRERR_NONE )
	//{
	//	AfxMessageBox( "Creating Name field failed.\n" );
	//	return FALSE;
	//}
	//  
	//OGRFieldDefn oField2("LBTG",OFTReal);
	//oField2.SetPrecision(3);// 浮点数
	//if( poLayer->CreateField( &oField2 ) != OGRERR_NONE )
	//{
	//	AfxMessageBox( "Creating Name field failed.\n" );
	//	return FALSE;
	//}
	//
	//OGRFieldDefn oField3("Number",OFTInteger);// 整型14
	//if( poLayer->CreateField( &oField3 ) != OGRERR_NONE )
	//{
	//	AfxMessageBox( "Creating Name field failed.\n" );
	//	return FALSE;
	//}

	////创建几何和Feature
	//OGRFeature *poFeature;
	//poFeature = new OGRFeature( poLayer->GetLayerDefn() );
	//poFeature->SetField( "GeoObjNum", strGeoObjNum );
	//poFeature->SetField( "LBTG", fLBTG );
	//poFeature->SetField( "Number", number );
	//OGRLineString *poLine = new OGRLineString();
	//poLine->setNumPoints(2);11 poLine->setPoint(0,startX,startY, 0.0);12 poLine->setPoint(1,endX,endY, 0.0);
	//poFeature->SetGeometryDirectly( poLine );
	//if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
	//{
	//	AfxMessageBox("Failed to create feature in shapefile.");
	//	return FALSE;
	//}

	//OGRDataSource::DestroyDataSource( poDS );
	//OGRCleanupAll();//资源清理
}

void appendStraightline2File(string straightlineFile, string imagePtfile, string groundPtfile, string imageFile)
{
	FILE *inpf;
	fopen_s(&inpf, straightlineFile.c_str(), "r+");

	FILE *outImagepf, *outGroundpf;
	fopen_s(&outImagepf, imagePtfile.c_str(), "a+");
	fopen_s(&outGroundpf, groundPtfile.c_str(), "a+");

	GdalRasterApp gdalApp;
	gdalApp.open(imageFile.c_str());

	int id;
	double ix1, iy1, gx1, gy1, gh1;
	double ix2, iy2, gx2, gy2, gh2;
	while (
		EOF != fscanf_s(inpf, "%d%lf%lf%lf%lf%lf", &id, &ix1, &iy1, &gx1, &gy1, &gh1)
		&& EOF != fscanf_s(inpf, "%d%lf%lf%lf%lf%lf", &id, &ix2, &iy2, &gx2, &gy2, &gh2)
		//EOF != fscanf_s(inpf, "%d%lf%lf%lf%lf%lf", &id, &gx1, &gy1, &ix1, &iy1, &gh1)
		//&& EOF != fscanf_s(inpf, "%d%lf%lf%lf%lf%lf", &id, &gx2, &gy2, &ix2, &iy2, &gh2)
		)
	{
		fprintf_s(outImagepf, "%d%03d\t%s\t%d\n", 2, id, "StraightLine", 2);
		fprintf_s(outImagepf, "%lf\t%lf\n", ix1, iy1);
		fprintf_s(outImagepf, "%lf\t%lf\n", ix2, iy2);

		fprintf_s(outGroundpf, "%d%03d\t%s\t%d\n", 2, id, "StraightLine", 2);
		fPoint world;
		gdalApp.linesample2world(fPoint(gx1, gy1), world);
		fprintf_s(outGroundpf, "%lf\t%lf\n", world.x, world.y);

		gdalApp.linesample2world(fPoint(gx2, gy2), world);
		fprintf_s(outGroundpf, "%lf\t%lf\n", world.x, world.y);

	}
	//StraightLine
	fclose(inpf);
	fclose(outImagepf);
	fclose(outGroundpf);
}

void jading()
{
	//wxString srcfile = "D:\\workspace\\testdata\\shp\\china\\china_100wan.shp";
	//string srcfile = "D:\\workspace\\Landsat\\areal features\\freeline_reference.shp";
	//wxString dstfile = "D:\\workspace\\testdata\\shp\\china\\parse.txt";
	//string dstfile = "D:\\workspace\\Landsat\\areal features\\areal_source.txt";
	//wxString srcfile = "C:\\Users\\Administrator\\Desktop\\china\\china_100wan.shp";
	//wxString dstfile = "C:\\Users\\Administrator\\Desktop\\china\\result";
	string outRefFile = "D:\\workspace\\jiading\\features\\reference.txt";
	string outSrcFile = "D:\\workspace\\jiading\\features\\source.txt";
	string fold = "D:\\workspace\\jiading\\features";
	vector<string> refFileList;
	//refFileList.push_back("D:\\workspace\\Landsat\\features\\polygon_reference.shp");
	refFileList.push_back("D:\\workspace\\jiading\\features\\freeline_reference.shp");
	//refFileList.push_back("D:\\workspace\\Landsat\\shapes\\lines_output.shp");

	vector<string> srcFileList;
	//srcFileList.push_back("D:\\workspace\\Landsat\\features\\polygon_source.shp");
	srcFileList.push_back("D:\\workspace\\jiading\\features\\freeline_source.shp");

	//Shp2Map(srcfile,dstfile);
	vector<shape_feature> refFeatureList;
	vector<shape_feature> srcFeatureList;
	for (vector<string>::const_iterator iter = refFileList.begin();
		iter != refFileList.end(); ++iter)
		Shp2Features(*iter, refFeatureList);
	for (vector<string>::const_iterator iter = srcFileList.begin();
		iter != srcFileList.end(); ++iter)
		Shp2Features(*iter, srcFeatureList);

	//for (vector<shape_feature>::iterator iter = refFeatureList.begin(); iter != refFeatureList.end(); iter++)
	//{
	//	for (size_t i = 0; i < iter->points.size(); i++)
	//	{
	//		//fPoint world(iter->points[i].x, iter->points[i].y);
	//		//fPoint linesample;
	//		//gdalApp.world2linesample(world, linesample);
	//		//iter->points[i].x = linesample.x;
	//		//iter->points[i].y = linesample.y;
	//		iter->points[i].y = -iter->points[i].y;
	//	}
	//}
	features2File(refFeatureList, outRefFile);

	string imageFile = "D:\\workspace\\jiading\\radarsat\\imagery_HH_LS_8_Lee.tif";
	//string imageFile = "D:\\workspace\\Landsat\\features\\edge\\source1.tif";
	GdalRasterApp gdalApp;
	gdalApp.open(imageFile.c_str());
	int Height = gdalApp.getDataset()->GetRasterYSize();
	for (vector<shape_feature>::iterator iter = srcFeatureList.begin(); iter != srcFeatureList.end(); iter++)
	{
		for (size_t i = 0; i < iter->points.size(); i++)
		{
			//fPoint world(iter->points[i].x, iter->points[i].y);
			//fPoint linesample;
			//gdalApp.world2linesample(world, linesample);
			//iter->points[i].x = linesample.x;
			//iter->points[i].y = linesample.y;
			iter->points[i].y = Height - iter->points[i].y;
		}
	}
	features2File(srcFeatureList, outSrcFile);

	//appendStraightline2File("D:\\workspace\\Landsat\\lines.txt", outSrcFile, outRefFile, "D:\\workspace\\Landsat\\features\\reference.tif");
}

void Landsat()
{
	string outRefFile = "D:\\workspace\\Landsat\\xinjiang\\features\\reference.txt";
	string outSrcFile = "D:\\workspace\\Landsat\\xinjiang\\features\\source.txt";
	string fold = "D:\\workspace\\Landsat\\features";
	vector<string> refFileList;
	//refFileList.push_back("D:\\workspace\\Landsat\\features\\polygon_reference.shp");
	//refFileList.push_back("D:\\workspace\\Landsat\\xinjiang\\features\\area_reference.shp");
	refFileList.push_back("D:\\workspace\\Landsat\\xinjiang\\features\\multiple_reference.shp");
	//refFileList.push_back("D:\\workspace\\Landsat\\xinjiang\\straightline_reference.shp");
	//refFileList.push_back("D:\\workspace\\Landsat\\shapes\\lines_output.shp");

	vector<string> srcFileList;
	//srcFileList.push_back("D:\\workspace\\Landsat\\features\\polygon_source.shp");
	//srcFileList.push_back("D:\\workspace\\Landsat\\xinjiang\\features\\area_source.shp");
	srcFileList.push_back("D:\\workspace\\Landsat\\xinjiang\\features\\multiple_source.shp");
	//srcFileList.push_back("D:\\workspace\\Landsat\\xinjiang\\straightline_source.shp");

	//Shp2Map(srcfile,dstfile);
	vector<shape_feature> refFeatureList;
	vector<shape_feature> srcFeatureList;
	for (vector<string>::const_iterator iter = refFileList.begin();
		iter != refFileList.end(); ++iter)
		Shp2Features(*iter, refFeatureList);
	for (vector<string>::const_iterator iter = srcFileList.begin();
		iter != srcFileList.end(); ++iter)
		Shp2Features(*iter, srcFeatureList);

	//for (vector<shape_feature>::iterator iter = refFeatureList.begin(); iter != refFeatureList.end(); iter++)
	//{
	//	for (size_t i = 0; i < iter->points.size(); i++)
	//	{
	//		//fPoint world(iter->points[i].x, iter->points[i].y);
	//		//fPoint linesample;
	//		//gdalApp.world2linesample(world, linesample);
	//		//iter->points[i].x = linesample.x;
	//		//iter->points[i].y = linesample.y;
	//		iter->points[i].y = -iter->points[i].y;
	//	}
	//}
	features2File(refFeatureList, outRefFile);

	//string imageFile = "D:\\workspace\\jiading\\radarsat\\imagery_HH_LS_8_Lee.tif";
	string imageFile = "D:\\workspace\\Landsat\\xinjiang\\features\\edge\\source1.tif";
	GdalRasterApp gdalApp;
	gdalApp.open(imageFile.c_str());
	int Height = gdalApp.getDataset()->GetRasterYSize();
	for (vector<shape_feature>::iterator iter = srcFeatureList.begin(); iter != srcFeatureList.end(); iter++)
	{
		for (size_t i = 0; i < iter->points.size(); i++)
		{
			//fPoint world(iter->points[i].x, iter->points[i].y);
			//fPoint linesample;
			//gdalApp.world2linesample(world, linesample);
			//iter->points[i].x = linesample.x;
			//iter->points[i].y = linesample.y;
			//iter->points[i].y = -iter->points[i].y;
			iter->points[i].y = Height  -  iter->points[i].y;
		}
	}
	features2File(srcFeatureList, outSrcFile);

	//appendStraightline2File("D:\\workspace\\Landsat\\lines.txt", outSrcFile, outRefFile, "D:\\workspace\\Landsat\\features\\reference.tif");
}



bool shape2points(string shpFile, vector<fPoint>& gptList)
{
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动

	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");//得到shp文件的处理器

	OGRDataSource* poDS = poDriver->Open(shpFile.c_str(), NULL);//打开文件

	//int iLayerCount = poDS->GetLayerCount();

	OGRLayer* poLayer = poDS->GetLayer(0);//获取shp图层

	OGRSpatialReference *inOSRS = poLayer->GetSpatialRef();
	OGRSpatialReference poLatLong;
	poLatLong.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *poTransform = OGRCreateCoordinateTransformation(inOSRS, &poLatLong);
	if (poTransform == NULL)
	{
		cout << "Projection Error!" << endl;
		system("Pause");
	}

	//读取几何和属性值
	OGRFeature *poFeature;
	poLayer->ResetReading();//确保是从该层的开头开始
	gptList.clear();
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
		int iField;
		for (iField = 0; iField < poFDefn->GetFieldCount(); iField++)
		{
			OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
			printf("%s=", poFieldDefn->GetNameRef());
			if (poFieldDefn->GetType() == OFTInteger)
			{
				printf("%d\n", poFeature->GetFieldAsInteger(iField));
			}
			else if (poFieldDefn->GetType() == OFTReal)
			{
				printf("%.3f\n", poFeature->GetFieldAsDouble(iField));
			}
			else if (poFieldDefn->GetType() == OFTString)
			{
				printf("%s\n", poFeature->GetFieldAsString(iField));
			}
			else
			{
				printf("%s\n", poFeature->GetFieldAsString(iField));
			}
		}

		OGRGeometry *poGeometry;
		poGeometry = poFeature->GetGeometryRef();
		OGRwkbGeometryType geoType = poGeometry->getGeometryType();
		if (poGeometry != NULL	&& wkbFlatten(poGeometry->getGeometryType()) == wkbPoint)//点文件
		{
			OGRPoint *poPoint = (OGRPoint *)poGeometry;
			printf("%.3f,%.3f\n", poPoint->getX(), poPoint->getY());
			double x = poPoint->getX();
			double y = poPoint->getY();
			if (!poTransform->Transform(1, &x, &y))
			{
				cout << "x="<< x << "\t" << "y=" << y << endl;
			}
			fPoint gpt(poPoint->getX(), poPoint->getY());
			gptList.push_back(gpt);
		}
		else if (poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPoint)//多点文件
		{
			unsigned char *pData;//unsigned char*类型指针，用于指向数据缓冲区         
			pData = new unsigned char[poGeometry->WkbSize()];

			printf("WkbSize:%d\n", poGeometry->WkbSize());

			poGeometry->exportToWkb(wkbNDR, pData);

			int pos = 0;
			char byteOrder;
			int wkbType;
			int num_wkbFeature;
			int nPoint;
			memcpy(&byteOrder, pData + pos, sizeof(byteOrder));
			pos += sizeof(byteOrder);
			memcpy(&wkbType, pData + pos, sizeof(wkbType));
			pos += sizeof(wkbType);
			memcpy(&nPoint, pData + pos, sizeof(nPoint));
			pos += sizeof(nPoint);
			printf("byteOrder=%d\nwkbType=%d\nPoint=%d\n", byteOrder, wkbType, nPoint);

			for (int iPoint = 0; iPoint < nPoint; iPoint++)
			{
				int WKBPoint_pos = pos;
				char WKBPoint_byteOrder;
				int WKBPoint_wkbType;
				int nPoint;
				memcpy(&byteOrder, pData + WKBPoint_pos, sizeof(WKBPoint_byteOrder));
				WKBPoint_pos += sizeof(WKBPoint_byteOrder);
				memcpy(&wkbType, pData + WKBPoint_pos, sizeof(WKBPoint_wkbType));
				WKBPoint_pos += sizeof(WKBPoint_wkbType);

				double x;
				double y;
				memcpy(&x, pData + WKBPoint_pos, sizeof(x));
				WKBPoint_pos += sizeof(x);
				memcpy(&y, pData + WKBPoint_pos, sizeof(y));
				WKBPoint_pos += sizeof(y);

				//double x;
				//double y;
				//memcpy(&x, pData + pos, sizeof(x));
				//pos += sizeof(x);
				//memcpy(&y, pData + pos, sizeof(y));
				//pos += sizeof(y);

				if (!poTransform->Transform(1, &x, &y))
				{
					cout << "x=" << x << "\t" << "y=" << y << endl;
				}
				fPoint gpt(x, y);
				gptList.push_back(gpt);
				printf("x=%lf,y=%lf\n", x, y);
				pos = WKBPoint_pos;
			}
			OGRFeature::DestroyFeature(poFeature);
		}
		else
		{
			printf("Error: the shape file is not a point layer.\n");
			return false;
		}
		printf("\n");
	}
	OGRDataSource::DestroyDataSource(poDS);
	OGRCleanupAll();//资源清理
	return true;
}

void main()
{
	//jading();
	//Landsat();
	////string shpFile = "D:\\shp\\cmr_point.shp";
	////vector<fPoint> gptList;
	//shape2points(shpFile,gptList);
	printf("OK\n");
}