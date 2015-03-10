#include "ShpApp.h"

bool Shp2Features(string srcfile, vector<shape_feature>& shape_feature_list, bool append_mode/* = true*/)
{
	if (!append_mode)
	{
		shape_feature_list.clear();
	}

	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动

	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");//得到shp文件的处理器
	
	OGRDataSource* poDS = poDriver->Open(srcfile.c_str(), NULL);//打开文件

	//int iLayerCount = poDS->GetLayerCount();

	OGRLayer* poLayer = poDS->GetLayer(0);//获取shp图层

	//读取几何和属性值
	OGRFeature *poFeature;
	poLayer->ResetReading();//确保是从该层的开头开始
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		shape_feature shpFeature;
		OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();

		OGRGeometry *poGeometry;
		poGeometry = poFeature->GetGeometryRef();
		OGRwkbGeometryType geoType = poGeometry->getGeometryType();

		if (poGeometry != NULL
			&& wkbFlatten(poGeometry->getGeometryType()) == wkbPoint)//点文件
		{
			shpFeature.wkbType = wkbPoint;
			shpFeature.nPoint = 1;
			shpFeature.type = FeatureType::featurePoint;
			OGRPoint *poPoint = (OGRPoint *)poGeometry;
			printf("%.3f,%.3f\n", poPoint->getX(), poPoint->getY());
			//fprintf(pf, "x=%lf,y=%lf\n",poPoint->getX(), poPoint->getY() );
			shpFeature.points.push_back(fPoint(poPoint->getX(), poPoint->getY()));
		}
		else if (poGeometry != NULL
			&& wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon)//多边形文件
		{
			unsigned char *pData = new unsigned char[poGeometry->WkbSize()];


			printf("WkbSize:%d\n", poGeometry->WkbSize());
			shpFeature.WkbSize = poGeometry->WkbSize();

			poGeometry->exportToWkb(wkbNDR, pData);
			int pos = 0;
			char byteOrder;
			int wkbType;
			int num_wkbFeature;
			memcpy(&byteOrder, pData + pos, sizeof(byteOrder));
			pos += sizeof(byteOrder);
			memcpy(&wkbType, pData + pos, sizeof(wkbType));
			pos += sizeof(wkbType);
			memcpy(&num_wkbFeature, pData + pos, sizeof(num_wkbFeature));
			pos += sizeof(num_wkbFeature);
			printf("byteOrder=%d\nwkbType=%d\nnum_wkbFeature=%d\n", byteOrder, wkbType, num_wkbFeature);
			shpFeature.byteOrder = byteOrder;
			shpFeature.wkbType = wkbType;
			shpFeature.num_wkbFeature = num_wkbFeature;
			shpFeature.type = FeatureType::featurePolygon;
			for (int iFeature = 0; iFeature < num_wkbFeature; iFeature++)
			{
				int nPoint;
				memcpy(&nPoint, pData + pos, sizeof(nPoint));
				pos += sizeof(nPoint);
				printf("nPoint=%d\n", nPoint);
				shpFeature.nPoint = nPoint;
				for (int iPoint = 0; iPoint < nPoint; iPoint++)
				{
					double x;
					double y;
					memcpy(&x, pData + pos, sizeof(x));
					pos += sizeof(x);
					memcpy(&y, pData + pos, sizeof(y));
					pos += sizeof(y);
					printf("x=%lf,y=%lf\n", x, y);
					//fprintf(pf, "x=%lf,y=%lf\n",x,y);
					shpFeature.points.push_back(fPoint(x, y));
				}
			}
			OGRFeature::DestroyFeature(poFeature);
		}
		else if (poGeometry != NULL
			&& wkbFlatten(poGeometry->getGeometryType()) == wkbLineString)//线文件
		{

			unsigned char *pData = new unsigned char[poGeometry->WkbSize()];


			printf("WkbSize:%d\n", poGeometry->WkbSize());
			shpFeature.WkbSize = poGeometry->WkbSize();

			poGeometry->exportToWkb(wkbNDR, pData);
			int pos = 0;
			char byteOrder;
			int wkbType;
			int nPoint;
			memcpy(&byteOrder, pData + pos, sizeof(byteOrder));
			pos += sizeof(byteOrder);
			memcpy(&wkbType, pData + pos, sizeof(wkbType));
			pos += sizeof(wkbType);
			memcpy(&nPoint, pData + pos, sizeof(nPoint));
			pos += sizeof(nPoint);
			printf("byteOrder=%d\nwkbType=%d\nPoint=%d\n", byteOrder, wkbType, nPoint);

			shpFeature.byteOrder = byteOrder;
			shpFeature.wkbType = wkbType;
			shpFeature.num_wkbFeature = 1;
			shpFeature.nPoint = nPoint;
			if (nPoint == 2)
			{
				shpFeature.type = FeatureType::featureStraightLine;
			}
			else
				shpFeature.type = FeatureType::featureFreeLine;
			for (int iPoint = 0; iPoint < nPoint; iPoint++)
			{
				double x;
				double y;
				memcpy(&x, pData + pos, sizeof(x));
				pos += sizeof(x);
				memcpy(&y, pData + pos, sizeof(y));
				pos += sizeof(y);
				printf("x=%lf,y=%lf\n", x, y);
				//fprintf(pf, "x=%lf,y=%lf\n",x,y);
				//fprintf(pf, "%lf\t%lf\n", x, y);
				shpFeature.points.push_back(fPoint(x, y));
			}
			OGRFeature::DestroyFeature(poFeature);
		}
		else
		{
			printf("no geometry\n");
		}
		printf("\n");
		shape_feature_list.push_back(shpFeature);
	}
	OGRDataSource::DestroyDataSource(poDS);
	OGRCleanupAll();//资源清理
	return TRUE;
}

void features2File(const vector<shape_feature>& featureList, string dstfile)
{
	FILE *pf;
	fopen_s(&pf, dstfile.c_str(), "w+");
	int nFeatures = (int)featureList.size();
	int counter[4] = {0};
	for (size_t i = 0; i < featureList.size(); i++)
	{
		fprintf_s(pf, "%d%03d\t%s\t%d\n", featureList[i].type, ++counter[featureList[i].type-1],
			FeatureTypeName[featureList[i].type].c_str(), featureList[i].nPoint);
		for (size_t j = 0; j < featureList[i].nPoint; j++)
		{
			fprintf_s(pf, "%lf\t%lf\n", featureList[i].points[j].x, featureList[i].points[j].y);
		}
	}
	fclose(pf);
}