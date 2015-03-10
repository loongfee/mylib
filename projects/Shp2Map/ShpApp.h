#ifndef SHP_APP_H
#define SHP_APP_H
#include <ogrsf_frmts.h>
#include <gdal.h>
#include <ogr_api.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>
#include <string>
#include <vector>
using namespace std;

enum FeatureType
{
	featurePoint = 1,
	featureStraightLine,
	featureFreeLine,
	featurePolygon
};
const string FeatureTypeName[] = {"Unknown",
	"Point", "StraightLine", "FreeLine", "Polygon"};

struct fPoint
{
	double x;
	double y;
	fPoint(){};
	fPoint(const double& ax, const double& ay)
	{
		x = ax; 
		y = ay;
	}
};

struct shape_feature{
	string Id;
	int WkbSize;
	int byteOrder;
	int wkbType;
	int num_wkbFeature;
	int nPoint;
	FeatureType type;
	std::vector<fPoint> points;
};

bool Shp2Features(string srcfile, vector<shape_feature>& shape_feature_list, bool append_mode = true);

void features2File(const vector<shape_feature>& featureList, string dstfile);


#endif