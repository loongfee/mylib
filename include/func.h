#ifndef func_HEADER
#define func_HEADER

#include <ossim/radi/ossimblockadjustment.h>
#include <ossim_plugin/radi/radiRpcModel.h>

#include "mprojectdefine.h"
#include <strUtil.h>

namespace mylib{
enum
{
	ERR_SOURCE,
	ERR_GCP,
	ERR_DESTINATION,
	ERR_OPTIMIZE,
	ERR_REPORT,
	SUCCESS,
};


struct ProjectionParameters{
	ossimString DatumName;
	double TrueOriginLongitude;
	double TrueOriginLatitude;
	double EastingFalse;
	double NorthingFalse;
	double PixelSizeX;
	double PixelSizeY;
	ProjectionParameters()
	{
		DatumName = "Î÷°²80";
		TrueOriginLongitude = 120;
		TrueOriginLatitude = 0.0;
		EastingFalse = 40500000.0;
		NorthingFalse = 0.0;
		PixelSizeX = 0.5;
		PixelSizeY = 0.5;
	};
};


//int OptimizeGcp(wxString sourcefile, wxString gcpfile, wxString destinationfile, wxString elevationpath);
//
//void AutoOptimize(wxArrayString sourceFiles, wxArrayString gcpFiles, wxArrayString destinationFiles, wxString elevationpath);

ossimProjection* newUtmView(const ossimGpt& centerGround,
						   const ossimDpt& metersPerPixel);


ossimTieGptSet *ImportGptFile(ossimFilename gcpFile);

void AddGpts(ossimFilename gptFile, vector<ossimBlockTieGpt> &gptList);
void AddCpts(ossimFilename cptFile, vector<ossimBlockTieGpt> &gptList);
void ToLatLon(ossimSensorModel* sensorModel, ossimTieGptSet* gptSet);
void ToLatLon(ossimMapProjection* proj, ossimTieGptSet* gptSet);
ossimTieGptSet* ReadGptFromFile(ossimFilename gcpFile);
void AppendResidue2File(ossimSensorModel* sensorModel, ossimTieGptSet* GptSet, fstream &fs);
void InitializePrj(MyProject& prj, ossimFilename sourceFile, ossimFilename demPath, ossimFilename projectionFile);

bool GetBlockElevations(MyProject prj[], vector< ossimBlockTieGpt >& gptList);

bool saveGptList(const char* GptFile, MyProject prj[], vector< ossimBlockTieGpt > gptList);


bool AppendTieAreas(vector<ossimTieFeature>& tieFeatureList, const vector<ossimDArea> &dptAreaList, const vector<ossimGArea> &gptAreaList);

bool AppendTieLines(vector<ossimTieFeature>& tieFeatureList, const vector<ossimDLine> &dptLineList, const vector<ossimGLine> &gptLineList);

bool AppendTieFeatures(vector<ossimTieFeature>& tieFeatureList, const vector<ossimDFeature> &dptFeatureList, const vector<ossimGFeature> &gptFeatureList);

bool AppendTieFeaturesFromTieGptSet(vector<ossimTieFeature>& tieFeatureList, ossimTieGptSet* pGptSet);

ossimRpcModel *createRpcModelFromPoints(ossimTieGptSet* gptSet, const ossimDpt& imageShift);


ossimRpcModel *createRpcModelFromProjection(const ossimDrect& imageBounds,
							 ossimProjection* imageProj,
							 ossim_uint32 xSamples,
							 ossim_uint32 ySamples,
							 bool shiftTo0Flag);

ossimTieGptSet* createTieGptSet(const ossimDrect& imageBounds,
									  const ossimSensorModel& proj,
									  ossim_uint32 xSamples = 8,
									  ossim_uint32 ySamples = 8,
									  bool latlon = false,
									  bool shiftTo0Flag = false);

void createTieGptSet(const ossimDrect& imageBounds,
					 const ossimSensorModel& proj,
					 double height,
					 ossimTieGptSet*& pTieGptSet,
					 ossim_uint32 xSamples = 8,
					 ossim_uint32 ySamples = 8,
					 bool latlon = false,
					 bool shiftTo0Flag = false);

void create3DGridPoints(const ossimDrect& imageBounds,
					 const ossimSensorModel& proj,
					 double height,
					 ossimTieGptSet*& pTieGptSet,
					 ossim_uint32 xSamples = 8,
					 ossim_uint32 ySamples = 8,
					 bool latlon = false,
					 bool shiftTo0Flag = false);

template<class T1, class T2>
bool ReadFeatures(string strFilename, vector < T1 > & featureList, ossimFeatureType featureType = ossimFeatureType::ossimUnknown)
{
	char strtmp[1024];
	std::vector<string> strList;
	string str;

	fstream os;
	os.open(strFilename.c_str(), ios_base::in);

	vector<string> delimiterList;
	delimiterList.push_back(" ");
	delimiterList.push_back("\t");
	while(os.getline(strtmp,1024))
	{
		str = strtmp;
		if(str.empty()) continue;
		strList.clear();
		splitString(str, delimiterList, strList);
		//SplitString(str, "	", strList, false);
		//if (strList.size()==1) {
		//	str = strList[0];
		//	strList.clear();
		//	SplitString(str, " ", strList, false);
		//}
		if (strList.size() != 3) continue;
		T1 tmpFeature;
		tmpFeature.strId = strList[0];
		int num = atoi(strList[2].c_str());

		if(0 == strList[1].compare("Point"))
		{
			if(ossimFeatureType::ossimUnknown != featureType && ossimFeatureType::ossimPointType != featureType ) continue;
			tmpFeature.m_featureType = ossimFeatureType::ossimPointType;
		}
		else if(0 == strList[1].compare("StraightLine"))
		{
			if(ossimFeatureType::ossimUnknown != featureType && ossimFeatureType::ossimStraightLineType != featureType ) continue;
			tmpFeature.m_featureType = ossimFeatureType::ossimStraightLineType;
		}
		else if(0 == strList[1].compare("FreeLine"))
		{
			if(ossimFeatureType::ossimUnknown != featureType && ossimFeatureType::ossimFreeLineType != featureType ) continue;
			tmpFeature.m_featureType = ossimFeatureType::ossimFreeLineType;
		}
		else if(0 == strList[1].compare("Polygon"))
		{
			if(ossimFeatureType::ossimUnknown != featureType && ossimFeatureType::ossimPolygonType != featureType ) continue;
			tmpFeature.m_featureType = ossimFeatureType::ossimPolygonType;
		}
		else continue;

		for(int i = 0;i < num;i++)
		{
			if(!os.getline(strtmp, 1024)) return false;
			string strPoint = strtmp;
			if(strPoint.empty()) continue;
			vector<string> strCorList;
			splitString(strPoint, delimiterList, strCorList);
			if(strCorList.size() > 3) continue;
			T2 pt(atof(strCorList[0].c_str()), atof(strCorList[1].c_str()));
			tmpFeature.m_Points.push_back(pt);
		}
		featureList.push_back(tmpFeature);
	}
	return true;
}

//bool ReadFeatures(string strFilename, vector < ossimGFeature > & featureList, ossimFeatureType featureType = ossimFeatureType::ossimUnknown);

bool FeatureReport(ossimFilename reportfile, ossimSensorModel* sensorModel, vector<ossimTieFeature> tieFeatureList);

NEWMAT::ColumnVector CalcResidue(ossimSensorModel* sensorModel,vector<ossimTieFeature> tieFeatureList);


ossimMapProjection* CreateProjection(ProjectionParameters proParam, ossimKeywordlist& MapProjection);

//bool projection2ll(ossimFilename inFile, ossimFilename outFile);

bool reprojectionPoints(ossimFilename inFile, ossimKeywordlist outKwl, ossimFilename outFile);
NEWMAT::ColumnVector CalcResidue(ossimSensorModel* sensorModel,ossimTieGptSet gptSet, bool bPixel = false, bool bLonLat = false);
bool OutputReport(ossimFilename reportfile, ossimSensorModel* sensorModel, ossimTieGptSet* ctrlSet, ossimTieGptSet* chkSet, bool bPixel = false, bool bLonLat = false);

bool readRPBFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct);
bool readRPCFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct);

bool readRPBFile(ossimFilename rpcFile, ossimplugins::radiRpcModel::rpcModelStruct& rpcStruct);
bool readRPCFile(ossimFilename rpcFile, ossimplugins::radiRpcModel::rpcModelStruct& rpcStruct);

ossimRpcModel* createRpcModelFromPoints(ossimFilename gcpFile);
}
#endif