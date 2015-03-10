#ifndef MyProject_HEADER
#define MyProject_HEADER
#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimTempFilename.h>
#include <ossim/base/ossimString.h>
//#include <ossim/base/ossimXmlNode.h>
#include <ossim/base/ossimDirectory.h>
#include <ossim/base/ossimDirectoryTree.h>
//#include <ossim/base/ossimXmlDocument.h>
#include <ossim/base/ossimGpt.h>
#include <ossim/base/ossimDpt.h>
#include <ossim/base/ossimKeywordNames.h>
#include <ossim/base/ossimTieGptSet.h> 
#include <ossim/base/ossimAdjustmentInfo.h>
#include <ossim/base/ossimIrect.h>
#include <ossim/base/ossimKeywordNames.h>
//#include <ossim/base/ossimStdOutProgress.h>
//#include <ossim/base/ossimLeastSquaresBilin.h>
#include <ossim/elevation/ossimElevManager.h>
#include <ossim/imaging/ossimImageHandler.h>
#include <ossim/imaging/ossimImageHandlerRegistry.h>
#include <ossim/imaging/ossimImageWriterFactoryRegistry.h>
#include <ossim/imaging/ossimImageFileWriter.h>
#include <ossim/imaging/ossimPolyCutter.h>
#include <ossim/imaging/ossimImageData.h>
#include <ossim/imaging/ossimBandSelector.h>
#include <ossim/imaging/ossimCacheTileSource.h>
#include <ossim/imaging/ossimImageRenderer.h>
#include <ossim/matrix/newmat.h>
#include <ossim/matrix/newmatrc.h>
#include <ossim/projection/ossimProjection.h>
#include <ossim/projection/ossimUtmProjection.h>
#include <ossim/projection/ossimTransMercatorProjection.h>
#include <ossim/projection/ossimMapProjectionFactory.h>
#include <ossim/projection/ossimProjectionFactoryRegistry.h>
#include <ossim/init/ossimInit.h>
//#include <ossim/projection/ossimIkonosRpcModel.h>
//#include <ossim/projection/ossimquickbirdrpcmodel.h>
#include <ossim/projection/ossimLandSatModel.h>
#include <ossim/projection/ossimspot5model.h>
//#include <ossim/radi/ossimTheosmodel.h>
#include <ossim_plugin/radi/ossimHj1Model.h>
#include <ossim/projection/ossimSensorModel.h>
#include <ossim/radi/ossimBlockAdjustment.h>
#include <ossim/projection/ossimRpcSolver.h>
#include "ossim/projection/ossimRpcModel.h"
#include <ossim/projection/ossimPolynomProjection.h>
#include <ossim/base/ossim2dLinearRegression.h>

#include <ossim/base/ossimEndian.h>
#include <ossim/base/ossimUnitTypeLut.h>
#include <ossim/base/ossimIoStream.h>
#include <ossim/base/ossimCommon.h>
#include <ossim/base/ossimRefPtr.h>
#include <ossim/base/ossimStreamFactoryRegistry.h>
#include <ossim/base/ossimFilename.h>
//#include <ossim/support_data/ossimFfL5.h>
//#include <ossim/support_data/ossimWorldviewRpcHeader.h>
#include <ossim/base/ossimStdOutProgress.h>
#include <ossim/imaging/ossimImageChain.h>
#include <ossim/util/ossimChipperUtil.h>

#include <string>
#include <vector>
#include <iostream>
using namespace std;

//#include "AlosApp.h"
#include <fileUtil.h>
#include <strUtil.h>
#include <gcpUtil.h>

namespace mylib{

enum SensorType
{
	modelLandsat5,
	modelSpot5,
	modelLandsat7,
	modelTheos,
	modelAlos,
	modelAlosAVNIR2_1B2,
	modelAlosPRISM_1B2,
	modelHJ1,
	UnknowMole
};

enum ModelType
{
	OrbitType,//代表卫星轨道模型被选中
	RPCType,//代表有理函数模型被选中
	PolynomialType,//代表多项式模型被选中
	AffineType,		//代表仿射变换模型
	CombinedAdjustmentType, //代表多景联合平差
	UnknowType
};

typedef struct point
{
	int indext;//序号
	double x;//待纠正坐标X
	double y;//待纠正坐标Y
	double x2;//标准坐标X
	double y2;//标准坐标Y
	double z;//DEM高程

}pointinfo;
struct r_cluster
{
	int indext;
	double r;
};

class myOutProgress: public ossimStdOutProgress
{
public:
	myOutProgress(ossim_uint32 precision = 0, bool flushStream=false)
		:ossimStdOutProgress(precision,flushStream)
	{};
	ossim_uint32 getPrecision() const {return thePrecision;};
};

class MyProject
{
public:
	MyProject();
	~MyProject();
public:
	ModelType   m_ModelType;
	ossimFilename m_ProjFileName;
	ossimFilename m_ImgFileNameUnc;     //待校正影像文件名
	ossimFilename m_ImgFileNameRef;     //参考影像文件名
	ossimFilename m_ImgFileNameout;    //校正结果影像文件名

	ossimFilename m_ImgFileNameUncHistogram;//1124
	ossimFilename m_ImgFileNameRefHistogram;//1124
	ossimString   m_SampleType;        //采样方法1125
	ossimString   m_FileOutType;      //输出文件类型1125
	ossimString   m_ScaleType;
	ossimString   m_PolyDegree;

	//ossimString m_SensorName;        
	SensorType m_SensorName;        //传感器名称 LANDSAT5 SPOT5 1120
	//ossimBandSelector * m_BandSelectUnc; //待校正影像波段选择
	//ossimBandSelector * m_BandSelectRef; //参考影像波段选择
	std::vector<ossim_uint32> m_RefBandList;//参考影像波段选择
	std::vector<ossim_uint32> m_UncBandList;//待校正影像波段选择
	std::vector<ossim_uint32> m_OutBandList;//输出影像波段选择

	ossimFilename  m_HomoPointsFilename;  //同名点文件路径
	ossimTieGptSet* m_CtrlGptSet;    //用于控制点选取界面中 控制点的保存
	ossimTieGptSet* m_ChkGptSet;     //用于控制点选取界面中 检查点的保存
	ossimKeywordlist geom;   //存储卫星优化后轨道参数

	ossimKeywordlist m_MapProjection;  //投影名称（包含投影参数和输出分辨率）
	ossimMapProjection * m_MapPar;       //投影名称（包含投影参数和输出分辨率）

	ossimTieGptSet* m_OptCtrlGptSet;   //用于控制点优化界面中 控制点的保存
	ossimTieGptSet* m_OptChkGptSet;    //用于控制点优化界面中 检查点的保存
	ossimTieGptSet* m_OptRejectGptSet;  //用于控制点优化界面中 剔除点的保存
	int starline;   ///以下四个变量表明子区大小，默认为全图
	int starpixel;
	int endpixel;
	int endline;
	ossimFilename  m_DemPath;            //DEM路径
	ossimElevManager* theMgr;
	ossimSensorModel* m_sensorModel;
public:
	void write(ossimFilename m_filename);
	void  read(ossimFilename m_filename);
	void ReadGcpAndProjection(ossimFilename strFilename,ossimTieGptSet* &m_gptset,ossimTieGptSet* &mcheck_gptset);
	void ReadGcpAndProjection(ossimFilename strFilename);
//1120
	void SavePointToFile(ossimFilename filenametosave,ossimTieGptSet* m_gptset,ossimTieGptSet* mcheck_gptset, bool printProjection=true);
	void SaveOptPointToFile(ossimFilename filenametosave,ossimTieGptSet* m_gptset);
//1120end
	SensorType getSensorType(ossimFilename imgFileName);
	ModelType   getModelType();

	bool CheckSenserModel(ossimSensorModel* &sensorModel,
		ossimTieGptSet* &allGptSet,
		ossimTieGptSet* &ctrlSet,
		ossimTieGptSet* &chkSet,
		ossimTieGptSet* &errSet,
		double meter_threshold,
		NEWMAT::ColumnVector &ctrlResidue,
		NEWMAT::ColumnVector &chkResidue,
		NEWMAT::ColumnVector &errResidue);

	bool DistributeOptimize(ossimTieGptSet* srcGptSet, ossimTieGptSet* &ctrlGptSet, ossimTieGptSet* &chkGtpSet,int nControl, int nCheck);

	void UpdateSensorModel(ossimTieGptSet tieGptSet,
		ossimSensorModel* &sensorModel,
		ossimKeywordlist& geom);
	void UpdateSensorModel(vector < ossimTieFeature > tieFeatureList,
		ossimSensorModel* &sensorModel,
		ossimKeywordlist& geom);
	int OptimizeGcp( ossimSensorModel* sensorModel,ossimTieGptSet* totalCtrlGptSet,
		ossimTieGptSet* &ctrlGptSet,ossimTieGptSet* &chkGptSet,ossimTieGptSet* &errGptSet,
		NEWMAT::ColumnVector &ctrlResidue,
		NEWMAT::ColumnVector &chkResidue,
		NEWMAT::ColumnVector &errResidue);
	bool GetElevations(ossimTieGptSet* &ctrlSet, double defaultElev = 0.0, bool forceDefault = false);
	bool GetElevations(vector < ossimBlockTieGpt > &tiePointList, double defaultElev = 0.0, bool forceDefault = false);
	bool InitiateSensorModel(ossimFilename sourcefile = "");
	bool Orthograph(ossimFilename outfile);
	NEWMAT::ColumnVector getResidue(ossimSensorModel *sensorModel, vector< ossimTieFeature > tieFeatureList, bool useImageObs = true /* 2010.1.18 loong*/);
	myOutProgress *getOutProgress() {return m_progress;};

	bool SaveFeaturetoShape(ossimFilename filenametosave,vector<ossimTieFeature> featureList,
		ossimFeatureType  featureType = ossimFeatureType::ossimPolygonType);
	bool UpdateFeatures(ossimSensorModel* sensorModel, vector<ossimTieFeature>& featureList);

	bool OutputReport1(ossimFilename reportfile, ossimSensorModel* sensorModel, ossimTieGptSet* ctrlSet, ossimTieGptSet* chkSet, bool bPixel = false);
	NEWMAT::ColumnVector CalcResidue1(ossimSensorModel* sensorModel,ossimTieGptSet gptSet, bool bPixel = false);
	bool CreateL5Projection(const ossimFilename &headerFile);
	void ReadImagePoints(ossimFilename strFilename, vector < ossimDFeature > & dptPointList);
	void ReadGroundPoints(ossimFilename strFilename, vector < ossimGFeature > & gptPointList);

	/*void ReadImageLines(ossimFilename strFilename, vector < ossimDFeature > & dptLineList);
	void ReadGroundLines(ossimFilename strFilename, vector < ossimGFeature > & gptLineList);

	void ReadImageAreas(ossimFilename strFilename, vector < ossimDFeature > & dptAreaList);
	void ReadGroundAreas(ossimFilename strFilename, vector < ossimGFeature > & gptAreaList);*/
	bool SavePoint2Shape(ossimFilename filenametosave,ossimTieGptSet* m_gptset);
	bool GetElevations(vector<ossimTieFeature> &tieFeatureList, double defaultElev = 0.0, bool forceDefault = false);
	bool OutputReport(ossimFilename reportfile, ossimSensorModel* sensorModel, ossimTieGptSet* ctrlSet, ossimTieGptSet* chkSet, bool bPixel=false);
	bool OutputReport(fstream& fs, ossimSensorModel* sensorModel, ossimTieGptSet* ctrlSet, ossimTieGptSet* chkSet);
	NEWMAT::ColumnVector CalcResidue(ossimSensorModel* sensorModel,ossimTieGptSet gptSet, bool bPixel=false);

protected: 
	bool WriteShpPoint(ossimFilename filenametosave,ossimTieGptSet* m_gptset,ossimTieGptSet* mcheck_gptset);
	bool readRPCFile(const char* rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct);
private:
	myOutProgress *m_progress;
}; 
}; // end of namespace mylib
//MyProject m_ProjectS;
#endif
