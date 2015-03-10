#pragma once

#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimTempFilename.h>
#include <ossim/base/ossimString.h>
#include <ossim/base/ossimDirectory.h>
#include <ossim/base/ossimGpt.h>
#include <ossim/base/ossimDpt.h>
#include <ossim/base/ossimKeywordNames.h>
#include <ossim/base/ossimTieGptSet.h> 
#include <ossim/base/ossimAdjustmentInfo.h>
#include <ossim/base/ossimIrect.h>
#include <ossim/base/ossimKeywordNames.h>
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
#include <ossim/projection/ossimLandSatModel.h>
#include <ossim/projection/ossimspot5model.h>
#include <ossim/projection/ossimSensorModel.h>
#include <ossim/projection/ossimRpcSolver.h>
#include "ossim/projection/ossimRpcModel.h"
#include <ossim/projection/ossimPolynomProjection.h>
#include <ossim/base/ossimEndian.h>
#include <ossim/base/ossimUnitTypeLut.h>
#include <ossim/base/ossimIoStream.h>
#include <ossim/base/ossimCommon.h>
#include <ossim/base/ossimRefPtr.h>
#include <ossim/base/ossimStreamFactoryRegistry.h>
#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimStdOutProgress.h>
#include <ossim/projection/ossimImageViewProjectionTransform.h>

#include <time.h>
//#include <qdir.h>
//#include <qstring.h>

#include <ogr_spatialref.h>
#include <gdal_priv.h>

#include <iostream>
#include <iterator>
#include <fstream>

#include "OrthProcessListener.h"
using namespace std;

class RpcModel
{
public:
	RpcModel(void);
	~RpcModel(){};

public:
	std::vector<ossim_uint32> getOutBandList() {return m_OutBandList;};

	ossimRpcModel::rpcModelStruct getRpcModelStruct() {return m_rpcStruct;};

	ossimMapProjection* getMapProjection() {return m_MapPar;};

	ossimKeywordlist getKeywordlist() {return m_MapProjection;};

	ossim_int32 getColumnsPerBand(){return m_Columns;};
	ossim_int32 getLinesPerBand(){return m_Lines;};

	bool getInitState() {return m_bInitState;};

	ossimString getSceneID(){return m_SceneID;};

	ossimString getPath(){return m_RSPPath;};
	ossimString getFrame(){return m_OSSIMrame;};
	ossimString getSceneCenterTime(){return m_SceneCenterTime;};
	ossimString getUtmZone(){return m_UTMZone;};
	//ossimDrect  getBoundUtm(){return m_theBoundUtm;};
	//ossimDrect  getBoundLatLon(){return m_theBoundLatLon;};
	ossimPolygon getBoundUtm(){return m_theBoundUtmPolygon;};
	ossimPolygon getBoundGeo(){return m_theBoundGeoPolygon;};
	ossim_int64	getPixelSize(){return m_PixelSize;};

	virtual int executeOrth(ossimFilename outFile, ossimFilename elevationPath = "") = 0;
	//QtRspfListener *getOutProgress() {return m_progress;};
	int GetProgress();
	bool m_IsFinished;
	virtual void setProjection(ossimMapProjection* MapPar) {m_MapPar = MapPar;};
	void setOutFileType(const char* type) {m_outFileType = ossimFilename(type);};
	void setOutSampleType(const char* type) {m_outSampleType = ossimFilename(type);};
	void setUseGcps(bool bUseGcps){m_bUseGcps = bUseGcps;};
	void setReplace(bool bReplace){m_bReplace = bReplace;};
	void SavePointToFile(ossimFilename filenametosave,ossimTieGptSet* m_gptset,ossimTieGptSet* mcheck_gptset);

	bool Orthograph(ossimFilename outfile);
	bool Orthograph_Tile(ossimFilename outfile, int tileSize);
	void UpdateSensorModel(ossimTieGptSet tieGptSet,
		ossimSensorModel* &sensorModel,
		ossimKeywordlist& geom);
	virtual bool GetElevations(ossimTieGptSet* &ctrlSet, double default_hgt = 0.0);
	bool Clone_Tile(ossimFilename outfile, int tileSize);
	NEWMAT::ColumnVector CalcResidue(ossimSensorModel* sensorModel,ossimTieGptSet gptSet);
	bool OutputReport(ossimFilename reportfile, ossimSensorModel* sensorModel, ossimTieGptSet* ctrlSet, ossimTieGptSet* chkSet);
	virtual void ReadGcpAndProjection(ossimFilename strFilename, ossimTieGptSet* &m_gptset, ossimTieGptSet* &mcheck_gptset);

	ossimTieGptSet*	m_ctrSet;
	ossimTieGptSet*	m_chkSet;
	ossimKeywordlist m_geom;   //存储卫星优化后轨道参数
	ossimFilename	m_outFileType;
	ossimFilename	m_outSampleType;
	ossimMapProjection* m_MapPar;
	ossimSensorModel* m_sensorModel;
	ossimElevManager* m_theElevManager;

protected:
	virtual bool init(const ossimFilename& AlosDir) = 0;
	virtual bool setProjection();
	virtual bool readRPCFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct) = 0;
	virtual bool readHeader(ossimFilename HDRfile) = 0;
	virtual ossimFilename GetOrthFileName(ossimFilename outPath) = 0;
	virtual ossimFilename GetGcpFileName() = 0;
	virtual ossimFilename GetReportFileName(ossimFilename outPath) = 0;
	//virtual void AppendGcpAndProjection(ossimFilename strFilename, ossimTieGptSet* &m_gptset, ossimTieGptSet* &mcheck_gptset);
protected:
	ossimString		m_SceneID;
	ossimString		m_ProductID;
	ossimString		m_Projection;
	ossim_float64	m_PixelSize;
	ossim_int32		m_Columns;
	ossim_int32		m_Lines;
	ossimString		m_UTMZone;
	ossimString		m_Datum;
	ossimString		m_EllipsoidModel;
	ossimString		m_SceneCenterTime;
	ossimDpt			m_theRefImgPt;
	ossimGpt			m_theRefGndPt;
	ossimDrect		m_theBoundUtm;
	ossimDrect		m_theBoundLatLon;
	//ossimDrect		m_theImageClipRect;
	ossimString		m_RSPPath;
	ossimString		m_OSSIMrame;
	ossimPolygon		m_theBoundUtmPolygon;
	ossimPolygon		m_theBoundGeoPolygon;

	ossimDpt			m_UtmUL;
	ossimDpt			m_UtmUR;
	ossimDpt			m_UtmLR;
	ossimDpt			m_UtmLL;

	ossimDpt			m_GeoUL;
	ossimDpt			m_GeoUR;
	ossimDpt			m_GeoLR;
	ossimDpt			m_GeoLL;

	ossimFilename m_FileOutType;
	//m_FileOutType = "pix";
	ossimFilename m_SampleType ;

	ossimRpcProjection m_rpcProjection;
	ossimRpcModel::rpcModelStruct m_rpcStruct;
	std::vector<ossim_uint32> m_OutBandList;
	ossimKeywordlist m_MapProjection;	//投影名称（包含投影参数和输出分辨率）
	ossimFilename	m_imageFileName;
	ossimFilename	m_resultName;
	ossimFilename	m_foldName;

	bool			m_bUseGcps;
	bool			m_bReplace;
	bool			m_bInitState;
	int				m_TileCount;
	int				m_CurrentTile;
	int				m_orthDone;

	OrthProcessListener *m_progress;
};