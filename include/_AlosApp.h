#ifndef ALOSAPP_HEADER
#define ALOSAPP_HEADER
#include <iostream>
#include <iterator>
#include <fstream>
using namespace std;

#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimTempFilename.h>
#include <ossim/base/ossimString.h>
#include <ossim/base/ossimDirectory.h>
#include <ossim/base/ossimDirectoryTree.h>
#include <ossim/base/ossimGpt.h>
#include <ossim/base/ossimDpt.h>
#include <ossim/base/ossimKeywordNames.h>
#include <ossim/base/ossimTieGptSet.h> 
#include <ossim/base/ossimAdjustmentInfo.h>
#include <ossim/base/ossimIrect.h>
#include <ossim/base/ossimKeywordNames.h>
#include <ossim/base/ossimCommon.h>
#include <ossim/base/ossimRefPtr.h>
#include <ossim/base/ossimStreamFactoryRegistry.h>
#include <ossim/base/ossimStdOutProgress.h>
#include <ossim/base/ossimUnitTypeLut.h>
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
#include <ossim/projection/ossimProjection.h>
#include <ossim/projection/ossimUtmProjection.h>
#include <ossim/projection/ossimTransMercatorProjection.h>
#include <ossim/projection/ossimMapProjectionFactory.h>
#include <ossim/projection/ossimProjectionFactoryRegistry.h>
#include <ossim/init/ossimInit.h>
#include <ossim/projection/ossimSensorModel.h>
#include <ossim/projection/ossimRpcSolver.h>
#include <ossim/projection/ossimRpcModel.h>


#include "strUtil.h"
#include "fileUtil.h"
#include <QFileInfo>
#include <QDir>
#define USE_QT 0

namespace mylib{

class AlosModel
{
public:
	AlosModel(){};
	~AlosModel() 
	{
		//delete m_MapPar;
		m_bInitState =false;
	};

	std::vector<ossim_uint32> getOutBandList() {return m_OutBandList;};

	ossimRpcModel::rpcModelStruct getRpcModelStruct() {return m_rpcStruct;};

	ossimMapProjection* getMapProjection() {return m_MapPar;};

	ossimKeywordlist getKeywordlist() {return m_MapProjection;};

	ossim_int32 getColumsPerBand(){return m_Columns;};
	ossim_int32 getLinesPerBand(){return m_Lines;};

	bool getInitState() {return m_bInitState;};

	ossimString getSceneID(){return m_SceneID;};

protected:
	virtual bool init(const ossimFilename& AlosDir) = 0;
	virtual bool setProjection();
	virtual bool readRPCFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct);
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
	ossimPolygon		m_theBoundGndPolygon;
	ossimDrect		m_theImageClipRect;

	ossimKeywordlist	m_MapProjection;
	ossimMapProjection* m_MapPar;
	ossimRpcProjection m_rpcProjection;
	ossimRpcModel::rpcModelStruct m_rpcStruct;
	std::vector<ossim_uint32> m_OutBandList;

	bool			m_bInitState;
};

class AlosPRISM: public AlosModel
{
public:
	AlosPRISM(const ossimFilename& AlosDir);

	ossimString		m_FileHDR;	//头文件
	ossimString		m_FileTIF;	//TIF文件
	ossimString		m_FileRPC;	//RPC文件

protected:
	virtual bool init(const ossimFilename& AlosDir);
};

class AlosAVNIR2: public AlosModel
{
public:
	AlosAVNIR2(const ossimFilename& AlosDir);

	ossimString				m_FileHDR;	//头文件
	//ossimString				m_FileTIF[4];	//TIF文件
	ossimString				m_FileTIF;	//TIF文件
	ossimString				m_FileRPC;	//RPC文件

protected:
	virtual bool init(const ossimFilename& AlosDir);
};
}; // end of namespace mylib
#endif