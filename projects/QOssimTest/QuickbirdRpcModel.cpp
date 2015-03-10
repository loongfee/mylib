#include "QuickbirdRpcModel.h"
#include <ossim/projection/ossimQuickbirdRpcModel.h>
#include <ossim/base/ossimKeywordlist.h>
#include <ossim/imaging/ossimQuickbirdTiffTileSource.h>
#include <ossim/imaging/ossimImageHandlerFactory.h>
#include <ossim/imaging/ossimQbTileFilesHandler.h>
#include <ossim/base/ossimRegExp.h>
#include "func.h"
QuickbirdRpcModel::QuickbirdRpcModel(void)
{
}

QuickbirdRpcModel::~QuickbirdRpcModel(void)
{
}

QuickbirdRpcModel::QuickbirdRpcModel(const ossimFilename& imgFileName)
{
	m_imageFileName = imgFileName;
	// 结果文件名
	m_resultName = imgFileName.fileNoExtension();
	vector<ossimString> strList = ossimString(m_resultName).split('-');
	if(3 == strList.size())
	{
		m_resultName = strList[0] + "-" + strList[2];
	}

	if(init(m_imageFileName))
	{
		m_bInitState = true;
	}
}
bool QuickbirdRpcModel::init(const ossimFilename& imgFileName)
{//Quick Bird
	ossimImageHandler* result = 0;
	result = new ossimQbTileFilesHandler();
	if(result->open(imgFileName))
	{
		ossimQuickbirdRpcModel *qbRpcModel = new ossimQuickbirdRpcModel((ossimQbTileFilesHandler*)result);
		qbRpcModel->saveRpcModelStruct(m_RpcStruct);
		m_sensorModel = new ossimQuickbirdRpcModel((ossimQbTileFilesHandler*)result);
		if(m_sensorModel)
		{
			return true;
		}
	}

	/*m_MapPar = PTR_CAST(ossimMapProjection, pProjection);
	m_sensorModel = PTR_CAST(ossimSensorModel, pProjection);
	if (m_MapPar && m_sensorModel) return true;*/
	return false;
}


bool QuickbirdRpcModel::readRPCFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct)
{
	return true;
}

bool QuickbirdRpcModel::readHeader(ossimFilename HDRfile)
{
	return true;
}

vector<ossimFilename> QuickbirdRpcModel::GetGcpFileNames()
{
	vector<ossimFilename> result;
	ossimDirectory(m_imageFileName.path()).findAllFilesThatMatch(result, "(gcp.txt|gcp.dat)$");
	return result;
}

void QuickbirdRpcModel::SetOutputBands()
{
	ossimInit::instance()->initialize();
	//ossimInit::instance()->loadPlugins(ossimFilename(getGdalPluginDir().toLatin1()));
	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(m_imageFileName);
	m_OutBandList.clear();
	//选择输出波段
	int nbands = handler->getNumberOfInputBands();
	vector<ossim_uint32> outBandList;
	if(1 == nbands)
	{
		m_OutBandList.push_back(0);
	}
	else if(3 == nbands)
	{
		//123
		m_OutBandList.push_back(0);
		m_OutBandList.push_back(1);
		m_OutBandList.push_back(2);
	}
	else if(6 == nbands)
	{
		//543
		m_OutBandList.push_back(4);
		m_OutBandList.push_back(3);
		m_OutBandList.push_back(2);
	}
	else if(7 == nbands)
	{
		//742
		m_OutBandList.push_back(6);
		m_OutBandList.push_back(3);
		m_OutBandList.push_back(1);
	}
}

int QuickbirdRpcModel::executeOrth(ossimFilename outPath, ossimFilename elevationPath/* = ""*/)
{
	return 0;
}