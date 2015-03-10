#pragma once
#include <ossim/projection/ossimQuickbirdRpcModel.h>
#include "RpcModel.h"

class QuickbirdRpcModel :
	public RpcModel
{
public:
	QuickbirdRpcModel(void);
	QuickbirdRpcModel(const ossimFilename& dataDir);
	~QuickbirdRpcModel();

	ossimString		m_FileHDR;	//头文件
	ossimString		m_FileTIF;	//TIF文件
	ossimString		m_FileRPC;	//RPC文件

	ossimRpcModel::rpcModelStruct	m_RpcStruct;

	virtual int executeOrth(ossimFilename outFile, ossimFilename elevationPath = "");

protected:
	virtual bool init(const ossimFilename& dataDir);
	virtual bool readRPCFile(ossimFilename rpcFile, ossimRpcModel::rpcModelStruct& rpcStruct);
	virtual bool readHeader(ossimFilename HDRfile);
	virtual ossimFilename GetOrthFileName(ossimFilename outPath){return "";};
	virtual ossimFilename GetGcpFileName(){return "";};
	virtual ossimFilename GetReportFileName(ossimFilename outPath){return "";};
	vector<ossimFilename> GetGcpFileNames();
	void SetOutputBands();
};