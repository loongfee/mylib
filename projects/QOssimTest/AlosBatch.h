#ifndef ALOSBATCH_HEADER
#define ALOSBATCH_HEADER
/************************************************************************/
/* Alos数据批量正射校正类                                               */
/************************************************************************/
#include <mprojectdefine.h>
//#include <AlosApp.h>
using namespace mylib;

class AlosBatch
{
public:
	enum AlosType{
		PRISM,
		AVNIR2,
	};
	struct AlosTask{
		ossimFilename inputPath;
		AlosType type; 
	};
public :
	AlosBatch(ossimFilename strInputPath, ossimFilename strOutputPath, int nType = 0);
	void setElevationPath(ossimFilename path){m_ElevationPath = path;};
	void setEgmFile(ossimFilename file){m_EgmFile = file;};
	bool Run();
	int GetFinishedCount() {return m_nFinished;};
	int GetTotalTaskCount() {return m_nTask;};

	bool Alos_PRISM_Rpc(ossimFilename AlosDir, ossimFilename elevationpath, ossimFilename EGMfile,
		ossimFilename outPath, ossimFilename logFile = "");
	bool Alos_AVNIR2_Rpc(ossimFilename AlosDir, ossimFilename elevationpath, ossimFilename EGMfile,
		ossimFilename outPath, ossimFilename logFile = "");

private:
	ossimFilename m_InputPath;
	ossimFilename m_OutputPath;

	ossimFilename m_ElevationPath;
	ossimFilename m_EgmFile;

	int m_Type;

	vector<AlosTask> m_AlosTask;

	int m_nTask;
	int m_nFinished;
	int m_nFailed;
	int m_nCurrentTask;

	bool addLogLine(ossimFilename strLogFile, ossimFilename strContet);
	bool setLogLine(ossimFilename strLogFile, int nLine, ossimFilename strContet);
};

#endif