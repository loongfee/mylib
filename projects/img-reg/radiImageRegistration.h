#ifndef radiImageRegistration_HEADER
#define radiImageRegistration_HEADER

#include <ossim/base/ossimString.h>
#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimOutputSource.h>
#include <ossim/base/ossimProcessInterface.h>
#include <ossim/base/ossimProcessProgressEvent.h>
#include <ossim/imaging/ossimFilterResampler.h>
#include <ossim/imaging/ossimImageChain.h>
#include <ossim/imaging/ossimImageHandler.h>
#include <ossim/imaging/ossimBandSelector.h>
#include <ossim/imaging/ossimImageRenderer.h>
#include <ossim/imaging/ossimCastTileSourceFilter.h>
#include <ossim/base/ossimTDpt.h>
#include <ossim_plugin/radi/radiBlockTieGpt.h>
#include <ossim_plugin/radi/radiBlockTieGptSet.h>
#include <vector>
//#include <mpi.h>
//#include <ossim/parallel/ossimMpi.h>

#include <ossim/elevation/ossimElevManager.h>

#include <ogrsf_frmts.h>
#include <gdal.h>
#include <ogr_api.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/xfeatures2d/nonfree.hpp"
//#include "opencv2/legacy/legacy.hpp"
//#include "opencv2/nonfree/nonfree.hpp"
#include <algorithm>
//#include "RANSAC/ransac/ransac.h"
//#include "RANSAC/estimators/Solver.h"
////#include "RANSAC/estimators/affineSolver.h"
////#include "RANSAC/estimators/affineError.h"
//#include "RANSAC/estimators/rigidSolver.h"
//#include "RANSAC/estimators/rigidError.h"
#include <GdalRasterApp.h>

#include <ez_ransac.h>

#include <levmar.h>
#include "mpi.h"
#include "functions.h"
using namespace mylib;
//using namespace groupsac;
using namespace std;
using namespace cv;
using namespace ossimplugins;

class ossimImageGeometry;
class ossimMapProjection;
class ossimListenerManager;

static char* strGoogle = "Google";
static char* strBing = "Bing";
static char* strMapquest = "Mapquest";
static char* strMapbox = "Mapbox";

//struct row_col{
//	double row_idx;
//	double col_idx;
//	row_col(double r, double c)
//	{
//		row_idx = r;
//		col_idx = c;
//	}
//};
//
//enum block_position{
//	left = 0,
//	right,
//	top,
//	bottom,
//	left_top,
//	left_bottom,
//	right_top,
//	right_bottom,
//	center,
//};
//
//struct image_block
//{
//	ossimIrect rect;
//	block_position position;
//};


#include <OpenThreads/Thread>
#include <OpenThreads/Mutex>
#include <OpenThreads/Barrier>
#include "radiImageRegistration.h"

#ifdef _WIN32
#include <process.h>
#define getpid() _getpid()
#else
#include <unistd.h>
#endif 

static OpenThreads::Barrier bar;
static int GLOBAL_NUM_THREADS;
static int totalBlocks;
static int finishedBlocks;

class radiImageRegistration;

//bool RowColCompare(const row_col& rc1, const row_col& rc2)
//{
//	return (fabs(rc1.row_idx)+fabs(rc1.col_idx)) < (fabs(rc2.row_idx)+fabs(rc2.col_idx));
//}

//static int row_col_step = 7;
//// multiple keywords
//static bool RowColCompare(row_col rc1, row_col rc2)
//{
//	return (fabs(rc1.row_idx)+fabs(rc1.col_idx)) < (fabs(rc2.row_idx)+fabs(rc2.col_idx));
//	//int fd1 = fabs(rc1.row_idx)+fabs(rc1.col_idx);
//	//int fd2 = fabs(rc2.row_idx)+fabs(rc2.col_idx);
//	//int d1 = (int)(fd1 + 0.5);
//	//int d2 = (int)(fd2 + 0.5);
//	//int m1 = d1%row_col_step;
//	//int m2 = d2 %row_col_step;
//	//if (m1 == m2)
//	//{
//	//	return fd1 < fd2;
//	//}
//	//return m1 < m2;
//}
////bool RowColCompare(row_col rc1, row_col rc2)
////{
////	//return (fabs(rc1.row_idx)+fabs(rc1.col_idx)) < (fabs(rc2.row_idx)+fabs(rc2.col_idx));
////	int fd1 = fabs(rc1.row_idx)+fabs(rc1.col_idx);
////	int fd2 = fabs(rc2.row_idx)+fabs(rc2.col_idx);
////	int d1 = (int)(fd1 + 0.5);
////	int d2 = (int)(fd2 + 0.5);
////	int m1 = d1%row_col_step;
////	int m2 = d2 %row_col_step;
////	if (m1 == m2)
////	{
////		return fd1 < fd2;
////	}
////	return m1 < m2;
////}
//static vector<row_col> assign_row_col_list(int tilerows, int tilecols, block_position position = block_position::center)
//{
//	vector<row_col> row_col_List;
//	double center_row = (tilerows-1) * 0.5;
//	double center_col = (tilecols-1) * 0.5;
//	double start_row = center_row;
//	double start_col = center_col;
//
//	//bool row_first = true;
//
//	//if (position == block_position::left
//	//	|| position == block_position::left_bottom
//	//	|| position == block_position::left_top)
//	//{
//	//	start_col = 0.0;
//	//	if (position == block_position::left)
//	//	{
//	//		row_first = false;
//	//	}
//	//}
//	//else if (position == block_position::right
//	//	|| position == block_position::right_bottom
//	//	|| position == block_position::right_top)
//	//{
//	//	start_col = tilecols;
//	//	if (position == block_position::right)
//	//	{
//	//		row_first = false;
//	//	}
//	//}
//	//else
//	//{
//	//	start_col = center_col;
//	//}
//
//	//if (position == block_position::top
//	//	|| position == block_position::left_top
//	//	|| position == block_position::right_top)
//	//{
//	//	start_row = 0;
//	//	if (position == block_position::top)
//	//	{
//	//		row_first = true;
//	//	}
//	//}
//	//else if (position == block_position::bottom
//	//	|| position == block_position::left_bottom
//	//	|| position == block_position::right_bottom)
//	//{
//	//	start_row = tilerows;
//	//	if (position == block_position::bottom)
//	//	{
//	//		row_first = true;
//	//	}
//	//}
//	//else
//	//{
//	//	start_row = center_row;
//	//}
//
//	if (position == block_position::left
//		|| position == block_position::left_bottom
//		|| position == block_position::left_top)
//	{
//		start_col = center_col;
//		start_col = 0.0;
//	}
//	else if (position == block_position::right
//		|| position == block_position::right_bottom
//		|| position == block_position::right_top)
//	{
//		start_col = tilecols;
//	}
//	else
//	{
//		start_col = center_col;
//		start_col = 0.0;
//	}
//
//	if (position == block_position::top
//		|| position == block_position::left_top
//		|| position == block_position::right_top)
//	{
//		start_row = center_row;
//		start_row = 0;
//	}
//	else if (position == block_position::bottom
//		|| position == block_position::left_bottom
//		|| position == block_position::right_bottom)
//	{
//		start_row = tilerows;
//	}
//	else
//	{
//		start_row = center_row;
//		start_row = 0;
//	}
//
//	//double weight = 2.0; 
//	//if (row_first)
//	//{
//	//	weight = tilecols*0.5;
//	//}
//	//else
//	//{
//	//	weight = tilerows*0.5;
//
//	//}
//	for (int i = 0; (i<tilerows); ++i)
//	{
//		for (int j = 0; (j<tilecols); ++j)
//		{
//			row_col_List.push_back(row_col(i - start_row, j - start_col));
//			//if (row_first)
//			//{
//			//	row_col_List.push_back(row_col((i - start_row)*weight, j - start_col));
//			//}
//			//else
//			//{
//			//	row_col_List.push_back(row_col((i - start_row)*weight, (j - start_col)*weight));
//
//			//}
//		}
//	}
//	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);
//	for (size_t i = 0; i < row_col_List.size(); i++)
//	{
//		//if (row_first)
//		//{
//		//	row_col_List[i].row_idx /= weight;
//		//}
//		//else
//		//{
//		//	row_col_List[i].col_idx /= weight;
//		//}
//		row_col_List[i].col_idx = int(row_col_List[i].col_idx + start_col);
//		row_col_List[i].row_idx = int(row_col_List[i].row_idx + start_row);
//		if (row_col_List[i].col_idx < 0)
//		{
//			row_col_List[i].col_idx = 0;
//		}
//		else if (row_col_List[i].col_idx >= tilecols)
//		{
//			row_col_List[i].col_idx = tilecols - 1;
//		}
//
//		if (row_col_List[i].row_idx < 0)
//		{
//			row_col_List[i].row_idx = 0;
//		}
//		else if (row_col_List[i].row_idx >= tilerows)
//		{
//			row_col_List[i].row_idx = tilerows - 1;
//		}
//	}
//	return row_col_List;
//}

class radiImageRegistration :
	public ossimOutputSource,
	public ossimProcessInterface
{
public:
   radiImageRegistration();
   radiImageRegistration(const radiImageRegistration& C);
   virtual ~radiImageRegistration();
   friend class radiMatchRectThread;
   enum point_type{
	   control = 0,
	   tie,
   };


   enum match_state{
	   success = 0,
	   slave_faild,
	   master_faild,
	   match_failed,
   };

   //accessors to parms
   inline void               setMaster(const ossimFilename& m) { theMaster=m; }
   inline const ossimFilename& getMaster()const { return theMaster; }
   inline void               setSlave(const ossimFilename& s) { theSlave=s; }
   inline const ossimFilename& getSlave()const { return theSlave; }
   inline void               setProjectionFile(const ossimFilename& s) { theProjectionFile=s; }
   inline const ossimFilename& getProjectionFile()const { return theProjectionFile; }
   //inline void               setMasterBand(ossim_uint32 b) { theMasterBand=b; }
   //inline ossim_uint32       getMasterBand()const { return theMasterBand; }
   //inline void               setSlaveBand(ossim_uint32 b) { theSlaveBand=b; }
   //inline ossim_uint32       getSlaveBand()const { return theSlaveBand; }
   inline void               setMasterBands(vector<ossim_uint32> bs) { theMasterBands = bs; }
   inline vector<ossim_uint32>       getMasterBands()const { return theMasterBands; }
   inline void               setSlaveBands(vector<ossim_uint32> bs) { theSlaveBands = bs; }
   inline vector<ossim_uint32>       getSlaveBands()const { return theSlaveBands; }
   inline void               setScaleRatio(const ossim_float64& r) { theScaleRatio=r; }
   inline ossim_float64      getScaleRatio()const { return theScaleRatio; }
   inline void               setSlaveAccuracy(const ossim_float64& a) { theSlaveAccuracy=a; }
   inline ossim_float64      getSlaveAccuracy()const { return theSlaveAccuracy; }
   inline void               setProjectionType(const ossimString& p) { theProjectionType=p; }
   inline const ossimString& getProjectionType()const { return theProjectionType; }
   inline void               setMasterPointProj(const ossimString& p) { theMasterPointProj=p; }
   inline const ossimString& getMasterPointProj()const { return theMasterPointProj; }
   inline void               setSlavePointProj(const ossimString& p) { theSlavePointProj=p; }
   inline const ossimString& getSlavePointProj()const { return theSlavePointProj; }
   inline void               setPointNumber(const ossim_uint32& n) { thePointNumber = n; }
   inline ossim_uint32       getPointNumber()const { return thePointNumber; }
   inline void               setTileSize(const ossim_float64& n) { theTileSize = n; }
   inline ossim_float64      getTileSize()const { return theTileSize; }
   inline void               setThreadNum(const ossim_uint32& nt) { theThreadNum = nt; }
   inline ossim_uint32       getThreadNum()const { return theThreadNum; }
   inline void               setSlaveId(const int& id) { theSlaveId=id; }
   inline const int& getSlaveId()const { return theSlaveId; }
   inline void               setMasterId(const int& id) { theMasterId=id; }
   inline const int& getMasterId()const { return theMasterId; }
   inline void               setPointType(const point_type& pt) { thePointType=pt; }
   inline const point_type& getPointType()const { return thePointType; }
   inline void               setUseGeographic(const bool& bGeographic) { theUseGeographic=bGeographic; }
   inline const bool& getUseGeographic()const { return theUseGeographic; }
   inline void               setDebug(const bool& bDebug) { theDebug=bDebug; }
   inline const bool& getDebug()const { return theDebug; }

   inline void               setSiftNfeatures(const int& n) { theSiftNfeatures=n; }
   inline int      getSiftNfeatures()const { return theSiftNfeatures; }   
   inline void               setSiftNOctaveLayers(const int& a) { theSiftNOctaveLayers=a; }
   inline int      getSiftNOctaveLayers()const { return theSiftNOctaveLayers; }
   inline ossim_float64       getSiftContrastThreshold()const { return theSiftContrastThreshold; }
   inline void               setSiftContrastThreshold(const ossim_float64& r) { theSiftContrastThreshold = r; }
   inline ossim_float64       getSiftEdgeThreshold()const { return theSiftEdgeThreshold; }
   inline void               setSiftEdgeThreshold(const ossim_float64& r) { theSiftEdgeThreshold = r; }
   inline ossim_float64       getSiftSigma()const { return theSiftSigma; }
   inline void               setSiftSigma(const ossim_float64& r) { theSiftSigma = r; }
   void setAreaOfInterest(const ossimIrect& rect);
   
   inline bool hasRun()const { return theHasRun; }

   // inherited methods
   virtual bool isOpen() const;
   virtual bool open();
   virtual void close();

   virtual bool  execute(); //also creates tie point file
   virtual       ossimObject* getObject()      { return this; }
   virtual const ossimObject* getObject()const { return this; }
   virtual       ossimObject* getObjectInterface() { return this; }
   
   virtual bool canConnectMyInputTo(ossim_int32 inputIndex,const ossimConnectableObject* object)const { return false; } //TBC : so far no input

   void setOutputName(const ossimString& filename);
   ossimIrect getMasterRect(ossimProjection* sProjection, const GdalRasterApp& mGdalApp,
	   ossimIrect sRect, double slaveAccuracy);
   ossimIrect getMasterRect(ossimProjection* sProjection, ossimProjection* mProjection,
	   ossimIrect sRect, double slaveAccuracy);
   ossimIrect getMasterRect(ossimIrect sRect, double slaveAccuracy);

   void appendTiePoints(const ossimFilename& filename);
   void appendControlPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection = NULL);   
   void writeTiePoints(const ossimFilename& filename);
   void writeControlPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection = NULL);
   void appendPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection = NULL);
   void writePoints(const ossimFilename& filename, ossimMapProjection* pMapProjection = NULL);
   void mergeTiePoints(radiBlockTieGptSet& totalTieGptSet, radiBlockTieGptSet* tieSet)
   {
	   totalTieGptSet.addImagesToList(tieSet->getImageList());
	   for (int i = 0; i < (int)tieSet->getTiePoints().size(); ++i)
	   {
		   totalTieGptSet.addTiePoint(tieSet->getTiePoints()[i]);
	   }
   }

   bool useOnline()
   {
	   if (0 == stricmp(theMaster.c_str(), strGoogle)
		   || 0 == stricmp(theMaster.c_str(), strBing)
		   || 0 == stricmp(theMaster.c_str(), strMapquest)
		   || 0 == stricmp(theMaster.c_str(), strMapbox))
	   {
		   return true;
	   }
	   return false;
   }

protected:  
	// loong
	void writeTiePoints(const radiBlockTieGptSet& tp);
	bool getMasterList(ossimFilename spatial_index_file, vector<ossimFilename>& masterList,
		ossimGpt ul, ossimGpt lr);
	inline bool getStoreFlag()const   { return theStoreFlag; }
	inline void setStoreFlag(bool sf) { theStoreFlag = sf; }
	bool getAllFeatures();
	bool getGridFeatures(const ossimIrect& rect);
	//int runMatch(const ossimIrect &srect, const ossimIrect &mrect, vector<ossimTDpt>& theTies, ossim_uint32 resLevel = 0);
	bool getGridFeaturesOnline(const image_block& block, radiBlockTieGptSet& tSet, const char* strRefSource);
	bool getGridFeaturesParallel(const ossimIrect& rect, void *pData);
	bool getGridFeaturesParallel(const image_block& rect, radiBlockTieGptSet& tSet);
	int runMatch(const cv::Mat& slaveMat, const cv::Mat& masterMat, ossimTDpt& tDpt, void *pData, bool bDebug = false);
	//int runMatchParallel1(const cv::Mat& slaveMat, const cv::Mat& masterMat, ossimTDpt& tDpt, void *pData, bool bDebug = false);
	int runMatchParallelNcc(const cv::Mat& slaveMat, const cv::Mat& masterMat, ossimTDpt& tDpt, void *pData, bool bDebug = false);
	bool createTileMat(const ossimRefPtr<ossimCastTileSourceFilter>& cast, const ossimIrect& rect, cv::Mat& outMat, ossim_uint32 resLevel = 0);
	ossimIpt slave2master(ossimProjection* slaveProjection,
		ossimProjection* masterProjection,
		ossimIpt slaveDpt);
	ossimDpt slave2master(ossimProjection* slaveProjection,
		ossimProjection* masterProjection,
		ossimDpt slaveDpt);
	ossimDpt slave2master(ossimDpt slaveDpt);
	ossimIpt slave2master(ossimIpt slaveDpt);
   ossimString         getRole() const;
   ossimImageHandler*  getProjectionHandler();
   ossimDrect imageRect2World(const ossimIrect& imageRect, ossimProjection* projection);
   ossimIrect worldRect2Image(const ossimDrect& imageRect, ossimProjection* projection);
   ossimDrect worldRectIntersection(const ossimDrect& r1, const ossimDrect& r2);
   
   ossimRefPtr<ossimImageGeometry> getOutputImageGeometry();

   ossimMapProjection* getOutputProjection();

   Eigen::VectorXd theGlobalAffineModel;
   
   bool buildRenderer(
      ossimImageChain* chain,
//         ossimImageSource* source, 
         ossimMapProjection* outProjection, 
         ossimImageRenderer* renderer,
         const ossimFilterResampler::ossimFilterResamplerType& stype =  ossimFilterResampler::ossimFilterResampler_CUBIC
         )const;

   ossimFilename     theLastMaster;
   bool              theStoreFlag;
#if OSSIM_HAS_MPI
   MPI_File          theFileStream;
#else
   fstream          theFileStream;
#endif
   ossimIrect        theAreaOfInterest;
   vector<ossimTDpt> theTiePoints;
   //ossim_uint32  theMasterBand;
   //ossim_uint32  theSlaveBand;
   vector<ossim_uint32>  theMasterBands;
   vector<ossim_uint32>  theSlaveBands;
   ossim_float64 theSlaveAccuracy;
   ossimString   theProjectionType;
   ossimString   theMasterPointProj;
   ossimString   theSlavePointProj;
   ossimProjection *theMasterProjection;
   ossimProjection *theSlaveProjection;
   ossim_float64 theScaleRatio;
   int   theSlaveId;
   int   theMasterId;
   point_type   thePointType;

   ossim_uint32  thePointNumber; 	// required number
   ossimFilename   theMaster;
   ossimFilename   theSlave;
   ossimFilename	  theFilename;
   ossimFilename   theProjectionFile;
   ossim_float64	  theTileSize;
   ossim_uint32	  theThreadNum;
   bool theUseGeographic;
   bool theDebug;
   int theSiftNfeatures;
   int theSiftNOctaveLayers;
   double theSiftContrastThreshold;
   double theSiftEdgeThreshold;
   double theSiftSigma;

   bool theHasRun; //to know whether execute has been run

   ossimRefPtr<ossimImageChain> theMChain;
   ossimRefPtr<ossimImageChain> theSChain;
   
   ossimRefPtr<ossimImageHandler>   handlerM;
   ossimRefPtr<ossimImageHandler>   handlerS;
   ossimRefPtr<ossimBandSelector>   theMasterBandSelector;
   ossimRefPtr<ossimBandSelector>   theSlaveBandSelector;
   ossimRefPtr<ossimImageRenderer>  rendererM;
   ossimRefPtr<ossimImageRenderer>  rendererS;
   std::vector<ossimRefPtr<ossimCastTileSourceFilter> > caster;
   radiBlockTieGptSet       theTset;

   //! Disallow operator=
   const radiImageRegistration& operator=(const radiImageRegistration& rhs) {return rhs;}

TYPE_DATA
};


////-----------------------------------------------------------------------------
//// Custom Thread Observer (w/finished count)
////
//class MyThreadObserver : public ThreadObserver {
//
//public:
//
//	MyThreadObserver() : ThreadObserver(), _finishedCount(0) {};
//
//	virtual ~MyThreadObserver() {};
//
//	void threadFinished(const int threadId) {
//
//		ThreadObserver::threadFinished(threadId);
//
//		++_finishedCount;
//	}
//
//	int getFinishedCount() {return _finishedCount;};
//
//
//private:
//
//	volatile int _finishedCount;
//
//};
//
//// check the working of OpenThreads::Thread::CurrentThread()
//static OpenThreads::Thread* CurrentChecker(){
//	return OpenThreads::Thread::CurrentThread();
//};
//
//class radiMatchRectThread: public OpenThreads::Thread, public ThreadReporter
//{
//public:
//	radiMatchRectThread(radiImageRegistration *pRegistration)
//		: OpenThreads::Thread(),
//		m_pRegistration(pRegistration){};
//	void setRect(const vector<ossimIrect>& rectList){m_rectList = rectList;};
//	virtual  ~radiMatchRectThread()
//	{
//		//m_pRegistration = NULL;
//	};
//	virtual void run()
//	{
//		if( CurrentChecker() != this)
//			notifyObserversMessage(getThreadId(), "Thread::CurrentThread is NOT working");
//		else
//			notifyObserversMessage(getThreadId(), "Thread::CurrentThread is working");
//
//		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.
//
//		for (int i = 0;i < (int)m_rectList.size();++i)
//		{
//			m_pRegistration->getGridFeaturesParallel(m_rectList[i], m_theTset);
//			_quitmutex.lock();
//			finishedBlocks++;
//			m_pRegistration->setPercentComplete((finishedBlocks)/(double)totalBlocks*100.0);
//			//m_Registration.setPercent((finishedBlocks)/(double)totalBlocks*100.0);
//			fflush( stdout );
//			_quitmutex.unlock();
//		}
//
//		notifyObserversMessage(getThreadId(), "Finished my work");
//
//		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.
//
//		//---------------------------------------------------------------------
//		// Now that we've done our work, wait for a sign that we should quit.
//		//
//		while (true) {
//
//			_quitmutex.lock();
//			if(_quitflag == true) break;
//			_quitmutex.unlock();
//
//			OpenThreads::Thread::YieldCurrentThread();
//		}
//
//
//		notifyObserversFinished(getThreadId());
//	};
//
//	void mergeTiePoints(radiBlockTieGpt& totalTieGptSet)
//	{
//		for (int i = 0;i < (int)m_theTset.getTiePoints().size();++i)
//		{
//			totalTieGptSet.addTiePoint(m_theTset.getTiePoints()[i]);
//		}
//	}
//	void quit() {
//		_quitmutex.lock();
//		_quitflag = true;
//		_quitmutex.unlock();
//	};
//	radiBlockTieGpt m_theTset;
//private:
//	bool getGridFeaturesParallel(const ossimIrect& rect);
//	radiImageRegistration *m_pRegistration;
//	vector<ossimIrect> m_rectList;
//	int *_dataPtr;
//	int _numElts;
//	volatile bool _quitflag;
//	OpenThreads::Mutex _quitmutex;
//};


class radiMatchRectThread: public OpenThreads::Thread
{
public:
	radiMatchRectThread(radiImageRegistration *pRegistration)
		: OpenThreads::Thread(),
		m_pRegistration(pRegistration){};
	void setBlocks(const vector<image_block>& blockList){ m_blockList = blockList; };
	virtual  ~radiMatchRectThread()
	{
		//m_pRegistration = NULL;
	};
	virtual void run()
	{
		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.
		int rank, nTasks;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

		for (int i = 0; i < (int)m_blockList.size(); ++i)
		{
			if (0 == stricmp(m_RefSource, strGoogle)
				|| 0 == stricmp(m_RefSource, strBing)
				|| 0 == stricmp(m_RefSource, strMapquest)
				|| 0 == stricmp(m_RefSource, strMapbox))
			{
				m_pRegistration->getGridFeaturesOnline(m_blockList[i], m_theTset, m_RefSource);
			}
			else
			{
				m_pRegistration->getGridFeaturesParallel(m_blockList[i], m_theTset);
			}
			_quitmutex.lock();
			finishedBlocks++;
			//m_pRegistration->setPercentComplete((finishedBlocks)/(double)totalBlocks*100.0);
			//int tatalFinishedBlocks = 0;
			//MPI_Reduce(&finishedBlocks, &tatalFinishedBlocks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			//m_Registration.setPercent((finishedBlocks)/(double)totalBlocks*100.0);
			//int p = (int)(tatalFinishedBlocks*100.0 / (double)totalBlocks + 0.5);
			//int p = (int)(finishedBlocks*100.0 / (double)totalBlocks + 0.5);
			//MPI_Allreduce(&finishedBlocks, &tatalFinishedBlocks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			//MPI_Reduce(&finishedBlocks, &tatalFinishedBlocks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			//int p = (int)(tatalFinishedBlocks*100.0 / (double)totalBlocks + 0.5);
			//int p = (int)(finishedBlocks*nTasks*100.0 / (double)m_totalBlocks + 0.5);
			int p = (int)(finishedBlocks*100.0 / (double)m_totalBlocks + 0.5);
			p = min(p, 100);
			if (0 == rank)
			{
				//MPI_Reduce(&finishedBlocks, &tatalFinishedBlocks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
				//m_Registration.setPercent((finishedBlocks)/(double)totalBlocks*100.0);
				//int p = (int)(finishedBlocks*100.0 / (double)totalBlocks + 0.5);
				printf("\r%d%%", p);
				fflush(stdout);
			}
			_quitmutex.unlock();
		}

		bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.

		//---------------------------------------------------------------------
		// Now that we've done our work, wait for a sign that we should quit.
		//
		while (true) {

			_quitmutex.lock();
			if(_quitflag == true) break;
			_quitmutex.unlock();

			OpenThreads::Thread::YieldCurrentThread();
		}
	};

	void mergeTiePoints(radiBlockTieGptSet& totalTieGptSet)
	{
		totalTieGptSet.addImagesToList(m_theTset.getImageList());
		for (int i = 0;i < (int)m_theTset.getTiePoints().size();++i)
		{
			totalTieGptSet.addTiePoint(m_theTset.getTiePoints()[i]);
		}
	}
	void quit() {
		_quitmutex.lock();
		_quitflag = true;
		_quitmutex.unlock();
	};
	radiBlockTieGptSet m_theTset;
	int m_totalBlocks;
	const char* m_RefSource;
private:
	bool getGridFeaturesParallel(const image_block& block);
	radiImageRegistration *m_pRegistration;
	vector<image_block> m_blockList;
	int *_dataPtr;
	int _numElts;
	volatile bool _quitflag;
	OpenThreads::Mutex _quitmutex;
};

#endif //radiImageRegistration_HEADER
