#include "radiMatchRect.h"

static int GLOBAL_NUM_THREADS;

//-----------------------------------------------------------------------------
// Custom Thread Observer (w/finished count)
//
class MyThreadObserver : public ThreadObserver {

public:

	MyThreadObserver() : ThreadObserver(), _finishedCount(0) {};

	virtual ~MyThreadObserver() {};

	void threadFinished(const int threadId) {

		ThreadObserver::threadFinished(threadId);

		++_finishedCount;
	}

	int getFinishedCount() {return _finishedCount;};


private:

	volatile int _finishedCount;

};

// check the working of OpenThreads::Thread::CurrentThread()
static OpenThreads::Thread* CurrentChecker(){
	return OpenThreads::Thread::CurrentThread();
};

radiMatchRect::radiMatchRect(radiImageRegistration* pRegistration)
	: OpenThreads::Thread(),
	m_pRegistration(pRegistration)
{

}

radiMatchRect::~radiMatchRect()
{
	m_pRegistration = NULL;
}

void radiMatchRect::run() {

	//if( CurrentChecker()!=this)
	//	notifyObserversMessage(getThreadId(), "Thread::CurrentThread is NOT working");
	//else
	//	notifyObserversMessage(getThreadId(), "Thread::CurrentThread is working");


	//bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.

	//char tmp[80];
	//sprintf(tmp, "StackSize: %d\n", static_cast<int>(getStackSize()));

	//notifyObserversStarted(getThreadId());
	//notifyObserversMessage(getThreadId(), "This is a thread message.");
	//notifyObserversMessage(getThreadId(), tmp);

	//register int i;
	//for (i=0; i<_numElts; ++i) {
	//	_dataPtr[i] = getThreadId();
	//}

	//notifyObserversMessage(getThreadId(), "Finished my work");

	//bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.

	////---------------------------------------------------------------------
	//// Now that we've done our work, wait for a sign that we should quit.
	////
	//while (true) {

	//	_quitmutex.lock();
	//	if(_quitflag == true) break;
	//	_quitmutex.unlock();

	//	OpenThreads::Thread::YieldCurrentThread();
	//}


	//notifyObserversFinished(getThreadId());

	if( CurrentChecker() != this)
		notifyObserversMessage(getThreadId(), "Thread::CurrentThread is NOT working");
	else
		notifyObserversMessage(getThreadId(), "Thread::CurrentThread is working");

	OpenThreads::Thread::
	bar.block(GLOBAL_NUM_THREADS);  // Sync the threads.

	for (int i = 0;i < (int)m_rectList.size();++i)
	{
		getGridFeaturesParallel(m_rectList[i]);
	}

	notifyObserversMessage(getThreadId(), "Finished my work");

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


	notifyObserversFinished(getThreadId());
}

bool radiMatchRect::getGridFeaturesParallel(const ossimIrect& rect)
{
	if (!m_pRegistration->theSlaveBandSelector)
	{
		ossimNotify(ossimNotifyLevel_WARN)
			<< "WARN ossimTieGenerator::scanForEdges():"
			<< "\nInput source is not a ossimImageChip.  Returning..." << std::endl;
		return false;
	}

	const ossim_int32 TILE_HEIGHT	= m_pRegistration->theTileSize;
	const ossim_int32 TILE_WIDTH	= m_pRegistration->theTileSize;
	const ossim_int32 START_LINE = rect.ul().y;
	const ossim_int32 STOP_LINE  = rect.lr().y;
	const ossim_int32 START_SAMP = rect.ul().x;
	const ossim_int32 STOP_SAMP  = rect.lr().x;
	int nWidth = STOP_SAMP-START_SAMP;
	int nHeight = STOP_LINE-START_LINE;

	// For percent complete status.
	ossim_int32 tilerows=(STOP_LINE-START_LINE+TILE_HEIGHT) / TILE_HEIGHT; //ceil : (stop-start+1+size-1)/size
	ossim_int32 tilecols=(STOP_SAMP-START_SAMP+TILE_WIDTH) / TILE_WIDTH;
	double total_tiles = ((double)tilerows)*tilecols;
	double tiles_processed = 0.0;

	// loop through all tiles
	// need to use a sequencer for parallelism in the future TBD
	ossim_int32 line=START_LINE;
	ossim_int32 i,j;

	vector<row_col> row_col_List;
	double center_row = tilerows * 0.5;
	double center_col = tilecols * 0.5;
	for (i=0;(i<tilerows);++i)
	{
		for (j=0;(j<tilecols);++j )
		{
			row_col_List.push_back(row_col((center_row-i-1),(center_col-j-1)));
		}
	}

	GdalRasterApp slaveApp;
	if (!slaveApp.open(m_pRegistration->getSlave().c_str()))
	{
		cerr<<"radiImageRegistration"<<"::execute can't open master image  "<< m_pRegistration->getSlave().c_str() <<endl;
		return false;
	}
	// select only one band (if multiple)
	ossim_uint32 sbc = slaveApp.nBand();

	//add a band selector
	ossim_uint32 sb = m_pRegistration->theSlaveBand;
	if (sb>=sbc) 
	{
		cerr<<"radiImageRegistration"<<"::execute Warning not enough bands in master, only "<< sbc <<endl;
		sb=0;
	}

	double sGSD = m_pRegistration->theSlaveProjection->getMetersPerPixel().x;
	//int ncore = omp_get_num_procs();//获取执行核的总数；  目前机器CPU的数量
	row_col_step = 5;
	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);	

	//GdalRasterApp slaveApp;
	//slaveApp.open(theSlave.c_str());
	bool bDebug = true;
	//if (theTset.size() == 20)
	//{
	//	bDebug = true;
	//}
	bool found = false;
	int N_PARALLEL = 1;//ncore*2;
	int search_size = (int)row_col_List.size();
	search_size = min(search_size, 1);
	for (int i=0;i < search_size && !found;)
		//for (int i=0;i < (int)row_col_List.size() && !found;)
	{
		int nParallel = min(N_PARALLEL, (int)row_col_List.size()-1-i);	// 保证末尾不越界

		std::vector<cv::Mat> slaveMatList;
		std::vector<ossimIrect> srectList;
		std::vector<bool>  slaveValidList(nParallel, false);
		for (int j = 0; j < nParallel;j++)
		{
			int icol = floor(row_col_List[i+j].col_idx+center_col+0.5);
			int irow = floor(row_col_List[i+j].row_idx+center_row+0.5);
			//int ii = (5 * i) % (int)row_col_List.size();
			//int icol = floor(row_col_List[ii+j].col_idx+center_col+0.5);
			//int irow = floor(row_col_List[ii+j].row_idx+center_row+0.5);
			ossim_int32 samp=START_SAMP+icol*TILE_WIDTH;
			ossim_int32 line=START_LINE+irow*TILE_HEIGHT;


			ossim_int32 BufHeight = TILE_HEIGHT;
			ossim_int32 BufWidth = TILE_WIDTH;
			//行末尾小块处理
			if (irow == tilerows-1)
			{
				BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;//得出当前块的宽度Bufsizex，高度Bufsizey
				BufHeight = min(BufHeight, TILE_HEIGHT);
			}
			//列末尾小块处理
			if (icol == tilecols-1)
			{
				BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;//得出当前块的宽度Bufsizex，高度Bufsizey
				BufWidth = min(BufWidth, TILE_WIDTH);
			}

			// slave
			cv::Mat slaveMat;
			ossimIrect srect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1));
			//if (slaveApp.getPrincipalRect2CvMatByte( srect, slaveMat))
			if (slaveApp.getRect2CvMatByte( srect, slaveMat, m_pRegistration->theSlaveBand))
			{
				slaveValidList[j] = true;
			}
			{
				//if (createTileMat(caster[0], srect, slaveMat, 0))
				//{
				//	slaveValidList[j] = true;
				//}
			}
			slaveMatList.push_back(slaveMat);
			srectList.push_back(srect);
		}
		for (int j = 0; j < nParallel;j++)
		{
			if (found || !slaveValidList[j])
			{
				continue;
			}
			ossimIpt sul = srectList[j].ul();
			ossimIpt slr = srectList[j].lr();
			ossimIpt delta_lr((ossim_int32)(ceil(m_pRegistration->theSlaveAccuracy)), (ossim_int32)(ceil(m_pRegistration->theSlaveAccuracy)) );

			vector<ossimTDpt> tp(1);
			// search in the reference library
			ossimGpt ul_latlon, lr_latlon;
			m_pRegistration->theSlaveProjection->lineSampleToWorld(sul-delta_lr, ul_latlon);
			m_pRegistration->theSlaveProjection->lineSampleToWorld(slr+delta_lr, lr_latlon);
			//ossimDpt tempDpt;
			//slaveApp.linesample2lonlat(sul-delta_lr, tempDpt);
			//ul_latlon = ossimGpt(tempDpt.y, tempDpt.x);
			//slaveApp.linesample2lonlat(slr+delta_lr, tempDpt);
			//lr_latlon = ossimGpt(tempDpt.y, tempDpt.x);
			vector<ossimFilename> masterFileList;
			if (m_pRegistration->theMaster.ext().upcase() == "SHP")
			{
				m_pRegistration->getMasterList(m_pRegistration->theMaster, masterFileList, ul_latlon, lr_latlon);
			}
			else{
				masterFileList.clear();
				masterFileList.push_back(m_pRegistration->theMaster);
			}

			//ossimRefPtr<ossimImageHandler> master_handler = NULL;
			//ossimProjection* master_projection = NULL;
			//ossimRefPtr<ossimBandSelector> master_bandselector = NULL;
			cv::Mat masterMat;
			//ossimRefPtr<ossimCastTileSourceFilter> master_caster = new ossimCastTileSourceFilter();
			//master_caster->setOutputScalarType(OSSIM_FLOAT64);
			//master_caster->setOutputScalarType(OSSIM_UCHAR);

			int match_result = m_pRegistration->match_state::success;
			for (int iFile = 0;iFile < (int)masterFileList.size();++iFile)
			{
				if (found || match_result == m_pRegistration->match_state::slave_faild)
				{
					continue;
				}

				ossimFilename lastMaster = masterFileList[iFile];
				m_pRegistration->theLastMaster = lastMaster;
				GdalRasterApp masterApp;
				if (!masterApp.open(lastMaster.c_str()))
				{
					cerr<<"radiImageRegistration"<<"::execute can't open master image  "<< lastMaster <<endl;
					continue;
				}
				// select only one band (if multiple)
				ossim_uint32 mbc = masterApp.nBand();

				//add a band selector
				ossim_uint32 mb = m_pRegistration->theMasterBand;
				if (mb>=mbc) 
				{
					//cerr<<"radiImageRegistration"<<"::execute Warning not enough bands in master, only "<< mbc <<endl;
					mb=0;
				}

				double gsd_scale;
				// master
				ossimIrect mrect;
				if (0 == strcmp(masterApp.getGetProjectionRef(), ""))
				{
					// 无投影
					mrect = ossimIrect(sul - delta_lr, slr + delta_lr);
					gsd_scale = 1.0;
				}
				else
				{
					double mGSD = masterApp.getGeoTransform()[1];
					gsd_scale = mGSD / sGSD;
					// master
					mrect = m_pRegistration->getMasterRect(m_pRegistration->theSlaveProjection, masterApp,
						srectList[j], m_pRegistration->theSlaveAccuracy);

				}
				if (!masterApp.getRect2CvMatByte( mrect, masterMat, mb, gsd_scale))
				{
					masterApp.close();
					continue;
				}
				match_result = m_pRegistration->runMatchParallel(slaveMatList[j], masterMat, tp[0], this, bDebug);
				if (match_result == m_pRegistration->match_state::slave_faild)
				{
					masterApp.close();
					continue;
				}
				if (m_pRegistration->theFilename != ossimFilename::NIL && m_pRegistration->match_state::success == match_result)
				{
					//ossimTieGpt tiePt;
					//tiePt.setImagePoint(ossimDpt(tp[0].getSlavePoint() + sul));
					//ossimGpt gpt;
					//masterApp.linesample2lonlat(tp[0].getMasterPoint() + mul, gpt);
					//tiePt.setGroundPoint(gpt);
					//theTset.addTiePoint(&tiePt);
					tp[0].setMasterPoint(tp[0].getMasterPoint() + mrect.ul());
					tp[0].setSlavePoint(tp[0].getSlavePoint() + sul);

					if (0 == strcmp(masterApp.getGetProjectionRef(), ""))
					{
						// 无投影
						ossimRefPtr<ossimTieGpt> tgi(new ossimTieGpt);
						ossimDpt dpt = tp[0].getMasterPoint();
						tgi->setGroundPoint(ossimGpt(dpt.y, dpt.x, 0.0));
						//set slave image position
						tgi->refImagePoint() = tp[0].getSlavePoint();
						//set score
						tgi->setScore(tp[0].score);

						//add to list
						m_pRegistration->theTset.addTiePoint(tgi);
					}
					else
					{
						// convert "Image to Image" to "Ground to Image" tie points    //TBC : use more generic tie points
						ossimRefPtr<ossimTieGpt> tgi(new ossimTieGpt);
						//set master ground pos
						ossimDpt lonlat;
						masterApp.linesample2lonlat(tp[0].getMasterPoint(), lonlat);
						double hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(ossimGpt(lonlat.y, lonlat.x));
						//ossimDpt dpt = outProj->forward(*tgi);
						tgi->setGroundPoint(ossimGpt(lonlat.y, lonlat.x, hgt));
						//set slave image position
						tgi->refImagePoint() = tp[0].getSlavePoint();
						//set score
						tgi->setScore(tp[0].score);
						//add to list
						m_pRegistration->theTset.addTiePoint(tgi);
					}
				}
				if (m_pRegistration->match_state::success == match_result)
				{
					masterApp.close();
					found = true;
					continue;
				}
				masterApp.close();
			}
		}
		slaveMatList.clear();
		slaveValidList.clear();
		srectList.clear();
		i += N_PARALLEL;
	}
	slaveApp.close();
	return found;
}