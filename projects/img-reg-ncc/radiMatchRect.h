#ifndef RADI_MATCH_RECT_H
#define RADI_MATCH_RECT_H
#include <OpenThreads/Thread>
#include <OpenThreads/Mutex>
#include <OpenThreads/Barrier>
#include "radiImageRegistration.h"'
#include "ThreadReporter.h"

#ifdef _WIN32
#include <process.h>
#define getpid() _getpid()
#else
#include <unistd.h>
#endif 

static OpenThreads::Barrier bar;

class radiMatchRect: public OpenThreads::Thread, public ThreadReporter
{
public:
	radiMatchRect(radiImageRegistration* pRegistration)
		: OpenThreads::Thread(),
		m_pRegistration(pRegistration){};
	void setRect(const vector<ossimIrect>& rectList){m_rectList = rectList;};
	virtual  ~radiMatchRect();
	virtual void run();
	void quit() {
		_quitmutex.lock();
		_quitflag = true;
		_quitmutex.unlock();
	};
private:
	bool getGridFeaturesParallel(const ossimIrect& rect);
	radiImageRegistration* m_pRegistration;
	vector<ossimIrect> m_rectList;
	int *_dataPtr;
	int _numElts;
	volatile bool _quitflag;
	OpenThreads::Mutex _quitmutex;
};

#endif