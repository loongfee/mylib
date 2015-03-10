#ifndef _ORTH_PROCESS_LISTENER_H_
#define _ORTH_PROCESS_LISTENER_H_
#include <ossim/base/ossimStdOutProgress.h>

class OrthProcessListener: public ossimProcessListener
{
public:
	ossimProcessProgressEvent  GaugeProgressEvent;
	double percentComplete;
public:
	OrthProcessListener();
	virtual void processEvent(ossimEvent& event);
	virtual void processProgressEvent(ossimProcessProgressEvent& event);
	double getPercentComplete();
};

#endif