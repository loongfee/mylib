#ifndef _STDOUTPROGRESS_H
#define _STDOUTPROGRESS_H
#include <ossim/base/ossimStdOutProgress.h>

class OrthProcessListener :
	public ossimProcessListener
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