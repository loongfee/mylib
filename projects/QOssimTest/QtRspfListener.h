#pragma once
#include <ossim/base/ossimStdOutProgress.h>

#include <QThread> 
class QtRspfListener : public QThread,
	public ossimProcessListener
{
	Q_OBJECT
signals:
	void updateProgress(double);
public:   
	QtRspfListener();   
	~QtRspfListener();
	ossimProcessProgressEvent  GaugeProgressEvent;
	virtual void processEvent(ossimEvent& event);
	virtual void processProgressEvent(ossimProcessProgressEvent& event);
};
