#include "QtRspfListener.h"

#include <QtGui>
#include <QApplication>

QtRspfListener::QtRspfListener()
{
}

QtRspfListener::~QtRspfListener()
{
}
void QtRspfListener::processEvent(ossimEvent& event)
{
	switch(event.getId())
	{
	case OSSIM_EVENT_PROCESS_PROGRESS_ID:
		{
			ossimProcessProgressEvent* eventCast = static_cast<ossimProcessProgressEvent*>(&event);
			processProgressEvent(*eventCast);
			break;
		}
	default:
		{
			ossimListener::processEvent(event);
			break;
		}
	}
}
void QtRspfListener::processProgressEvent(ossimProcessProgressEvent& event)
{
	emit updateProgress(event.getPercentComplete());
}