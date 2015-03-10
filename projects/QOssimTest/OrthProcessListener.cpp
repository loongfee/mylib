#include "OrthProcessListener.h"


OrthProcessListener::OrthProcessListener():ossimProcessListener()
{}
void OrthProcessListener::processEvent(ossimEvent& event)
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
void OrthProcessListener::processProgressEvent(ossimProcessProgressEvent& event)
{
	percentComplete = event.getPercentComplete();
}

double OrthProcessListener::getPercentComplete()
{
	return percentComplete;
}