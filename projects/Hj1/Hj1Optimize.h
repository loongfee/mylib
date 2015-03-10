#pragma once
#include <mprojectdefine.h>
#include <ossim/base/ossimPreferences.h>
using namespace std;
using namespace mylib;
using namespace ossimplugins;

namespace Hj1Optimize{
	struct HJ1OptStruct
	{
		ossimHj1Model* hj1Model;
		ossimTieGptSet* gcpSet;

		HJ1OptStruct()
			:hj1Model(NULL),
			gcpSet(NULL)
		{

		}

		HJ1OptStruct(ossimHj1Model* ahj1Model, ossimTieGptSet* agcpSet)
			:hj1Model(ahj1Model)
		{
			gcpSet = new ossimTieGptSet;
			*gcpSet = *agcpSet;
		}
	};
}