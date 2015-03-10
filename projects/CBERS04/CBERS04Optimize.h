#pragma once
#include <mprojectdefine.h>
#include <ossim/base/ossimPreferences.h>
#include <ossim_plugin/radi/radiCbers04Model.h>
using namespace std;
using namespace mylib;
using namespace ossimplugins;

namespace Cbers04Optimize{
	struct Cbers04OptStruct
	{
		radiCbers04Model* cbers04Model;
		ossimTieGptSet* gcpSet;

		Cbers04OptStruct()
			:cbers04Model(NULL),
			gcpSet(NULL)
		{

		}

		Cbers04OptStruct(radiCbers04Model* acbers04Model, ossimTieGptSet* agcpSet)
			:cbers04Model(acbers04Model)
		{
			gcpSet = new ossimTieGptSet;
			*gcpSet = *agcpSet;
		}
	};
}