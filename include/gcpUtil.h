#pragma once

#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimTempFilename.h>
#include <ossim/base/ossimString.h>
//#include <ossim/base/ossimXmlNode.h>
#include <ossim/base/ossimDirectory.h>
#include <ossim/base/ossimDirectoryTree.h>
//#include <ossim/base/ossimXmlDocument.h>
#include <ossim/base/ossimGpt.h>
#include <ossim/base/ossimDpt.h>
#include <ossim/base/ossimKeywordNames.h>
#include <ossim/base/ossimTieGptSet.h> 
#include <ossim/elevation/ossimElevManager.h>
#include <ossim/imaging/ossimImageHandler.h>
#include <ossim/imaging/ossimImageHandlerRegistry.h>
#include <ossim/projection/ossimProjection.h>
#include <ossim/projection/ossimUtmProjection.h>
#include <ossim/projection/ossimTransMercatorProjection.h>
#include <ossim/projection/ossimMapProjectionFactory.h>
#include <ossim/projection/ossimProjectionFactoryRegistry.h>
#include <ossim/init/ossimInit.h>
#include <ossim/base/ossimIoStream.h>
#include <ossim/base/ossimCommon.h>
#include <ossim/base/ossimRefPtr.h>
#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimRtti.h>

#include <ossim/base/ossimXmlDocument.h>

#include <string>
#include <vector>
#include <iostream>
#include <direct.h>
using namespace std;

#include "strUtil.h"

namespace mylib{

void readTextGcpFile(ossimFilename strFilename,
				 ossimTieGptSet* &gcpSet,
				 ossimTieGptSet* &chkSet,
				 ossimKeywordlist* prjKwl = NULL);

void readXmlGcpFile(ossimFilename strFilename,
						ossimTieGptSet* &gcpSet,
						ossimTieGptSet* &chkSet,
						ossimKeywordlist* prjKwl = NULL);

void readGcpFile(ossimFilename strFilename,
						ossimTieGptSet* &gcpSet,
						ossimTieGptSet* &chkSet,
						ossimKeywordlist* prjKwl = NULL);

void saveGcpFile(ossimFilename filenametosave,
					 ossimTieGptSet* gcpSet, 
					 ossimTieGptSet* chkSet = NULL, 
					 ossimKeywordlist* prjKwl = NULL,
					 bool extern_file = true);

bool get_elevation(ossimTieGptSet* &ctrlSet, ossimKeywordlist prjKwl, double defaultElev = 0.0);
bool get_elevation(ossimTieGptSet* &ctrlSet, double defaultElev = 0.0);

void get_elevation(ossimFilename infile, ossimFilename outfile, ossimFilename elevationPath = "", double defaultElev = 0.0);
ossimTieGptSet* datum_shift(ossimTieGptSet* inGcpSet, const char* outDatumString);
void batch_get_elevation(ossimFilename inDir, ossimFilename outDir, 
						 ossimFilename filter="*.txt", 
						 ossimFilename elevationPath="",
						 double defaultElev = 0.0);


ossimTieGptSet* projection2ll(ossimTieGptSet* inGcpSet, ossimKeywordlist prjKwl);

bool projection2ll(ossimFilename inFile, ossimFilename outFile, bool bGetElevation = false);

ossimTieGptSet* reprojectionPoints(ossimTieGptSet* inGcpSet, ossimKeywordlist inPrjKwl, ossimKeywordlist outPrjKwl);

bool reprojectionPoints(ossimFilename inFile, ossimKeywordlist outPrjKwl, ossimFilename outFile);

}; // end of namespace mylib