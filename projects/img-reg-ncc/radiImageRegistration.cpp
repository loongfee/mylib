// radiImageRegistration.cpp
#include <iostream>

#include "radiImageRegistration.h"
#include <ossim/imaging/ossimImageSource.h>
#include <ossim/imaging/ossimImageHandler.h>
#include <ossim/imaging/ossimCastTileSourceFilter.h>
#include <ossim/imaging/ossimImageRenderer.h>
#include <ossim/imaging/ossimCacheTileSource.h>
#include <ossim/imaging/ossimImageChain.h>
#include <ossim/imaging/ossimBandSelector.h>
#include <ossim/imaging/ossimImageHandlerRegistry.h>
#include <ossim/projection/ossimMapProjection.h>
#include <ossim/projection/ossimProjectionFactoryRegistry.h>
#include <ossim/projection/ossimImageViewProjectionTransform.h>
#include <ossim/base/ossimTrace.h>
#include <ossim/base/ossimXmlDocument.h>
#include <ossim/base/ossimFilenameProperty.h>
#include <ossim/base/ossimStringProperty.h>
#include <ossim/imaging/ossimImageGeometry.h>
#include <ossim/projection/ossimProjection.h>
#include <ossim/base/ossimXmlDocument.h>

extern "C"{
#include <vl/generic.h>
#include <vl/stringop.h>
#include <vl/pgm.h>
#include <vl/sift.h>
#include <vl/getopt_long.h>
};

#include <ctime>
#include <iostream>

#include <QString>
#include <QDir>
#include <QXmlStreamWriter>
#include <QMapIterator>
using namespace std;

//#include "radiMatchRect.h"

//#pragma comment(lib, "opencv_core248.lib")
//#pragma comment(lib, "opencv_features2d248.lib")
//#pragma comment(lib, "opencv_flann248.lib")
//#pragma comment(lib, "opencv_highgui248.lib")
//#pragma comment(lib, "opencv_nonfree248.lib")
#pragma comment(lib, "vl.lib")

RTTI_DEF2(radiImageRegistration, "radiImageRegistration", ossimOutputSource, ossimProcessInterface);

double radiImageRegistration::theLMS[6 * 9] = {
	-1.1111111111111116e-001, 2.2222222222222210e-001, -1.1111111111111116e-001, 2.2222222222222210e-001, 5.5555555555555536e-001, 2.2222222222222210e-001, -1.1111111111111116e-001, 2.2222222222222210e-001, -1.1111111111111116e-001,
	-1.6666666666666666e-001, 0.0000000000000000e+000, 1.6666666666666666e-001, -1.6666666666666666e-001, 0.0000000000000000e+000, 1.6666666666666666e-001, -1.6666666666666666e-001, 0.0000000000000000e+000, 1.6666666666666666e-001,
	-1.6666666666666666e-001, -1.6666666666666666e-001, -1.6666666666666666e-001, 0.0000000000000000e+000, 0.0000000000000000e+000, 0.0000000000000000e+000, 1.6666666666666666e-001, 1.6666666666666666e-001, 1.6666666666666666e-001,
	2.5000000000000000e-001, 0.0000000000000000e+000, -2.5000000000000000e-001, 0.0000000000000000e+000, 0.0000000000000000e+000, 0.0000000000000000e+000, -2.5000000000000000e-001, 0.0000000000000000e+000, 2.5000000000000000e-001,
	1.6666666666666669e-001, -3.3333333333333331e-001, 1.6666666666666669e-001, 1.6666666666666674e-001, -3.3333333333333326e-001, 1.6666666666666674e-001, 1.6666666666666669e-001, -3.3333333333333331e-001, 1.6666666666666669e-001,
	1.6666666666666669e-001, 1.6666666666666674e-001, 1.6666666666666669e-001, -3.3333333333333331e-001, -3.3333333333333326e-001, -3.3333333333333331e-001, 1.6666666666666669e-001, 1.6666666666666674e-001, 1.6666666666666669e-001
};

radiImageRegistration::radiImageRegistration()
	: ossimOutputSource(NULL, // owner
	0,
	0,
	true,
	true),
 ossimProcessInterface(),
   theMaster(),
   theSlave(),
   theMasterBand(0),
   theSlaveBand(0),
   theScaleRatio(1.0),
   theSlaveAccuracy(50.0),
   thePointNumber(25),
   theSiftNfeatures(0),
   theSiftNOctaveLayers(3),
   theSiftContrastThreshold(0.04),
   theSiftEdgeThreshold(10),
   theSiftSigma(1.6),
   theProjectionType("M"),
   theMasterPointProj("G"),
   theSlavePointProj("I"),
   theHasRun(false),
   handlerM(NULL),
   handlerS(NULL),
   theMasterBandSelector(NULL),
   theSlaveBandSelector(NULL),
   theTileSize(256),
   theSlaveId(1),
   theMasterId(2),
   thePointType(point_type::control),
   theUseGeographic(false),
   theDebug(true),
   theThreadNum(0),
   theTset(),
   theMinNCC(0.75),
   theStoreFlag(false)
{
   //default output name : XML
   //1.add all OSSIM blocks
   // ingredients :
   // 2x ossimImageSource (for master & slave)
   // 2x ossimImagerRenderer
   // 1x ossimHarrisCorners
   // 2x ossimCastTileFilter (to get floating point
   // 1x ossimChipMatch
   // 1x ossimTieGenerator

   // don't create sources (factories will do it)
   caster.push_back(new ossimCastTileSourceFilter());
   caster.push_back(new ossimCastTileSourceFilter());  
   cornerDetector      = new ossimHarrisCorners();
#if OSSIM_HAS_MPI
   theFileStream = NULL;
#endif
}

radiImageRegistration::radiImageRegistration(const radiImageRegistration& a)
	: ossimOutputSource(NULL, // owner
	0,
	0,
	true,
	true),
	ossimProcessInterface(),
	theMaster(a.theMaster),
	theSlave(a.theSlave),
	theMasterBand(a.theMasterBand),
	theSlaveBand(a.theSlaveBand),
	theScaleRatio(a.theScaleRatio),
	theSlaveAccuracy(a.theSlaveAccuracy),
	thePointNumber(a.thePointNumber),
	theSiftNfeatures(a.theSiftNfeatures),
	theSiftNOctaveLayers(a.theSiftNOctaveLayers),
	theSiftContrastThreshold(a.theSiftContrastThreshold),
	theSiftEdgeThreshold(a.theSiftEdgeThreshold),
	theSiftSigma(a.theSiftSigma),
	theProjectionType(a.theProjectionType),
	theMasterPointProj(a.theMasterPointProj),
	theSlavePointProj(a.theSlavePointProj),
	theHasRun(a.theHasRun),
	handlerM(a.handlerM),
	handlerS(a.handlerS),
	theMasterBandSelector(a.theMasterBandSelector),
	theSlaveBandSelector(a.theSlaveBandSelector),
	theTileSize(a.theTileSize),
	theSlaveId(a.theSlaveId),
	theMasterId(a.theMasterId),
	thePointType(a.thePointType),
	theThreadNum(a.theThreadNum),
	theSlaveProjection(a.theSlaveProjection),
	theTset(),
	theMinNCC(0.75),
	theUseGeographic(a.theUseGeographic),
	theDebug(a.theDebug),
	theStoreFlag(a.theStoreFlag)
{
}

radiImageRegistration::~radiImageRegistration()
{
   //TBC : delete handlers created by factory?
   
   if(caster.size())
   {
      caster[0]->disconnect();
      caster[1]->disconnect();
      caster.clear();
   }
   if(theMChain.valid())
   {
      theMChain->disconnect();
   }
   if(theSChain.valid())
   {
      theSChain->disconnect();
   }
   theMChain = 0;
   theSChain = 0;

}

ossimString radiImageRegistration::getRole() const
{
   ossimString role = "unknown";
   
   //use slave or master projection
   if (theProjectionType == "S")
   {
      role="slave";
   }
   else if (theProjectionType == "M")
   {
      role="master";
   }
   else
   {
      cerr<<"radiImageRegistration::getRole unknown output projection, need to supply it"<<endl;
   }

   return role;
}

ossimImageHandler*  radiImageRegistration::getProjectionHandler()
{
   //use slave or master projection
   ossimImageHandler* projHandler = 0;
   if (theProjectionType == "S")
   {
      projHandler = handlerS.get();
   }
   else if (theProjectionType == "M")
   {
	   if (!handlerM)
	   {
		   handlerM = ossimImageHandlerRegistry::instance()->open(theLastMaster);
		   if (!handlerM)
		   {
			   cerr<<"radiImageRegistration"<<"::execute can't create handler for slave image  "<< theSlave <<endl;
			   return false;
		   }
		   if(theMChain.valid())
		   {
			   theMChain->disconnect();
		   }
		   theMChain = new ossimImageChain;
		   theMChain->add(handlerM.get());
	   }
      projHandler = handlerM.get();
   }
   else
   {
      cerr<<"radiImageRegistration::getProjectionHandler cannot get handler for " << getRole() <<endl;
   }
   return projHandler;
}

ossimRefPtr<ossimImageGeometry> radiImageRegistration::getOutputImageGeometry()
{
   ossimRefPtr<ossimImageGeometry> geom = 0;
   ossimKeywordlist prjKwl;
   if (!theProjectionFile.empty())
   {
	   prjKwl.addFile(theProjectionFile);
	   geom->loadState(prjKwl);
   }
   else
   {
	   ossimImageHandler* projHandler = getProjectionHandler();
	   if(projHandler)
	   {
		   geom = projHandler->getImageGeometry();
	   }
   }
   return geom;
}

//getOutputProjection() - define output projection
// according to projType
ossimMapProjection* radiImageRegistration::getOutputProjection()
{
   ossimMapProjection* mop = 0;

   ossimRefPtr<ossimImageGeometry> geom = getOutputImageGeometry();
   if( geom.valid() )
   {
      if ( geom->getProjection() )
      {
         mop = PTR_CAST(ossimMapProjection, geom->getProjection());
         if( !mop )
         {
            ossimDpt mpp = geom->getMetersPerPixel();
            ossimProjection* outProjection =
               ossimProjectionFactoryRegistry::instance()->
               createProjection(ossimString("ossimEquDistCylProjection"));
            mop = PTR_CAST(ossimMapProjection, outProjection);
         }
         
         if(mop)
         {
            mop->setDatum(ossimDatumFactory::instance()->wgs84());
            mop->update();

            // apply user scale factor (resize)
            // then hopefully overviews can be used
            if ( (theScaleRatio != 1.0) && (theScaleRatio>0) )
            {
               cout << "applying scale ratio on " << getRole() <<
                  ": "<<theScaleRatio << endl; //TBR?

               mop->applyScale(ossimDpt(1.0/theScaleRatio,1.0/theScaleRatio),
                               false);
            }
         }
      }
      else
      {
         cerr << "radiImageRegistration::getOutputProjection cannot create projection from " << getRole() <<" geometry." <<endl;
      }
   }
   else
   {
      cerr << "radiImageRegistration::getOutputProjection cannot get "
           <<getRole() << " geometry." << endl;
   }

   return mop;
}

// buildRenerer() - builds renderer for an imageSource
// accounts for :
//  -scale factor
//  -required projection
bool radiImageRegistration::buildRenderer(
   ossimImageChain* chain,
   ossimMapProjection* outProjection, 
   ossimImageRenderer* renderer,
   const ossimFilterResampler::ossimFilterResamplerType& stype ) const
{
   if(chain)
   {
      chain->add(new ossimCacheTileSource);
      ossimRefPtr<ossimImageGeometry> geom = chain->getImageGeometry();
      if(geom.valid()&&geom->getProjection())
      {       
         ossimImageViewProjectionTransform* transform = new ossimImageViewProjectionTransform;
         transform->setImageGeometry(geom.get());
         transform->setViewGeometry(new ossimImageGeometry(0, outProjection));
         renderer->setImageViewTransform(transform);
         renderer->getResampler()->setFilterType(stype);
         chain->add(renderer);
         chain->add(new ossimCacheTileSource);
      }
      else
      {
         cerr<<"radiImageRegistration"<<"::buildRenderer cannot get projection from master/slave"<<endl;
         return false;
      }      
   }
   else
   {
      cerr<<"radiImageRegistration"<<"::buildRenderer NULL source"<<endl;
      return false;
   }
   return true;
}


bool kpt_compare(const cv::KeyPoint& d1, const cv::KeyPoint& d2)
{
	return d1.response > d2.response;
}

void normalization(const cv::Mat& inMat, cv::Mat& outMat)
{
	cv::Scalar mean_value;// = cv::mean(inMat);
	cv::Mat stdDevMat;
	cv::meanStdDev(inMat, mean_value, stdDevMat);
	cv::divide(inMat-mean_value, stdDevMat, outMat);
}

void findGoodMatches(vector< vector< DMatch >  > all_matches_2, vector< DMatch >& good_matches, float nndrRatio = 0.80f)
{
	good_matches.clear();
	//for (int i = 0; i < (int)matches.size(); i++)
	//{
	//	good_matches.push_back(matches[i][0]);
	//}
	good_matches.reserve(all_matches_2.size());

	for (size_t i = 0; i < all_matches_2.size(); ++i)
	{ 
		if (all_matches_2[i].size() < 2)
			continue;

		const DMatch &m1 = all_matches_2[i][0];
		const DMatch &m2 = all_matches_2[i][1];

		if(m1.distance <= nndrRatio * m2.distance)
			good_matches.push_back(m1);     
	}
}


// loong
ossimIpt radiImageRegistration::slave2master(ossimProjection* slaveProjection,
									ossimProjection* masterProjection,
									ossimIpt slaveDpt)
{
	ossimDpt masterDpt;
	ossimGpt gpt;
	slaveProjection->lineSampleToWorld(slaveDpt, gpt);
	masterProjection->worldToLineSample(gpt, masterDpt);
	return ossimIpt(masterDpt.x, masterDpt.y);
}

ossimDpt radiImageRegistration::slave2master(ossimProjection* slaveProjection,
									ossimProjection* masterProjection,
									ossimDpt slaveDpt)
{
	ossimDpt masterDpt;
	ossimGpt gpt;
	slaveProjection->lineSampleToWorld(slaveDpt, gpt);
	masterProjection->worldToLineSample(gpt, masterDpt);
	return masterDpt;
}

ossimDpt radiImageRegistration::slave2master(ossimDpt slaveDpt)
{
	ossimDpt masterDpt;
	ossimGpt gpt;
	theSlaveProjection->lineSampleToWorld(slaveDpt, gpt);
	theMasterProjection->worldToLineSample(gpt, masterDpt);
	return masterDpt;
}

ossimIpt radiImageRegistration::slave2master(ossimIpt slaveDpt)
{
	ossimDpt masterDpt;
	ossimGpt gpt;
	theSlaveProjection->lineSampleToWorld(slaveDpt, gpt);
	theMasterProjection->worldToLineSample(gpt, masterDpt);
	return ossimIpt(masterDpt.x, masterDpt.y);
}

ossimDrect radiImageRegistration::imageRect2World(const ossimIrect& imageRect, ossimProjection* projection)
{
	ossimGpt p[4];
	if (point_type::control == thePointType)
	{
		//p[0] = projection->inverse(imageRect.ul());
		//p[1] = projection->inverse(imageRect.ur());
		//p[2] = projection->inverse(imageRect.lr());
		//p[3] = projection->inverse(imageRect.ll());
		projection->lineSampleToWorld(imageRect.ul(), p[0]);
		projection->lineSampleToWorld(imageRect.ur(), p[1]);
		projection->lineSampleToWorld(imageRect.lr(), p[2]);
		projection->lineSampleToWorld(imageRect.ll(), p[3]);
	}
	else
	{
		projection->lineSampleHeightToWorld(imageRect.ul(), 0.0, p[0]);
		projection->lineSampleHeightToWorld(imageRect.ur(), 0.0, p[1]);
		projection->lineSampleHeightToWorld(imageRect.lr(), 0.0, p[2]);
		projection->lineSampleHeightToWorld(imageRect.ll(), 0.0, p[3]);
	}

	double xmin = p[0].lon;
	double xmax = p[0].lon;
	double ymin = p[0].lat;
	double ymax = p[0].lat;

	for (int i=1;i<4;++i)
	{
		if (xmin > p[i].lon)
		{
			xmin = p[i].lon;
		}
		if (xmax < p[i].lon)
		{
			xmax = p[i].lon;
		}
		if (ymin > p[i].lat)
		{
			ymin = p[i].lat;
		}
		if (ymax < p[i].lat)
		{
			ymax = p[i].lat;
		}
	}

	return ossimDrect(xmin, ymin, xmax, ymax);
}

ossimIrect radiImageRegistration::worldRect2Image(const ossimDrect& wroldRect, ossimProjection* projection)
{
	ossimDpt dpt;
	ossimDpt p[4];
	dpt = wroldRect.ul();
	projection->worldToLineSample(ossimGpt(dpt.y, dpt.x, 0.0), p[0]);
	dpt = wroldRect.ur();
	projection->worldToLineSample(ossimGpt(dpt.y, dpt.x, 0.0), p[1]);
	dpt = wroldRect.lr();
	projection->worldToLineSample(ossimGpt(dpt.y, dpt.x, 0.0), p[2]);
	dpt = wroldRect.ll();
	projection->worldToLineSample(ossimGpt(dpt.y, dpt.x, 0.0), p[3]);

	int xmin = p[0].x;
	int xmax = p[0].x;
	int ymin = p[0].y;
	int ymax = p[0].y;

	for (int i=1;i<4;++i)
	{
		if (xmin > p[i].x)
		{
			xmin = p[i].x;
		}
		if (xmax < p[i].x)
		{
			xmax = p[i].x;
		}
		if (ymin > p[i].y)
		{
			ymin = p[i].y;
		}
		if (ymax < p[i].y)
		{
			ymax = p[i].y;
		}
	}

	return ossimIrect(xmin, ymin, xmax, ymax);
}

ossimDrect radiImageRegistration::worldRectIntersection(const ossimDrect& r1, const ossimDrect& r2)
{
	ossimDrect result;
	result.makeNan();
	if(r1.hasNans() || r2.hasNans())
	{

		return result;
	}


	double x0 = max(r1.ul().x, r2.ul().x);
	double x1 = min(r1.lr().x, r2.lr().x);
	double y0, y1;

	y0 = max(r1.ul().y, r2.ul().y);
	y1 = min(r1.lr().y, r2.lr().y);

	if( (x1 < x0) || (y1 < y0) )
		return result;
	else
		result = ossimDrect(x0, y0, x1, y1);
	return result;
}

bool radiImageRegistration::createTileMat(const ossimRefPtr<ossimCastTileSourceFilter>& cast, const ossimIrect& rect, cv::Mat& outMat, ossim_uint32 resLevel)
{
	//get Inputs
	ossimImageSource* imageSource = PTR_CAST(ossimImageSource, cast.get());
	if (!imageSource)
	{
		return false;
	}
	ossimRefPtr<ossimImageData> imageData = NULL;
	imageData = imageSource->getTile(rect, resLevel); //same resLevel?? TBC
	if(!imageData.valid() || imageData->getDataObjectStatus() == OSSIM_EMPTY
		|| imageData->getDataObjectStatus() == OSSIM_PARTIAL)
	{
		imageSource = NULL;
		return false;
	}

	//outMat = cv::Mat(cv::Size(imageData->getWidth(), imageData->getHeight()), CV_64FC1);
	//outMat.data = static_cast<uchar*>(imageData->getBuf(0));
	//outMat.convertTo(outMat, CV_8UC1);
	outMat = cv::Mat(cv::Size(rect.width(), rect.height()), CV_8UC1);
	outMat.data = static_cast<uchar*>(imageData->getBuf(0));
	//memcpy(outMat.data, static_cast<uchar*>(imageData->getBuf(0)), imageData->getWidth()*imageData->getHeight());
	//imageData->unloadTile((void*)outMat.data, rect, ossimInterleaveType::OSSIM_INTERLEAVE_UNKNOWN);
	imageSource = NULL;
	return true;
}

bool radiImageRegistration::getMasterList(ossimFilename spatial_index_file, vector<ossimFilename>& masterList,
				   ossimGpt ul, ossimGpt lr)
{
	masterList.clear();
	GDALAllRegister();
	OGRRegisterAll();//注册所有的文件格式驱动

	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");//得到shp文件的处理器

	OGRDataSource* poDS = poDriver->Open(spatial_index_file.c_str(), NULL);//打开文件
	
	OGRLayer* poLayer = poDS->GetLayer(0);//获取shp图层

	poLayer->SetSpatialFilterRect(ul.lon, ul.lat, lr.lon, lr.lat);
	//poLayer->SetSpatialFilterRect(ul.lat, ul.lon, lr.lat, lr.lon);

	OGRFeature *feature;
	while (feature = poLayer->GetNextFeature())
	{
		masterList.push_back(feature->GetFieldAsString("path"));
	}
	return true;
}

bool radiImageRegistration::execute()
{
	bool result=true;
	// check point type
	//ossimRefPtr<ossimImageHandler> mHandler;
	//ossimRefPtr<ossimProjection> mProjection;
	if (theMaster.ext().upcase() != "SHP")
	{ 
		handlerM = ossimImageHandlerRegistry::instance()->open(theMaster);
		if (NULL == handlerM.get())
		{
			cerr<<"radiImageRegistration"<<"::execute can't open master image  "<< theMaster <<endl;
			return false;
		}
		if(thePointType == point_type::control &&
			(NULL == handlerM->getImageGeometry().get() || NULL == (theMasterProjection = handlerM->getImageGeometry()->getProjection())))
		{
			// 无投影
			cerr<<"radiImageRegistration::execute can't get control points as"
				"the master image does not contain projection information"<<endl;
			return false;
		}
	}

	std::vector< std::pair<int, ossimFilename> > imageList;
	imageList.push_back(std::pair<int, ossimFilename>(theSlaveId, theSlave));
	imageList.push_back(std::pair<int, ossimFilename>(theMasterId, theMaster));
	theTset.clearImageList();
	theTset.addImagesToList(imageList);

	if(theSChain.valid())
	{
		theSChain->disconnect();
	}
	theSChain = new ossimImageChain;

	handlerS = ossimImageHandlerRegistry::instance()->open(theSlave);
	if (!handlerS)
	{
		cerr<<"radiImageRegistration"<<"::execute can't create handler for slave image  "<< theSlave <<endl;
		return false;
	}
	theSChain->add(handlerS.get());

	ossim_uint32 sbc = handlerS->getNumberOfOutputBands();
	//add a band selector
	ossim_uint32 sb = getSlaveBand();
	if (sb>=sbc) 
	{
		cerr<<"radiImageRegistration"<<"::execute Warning not enough bands in slave, only "<< sbc <<endl;
		sb=0;
	}
#if OSSIM_HAS_MPI
	if (ossimMpi::instance()->getRank() == 0)
#endif
	{
		cout<<"Using band "<<sb+1<<" for slave"<<endl; //TBR
	}
	theSlaveBandSelector = new ossimBandSelector;
	theSlaveBandSelector->connectMyInputTo(0, handlerS.get());
	theSlaveBandSelector->setOutputBandList(vector<ossim_uint32>(1,sb));
	//      theSlaveSource = theSlaveBandSelector;
	theSChain->add(theSlaveBandSelector.get());

	//init casters
	//caster[0]->setOutputScalarType(OSSIM_UCHAR);
	//caster[1]->setOutputScalarType(OSSIM_UCHAR);
	caster[0]->setOutputScalarType(OSSIM_FLOAT64);
	caster[1]->setOutputScalarType(OSSIM_FLOAT64);

	//init gen
	setStoreFlag(true); //also keep all tie points in memory

	//TBD : set area of interest to buffer around slave?

	// -- 3 -- tie blocks, from sources to outputs

	caster[0]->connectMyInputTo(0, theSlaveBandSelector.get());


	cornerDetector->setK(0.05); //hardcoded
	cornerDetector->setGaussStd(7 / 2.0); //TBC : hardcoded ratio
	cornerDetector->setDensity(0.002);

	cornerDetector->connectMyInputTo(0, theSChain.get());
	//////////////////////////////////////////////////////////////////////////
	//ossimKeywordlist kwl;
	//handlerS->getImageGeometry()->getProjection()->saveState(kwl);
	//theSlaveProjection = ossimProjectionFactoryRegistry::instance()->createProjection(kwl);

	theSlaveProjection = handlerS->getImageGeometry()->getProjection();
	//////////////////////////////////////////////////////////////////////////
	if(theAreaOfInterest.hasNans() || theAreaOfInterest.width() <= 1 || theAreaOfInterest.height() <= 1)
	{
		theAreaOfInterest = handlerS->getBoundingRect(0);
	}
	else
	{
		theAreaOfInterest = theAreaOfInterest.clipToRect(handlerS->getBoundingRect(0));
	}
	if (theMaster.ext().upcase() != "SHP" &&
		!(NULL == handlerM->getImageGeometry().get() || NULL == (theMasterProjection = handlerM->getImageGeometry()->getProjection())))
	{
		ossimIrect mBoundary = handlerM->getBoundingRect(0);
		ossimDrect mWorldRect = imageRect2World(mBoundary, theMasterProjection);
		ossimDrect sWorldRect = imageRect2World(theAreaOfInterest, theSlaveProjection);
		//sWorldRect = sWorldRect.clipToRect(mWorldRect);
		sWorldRect = worldRectIntersection(sWorldRect, mWorldRect);
		//theAreaOfInterest = theAreaOfInterest.clipToRect(worldRect2Image(sWorldRect, theSlaveProjection));
		if (sWorldRect.hasNans())
		{
			theAreaOfInterest.makeNan();
		}
		else
		{
			ossimGpt gpt;
			ossimDpt mul, mlr;
			if (point_type::control == thePointType)
			{
				theMasterProjection->lineSampleToWorld(mBoundary.ul(), gpt);
				theSlaveProjection->worldToLineSample(gpt, mul);
				theMasterProjection->lineSampleToWorld(mBoundary.lr(), gpt);
				theSlaveProjection->worldToLineSample(gpt, mlr);
			}
			else
			{
				theMasterProjection->lineSampleHeightToWorld(mBoundary.ul(), 0.0, gpt);
				theSlaveProjection->worldToLineSample(gpt, mul);
				theMasterProjection->lineSampleHeightToWorld(mBoundary.lr(), 0.0, gpt);
				theSlaveProjection->worldToLineSample(gpt, mlr);
			}
			//mul = mul + ossimDpt(-theSlaveAccuracy*2, -theSlaveAccuracy*2);
			//mlr = mlr + ossimDpt(theSlaveAccuracy*2, theSlaveAccuracy*2);
			ossimIrect mBoundary2s(mul.x, mul.y, mlr.x, mlr.y);
			theAreaOfInterest = theAreaOfInterest.clipToRect(mBoundary2s);
		}
	}
	if(theAreaOfInterest.hasNans() || theAreaOfInterest.width() <= 1 || theAreaOfInterest.height() <= 1)
	{
		cerr<<"no overlap region is found."<<endl;
		return false;
	}

	//open();

	theTset.clearTiePoints();
	// -- 4 -- run
	result = getAllFeatures();

	cout<<theTset.getTiePoints().size()<<" tie points are found."<<endl;

//
//#if OSSIM_HAS_MPI
//	int myid = ossimMpi::instance()->getRank();
//	char buf[1024];
//	sprintf_s(buf, "%s_%d.%s",
//		theFilename.fileNoExtension(),
//		myid, theFilename.ext());
//	ossimFilename tempFilename =  ossimFilename(buf);
//	fstream ofs;
//	ofs.open(tempFilename.c_str(), ios_base::out);
//	ossimMapProjection* outProjection = radiImageRegistration::getOutputProjection();
//	vector< ossimRefPtr<ossimTieGpt> >::const_iterator it;
//	int icount = 0;
//	for (it = theTset.getTiePoints().begin();it!=theTset.getTiePoints().end();++it)
//	{
//		char buf[1024];
//		double hgt = ossimElevManager::instance()->getHeightAboveEllipsoid((*it)->getGroundPoint());
//		ossimDpt dpt = (*it)->getImagePoint();
//		ossimGpt gpt = (*it)->getGroundPoint();
//		if (!outProjection->isGeographic())
//		{
//			ossimDpt d1 = outProjection->forward(gpt);
//			gpt = ossimGpt(d1.x, d1.y, hgt);
//			sprintf_s(buf, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n\0", 
//				icount+=1,
//				dpt.x, 
//				dpt.y,
//				gpt.lat,
//				gpt.lon,
//				hgt);
//			int nLength = strlen(buf);
//			ofs<<buf;
//		}
//		else
//		{
//			sprintf_s(buf, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n\0", 
//				icount+=1,
//				dpt.x, 
//				dpt.y,
//				gpt.lon,
//				gpt.lat,
//				hgt);
//			int nLength = strlen(buf);
//			ofs<<buf;
//		}
//	}
//	ofs.close();
//	MPI_Barrier(MPI_COMM_WORLD);
//	if (ossimMpi::instance()->getRank() == 0)
//	{
//		fstream outFs;
//		outFs.open(theFilename.c_str(), ios_base::out);
//		for (int i = 0;i < ossimMpi::instance()->getNumberOfProcessors();++i)
//		{
//			char buf_[1024];
//			sprintf_s(buf_, "%s_%d.%s",
//			theFilename.fileNoExtension(),
//			i, theFilename.ext());
//			ossimFilename tmpfile =  ossimFilename(buf_);
//			fstream inFs;
//			inFs.open(tmpfile.c_str(), ios_base::in);
//			char charBuf[1024];
//			while (inFs.getline(charBuf, 1024))
//			{
//				outFs<<charBuf<<endl;
//			}
//			inFs.close();
//			DeleteFileA(tmpfile);
//		}
//		outFs.close();
//
//		fstream geomFs;
//		geomFs.open(theFilename.setExtension("geom").c_str(), ios_base::out);
//		ossimKeywordlist prjKwl;
//		getOutputProjection()->saveState(prjKwl);
//		geomFs<<prjKwl;
//		geomFs.close();
//	}
//#else
//	fstream ofs;
//	ofs.open(theFilename.c_str(), ios_base::out);
//	ossimMapProjection* outProjection = radiImageRegistration::getOutputProjection();
//	vector< ossimRefPtr<ossimTieGpt> >::const_iterator it;
//	int icount = 0;
//	for (it = theTset.getTiePoints().begin();it!=theTset.getTiePoints().end();++it)
//	{
//		char buf[1024];
//		//double hgt = ossimElevManager::instance()->getHeightAboveEllipsoid((*it)->getGroundPoint());
//		double hgt = (*it)->getGroundPoint().hgt;
//		ossimDpt dpt = (*it)->getImagePoint();
//		ossimGpt gpt = (*it)->getGroundPoint();
//		if(!outProjection || outProjection->isGeographic())
//		{
//			sprintf_s(buf, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n\0", 
//				icount+=1,
//				dpt.x, 
//				dpt.y,
//				gpt.lon,
//				gpt.lat,
//				hgt);
//			int nLength = strlen(buf);
//			ofs<<buf;
//		}
//		else
//		{
//			ossimDpt d1 = outProjection->forward(gpt);
//			gpt = ossimGpt(d1.x, d1.y, hgt);
//			sprintf_s(buf, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n\0", 
//				icount+=1,
//				dpt.x, 
//				dpt.y,
//				gpt.lat,
//				gpt.lon,
//				hgt);
//			int nLength = strlen(buf);
//			ofs<<buf;
//		}
//	}
//	ofs.close();
//
//	fstream geomFs;
//	geomFs.open(theFilename.setExtension("geom").c_str(), ios_base::out);
//	ossimKeywordlist prjKwl;
//	getOutputProjection()->saveState(prjKwl);
//	geomFs<<prjKwl;
//	geomFs.close();
//#endif

	//close();

	theHasRun = true;
	return true;
}

void radiImageRegistration::appendTiePoints(const ossimFilename& filename)
{
	int nPoints = int(theTset.getTiePoints().size());
	if (nPoints < 1)
		return;
	if(0 == strcmp(filename.ext().upcase().c_str(), "TXT"))
	{
		fstream ofs;
		ofs.open(theFilename.c_str(), ios_base::app);
		ofs<<"# point type: tie point"<<endl;
		ofs<<"# slave image id: "<<theSlaveId<<endl;
		ofs<<"# master image id: "<<theMasterId<<endl;
		ofs<<"# slave image path: "<<theSlave<<endl;
		ofs<<"# master image path: "<<theMaster<<endl;
		ofs<<"# point count: "<<theTset.getTiePoints().size()<<endl;
		vector< ossimRefPtr<radiBlockTieGpt> >::const_iterator it;
		int icount = 0;
		for (it = theTset.getTiePoints().begin();it!=theTset.getTiePoints().end();++it)
		{
			char buf[1024];
			//double hgt = ossimElevManager::instance()->getHeightAboveEllipsoid((*it)->getGroundPoint());
			double hgt = (*it)->getGroundPoint().hgt;
			ossimDpt dpt = (*it)->getImagePoint();
			ossimGpt gpt = (*it)->getGroundPoint();
			sprintf_s(buf, "%d\t%lf\t%lf\t%lf\t%lf\n\0", 
				icount+=1,
				dpt.x, 
				dpt.y,
				gpt.lon,
				gpt.lat);
			int nLength = strlen(buf);
			ofs<<buf;
		}
		ofs.close();
	}
	else if (0 == strcmp(filename.ext().upcase().c_str(), "XML"))
	{
		ossimXmlDocument oldGmlDoc(filename);
		ossimRefPtr<ossimXmlNode> root_node = oldGmlDoc.getRoot();
		theTset.appendAsGmlNode(root_node);
		ossimXmlDocument newGmlDoc;
		newGmlDoc.initRoot(root_node); //need namespaces etc...
		newGmlDoc.write(filename);
	}
}

void radiImageRegistration::appendControlPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
{
	int nPoints = int(theTset.getTiePoints().size());
	if (nPoints < 1)
		return;
	if(0 == strcmp(filename.ext().upcase().c_str(), "TXT"))
	{
		fstream ofs;
		ofs.open(theFilename.c_str(), ios_base::app);
		ofs<<"# point type: control point"<<endl;
		ofs<<"# slave image id: "<<theSlaveId<<endl;
		ofs<<"# master image id: "<<theMasterId<<endl;
		ofs<<"# slave image path: "<<theSlave<<endl;
		ofs<<"# master image path: "<<theMaster<<endl;
		ofs<<"# point count: "<<theTset.getTiePoints().size()<<endl;
		vector< ossimRefPtr<radiBlockTieGpt> >::const_iterator it;
		int icount = 0;
		if(NULL == pMapProjection)
		{
			if (theUseGeographic)
			{
				pMapProjection = PTR_CAST(ossimMapProjection, 
					ossimProjectionFactoryRegistry::instance()->
					createProjection(ossimString("ossimEquDistCylProjection")));
			}
			else
			{
				pMapProjection = getOutputProjection();
			}
		}
		for (it = theTset.getTiePoints().begin();it!=theTset.getTiePoints().end();++it)
		{
			char buf[1024];
			double hgt = (*it)->getGroundPoint().hgt;
			ossimDpt dpt = (*it)->getImagePoint();
			ossimGpt gpt = (*it)->getGroundPoint();
			if(!pMapProjection || pMapProjection->isGeographic())
			{
				sprintf_s(buf, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n\0", 
					icount+=1,
					dpt.x, 
					dpt.y,
					gpt.lon,
					gpt.lat,
					hgt);
				int nLength = strlen(buf);
				ofs<<buf;
			}
			else
			{
				ossimDpt d1 = pMapProjection->forward(gpt);
				gpt = ossimGpt(d1.x, d1.y, hgt);
				sprintf_s(buf, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n\0", 
					icount+=1,
					dpt.x, 
					dpt.y,
					gpt.lat,
					gpt.lon,
					hgt);
				int nLength = strlen(buf);
				ofs<<buf;
			}
		}
		ofs.close();

		if (pMapProjection)
		{
			fstream geomFs;
			geomFs.open(theFilename.setExtension("geom").c_str(), ios_base::out);
			ossimKeywordlist prjKwl;
			pMapProjection->saveState(prjKwl);
			geomFs<<prjKwl;
			geomFs.close();
		}
	}
	else if (0 == strcmp(filename.ext().upcase().c_str(), "XML"))
	{
		ossimXmlDocument oldGmlDoc(filename);
		ossimRefPtr<ossimXmlNode> root_node = oldGmlDoc.getRoot();
		theTset.appendAsGmlNode(root_node);
		ossimXmlDocument newGmlDoc;
		newGmlDoc.initRoot(root_node); //need namespaces etc...
		newGmlDoc.write(filename);
	}
}


void radiImageRegistration::writeTiePoints(const ossimFilename& filename)
{
	fstream ofs;
	ofs.open(theFilename.c_str(), ios_base::out);
	ofs.close();
	appendTiePoints(filename);
}
void radiImageRegistration::writeControlPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
{
	fstream ofs;
	ofs.open(theFilename.c_str(), ios_base::out);
	ofs.close();
	appendControlPoints(filename, pMapProjection);
}

void radiImageRegistration::appendPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
{
	if (point_type::tie == thePointType)
	{
		appendTiePoints(filename);
	}
	else
	{
		appendControlPoints(filename, pMapProjection);
	}
}

void radiImageRegistration::writePoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
{
	if (point_type::tie == thePointType)
	{
		writeTiePoints(filename);
	}
	else
	{
		writeControlPoints(filename, pMapProjection);
	}
}

//int radiImageRegistration::runMatchParallel(const cv::Mat& slaveMat, const cv::Mat& masterMat, ossimTDpt& tDpt, void* pData, bool bDebug)
//{
//	//// Declare Ipoints and other stuff
//	//IpPairVec matches;
//	//IpVec ipts_slave;
//	//IpVec ipts_master;
//
//	//// Detect and describe interest points in the image
//	//surfDetDes(&IplImage(slaveMat), ipts_slave, false, 5, 4, 2, 0.00014f); 
//	//surfDetDes(&IplImage(masterMat), ipts_master, false, 5, 4, 2, 0.00014f); 
//
//	//// Fill match vector
//	//getMatches(ipts_slave,ipts_master,matches);
//
//	//if (matches.size() < 4)
//	//{
//	//	//handlerM->close();
//	//	return match_state::match_failed;
//	//}
//
//	////-- Create input data
//	//Eigen::MatrixXd dataPoints((int)matches.size(), 4);
//	//for(unsigned int i = 0;i < matches.size();++i)
//	//{
//	//	dataPoints(i, 0) = matches[i].first.x;
//	//	dataPoints(i, 1) = matches[i].first.y;
//	//	dataPoints(i, 2) = matches[i].second.x;
//	//	dataPoints(i, 3) = matches[i].second.y;
//	//}
//
//	//// RANSAC detect outliers
//	//auto_ptr< estimators::Solver<Eigen::MatrixXd, Eigen::VectorXd> > ptrSolver(
//	//	new estimators::affineSolver<Eigen::MatrixXd,Eigen::VectorXd>);
//	//vector<int> inliers;
//	////for (int i = 0; i < (int)good_matches.size(); i++) inliers.push_back(i);
//	//vector<Eigen::VectorXd> models;
//
//	//ransac::Ransac_Handler ransac_fun_Handler;
//	//bool result = ransac::Ransac_RobustEstimator
//	//	(
//	//	dataPoints, // the input data
//	//	estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::extractor, // How select sampled point from indices
//	//	dataPoints.rows(),  // the number of putatives data
//	//	*(ptrSolver.get()),  // compute the underlying model given a sample set
//	//	estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::defaultEvaluator,  // the function to evaluate a given model
//	//	//Ransac Object that contain function:
//	//	// CandidatesSelector, Sampler and TerminationFunction
//	//	ransac_fun_Handler, // the basic ransac object
//	//	1000,  // the maximum rounds for RANSAC routine
//	//	inliers, // inliers to the final solution
//	//	models, // models array that fit input data
//	//	0.95 // the confidence want to achieve at the end
//	//	);
//	//if (inliers.size() < 6)
//	//{
//	//	//handlerM->close();
//	//	return match_state::match_failed;
//	//}
//
//	//if (fabs(1.0 - fabs(models[0][1] * models[0][5] - models[0][2] * models[0][4])) > 0.5)
//	//{
//	//	//handlerM->close();
//	//	return match_state::match_failed;
//	//}
//
//	//// maybe a lot of correspondences are found, but we need only one correspondence for a tile
//	//ossimIpt sc = ossimIpt(matches[inliers[0]].first.x, matches[inliers[0]].first.y);
//	//ossimIpt mc = ossimIpt(matches[inliers[0]].second.x, matches[inliers[0]].second.y);
//
//	//tDpt = ossimTDpt( mc, sc, 0.0 );
//
//	//return match_state::success;
////////////////////////////////////////////////////////////////////////////
//
//	radiImageRegistration* pThis = (radiImageRegistration*)pData;
//	// detect corners
//	cv::initModule_nonfree();
//	std::vector<KeyPoint> skeypoints, mkeypoints;
//	Mat sdescriptors, mdescriptors;
//	SIFT detector(pThis->theSiftNfeatures, pThis->theSiftNOctaveLayers, 
//		pThis->theSiftContrastThreshold, pThis->theSiftEdgeThreshold, pThis->theSiftSigma);
//
//	if (bDebug)
//	{
//		cv::imwrite("slave.png", slaveMat);
//		cv::imwrite("master.png", masterMat);
//	}
//
//	// detect
//	detector.detect( slaveMat, skeypoints );
//	if (skeypoints.size() < 10 )
//	{
//		return match_state::slave_faild;
//
//	}
//	// detect
//	detector.detect( masterMat, mkeypoints );
//	if(mkeypoints.size() < 10)
//	{
//		//handlerM->close();
//		return match_state::master_faild;
//	}
//
//	// extract
//	cv::SiftDescriptorExtractor extractor;
//	extractor.compute( slaveMat, skeypoints, sdescriptors );
//	extractor.compute( masterMat, mkeypoints, mdescriptors );
//
//
//	BFMatcher matcher(NORM_L1, false); 
//	vector< vector< DMatch >  > matches;
//	matcher.knnMatch( sdescriptors, mdescriptors, matches, 2 );
//	
//	vector< DMatch > good_matches;
//	findGoodMatches(matches, good_matches, 0.70f);
//
//	if (good_matches.size() < 4)
//	{
//		//handlerM->close();
//		return match_state::match_failed;
//	}
//
//	//-- Create input data
//	Eigen::MatrixXd dataPoints((int)good_matches.size(), 4);
//	for(unsigned int i = 0;i < good_matches.size();++i)
//	{
//		dataPoints(i, 0) = skeypoints[good_matches[i].queryIdx].pt.x;
//		dataPoints(i, 1) = skeypoints[good_matches[i].queryIdx].pt.y;
//		dataPoints(i, 2) = mkeypoints[good_matches[i].trainIdx].pt.x;
//		dataPoints(i, 3) = mkeypoints[good_matches[i].trainIdx].pt.y;
//	}
//
//	// RANSAC detect outliers
//	auto_ptr< estimators::Solver<Eigen::MatrixXd, Eigen::VectorXd> > ptrSolver(
//		new estimators::affineSolver<Eigen::MatrixXd,Eigen::VectorXd>);
//	vector<int> inliers;
//	//for (int i = 0; i < (int)good_matches.size(); i++) inliers.push_back(i);
//	vector<Eigen::VectorXd> models;
//
//	ransac::Ransac_Handler ransac_fun_Handler;
//	bool result = ransac::Ransac_RobustEstimator
//		(
//		dataPoints, // the input data
//		estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::extractor, // How select sampled point from indices
//		dataPoints.rows(),  // the number of putatives data
//		*(ptrSolver.get()),  // compute the underlying model given a sample set
//		estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::defaultEvaluator,  // the function to evaluate a given model
//		//Ransac Object that contain function:
//		// CandidatesSelector, Sampler and TerminationFunction
//		ransac_fun_Handler, // the basic ransac object
//		1000,  // the maximum rounds for RANSAC routine
//		inliers, // inliers to the final solution
//		models, // models array that fit input data
//		0.95 // the confidence want to achieve at the end
//		);
//	if (inliers.size() < 6)
//	{
//		//handlerM->close();
//		return match_state::match_failed;
//	}
//
//	//if (fabs(models[0][1] * models[0][5] - models[0][2] * models[0][4]) < 0.5)
//	//{
//	//	//handlerM->close();
//	//	return match_state::match_failed;
//	//}
//	if (fabs(1.0 - fabs(models[0][1] * models[0][5] - models[0][2] * models[0][4])) > 0.5)
//	{
//		//handlerM->close();
//		return match_state::match_failed;
//	}
//
//	// maybe a lot of correspondences are found, but we need only one correspondence for a tile
//	ossimIpt mc = ossimIpt(mkeypoints[good_matches[inliers[0]].trainIdx].pt.x, mkeypoints[good_matches[inliers[0]].trainIdx].pt.y);
//	ossimIpt sc = ossimIpt(skeypoints[good_matches[inliers[0]].queryIdx].pt.x, skeypoints[good_matches[inliers[0]].queryIdx].pt.y);
//
//	tDpt = ossimTDpt( mc, sc, good_matches[inliers[0]].distance );
//
//
//	if (bDebug)
//	{
//		//cout<<models[0]<<endl;
//		vector< DMatch > final_matches;
//		for (int i = 0;i < (int)inliers.size();++i)
//		{
//			final_matches.push_back(good_matches[inliers[i]]);
//		}
//		// Draw matches
//		cv::Mat imgMatch;
//		drawMatches(slaveMat, skeypoints, masterMat, mkeypoints, final_matches, imgMatch);
//		cv::imwrite("result.png", imgMatch);
//	}
//
//	//handlerM->close();
//	return match_state::success;
//}

void radiImageRegistration::VLFeatSift(const cv::Mat& inMat, vector<KeyPoint>& kpts, cv::Mat& descriptors)
{
	int noctaves=2,nlevels=4,o_min=0;
	// noctaves=(int)(log(min)/log(2));
	vl_sift_pix *ImageData=new vl_sift_pix[inMat.rows * inMat.cols];
	unsigned char *Pixel;
	for (int i=0;i<inMat.rows;i++)
	{
		for (int j=0;j<inMat.cols;j++)
		{
			Pixel=(unsigned char*)(inMat.data+i*inMat.cols+j);
			ImageData[i*inMat.cols+j]=*(Pixel);
		}
	}
	VlSiftFilt *SiftFilt=NULL;
	SiftFilt=vl_sift_new(inMat.cols, inMat.rows, noctaves, nlevels, o_min);
	double edge_thresh = 30 ;  //-1 will use the default (as in matlab)
	double peak_thresh = 0.001 ;
	double norm_thresh = -1 ;
	double magnif      = -1 ;
	double window_size = -1 ;
	if (peak_thresh >= 0) vl_sift_set_peak_thresh (SiftFilt, peak_thresh) ;
	if (edge_thresh >= 0) vl_sift_set_edge_thresh (SiftFilt, edge_thresh) ;
	if (norm_thresh >= 0) vl_sift_set_norm_thresh (SiftFilt, norm_thresh) ;
	if (magnif      >= 0) vl_sift_set_magnif      (SiftFilt, magnif) ;
	if (window_size >= 0) vl_sift_set_window_size (SiftFilt, window_size);
	int nKeyPoint=0;
	int idx=0;

	//descriptor = cv::Mat(cv::Size(128, 1), CV_32FC1);
	if (vl_sift_process_first_octave(SiftFilt, ImageData)!=VL_ERR_EOF)
	{
		while (TRUE)
		{
			//计算每组中的关键点
			vl_sift_detect(SiftFilt);
			//遍历并绘制每个点			
			nKeyPoint += SiftFilt->nkeys;
			VlSiftKeypoint *pKeyPoint=SiftFilt->keys;
			for (int i=0;i<SiftFilt->nkeys;i++)
			{
				VlSiftKeypoint TemptKeyPoint=*pKeyPoint;
				pKeyPoint++;
				//cv::KeyPoint kpt(float x, float y, float _size, float _angle=-1,
				//	float _response=0, int _octave=0, int _class_id=-1);
				cv::KeyPoint kpt(TemptKeyPoint.x, TemptKeyPoint.y, TemptKeyPoint.sigma/2);
				//kpts.push_back(kpt);

				//cvDrawCircle(Image, cvPoint(TemptKeyPoint.x,TemptKeyPoint.y),TemptKeyPoint.sigma/2,CV_RGB(255,0,0));
				idx++;
				//计算并遍历每个点的方向
				double angles[4];
				int angleCount=vl_sift_calc_keypoint_orientations(SiftFilt,angles,&TemptKeyPoint);
				for (int j=0;j<angleCount;j++)
				{
					double TemptAngle=angles[j];
					//printf("%d: %f\n",j,TemptAngle);
					//计算每个方向的描述
					cv::Mat aDescriptor = cv::Mat(cv::Size(128, 1), CV_32FC1);
					//float *Descriptors=new float[128];
					vl_sift_calc_keypoint_descriptor(SiftFilt, (float*)aDescriptor.data, &TemptKeyPoint, TemptAngle);
					descriptors.push_back(aDescriptor);
					kpts.push_back(kpt);
					//int k=0;
					//while (k<128)
					//{
					//	printf("%d: %f",k,Descriptors[k]);
					//	printf("; ");
					//	k++;
					//}

					//printf("\n");
					//delete []Descriptors;
					//Descriptors=NULL;
				}

				//cv::Mat aDescriptor = cv::Mat(cv::Size(128, 1), CV_32FC1);
				//vl_sift_calc_keypoint_descriptor(SiftFilt, (float*)aDescriptor.data, &TemptKeyPoint, angles[0]);
				//descriptors.push_back(aDescriptor);
			}
			//下一阶
			if (vl_sift_process_next_octave(SiftFilt)==VL_ERR_EOF)
			{
				break;
			}
			//free(pKeyPoint);
			nKeyPoint = NULL;
		}
	}
	vl_sift_delete(SiftFilt);
	delete []ImageData;
	ImageData=NULL;
}

static void _prepareImgAndDrawKeylines( const cv::Mat& img1, 
									   const cv::Mat& img2, 
									   ossimTDpt tpt,
									   cv::Mat& outImg,
									   const cv::Scalar& singlePointColor)
{
	Size size( img1.cols + img2.cols, MAX(img1.rows, img2.rows) );

	outImg.create( size, CV_MAKETYPE(img1.depth(), 3) );
	cv::Mat outImg1 = outImg( Rect(0, 0, img1.cols, img1.rows) );
	cv::Mat outImg2 = outImg( Rect(img1.cols, 0, img2.cols, img2.rows) );
	if( img1.type() == CV_8U )
		cvtColor( img1, outImg1, CV_GRAY2BGR );
	else
		img1.copyTo( outImg1 );

	if( img2.type() == CV_8U )
		cvtColor( img2, outImg2, CV_GRAY2BGR );
	else
		img2.copyTo( outImg2 );
	// draw keypoints

	ossimDpt slavePt = tpt.getSlavePoint();
	ossimDpt masterPt = tpt.getMasterPoint();
	int semiCrossWidth = 5;
	cv::line(outImg1, cv::Point(slavePt.x - semiCrossWidth, slavePt.y), cv::Point(slavePt.x + semiCrossWidth, slavePt.y),
		singlePointColor, 1);
	cv::line(outImg1, cv::Point(slavePt.x, slavePt.y - semiCrossWidth), cv::Point(slavePt.x, slavePt.y + semiCrossWidth),
		singlePointColor, 1);

	cv::line(outImg2, cv::Point(masterPt.x - semiCrossWidth, masterPt.y), cv::Point(masterPt.x + semiCrossWidth, masterPt.y),
		singlePointColor, 1);
	cv::line(outImg2, cv::Point(masterPt.x, masterPt.y - semiCrossWidth), cv::Point(masterPt.x, masterPt.y + semiCrossWidth),
		singlePointColor, 1);
}

static void drawTogether(const cv::Mat& img1,
	const cv::Mat& img2,
	cv::Mat& outImg)
{
	cv::Mat outImg1;
	cv::Mat outImg2;
	Size size(img1.cols + img2.cols, MAX(img1.rows, img2.rows));

	outImg.create(size, CV_MAKETYPE(img1.depth(), 3));
	outImg1 = outImg(Rect(0, 0, img1.cols, img1.rows));
	outImg2 = outImg(Rect(img1.cols, 0, img2.cols, img2.rows));

	if (img1.type() == CV_8U)
		cvtColor(img1, outImg1, CV_GRAY2BGR);
	else
		img1.copyTo(outImg1);

	if (img2.type() == CV_8U)
		cvtColor(img2, outImg2, CV_GRAY2BGR);
	else
		img2.copyTo(outImg2);
}

bool radiImageRegistration::runMatchNcc(ossimIrect srect, double gsd_scale, ossimTDpt& tDpt, void* pData, bool bDebug)
{
	radiImageRegistration* pThis = (radiImageRegistration*)pData;
	int resLevel = 0;
	theMinNCC = 0.8;
	int sRadius = 11;

	//if (!cornerDetector)
	//{
	//	return false;
	//}

	long w = srect.width();
	long h = srect.height();

	////get corner data tile (same size as inner tile)
	//_cornermutex.lock();
	//ossimRefPtr<ossimImageData> cornerData = cornerDetector->getTile(srect, resLevel);
	//_cornermutex.unlock();
	//if (!cornerData.valid() || !isSourceEnabled())
	//{
	//	return false;
	//}


	GdalRasterApp slaveApp;
	if (!slaveApp.open(getSlave().c_str()))
	{
		cerr << "radiImageRegistration" << "::execute can't open slave image  " << getSlave().c_str() << endl;
		return false;
	}
	// select only one band (if multiple)
	ossim_uint32 sbc = slaveApp.nBand();
	//add a band selector
	ossim_uint32 sb = theSlaveBand;
	if (sb >= sbc)
	{
		cerr << "radiImageRegistration" << "::execute Warning not enough bands in slave, only " << sbc << endl;
		sb = 0;
	}

	cv::Mat slaveMatAll;
	slaveApp.getRect2CvMatByte(srect, slaveMatAll, theSlaveBand, 1.0);
	if (countNonZero(slaveMatAll) < 1)
	{
		return false;
	}


	GdalRasterApp masterApp;
	if (!masterApp.open(getMaster().c_str()))
	{
		cerr << "radiImageRegistration" << "::execute can't open master image  " << getMaster().c_str() << endl;
		return false;
	}
	// select only one band (if multiple)
	ossim_uint32 mbc = masterApp.nBand();
	//add a band selector
	ossim_uint32 mb = theMasterBand;
	if (mb >= mbc)
	{
		cerr << "radiImageRegistration" << "::execute Warning not enough bands in master, only " << mbc << endl;
		mb = 0;
	}

	ossimIpt delta_sr(sRadius, sRadius);
	ossimIpt delta_mr(sRadius + (ossim_int32)(ceil(theSlaveAccuracy)),
		sRadius + (ossim_int32)(ceil(theSlaveAccuracy)));
	delta_mr = delta_mr * gsd_scale;

	ossimNCC_FFTW* theNCCengine = NULL;
	//TBD: use pixel size in meters to change delta_lr according to zoom

	vector<ossimTDpt> theTies;

	// 改进的harris角点检测方法
	std::vector<cv::Point> corners;
	cv::goodFeaturesToTrack(slaveMatAll, corners,
		200,
		//角点最大数目
		0.01,
		// 质量等级，这里是0.01*max（min（e1，e2）），e1，e2是harris矩阵的特征值
		10);

	for (size_t ip = 0; ip < corners.size(); ip++)
	{
		int i = corners[ip].x;
		int j = corners[ip].y;

		//circle(slaveMatAll, Point(i, j), 5, Scalar(0), 2, 8, 0);
		//get slave data for specified center + radius
		// radius doesn't change with resLevel
		ossimIpt delta_sc(i, j);
		ossimIrect srect_tile(srect.ul() + delta_sc - delta_sr, srect.ul() + delta_sc + delta_sr); //square, size 2*radius+1 pixels

		cv::Mat slaveMat;
		slaveApp.getRect2CvMatByte(srect_tile, slaveMat, sb, 1.0);


		ossimIpt sul = srect_tile.ul();
		ossimIpt slr = srect_tile.lr();
		//ossimIpt delta_lr((ossim_int32)(ceil(theSlaveAccuracy)), (ossim_int32)(ceil(theSlaveAccuracy)));

		ossimDpt sc_tile;
		ossimDpt delta_mc;
		srect_tile.getCenter(sc_tile);
		ossimGpt gpt;
		theSlaveProjection->lineSampleToWorld(sc_tile, gpt);
		theMasterProjection->worldToLineSample(gpt, delta_mc);
		ossimIrect mrect_tile(delta_mc - delta_mr, delta_mc + delta_mr); //square, size 2*radius+1 pixels

		cv::Mat masterMat;
		masterApp.getRect2CvMatByte(mrect_tile, masterMat, mb, gsd_scale);
		if (countNonZero(masterMat) < 1)
		{
			return false;
		}


		//find normalized cross-correlation maximum
		//TBD: assuming floating point input
		double dx = 0.0;
		double dy = 0.0;
		double ncor = 0.0;

		getMaxCorrelation(theNCCengine, slaveMat, masterMat, &dx, &dy, &ncor);

		//filter on NCC value
		if (ncor >= theMinNCC)
		{
			//create tie point & store

			theTies.push_back(ossimTDpt(ossimDpt(dx, dy) / gsd_scale + delta_mc, srect.ul() + delta_sc, ncor));
			tDpt = ossimTDpt(ossimDpt(dx, dy) / gsd_scale + delta_mr, delta_sc, ncor);

			if (bDebug)
			{
				cv::Mat outMat;
				_prepareImgAndDrawKeylines(slaveMat,
					masterMat,
					ossimTDpt(ossimDpt(dx, dy) / gsd_scale + delta_mr, delta_sr, ncor),
					outMat,
					cv::Scalar(0, 0, 255));
				cv::imwrite("matched_.png", outMat);
				int mmm = 0;

				_prepareImgAndDrawKeylines(slaveMatAll,
					masterMat,
					tDpt,
					outMat,
					cv::Scalar(0, 0, 255));
				cv::imwrite("matched.png", outMat);
			}
			return true;
		}
	}
	return false;

	///// Detector parameters
	//int blockSize = 2;
	//int apertureSize = 3;
	//double k = 0.04;
	//int thresh = 200;
	///// Detecting corners
	//cv::Mat cor;
	//cornerHarris(slaveMatAll, cor, blockSize, apertureSize, k, BORDER_DEFAULT);
	///// Normalizing
	//normalize(cor, cor, 0, 255, NORM_MINMAX, CV_32FC1, Mat());
	//convertScaleAbs(cor, cor);

	////loop on corners (<>NULL & >=2 TBC)
	////ossim_uint32 coff = 0; //offset (speedup)
	////ossim_int32 ci = 0;
	////ossim_int32 cj = 0;
	////chip image radii
	//for (int j = 0; j < cor.rows; j++)
	//{
	//	for (int i = 0; i < cor.cols; i++)
	//	{
	//		if ((int)cor.at<float>(j, i) > thresh)
	//		{
	//			circle(slaveMatAll, Point(i, j), 5, Scalar(0), 2, 8, 0);
	//			//get slave data for specified center + radius
	//			// radius doesn't change with resLevel
	//			ossimIpt delta_sc(i, j);
	//			ossimIrect srect_tile(srect.ul() + delta_sc - delta_sr, srect.ul() + delta_sc + delta_sr); //square, size 2*radius+1 pixels

	//			cv::Mat slaveMat;
	//			slaveApp.getRect2CvMatByte(srect_tile, slaveMat, sb);


	//			ossimIpt sul = srect_tile.ul();
	//			ossimIpt slr = srect_tile.lr();
	//			//ossimIpt delta_lr((ossim_int32)(ceil(theSlaveAccuracy)), (ossim_int32)(ceil(theSlaveAccuracy)));

	//			ossimDpt sc_tile;
	//			ossimDpt delta_mc;
	//			srect_tile.getCenter(sc_tile);
	//			ossimGpt gpt;
	//			theSlaveProjection->lineSampleToWorld(sc_tile, gpt);
	//			theMasterProjection->worldToLineSample(gpt, delta_mc);
	//			ossimIrect mrect_tile(delta_mc - delta_mr, delta_mc + delta_mr); //square, size 2*radius+1 pixels

	//			cv::Mat masterMat;
	//			masterApp.getRect2CvMatByte(mrect_tile, masterMat, mb, gsd_scale);
	//			if (countNonZero(masterMat) < 1)
	//			{
	//				return false;
	//			}


	//			//find normalized cross-correlation maximum
	//			//TBD: assuming floating point input
	//			double dx = 0.0;
	//			double dy = 0.0;
	//			double ncor = 0.0;

	//			getMaxCorrelation(theNCCengine, slaveMat, masterMat, &dx, &dy, &ncor);

	//			//filter on NCC value
	//			if (ncor >= theMinNCC)
	//			{
	//				//create tie point & store

	//				theTies.push_back(ossimTDpt(ossimDpt(dx, dy) / gsd_scale + mrect_tile.ul(), srect.ul() + delta_mc, ncor));
	//			}

	//		}
	//	}
	//}



	if (bDebug)
	{
		cv::imwrite("slave_corner.png", slaveMatAll);
	}


	//if ((cornerData->getDataObjectStatus() != OSSIM_NULL) && (cornerData->getDataObjectStatus() != OSSIM_EMPTY))
	//{
	//	//loop on corners (<>NULL & >=2 TBC)
	//	ossim_uint32 coff = 0; //offset (speedup)
	//	ossim_int32 ci = 0;
	//	ossim_int32 cj = 0;
	//	//chip image radii
	//	ossimIpt delta_sr(sRadius, sRadius);
	//	ossimIpt delta_mr(sRadius + (ossim_int32)(ceil(theSlaveAccuracy)),
	//		sRadius + (ossim_int32)(ceil(theSlaveAccuracy)));
	//	delta_mr = delta_mr * gsd_scale;

	//	for (cj = 0; cj<h; ++cj) //rows
	//	{
	//		for (ci = 0; ci<w; ++ci) //cols
	//		{
	//			if (!cornerData->isNull(coff, 0))
	//			{
	//				//get slave data for specified center + radius
	//				// radius doesn't change with resLevel
	//				ossimIpt delta_sc(ci, cj);
	//				ossimIrect srect_tile(srect.ul() + delta_sc - delta_sr, srect.ul() + delta_sc + delta_sr); //square, size 2*radius+1 pixels

	//				cv::Mat slaveMat;
	//				slaveApp.getRect2CvMatByte(srect_tile, slaveMat, sb);


	//				ossimIpt sul = srect_tile.ul();
	//				ossimIpt slr = srect_tile.lr();
	//				//ossimIpt delta_lr((ossim_int32)(ceil(theSlaveAccuracy)), (ossim_int32)(ceil(theSlaveAccuracy)));

	//				ossimDpt sc_tile;
	//				ossimDpt delta_mc;
	//				srect_tile.getCenter(sc_tile);
	//				ossimGpt gpt;
	//				theSlaveProjection->lineSampleToWorld(sc_tile, gpt);
	//				theMasterProjection->worldToLineSample(gpt, delta_mc);
	//				ossimIrect mrect_tile(delta_mc - delta_mr, delta_mc + delta_mr); //square, size 2*radius+1 pixels

	//				cv::Mat masterMat;
	//				masterApp.getRect2CvMatByte(mrect_tile, masterMat, mb, gsd_scale);


	//				//find normalized cross-correlation maximum
	//				//TBD: assuming floating point input
	//				double dx = 0.0;
	//				double dy = 0.0;
	//				double ncor = 0.0;

	//				getMaxCorrelation(slaveMat, masterMat, &dx, &dy, &ncor);

	//				//filter on NCC value
	//				if (ncor >= theMinNCC)
	//				{
	//					//create tie point & store
	//					
	//					theTies.push_back(ossimTDpt(ossimDpt(dx, dy) / gsd_scale + mrect_tile.ul(), srect.ul() + delta_mc, ncor));
	//				}
	//			}
	//			++coff;
	//		}
	//	}
	//}

	//no need to erase theTile (automatic)
	if (theNCCengine != NULL)
	{
		delete theNCCengine;
		theNCCengine = NULL;
	}

	if (theTies.size() < 4)
	{
		//handlerM->close();
		return false;
	}

	////-- Create input data
	//Eigen::MatrixXd dataPoints((int)theTies.size(), 4);
	//for (unsigned int i = 0; i < theTies.size(); ++i)
	//{
	//	dataPoints(i, 0) = theTies[i].getSlavePoint().x;
	//	dataPoints(i, 1) = theTies[i].getSlavePoint().y;
	//	dataPoints(i, 2) = theTies[i].getMasterPoint().x;
	//	dataPoints(i, 3) = theTies[i].getMasterPoint().y;
	//}
	//// RANSAC detect outliers
	//auto_ptr< estimators::Solver<Eigen::MatrixXd, Eigen::VectorXd> > ptrSolver(
	//	new estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>);
	//vector<int> inliers;
	////for (int i = 0; i < (int)good_matches.size(); i++) inliers.push_back(i);
	//vector<Eigen::VectorXd> models;

	//ransac::Ransac_Handler ransac_fun_Handler;
	//bool result = ransac::Ransac_RobustEstimator
	//	(
	//	dataPoints, // the input data
	//	estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::extractor, // How select sampled point from indices
	//	dataPoints.rows(),  // the number of putatives data
	//	*(ptrSolver.get()),  // compute the underlying model given a sample set
	//	estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::defaultEvaluator,  // the function to evaluate a given model
	//	//Ransac Object that contain function:
	//	// CandidatesSelector, Sampler and TerminationFunction
	//	ransac_fun_Handler, // the basic ransac object
	//	1000,  // the maximum rounds for RANSAC routine
	//	inliers, // inliers to the final solution
	//	models, // models array that fit input data
	//	0.95 // the confidence want to achieve at the end
	//	);
	//if (inliers.size() < 4)
	//{
	//	//handlerM->close();
	//	return false;
	//}
	//double delta_energy = 1.0 - fabs(models[0][1] * models[0][5] - models[0][2] * models[0][4]);
	//if (fabs(delta_energy) > 0.5)
	//{
	//	//handlerM->close();
	//	return false;
	//}

	//double max_score = -99999.0;
	//for (int i = 1; i < (int)inliers.size(); ++i)
	//{
	//	if (theTies[inliers[i]].score > max_score)
	//	{
	//		tDpt = theTies[inliers[i]];
	//		max_score = tDpt.score;
	//	}
	//}

	//if (bDebug)
	//{
	//	cv::Mat outMat;
	//	_prepareImgAndDrawKeylines(slaveMatAll,
	//		masterMat,
	//		tDpt,
	//		outMat,
	//		cv::Scalar(0, 0, 255));
	//	cv::imwrite("matched.png", outMat);
	//}

	return true;

	//if (theTies.size() > 0)
	//{
	//	for (size_t i = 0; i < theTies.size(); i++)
	//	{
	//		if (theTies[i].score > max_score)
	//		{
	//			tDpt = theTies[i];
	//			max_score = tDpt.score;
	//		}
	//	}
	//	return true;
	//}
	//return false;
}


//void
//radiImageRegistration::getMaxCorrelation(ossimRefPtr<ossimImageData> Mchip, ossimRefPtr<ossimImageData> Schip,
//double* pdispx, double* pdispy, double* pcor)
//{
//	//use FFTW 3.0.1
//	//assume displacement between center of master to center of slave buffer
//	// Mchip must smaller than Schip (Schip incorporates error buffer)
//
//	ossim_uint32 sx = Schip->getWidth();
//	ossim_uint32 sy = Schip->getHeight();
//	ossim_uint32 mx = Mchip->getWidth();
//	ossim_uint32 my = Mchip->getHeight();
//	//cout<<"mx="<<mx<<" my="<<my<<" sx="<<sx<<" sy="<<sy<<endl; //TBR
//
//	int cx = sx + mx - 1;
//	int cy = sy + my - 1;
//
//	if (theNCCengine != NULL)
//	{
//		//check correlation size
//		if (!theNCCengine->sameDims(cy, cx))
//		{
//			//re build NCC engine //TBD : use wisdom
//			delete theNCCengine;
//			theNCCengine = NULL;
//		}
//	}
//	if (theNCCengine == NULL)
//	{
//		//build a new NCC engine //TBD : use wisdom
//		theNCCengine = new ossimNCC_FFTW(cy, cx);
//	}
//
//	theNCCengine->ingestMaster(my, mx, Mchip->getDoubleBuf());
//	theNCCengine->ingestSlave(sy, sx, Schip->getDoubleBuf());
//
//	if (!theNCCengine->calculateNCC())
//	{
//		// TBD err mngt
//		if (pcor) *pcor = 0.0;
//		if (pdispx) *pdispx = 0.0;
//		if (pdispy) *pdispy = 0.0;
//		cout << "Error in NCC calculation" << endl;
//		return;
//	}
//	int mj = theNCCengine->getMaxCorrX();
//	int mi = theNCCengine->getMaxCorrY();
//	double bestcorr = theNCCengine->getMaxCorr();
//	int oj = (cx - 1) / 2;//we know that cx and cy are odd!!
//	int oi = (cy - 1) / 2;
//	int deltaj = (sx - mx) / 2; //we know that sx-mx is even
//	int deltai = (sy - my) / 2;
//
//	//original best shift (integer shift for for max value)
//	double dmcx = mj - oj;
//	double dmcy = mi - oi;
//
//	//find maximum, sub-pixel precision
//	//use least-square fit on 2nd order polynomial
//	if ((mj > oj - deltaj) && (mj < oj + deltaj) && (mi > oi - deltai) && (mi < oi + deltai))
//	{
//		//then there's a 3x3 neighborhood we can use to get better precision
//		vector<double> p2c(6); //2nd order x y polynomial coefficents (see theLMS comments)
//		vector<double>::iterator it = p2c.begin();
//		double* pm = theLMS;
//		const ossimNCC_FFTW::cMatrix& corrmat = theNCCengine->getNcc();
//		//matrix product with values of 3x3 neighborhood
//		for (int k = 0; k<6; ++k)
//		{
//			*it = 0.0;
//			for (int i = -1; i <= 1; ++i)
//			{
//				for (int j = -1; j <= 1; ++j)
//				{
//					*it += *(pm++) * corrmat(mi + i, mj + j);
//				}
//			}
//			++it;
//		}
//		//check convexity (det>0) + downwards orientation (trace<0)
//		double trace = p2c[4] + p2c[5];
//		if (trace<-1e-13) //TBC : -epsilon
//		{
//			double det = p2c[4] * p2c[5] - 0.25*p2c[3] * p2c[3];
//			if (det>1e-13) //TBC : epsilon
//			{
//				//ok : convex + downwards
//				//find maximum position
//				double optx = (p2c[3] * p2c[2] - 2.0 * p2c[5] * p2c[1]) / det * 0.25;
//				double opty = (p2c[3] * p2c[1] - 2.0 * p2c[4] * p2c[2]) / det * 0.25;
//				//limit new position to center pixel square
//				//TBD : need to find better model for NCC subpixel
//				if ((fabs(optx) <= 0.501) && (fabs(opty) <= 0.501))
//				{
//					dmcx += optx;
//					dmcy += opty;
//					//change corelation max value (dangerous) : TBD ? TBC
//				}
//			}
//		}
//	}
//
//	//give results back  
//	if (pcor)   *pcor = bestcorr;
//	if (pdispx) *pdispx = dmcx;
//	if (pdispy) *pdispy = dmcy;
//}

void toDouble(const cv::Mat& mat, double** data)
{
	ossim_uint32 sx = mat.cols;
	ossim_uint32 sy = mat.rows;
	*data = new double[sx*sy];
	for (size_t i = 0; i < sx*sy; i++)
	{
		(*data)[i] = (double)mat.data[i];
	}	
}

void
radiImageRegistration::getMaxCorrelation(
ossimNCC_FFTW* &theNCCengine, cv::Mat masterMat, cv::Mat slaveMat,
double* pdispx, double* pdispy, double* pcor)
{
	//use FFTW 3.0.1
	//assume displacement between center of master to center of slave buffer
	// Mchip must smaller than Schip (Schip incorporates error buffer)

	ossim_uint32 sx = slaveMat.cols;
	ossim_uint32 sy = slaveMat.rows;
	ossim_uint32 mx = masterMat.cols;
	ossim_uint32 my = masterMat.rows;
	//cout<<"mx="<<mx<<" my="<<my<<" sx="<<sx<<" sy="<<sy<<endl; //TBR

	int cx = sx + mx - 1;
	int cy = sy + my - 1;

	if (theNCCengine != NULL)
	{
		//check correlation size
		if (!theNCCengine->sameDims(cy, cx))
		{
			//re build NCC engine //TBD : use wisdom
			delete theNCCengine;
			theNCCengine = NULL;
		}
	}
	if (theNCCengine == NULL)
	{
		//build a new NCC engine //TBD : use wisdom
		theNCCengine = new ossimNCC_FFTW(cy, cx);
	}

	double* mData;
	double* sData;
	toDouble(masterMat, &mData);
	toDouble(slaveMat, &sData);
	theNCCengine->ingestMaster(my, mx, mData);
	theNCCengine->ingestSlave(sy, sx, sData);

	if (!theNCCengine->calculateNCC())
	{
		// TBD err mngt
		if (pcor) *pcor = 0.0;
		if (pdispx) *pdispx = 0.0;
		if (pdispy) *pdispy = 0.0;
		cout << "Error in NCC calculation" << endl;
		delete[] mData;
		mData = NULL;
		delete[] sData;
		sData = NULL;
		return;
	}
	int mj = theNCCengine->getMaxCorrX();
	int mi = theNCCengine->getMaxCorrY();
	double bestcorr = theNCCengine->getMaxCorr();
	int oj = (cx - 1) / 2;//we know that cx and cy are odd!!
	int oi = (cy - 1) / 2;
	int deltaj = (sx - mx) / 2; //we know that sx-mx is even
	int deltai = (sy - my) / 2;

	//original best shift (integer shift for for max value)
	double dmcx = mj - oj;
	double dmcy = mi - oi;

	//find maximum, sub-pixel precision
	//use least-square fit on 2nd order polynomial
	if ((mj > oj - deltaj) && (mj < oj + deltaj) && (mi > oi - deltai) && (mi < oi + deltai))
	{
		//then there's a 3x3 neighborhood we can use to get better precision
		vector<double> p2c(6); //2nd order x y polynomial coefficents (see theLMS comments)
		vector<double>::iterator it = p2c.begin();
		double* pm = theLMS;
		const ossimNCC_FFTW::cMatrix& corrmat = theNCCengine->getNcc();
		//matrix product with values of 3x3 neighborhood
		for (int k = 0; k<6; ++k)
		{
			*it = 0.0;
			for (int i = -1; i <= 1; ++i)
			{
				for (int j = -1; j <= 1; ++j)
				{
					*it += *(pm++) * corrmat(mi + i, mj + j);
				}
			}
			++it;
		}
		//check convexity (det>0) + downwards orientation (trace<0)
		double trace = p2c[4] + p2c[5];
		if (trace<-1e-13) //TBC : -epsilon
		{
			double det = p2c[4] * p2c[5] - 0.25*p2c[3] * p2c[3];
			if (det>1e-13) //TBC : epsilon
			{
				//ok : convex + downwards
				//find maximum position
				double optx = (p2c[3] * p2c[2] - 2.0 * p2c[5] * p2c[1]) / det * 0.25;
				double opty = (p2c[3] * p2c[1] - 2.0 * p2c[4] * p2c[2]) / det * 0.25;
				//limit new position to center pixel square
				//TBD : need to find better model for NCC subpixel
				if ((fabs(optx) <= 0.501) && (fabs(opty) <= 0.501))
				{
					dmcx += optx;
					dmcy += opty;
					//change corelation max value (dangerous) : TBD ? TBC
				}
			}
		}
	}

	//give results back  
	if (pcor)   *pcor = bestcorr;
	if (pdispx) *pdispx = dmcx;
	if (pdispy) *pdispy = dmcy;

	delete[] mData;
	mData = NULL;
	delete[] sData;
	sData = NULL;
}

ossimIrect radiImageRegistration::getMasterRect(ossimProjection* sProjection, const GdalRasterApp& mGdalApp,
												ossimIrect sRect, double slaveAccuracy)
{
	//ossimIpt delta((ossim_int32)(ceil(slaveAccuracy)), (ossim_int32)(ceil(slaveAccuracy)) );
	ossimIpt mul, mlr;
	ossimGpt gpt;
	ossimIpt p[4];
	sProjection->lineSampleToWorld(sRect.ul() + ossimIpt(-slaveAccuracy,-slaveAccuracy), gpt);
	mGdalApp.lonlat2linesample(ossimDpt(gpt.lon, gpt.lat), p[0]);

	sProjection->lineSampleToWorld(sRect.ur() + ossimIpt(slaveAccuracy,-slaveAccuracy), gpt);
	mGdalApp.lonlat2linesample(ossimDpt(gpt.lon, gpt.lat), p[1]);

	sProjection->lineSampleToWorld(sRect.lr() + ossimIpt(slaveAccuracy,slaveAccuracy), gpt);
	mGdalApp.lonlat2linesample(ossimDpt(gpt.lon, gpt.lat), p[2]);

	sProjection->lineSampleToWorld(sRect.ll() + ossimIpt(-slaveAccuracy, slaveAccuracy), gpt);
	mGdalApp.lonlat2linesample(ossimDpt(gpt.lon, gpt.lat), p[3]);

	int xmin = p[0].x;
	int xmax = p[0].x;
	int ymin = p[0].y;
	int ymax = p[0].y;

	for (int i=1;i<4;++i)
	{
		if (xmin > p[i].x)
		{
			xmin = p[i].x;
		}
		if (xmax < p[i].x)
		{
			xmax = p[i].x;
		}
		if (ymin > p[i].y)
		{
			ymin = p[i].y;
		}
		if (ymax < p[i].y)
		{
			ymax = p[i].y;
		}
	}
	return ossimIrect(xmin, ymin, xmax, ymax);
}

ossimIrect radiImageRegistration::getMasterRect(ossimProjection* sProjection, ossimProjection* mProjection,
												ossimIrect sRect, double slaveAccuracy)
{
	//ossimIpt sul = sRect.ul() + ossimIpt(-slaveAccuracy, -slaveAccuracy);
	//ossimIpt sll = sRect.ll() + ossimIpt(-slaveAccuracy, slaveAccuracy);

	ossimIpt sul = sRect.ul() + ossimIpt(-slaveAccuracy, -slaveAccuracy);
	ossimIpt sur = sRect.ur() + ossimIpt(slaveAccuracy, -slaveAccuracy);
	ossimIpt slr = sRect.lr() + ossimIpt(slaveAccuracy, slaveAccuracy);
	ossimIpt sll = sRect.ll() + ossimIpt(-slaveAccuracy, slaveAccuracy);
	ossimDpt mul, mlr;
	ossimGpt gpt;
	if (point_type::control == thePointType)
	{
		//sProjection->lineSampleToWorld(sul, gpt);
		//mProjection->worldToLineSample(gpt, mul);
		//sProjection->lineSampleToWorld(slr, gpt);
		//mProjection->worldToLineSample(gpt, mlr);
		//ossimDpt p1, p2, p3, p4;
		ossimDpt p[4];
		sProjection->lineSampleToWorld(sul, gpt);
		mProjection->worldToLineSample(gpt, p[0]);
		sProjection->lineSampleToWorld(sur, gpt);
		mProjection->worldToLineSample(gpt, p[1]);
		sProjection->lineSampleToWorld(slr, gpt);
		mProjection->worldToLineSample(gpt, p[2]);
		sProjection->lineSampleToWorld(sll, gpt);
		mProjection->worldToLineSample(gpt, p[3]);

		int xmin = p[0].x;
		int xmax = p[0].x;
		int ymin = p[0].y;
		int ymax = p[0].y;

		for (int i = 1; i<4; ++i)
		{
			if (xmin > p[i].x)
			{
				xmin = p[i].x;
			}
			if (xmax < p[i].x)
			{
				xmax = p[i].x;
			}
			if (ymin > p[i].y)
			{
				ymin = p[i].y;
			}
			if (ymax < p[i].y)
			{
				ymax = p[i].y;
			}
		}
		return ossimDrect(xmin, ymin, xmax, ymax);
	}
	else
	{
		sProjection->lineSampleHeightToWorld(sul, 0.0, gpt);
		mProjection->worldToLineSample(gpt, mul);
		sProjection->lineSampleHeightToWorld(slr, 0.0, gpt);
		mProjection->worldToLineSample(gpt, mlr);
	}
	return ossimIrect(mul.x, mul.y, mlr.x, mlr.y);
}

bool radiImageRegistration::getGridFeatures(const ossimIrect& rect, radiBlockTieGptSet& tSet)
{

	//cout<<rect.ul()<<"\t"<<rect.lr()<<endl;
	if (!theSlaveBandSelector)
	{
		ossimNotify(ossimNotifyLevel_WARN)
			<< "WARN ossimTieGenerator::scanForEdges():"
			<< "\nInput source is not a ossimImageChip.  Returning..." << std::endl;
		return false;
	}

	const ossim_int32 TILE_HEIGHT	= min(theTileSize, rect.height());
	const ossim_int32 TILE_WIDTH	= min(theTileSize, rect.width());
	const ossim_int32 START_LINE = rect.ul().y;
	const ossim_int32 STOP_LINE  = rect.lr().y;
	const ossim_int32 START_SAMP = rect.ul().x;
	const ossim_int32 STOP_SAMP  = rect.lr().x;
	int nWidth = STOP_SAMP-START_SAMP;
	int nHeight = STOP_LINE-START_LINE;

	//// For percent complete status.
	//ossim_int32 tilerows = ossim_int32((STOP_LINE-START_LINE+TILE_HEIGHT) / TILE_HEIGHT + 0.5); //ceil : (stop-start+1+size-1)/size
	//ossim_int32 tilecols = ossim_int32((STOP_SAMP-START_SAMP+TILE_WIDTH) / TILE_WIDTH + 0.5);
	ossim_int32 tilerows = ossim_int32((STOP_LINE-START_LINE+1) / (double)TILE_HEIGHT + 0.5); //ceil : (stop-start+1+size-1)/size
	ossim_int32 tilecols = ossim_int32((STOP_SAMP-START_SAMP+1) / (double)TILE_WIDTH + 0.5);
	//double total_tiles = ((double)tilerows)*tilecols;
	//double tiles_processed = 0.0;

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
			//row_col_List.push_back(row_col((center_row - i - 1), (center_col - j - 1)));
			row_col_List.push_back(row_col(i, j));
		}
	}

	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);

	bool found = false;
	int search_size = (int)row_col_List.size();

	for (int i = 0; i < search_size; ++i)
	{
		int icol = int(row_col_List[i].col_idx);
		int irow = int(row_col_List[i].row_idx);
		ossim_int32 samp = START_SAMP + icol*TILE_WIDTH;
		ossim_int32 line = START_LINE + irow*TILE_HEIGHT;

		ossim_int32 BufHeight = min(TILE_HEIGHT, nHeight);
		ossim_int32 BufWidth = min(TILE_WIDTH, nWidth);
		//行末尾小块处理
		if (irow == tilerows - 1)
		{
			BufHeight = nHeight - (tilerows - 1) * TILE_HEIGHT;
			BufHeight = min(BufHeight, TILE_HEIGHT);
		}
		//列末尾小块处理
		if (icol == tilecols - 1)
		{
			BufWidth = nWidth - (tilecols - 1) * TILE_WIDTH;
			BufWidth = min(BufWidth, TILE_WIDTH);
		}

		double gsd_scale;

		double sGSD = 1.0;// = min(theSlaveProjection->getMetersPerPixel().x, theSlaveProjection->getMetersPerPixel().y);
		if (NULL != theSlaveProjection)
		{
			sGSD = min(theSlaveProjection->getMetersPerPixel().x, theSlaveProjection->getMetersPerPixel().y);
		}

		// master
		ossimIrect mrect;
		//if (0 == strcmp(masterApp.getGetProjectionRef(), ""))
		if (NULL == handlerM->getImageGeometry().get() || NULL == (theMasterProjection = handlerM->getImageGeometry()->getProjection()))
		{
			// 无投影
			gsd_scale = 1.0;
		}
		else
		{
			double mGSD = theMasterProjection->getMetersPerPixel().x;
			//double mGSD = masterApp.getGeoTransform()[1];
			gsd_scale = mGSD / sGSD;
		}

		ossimTDpt tp;
		// slave
		cv::Mat slaveMat;
		ossimIrect srect(ossimIpt(samp, line), ossimIpt(samp + BufWidth - 1, line + BufHeight - 1));


		if (true == runMatchNcc(srect, gsd_scale, tp, this, theDebug))
		{
			if (NULL == handlerM->getImageGeometry().get()
				|| NULL == (theMasterProjection = handlerM->getImageGeometry()->getProjection())
				|| thePointType == point_type::tie)
				//if (0 == strcmp(masterApp.getGetProjectionRef(), "")
				//	|| thePointType == point_type::tie)
			{
				tp.setMasterPoint(tp.getMasterPoint());
				tp.setSlavePoint(tp.getSlavePoint());
				// 无投影
				ossimRefPtr<radiBlockTieGpt> tgi(new radiBlockTieGpt);
				ossimDpt dpt = tp.getMasterPoint();
				tgi->setGroundPoint(ossimGpt(dpt.y, dpt.x, 0.0));
				//set slave image position
				tgi->refImagePoint() = tp.getSlavePoint();
				tgi->setSlaveId(theSlaveId);
				tgi->setMasterId(theMasterId);
				if (thePointType == point_type::tie)
				{
					tgi->m_DptList.push_back(pair<int, ossimDpt>(theSlaveId, tp.getSlavePoint()));
					tgi->m_DptList.push_back(pair<int, ossimDpt>(theMasterId, tp.getMasterPoint()));
					tgi->setPointType(radiBlockTieGpt::unknown_tie_image_points);
				}
				else
				{
					tgi->setPointType(radiBlockTieGpt::known_ground_control_points);
				}
				//set score
				tgi->setScore(tp.score);
				//add to list
				tSet.addTiePoint(tgi);
			}
			else
			{
				tp.setMasterPoint(tp.getMasterPoint());
				tp.setSlavePoint(tp.getSlavePoint());
				// convert "Image to Image" to "Ground to Image" tie points    //TBC : use more generic tie points
				ossimRefPtr<radiBlockTieGpt> tgi(new radiBlockTieGpt);
				ossimGpt ll;
				theMasterProjection->lineSampleToWorld(tp.getMasterPoint(), ll);
				ll.hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(ll);
				tgi->setGroundPoint(ll);
				//set slave image position
				tgi->refImagePoint() = tp.getSlavePoint();
				tgi->setSlaveId(theSlaveId);
				tgi->setMasterId(theMasterId);
				//if (thePointType == point_type::tie)
				//{
				//	tgi->setPointType(radiBlockTieGpt::unknown_tie_image_points);
				//}
				//else
				//{
				tgi->setPointType(radiBlockTieGpt::known_ground_control_points);
				//}
				//set score
				tgi->setScore(tp.score);
				//add to list
				tSet.addTiePoint(tgi);
			}
			break;
		}
	}
	return false;
}

bool radiImageRegistration::getAllFeatures()
{
	static const char MODULE[] = "ossimTieGenerator::getAllFeatures";
	
	if (!theSlaveBandSelector)
	{
		ossimNotify(ossimNotifyLevel_WARN)
			<< "WARN ossimTieGenerator::scanForEdges():"
			<< "\nInput source is not a ossimImageChip.  Returning..." << std::endl;
		return false;
	}

	int nWidth = theAreaOfInterest.width();
	int nHeight =theAreaOfInterest.height();
	int nPointRequired = thePointNumber;
	// Some constants needed throughout...
	const ossim_int32 START_LINE = theAreaOfInterest.ul().y;
	const ossim_int32 STOP_LINE  = theAreaOfInterest.lr().y;
	const ossim_int32 START_SAMP = theAreaOfInterest.ul().x;
	const ossim_int32 STOP_SAMP  = theAreaOfInterest.lr().x;

	// For percent complete status.
	ossim_int32 tilerows = ceil(sqrt(nPointRequired * nHeight / (double)nWidth ));
	ossim_int32 tilecols = ceil(sqrt(nPointRequired * nWidth / (double)nHeight ));
	const ossim_int32 TILE_HEIGHT = int((nHeight - theTileSize*2) / (tilerows - 1));
	const ossim_int32 TILE_WIDTH = int((nWidth - theTileSize*2) / (tilecols - 1));
	//const ossim_int32 TILE_HEIGHT = ceil(nHeight / (double)tilerows - 0.5);
	//const ossim_int32 TILE_WIDTH = ceil(nWidth / (double)tilecols - 0.5);
	double total_tiles = ((double)tilerows)*tilecols;
	double total_tiles_processed;
	double tiles_processed = 0.0;
	// Set the status message to be "scanning source for edges..."
#if OSSIM_HAS_MPI
	if (ossimMpi::instance()->getRank() == 0)
#endif
	{
		ossimNotify(ossimNotifyLevel_INFO) << "Getting tie points..." << std::endl;
	}


	// loop through all tiles
	// need to use a sequencer for parallelism in the future TBD
	theTiePoints.clear();

	vector<row_col> row_col_List;
	for (int i=0;i<tilerows;++i)
	{
		for (int j=0;j<tilecols;++j )
		{
			row_col_List.push_back(row_col(i,j));
		}
	}

	// Start off with a percent complete at 0...
	//setPercentComplete(0.0);
#if OSSIM_HAS_MPI
	MPI_Bcast(&tiles_processed, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int myid = ossimMpi::instance()->getRank();
	int numprocs = ossimMpi::instance()->getNumberOfProcessors();
#else
	int myid = 0;
	int numprocs = 1;
#endif
	ossimRefPtr<ossimCastTileSourceFilter> slaveCaster = caster[0];
	//for (int i = myid;i < (int)row_col_List.size(); i += numprocs)
	//{
	//	int irow = (int)row_col_List[i].row_idx;

	//	ossim_int32 line=START_LINE+irow*TILE_HEIGHT;
	//	ossim_int32 BufHeight = TILE_HEIGHT;
	//	//列末尾小块处理
	//	if (irow == tilerows-1)
	//	{
	//		BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;
	//		BufHeight = min(BufHeight, TILE_HEIGHT);
	//	}

	//	int icol = (int)row_col_List[i].col_idx;
	//	ossim_int32 samp=START_SAMP+icol*TILE_WIDTH;
	//	ossim_int32 BufWidth = TILE_WIDTH;
	//	//列末尾小块处理
	//	if (icol == tilecols-1)
	//	{
	//		BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;
	//		BufWidth = min(BufWidth, TILE_WIDTH);
	//	}
	//	// Get the tie points
	//	getGridFeaturesParallel(ossimIrect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1)));

	//	// Set the percent complete.
	//	tiles_processed += 1.0;
	//	//MPI_Reduce(&tiles_processed, 
	//	//	&total_tiles_processed, 
	//	//	1, 
	//	//	MPI_DOUBLE, 
	//	//	MPI_SUM, 
	//	//	0, 
	//	//	MPI_COMM_WORLD);
	//	if (myid == 0)
	//	{
	//		setPercentComplete((i+1)/total_tiles*100.0);
	//	}
	//	//if (i+numprocs > (int)row_col_List.size())
	//	//{
	//	//	MPI_Barrier(MPI_COMM_WORLD);
	//	//}
	//	
	//	fflush( stdout );
	//}
	totalBlocks = (int)row_col_List.size();
	if (theThreadNum == 0)
	{
		theThreadNum = OpenThreads::GetNumberOfProcessors() * 2;
	}
	//if(num_threads > totalBlocks) num_threads = totalBlocks;
	//num_threads = 1;
	cout<<"using "<<theThreadNum<<" threads..."<<endl;
	GLOBAL_NUM_THREADS = theThreadNum + 1;
	std::vector<radiMatchRectThread *> threads(theThreadNum);
	OpenThreads::Thread::SetConcurrency(theThreadNum);
	OpenThreads::Thread::Init();
	//totalBlocks = 20; // for debug
	for(int i=0; i<theThreadNum; ++i) {
		threads[i] = new radiMatchRectThread(this);
		vector<ossimIrect> rectList;
		for (int j = i;j < totalBlocks; j += theThreadNum)
		{
			int irow = (int)row_col_List[j].row_idx;

			ossim_int32 line=START_LINE+irow*TILE_HEIGHT;
			ossim_int32 BufHeight = TILE_HEIGHT;
			ossim_int32 BufWidth = TILE_WIDTH;

			//ossim_int32 BufHeight = min(TILE_HEIGHT, nHeight);
			//ossim_int32 BufWidth = min(TILE_WIDTH, nWidth);
			//列末尾小块处理
			if (irow == tilerows-1)
			{
				BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;
				BufHeight = min(BufHeight, TILE_HEIGHT);
			}
			//cout<<BufWidth<<"\t"<<BufHeight<<endl;

			int icol = (int)row_col_List[j].col_idx;
			ossim_int32 samp=START_SAMP+icol*TILE_WIDTH;
			//列末尾小块处理
			if (icol == tilecols-1)
			{
				BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;
				BufWidth = min(BufWidth, TILE_WIDTH);
			}
			rectList.push_back(ossimIrect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1)));
		}
		threads[i]->setRect(rectList);
		threads[i]->start();
	}
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished

	
	//for (int i = myid;i < (int)row_col_List.size(); i += numprocs)
	//{
	//	int irow = (int)row_col_List[i].row_idx;

	//	ossim_int32 line=START_LINE+irow*TILE_HEIGHT;
	//	ossim_int32 BufHeight = TILE_HEIGHT;
	//	//列末尾小块处理
	//	if (irow == tilerows-1)
	//	{
	//		BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;
	//		BufHeight = min(BufHeight, TILE_HEIGHT);
	//	}

	//	int icol = (int)row_col_List[i].col_idx;
	//	ossim_int32 samp=START_SAMP+icol*TILE_WIDTH;
	//	ossim_int32 BufWidth = TILE_WIDTH;
	//	//列末尾小块处理
	//	if (icol == tilecols-1)
	//	{
	//		BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;
	//		BufWidth = min(BufWidth, TILE_WIDTH);
	//	}
	//	// Get the tie points
	//	getGridFeaturesParallel(ossimIrect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1)));

	//	// Set the percent complete.
	//	tiles_processed += 1.0;
	//	//MPI_Reduce(&tiles_processed, 
	//	//	&total_tiles_processed, 
	//	//	1, 
	//	//	MPI_DOUBLE, 
	//	//	MPI_SUM, 
	//	//	0, 
	//	//	MPI_COMM_WORLD);
	//	if (myid == 0)
	//	{
	//		setPercentComplete((i+1)/total_tiles*100.0);
	//	}
	//	//if (i+numprocs > (int)row_col_List.size())
	//	//{
	//	//	MPI_Barrier(MPI_COMM_WORLD);
	//	//}

	//	fflush( stdout );
	//}

#if OSSIM_HAS_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	//writeTiePoints(theTset);
//	ossim_int32 i,j;
//#pragma omp parallel for
//	//for (i=0;(i<tilerows)&&!needsAborting();++i)
//	for (i=0;i<tilerows;++i)
//	{
//		ossim_int32 line=START_LINE+i*TILE_HEIGHT;
//		ossim_int32 BufHeight = TILE_HEIGHT;
//		//列末尾小块处理
//		if (i == tilerows-1)
//		{
//			BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;//得出当前块的宽度Bufsizex，高度Bufsizey
//			BufHeight = min(BufHeight, TILE_HEIGHT);
//		}
//		//for (j=0;(j<tilecols)&&!needsAborting();++j )
//		for (j=0;j<tilecols;++j )
//		{
//			ossim_int32 samp=START_SAMP+j*TILE_WIDTH;
//			ossim_int32 BufWidth = TILE_WIDTH;
//			//列末尾小块处理
//			if (j == tilecols-1)
//			{
//				BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;//得出当前块的宽度Bufsizex，高度Bufsizey
//				BufWidth = min(BufWidth, TILE_WIDTH);
//			}
//			// Get the tie points
//
//#pragma omp critical
//			{
//				getGridFeaturesParallel(ossimIrect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1)));
//			}
//
//			// Set the percent complete.
//			tiles_processed += 1.0;
//			setPercentComplete(tiles_processed/total_tiles*100.0);
//
//		}
	//	}
	if (myid == 0)
	{
		//setPercentComplete(100.0);
	}

	// mergeTiePoints
	theTset.clearTiePoints();
	for (int i = 0;i < theThreadNum;++i)
	{
		//cout<<"thread "<<i+1<<": "<<threads[i]->m_theTset.getTiePoints().size()<<endl;
		threads[i]->mergeTiePoints(theTset);
	}
	return true;
}

void radiImageRegistration::writeTiePoints(const radiBlockTieGptSet& tp)
{
	ossimMapProjection* outProjection = radiImageRegistration::getOutputProjection();
	vector< ossimRefPtr<radiBlockTieGpt> >::const_iterator it;
	int icount = 0;
	for (it = tp.getTiePoints().begin();it!=tp.getTiePoints().end();++it)
	{
		char buf[1024];
		double hgt = ossimElevManager::instance()->getHeightAboveEllipsoid((*it)->getGroundPoint());
		ossimDpt dpt = (*it)->getImagePoint();
		ossimGpt gpt = (*it)->getGroundPoint();
		if (!outProjection->isGeographic())
		{
			ossimDpt d1 = outProjection->forward(gpt);
			gpt = ossimGpt(d1.x, d1.y, hgt);
		}
		sprintf_s(buf, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n\0", 
			icount+=1,
			dpt.x, 
			dpt.y,
			gpt.lon,
			gpt.lat,
			hgt);
#if OSSIM_HAS_MPI
		int nLength = strlen(buf);
		MPIO_Request request;
		//MPI_File_iwrite_shared(theFileStream,
		//	buf,
		//	nLength,
		//	MPI_CHAR,
		//	&request);
		MPI_Status status;
		MPI_File_write_shared(theFileStream,
			buf,
			nLength,
			MPI_CHAR,
			&status);
#else
		theFileStream<<buf<<endl;
#endif
	}
}

void radiImageRegistration::setOutputName(const ossimString& filename)
{
	ossimOutputSource::setOutputName(filename);

	if (isOpen()) close();

	if (filename != "")
	{
		theFilename = filename;
	}
}

void radiImageRegistration::setAreaOfInterest(const ossimIrect& rect)
{
	theAreaOfInterest = rect;
}

bool radiImageRegistration::isOpen()const
{
#if OSSIM_HAS_MPI
	return (theFileStream != NULL);
#else
	return const_cast<fstream*>(&theFileStream)->is_open();
#endif
}

bool radiImageRegistration::open()
{
	if(isOpen())
	{
		close();
	}

	if (theFilename == ossimFilename::NIL)
	{
		return false;
	}
#if OSSIM_HAS_MPI
	MPI_File_delete((char*)theFilename.c_str(), MPI_INFO_NULL);
	MPI_File_open(MPI_COMM_WORLD,
		(char*)theFilename.c_str(),
		MPI_MODE_CREATE |  MPI_MODE_WRONLY,
		MPI_INFO_NULL,
		&theFileStream);
#else
	theFileStream.open(theFilename.c_str());
#endif
	
	return theFileStream.is_open();
}

void radiImageRegistration::close()
{
	//if (isOpen()) theFileStream.close();
	if (isOpen())
	{
#if OSSIM_HAS_MPI
		MPI_File_close(&theFileStream);
		theFileStream = NULL;
		if (ossimMpi::instance()->getRank() == 0)
#else
		theFileStream.close();
#endif
		{
			fstream ofs;
			ofs.open(theFilename.setExtension("geom").c_str(), ios_base::out);
			ossimKeywordlist prjKwl;
			getOutputProjection()->saveState(prjKwl);
			ofs<<prjKwl;
			ofs.close();
		}
	}
}
