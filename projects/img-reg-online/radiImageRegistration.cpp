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
#include <boost/filesystem.hpp>
#include <curl/curl.h>
using namespace std;

static int matched_counter = 0;
//#include "radiMatchRect.h"

//#pragma comment(lib, "opencv_core248.lib")
//#pragma comment(lib, "opencv_features2d248.lib")
//#pragma comment(lib, "opencv_flann248.lib")
//#pragma comment(lib, "opencv_highgui248.lib")
//#pragma comment(lib, "opencv_nonfree248.lib")
#pragma comment(lib, "vl.lib")

RTTI_DEF2(radiImageRegistration, "radiImageRegistration", ossimOutputSource, ossimProcessInterface);

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
   theUseGeographic(false),
   theDebug(true),
   theThreadNum(0),
   theTset(),
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
	theThreadNum(a.theThreadNum),
	theSlaveProjection(a.theSlaveProjection),
	theTset(),
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

ossimDpt gsd_degree(int nZoolLevel = 0)
{
	//Map resolution = 156543.04 meters / pixel * cos(latitude) / (2 ^ zoomlevel)
	double x = 1.40625 / (1 << nZoolLevel);
	double y = 0.703125 / (1 << nZoolLevel);
	return ossimDpt(x,y);
}

ossimDpt gsd_meter(double lat = 0, int nZoolLevel = 0)
{
	//Map resolution = 156543.04 meters / pixel * cos(latitude) / (2 ^ zoomlevel)
	double d = 156543.033900000 * cos(lat*0.0174532925) / (1 << nZoolLevel);
	return ossimDpt(d, d);
}

int getZoomLevel(double res, double lat = 0)
{
	//return (int)(log(156543.033900000 * cos(lat*0.0174532925) / res) / log(2) + 0.5);
	//return (int)(log(156543.033900000 * cos(lat*0.0174532925) / res) / log(2));
	return ceil(log(156543.033900000 * cos(lat*0.0174532925) / res) / log(2.0));
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
	projection->lineSampleHeightToWorld(imageRect.ul(), 0.0, p[0]);
	projection->lineSampleHeightToWorld(imageRect.ur(), 0.0, p[1]);
	projection->lineSampleHeightToWorld(imageRect.lr(), 0.0, p[2]);
	projection->lineSampleHeightToWorld(imageRect.ll(), 0.0, p[3]);

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
	int numtasks, taskid, len;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	bool result=true;

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
	if (NULL == handlerS->getImageGeometry().get() || NULL == handlerS->getImageGeometry()->getProjection())
	{
		// 无投影
		cerr << "radiImageRegistration::execute can't get control points as"
			"the slave image does not contain geomatical information" << endl;
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

	if (0 == taskid)
	{
		cout<<"Using band "<<sb+1<<" for slave"<<endl; //TBR
	}
	theSlaveBandSelector = new ossimBandSelector;
	theSlaveBandSelector->connectMyInputTo(0, handlerS.get());
	theSlaveBandSelector->setOutputBandList(vector<ossim_uint32>(1,sb));
	//      theSlaveSource = theSlaveBandSelector;
	theSChain->add(theSlaveBandSelector.get());

	//init casters
	caster[0]->setOutputScalarType(OSSIM_UCHAR);
	caster[1]->setOutputScalarType(OSSIM_UCHAR);

	//init gen
	setStoreFlag(true); //also keep all tie points in memory

	//TBD : set area of interest to buffer around slave?

	// -- 3 -- tie blocks, from sources to outputs

	caster[0]->connectMyInputTo(0, theSlaveBandSelector.get());
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
	//if (theMaster.ext().upcase() != "SHP" &&
	//	!(NULL == mHandler->getImageGeometry().get() || NULL == (mProjection = mHandler->getImageGeometry()->getProjection()).get()))
	//{
	//	ossimIrect mBoundary = mHandler->getBoundingRect(0);
	//	ossimDrect mWorldRect = imageRect2World(mBoundary, mProjection.get());
	//	ossimDrect sWorldRect = imageRect2World(theAreaOfInterest, theSlaveProjection);
	//	//sWorldRect = sWorldRect.clipToRect(mWorldRect);
	//	sWorldRect = worldRectIntersection(sWorldRect, mWorldRect);
	//	//theAreaOfInterest = theAreaOfInterest.clipToRect(worldRect2Image(sWorldRect, theSlaveProjection));
	//	if (sWorldRect.hasNans())
	//	{
	//		theAreaOfInterest.makeNan();
	//	}
	//	else
	//	{
	//		ossimGpt gpt;
	//		ossimDpt mul, mlr;
	//		if (point_type::control == thePointType)
	//		{
	//			mProjection->lineSampleToWorld(mBoundary.ul(), gpt);
	//			theSlaveProjection->worldToLineSample(gpt, mul);
	//			mProjection->lineSampleToWorld(mBoundary.lr(), gpt);
	//			theSlaveProjection->worldToLineSample(gpt, mlr);
	//		}
	//		else
	//		{
	//			mProjection->lineSampleHeightToWorld(mBoundary.ul(), 0.0, gpt);
	//			theSlaveProjection->worldToLineSample(gpt, mul);
	//			mProjection->lineSampleHeightToWorld(mBoundary.lr(), 0.0, gpt);
	//			theSlaveProjection->worldToLineSample(gpt, mlr);
	//		}
	//		//mul = mul + ossimDpt(-theSlaveAccuracy*2, -theSlaveAccuracy*2);
	//		//mlr = mlr + ossimDpt(theSlaveAccuracy*2, theSlaveAccuracy*2);
	//		ossimIrect mBoundary2s(mul.x, mul.y, mlr.x, mlr.y);
	//		theAreaOfInterest = theAreaOfInterest.clipToRect(mBoundary2s);
	//	}
	//}
	if(theAreaOfInterest.hasNans() || theAreaOfInterest.width() <= 1 || theAreaOfInterest.height() <= 1)
	{
		cerr<<"no overlap region is found."<<endl;
		return false;
	}

	//open();

	theTset.clearTiePoints();
	// -- 4 -- run
	result = getAllFeatures();

	if (0 == taskid)
	{
		cout << theTset.getTiePoints().size() << " tie points are found." << endl;
	}

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
		ofs.open(filename.c_str(), ios_base::app);
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
		ofs.open(filename.c_str(), ios_base::app);
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
	ofs.open(filename.c_str(), ios_base::out);
	ofs.close();
	appendTiePoints(filename);
}
void radiImageRegistration::writeControlPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
{
	fstream ofs;
	ofs.open(filename.c_str(), ios_base::out);
	ofs.close();
	appendControlPoints(filename, pMapProjection);
}

void radiImageRegistration::appendPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
{
	appendControlPoints(filename, pMapProjection);
}

void radiImageRegistration::writePoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
{
	writeControlPoints(filename, pMapProjection);
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
	//int noctaves = 2, nlevels = 4, o_min = 0;
	int noctaves = 4, nlevels = 4, o_min = 0;
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
	//double edge_thresh = 10 ;  //-1 will use the default (as in matlab)
	//double peak_thresh = 0.04;
	double edge_thresh = 10;  //-1 will use the default (as in matlab)
	double peak_thresh = 0.05;
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
					kpt.angle = TemptAngle;
					kpt.response = TemptKeyPoint.contrast;
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


float get_sub_pix(cv::Mat const &img, cv::Point2f const &pt)
{
	int x = static_cast<int>(pt.x);
	int y = static_cast<int>(pt.y);

	int x0 = cv::borderInterpolate(x, img.cols, cv::BORDER_REPLICATE);
	int x1 = cv::borderInterpolate(x + 1, img.cols, cv::BORDER_REPLICATE);
	int y0 = cv::borderInterpolate(y, img.rows, cv::BORDER_REPLICATE);
	int y1 = cv::borderInterpolate(y + 1, img.rows, cv::BORDER_REPLICATE);

	float a = pt.x - x;
	float c = pt.y - y;

	float x1_interpolate = (img.at<uchar>(y0, x0) * (1.0 - a)
		+ img.at<uchar>(y0, x1) * a);
	float x2_interpolate = (img.at<uchar>(y1, x0) * (1.0 - a)
		+ img.at<uchar>(y1, x1) * a);
	float target = x1_interpolate * (1.0 - c) + x2_interpolate * c;

	return target;
}

struct LSM_STRUCT
{
	cv::Mat templateMat;
	cv::Mat searchMat;
};

void funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	// a1 a2 a3 b1 b2 b3 k1 k2
	// x' = a1*x + a2*y + a3
	// y' = b1*x + b2*y + b3
	// f(x,y) = k1*g(x',y') + k2
	LSM_STRUCT *pThis = (LSM_STRUCT*)adata;
	int nX = pThis->templateMat.cols;
	int nY = pThis->templateMat.rows;

	assert(nequation == (nX*nY));

	int pos = 0;
	int i;
	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;

	for (size_t ix = 0; ix < nX; ix++)
	{
		for (size_t iy = 0; iy < nY; iy++)
		{
			float svalue = get_sub_pix(pThis->templateMat, cv::Point2f(ix, iy));
			double mx = param[0] * ix + param[1] * iy + param[2];
			double my = param[3] * ix + param[4] * iy + param[5];
			float mvalue = get_sub_pix(pThis->searchMat, cv::Point2f(mx, my));
			hx[pos++] = (svalue - param[6] * mvalue - param[7]);
		}
	}
}

void jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata)
{
	// a1 a2 a3 b1 b2 b3 k1 k2
	// x' = a1*x + a2*y + a3
	// y' = b1*x + b2*y + b3
	// f(x,y) = k1*g(x',y') + k2
	LSM_STRUCT *pThis = (LSM_STRUCT*)adata;
	int nX = pThis->templateMat.cols;
	int nY = pThis->templateMat.rows;

	assert(nequation == (nX*nY));

	int pos = 0;
	int i;

	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;
	int c = 0;
	for (size_t ix = 0; ix < nX; ix++)
	{
		for (size_t iy = 0; iy < nY; iy++)
		{
			float svalue = get_sub_pix(pThis->templateMat, cv::Point2f(ix, iy));
			for (int p = 0; p<6; ++p)
			{
				// a1 a2 a3 b1 b2 b3
				double middle = param[p];
				param[p] = middle + pstep_scale;
				double mx = param[0] * ix + param[1] * iy + param[2];
				double my = param[3] * ix + param[4] * iy + param[5];
				float mvalue1 = get_sub_pix(pThis->searchMat, cv::Point2f(mx, my));
				float res1 = svalue - param[6] * mvalue1 - param[7];

				param[p] = middle - pstep_scale;
				mx = param[0] * ix + param[1] * iy + param[2];
				my = param[3] * ix + param[4] * iy + param[5];
				float mvalue2 = get_sub_pix(pThis->searchMat, cv::Point2f(mx, my));
				float res2 = svalue - param[6] * mvalue2 - param[7];

				j[c++] = (res1 - res2)*den;
				param[p] = middle;
			}
			double mx = param[0] * ix + param[1] * iy + param[2];
			double my = param[3] * ix + param[4] * iy + param[5];
			float mvalue = get_sub_pix(pThis->searchMat, cv::Point2f(mx, my));

			// k1
			j[c++] = -mvalue;
			// k2
			j[c++] = -1.0;
		}
	}
}

void leastSquareMatching(cv::Mat templateMat, cv::Mat searchMat, double lsm_model[8])
{
	//cv::Mat img_G0;
	//cv::Mat img_G1;

	//cv::GaussianBlur(templateMat, img_G0, Size(3, 3), 1.6);
	////cv::GaussianBlur(img_G0, img_G1, Size(3, 3), 0);
	////templateMat = img_G0 - img_G1;
	//templateMat = cv::abs(templateMat - img_G0);
	////normalize(templateMat, templateMat, 255, 0, CV_MINMAX);

	//cv::GaussianBlur(searchMat, img_G0, Size(3, 3), 1.6);
	////cv::GaussianBlur(img_G0, img_G1, Size(3, 3), 0);
	////searchMat = img_G0 - img_G1;
	//searchMat = cv::abs(searchMat - img_G0);
	////normalize(searchMat, searchMat, 255, 0, CV_MINMAX);

	int nX = templateMat.cols;
	int nY = templateMat.rows;

	int nparam = 8;

	arma::vec outParameters(nparam);
	double *x = new double[nX*nY];
	for (int i = 0; i<nX*nY; i++) x[i] = 0.0;

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0] = LM_INIT_MU; opts[1] = 1E-15; opts[2] = 1E-15; opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	//int ret = dlevmar_dif(funcErrorEquation, p, x, nparam, nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
	LSM_STRUCT lsm_data;
	lsm_data.searchMat = searchMat;
	lsm_data.templateMat = templateMat;
	int ret = dlevmar_der(funcErrorEquation, jacErrorEquation, lsm_model, x, nparam, nX*nY, 1000, opts, info, NULL, NULL, &lsm_data); // with analytic Jacobian
	delete[]x;
	x = NULL;
	int ix = nX / 2;
	int iy = nY / 2;
}

int radiImageRegistration::runMatchParallel(const cv::Mat& slaveMat, const cv::Mat& masterMat, ossimTDpt& tDpt, void* pData, bool bDebug)
{
	radiImageRegistration* pThis = (radiImageRegistration*)pData;
	// detect corners
	//cv::initModule_nonfree();
	cv::initModule_features2d();
	std::vector<KeyPoint> skeypoints, mkeypoints;
	Mat sdescriptors, mdescriptors;
	//SIFT detector(pThis->theSiftNfeatures, pThis->theSiftNOctaveLayers, 
	//	pThis->theSiftContrastThreshold, pThis->theSiftEdgeThreshold, pThis->theSiftSigma);
	//if (bDebug)
	//{
	//	//int pid = getpid();
	//	//string strId = string(QString::number(pid).toLatin1());
	//	//cv::imwrite("slave_"+strId+".png", slaveMat);
	//	//cv::imwrite("master_"+strId+".png", masterMat);
	//	cv::imwrite("slave.png", slaveMat);
	//	cv::imwrite("master.png", masterMat);
	//}
	if (bDebug)
	{
		cv::imwrite("slave.png", slaveMat);
		cv::imwrite("master.png", masterMat);

		cv::Mat img_outImage;
		drawTogether(slaveMat, masterMat, img_outImage);
		cv::imwrite("image.png", img_outImage);
	}

	//// detect
	////detector.detect( slaveMat, skeypoints );
	//if (skeypoints.size() < 10 )
	//{
	//	return match_state::slave_faild;

	//}
	//// detect
	//detector.detect( masterMat, mkeypoints );
	//if(mkeypoints.size() < 10)
	//{
	//	//handlerM->close();
	//	return match_state::master_faild;
	//}

	//// extract
	//cv::SiftDescriptorExtractor extractor;
	//extractor.compute( slaveMat, skeypoints, sdescriptors );
	//extractor.compute( masterMat, mkeypoints, mdescriptors );


	VLFeatSift(slaveMat, skeypoints, sdescriptors);
	if (skeypoints.size() < 10 )
	{
		return match_state::slave_faild;

	}
	VLFeatSift(masterMat, mkeypoints, mdescriptors);
	if(mkeypoints.size() < 10)
	{
		//handlerM->close();
		return match_state::master_faild;
	}

	BFMatcher matcher(NORM_L1, false); 
	vector< vector< DMatch >  > matches;
	matcher.knnMatch( sdescriptors, mdescriptors, matches, 2 );
	// inverse mathcing
	vector< vector< DMatch >  > matches2;
	matcher.knnMatch(mdescriptors, sdescriptors, matches2, 2);

	// "cross-matching" and "first and second minimum distances ratio test"
	vector< DMatch > good_matches;
	for (size_t i = 0; i < matches.size(); i++)
	{
		if (matches[i][0].distance / (matches[i][1].distance+FLT_EPSILON) < 0.6)
		{
			good_matches.push_back(matches[i][0]);
			continue;
		}
		int queryIdx = matches[i][0].queryIdx;
		int trainIdx = matches[i][0].trainIdx;
		for (size_t j = 0; j < matches2.size(); j++)
		{
			int queryIdx2 = matches2[j][0].trainIdx;
			int trainIdx2 = matches2[j][0].queryIdx;

			if (queryIdx == queryIdx2 && trainIdx == trainIdx2)
			{
				good_matches.push_back(matches[i][0]);
				break;
			}
		}
	}

	////vector< DMatch > good_matches;
	////findGoodMatches(matches, good_matches, 0.75f);

	//float std_scale_diff_threshold = 0.2f;
	//for (size_t i = 0; i < good_matches.size();)
	//{
	//	float s1 = mkeypoints[good_matches[i].trainIdx].size;
	//	float s2 = skeypoints[good_matches[i].queryIdx].size;
	//	float std_scale_diff = (s1 - s2) / sqrt(s1*s1 + s2*s2 + FLT_EPSILON);
	//	if (std_scale_diff > std_scale_diff_threshold)
	//	{
	//		good_matches.erase(good_matches.begin() + i);
	//		continue;
	//	}
	//	i++;
	//}
	
	if (good_matches.size() < 3)
	{
		//handlerM->close();
		return match_state::match_failed;
	}

	//-- Create input data
	Eigen::MatrixXd dataPoints((int)good_matches.size(), 4);
	for(unsigned int i = 0;i < good_matches.size();++i)
	{
		dataPoints(i, 0) = skeypoints[good_matches[i].queryIdx].pt.x;
		dataPoints(i, 1) = skeypoints[good_matches[i].queryIdx].pt.y;
		dataPoints(i, 2) = mkeypoints[good_matches[i].trainIdx].pt.x;
		dataPoints(i, 3) = mkeypoints[good_matches[i].trainIdx].pt.y;
	}
	
	// RANSAC detect outliers
	auto_ptr< estimators::Solver<Eigen::MatrixXd, Eigen::VectorXd> > ptrSolver(
		new estimators::rigidSolver<Eigen::MatrixXd, Eigen::VectorXd>);
	vector<int> inliers;
	//for (int i = 0; i < (int)good_matches.size(); i++) inliers.push_back(i);
	vector<Eigen::VectorXd> models;

	ransac::Ransac_Handler ransac_fun_Handler;
	bool result = ransac::Ransac_RobustEstimator
		(
		dataPoints, // the input data
		estimators::rigidSolver<Eigen::MatrixXd, Eigen::VectorXd>::extractor, // How select sampled point from indices
		dataPoints.rows(),  // the number of putatives data
		*(ptrSolver.get()),  // compute the underlying model given a sample set
		estimators::rigidSolver<Eigen::MatrixXd, Eigen::VectorXd>::defaultEvaluator,  // the function to evaluate a given model
		//Ransac Object that contain function:
		// CandidatesSelector, Sampler and TerminationFunction
		ransac_fun_Handler, // the basic ransac object
		1000,  // the maximum rounds for RANSAC routine
		inliers, // inliers to the final solution
		models, // models array that fit input data
		0.9//0.95 // the confidence want to achieve at the end
		);
	double s = sqrt(models[0][0] * models[0][0] + models[0][1] * models[0][1]);
	double angle = atan(models[0][1] / models[0][0]);

	double angle_threshold = 0.5; // rad -- 28.6478897565deg
	
	if (fabs(angle) > PI*0.25)
	{
		return match_state::match_failed;
	}
	for (size_t i = 0; i < inliers.size(); i++)
	{
		float a1 = mkeypoints[good_matches[inliers[i]].trainIdx].angle;
		float a2 = skeypoints[good_matches[inliers[i]].queryIdx].angle;
		double angle_diff = a1 - a2;
		if (fabs(angle_diff + angle) > angle_threshold)
		{
			//return match_state::match_failed;
			inliers.erase(inliers.begin() + i);
			i--;
		}
	}
	if (inliers.size() < 3)
	{
		//handlerM->close();
		return match_state::match_failed;
	}

	if (fabs(s - 1) > 0.5)
	{
		//std::vector<double> scale_ratioList;
		//std::vector<double> scale_diffList;
		//std::vector<double> angle_diffList;
		//for (size_t i = 0; i < inliers.size(); i++)
		//{
		//	float s1 = mkeypoints[good_matches[inliers[i]].trainIdx].size;
		//	float s2 = skeypoints[good_matches[inliers[i]].queryIdx].size;
		//	scale_ratioList.push_back(s1 / (s2 + FLT_EPSILON));
		//	scale_diffList.push_back((s1 - s2) / sqrt(s1*s1 + s2*s2 + FLT_EPSILON));

		//	float a1 = mkeypoints[good_matches[inliers[i]].trainIdx].angle;
		//	float a2 = skeypoints[good_matches[inliers[i]].queryIdx].angle;
		//	angle_diffList.push_back(a1 - a2);
		//}
		if (bDebug)
		{
			//cout<<models[0]<<endl;
			vector< DMatch > final_matches;
			for (int i = 0; i < (int)inliers.size(); ++i)
			{
				final_matches.push_back(good_matches[inliers[i]]);
			}
			// Draw matches
			cv::Mat imgMatch;
			drawMatches(slaveMat, skeypoints, masterMat, mkeypoints, final_matches, imgMatch);
			cv::imwrite("result.png", imgMatch);

			cv::Mat outMat;
			_prepareImgAndDrawKeylines(slaveMat,
				masterMat,
				tDpt,
				outMat,
				cv::Scalar(0, 0, 255));
			cv::imwrite("matched.png", outMat);
		}
		return match_state::match_failed;
	}

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
	//if (inliers.size() < 3)
	//{
	//	//handlerM->close();
	//	return match_state::match_failed;
	//}

	////if (fabs(models[0][1] * models[0][5] - models[0][2] * models[0][4]) < 0.5)
	////{
	////	//handlerM->close();
	////	return match_state::match_failed;
	////}
	//double delta_energy = 1.0 - fabs(models[0][1] * models[0][5] - models[0][2] * models[0][4]);
	//if (fabs(delta_energy) > 0.2)
	//{
	//	std::vector<double> scale_ratioList;
	//	std::vector<double> scale_diffList;
	//	for (size_t i = 0; i < inliers.size(); i++)
	//	{
	//		float s1 = mkeypoints[good_matches[inliers[i]].trainIdx].size;
	//		float s2 = skeypoints[good_matches[inliers[i]].queryIdx].size;
	//		scale_ratioList.push_back(s1 / (s2 + FLT_EPSILON));
	//		scale_diffList.push_back((s1 - s2) / sqrt(s1*s1 + s2*s2 + FLT_EPSILON));
	//	}
	//	//handlerM->close();
	//	return match_state::match_failed;
	//}


	int slaveTemplateSize = 9;	// should be odd
	int semiTemplateSize = slaveTemplateSize / 2; // masterTemplateSize is odd
	int searchSize = 15;	// should be odd
	int semisearchSize = searchSize / 2; // masterTemplateSize is odd
	// maybe a lot of correspondences are found, but we need only one correspondence for a tile
	// whose contrast is the largest?
	float maxId = inliers[0];
	float maxContrast = mkeypoints[good_matches[inliers[0]].trainIdx].response + skeypoints[good_matches[inliers[0]].queryIdx].response;
	for (int i = 1; i < (int)inliers.size(); ++i)
	{
		double c = mkeypoints[good_matches[inliers[i]].trainIdx].response + skeypoints[good_matches[inliers[i]].queryIdx].response;

		int ix = (int)(mkeypoints[good_matches[inliers[i]].trainIdx].pt.x + 0.5);
		int iy = (int)(mkeypoints[good_matches[inliers[i]].trainIdx].pt.y + 0.5);
		int isx = (int)(skeypoints[good_matches[inliers[i]].queryIdx].pt.x + 0.5);
		int isy = (int)(skeypoints[good_matches[inliers[i]].queryIdx].pt.y + 0.5);
		if (maxContrast < c)
		{
			if (ix < semisearchSize || iy < semisearchSize || ix >= masterMat.cols - semisearchSize || iy >= masterMat.rows - semisearchSize)
			{
				continue;
			}
			if (isx < semiTemplateSize || isy < semiTemplateSize || isx >= slaveMat.cols - semiTemplateSize || isy >= slaveMat.rows - semiTemplateSize)
			{
				continue;
			}

			maxContrast = c;
			maxId = inliers[i];
		}
	}

	int imx = (int)(mkeypoints[good_matches[maxId].trainIdx].pt.x + 0.5);
	int imy = (int)(mkeypoints[good_matches[maxId].trainIdx].pt.y + 0.5);
	int isx = (int)(skeypoints[good_matches[maxId].queryIdx].pt.x + 0.5);
	int isy = (int)(skeypoints[good_matches[maxId].queryIdx].pt.y + 0.5);

	if (imx < semisearchSize || imy < semisearchSize || imx >= masterMat.cols - semisearchSize || imy >= masterMat.rows - semisearchSize)
	{
		semisearchSize = min(min(abs(imx), abs(imy)),
			min(abs(imx - (masterMat.cols - 1)), abs(imy - (masterMat.rows - 1))));
		searchSize = semisearchSize * 2 + 1;
	}
	if (isx < semiTemplateSize || isy < semiTemplateSize || isx >= slaveMat.cols - semiTemplateSize || isy >= slaveMat.rows - semiTemplateSize)
	{
		semiTemplateSize = min(min(abs(isx), abs(isy)),
			min(abs(isx - (slaveMat.cols - 1)), abs(isy - (slaveMat.rows - 1))));

		slaveTemplateSize = semiTemplateSize * 2 + 1;
	}
	// matching precisely (subpixel)
	cv::Mat mSearchMat = masterMat(cv::Rect(imx - semisearchSize, imy - semisearchSize, searchSize, searchSize));
	cv::Mat sTemplateMat = slaveMat(cv::Rect(isx - semiTemplateSize, isy - semiTemplateSize, slaveTemplateSize, slaveTemplateSize));
	double lsm_model[8];
	// a1 a2 a3 b1 b2 b3 k1 k2
	// x' = a1*x + a2*y + a3
	// y' = b1*x + b2*y + b3
	// f(x,y) = k1*g(x',y') + k2
	//lsm_model[0] = 1.0;
	//lsm_model[1] = 0.0;
	//lsm_model[2] = 0.0 - semiTemplateSize + semisearchSize;
	//lsm_model[3] = 0.0;
	//lsm_model[4] = 1.0;
	//lsm_model[5] = 0.0 - semiTemplateSize + semisearchSize;
	//lsm_model[6] = 1.0;
	//lsm_model[7] = 0.0;
	lsm_model[0] = models[0][0];
	lsm_model[1] = models[0][1];
	lsm_model[2] = models[0][2];
	lsm_model[3] = -models[0][1];
	lsm_model[4] = models[0][0];
	lsm_model[5] = models[0][3];
	lsm_model[6] = 1.0;
	lsm_model[7] = 0.0;
	lsm_model[2] += lsm_model[0] * (isx - semiTemplateSize) + lsm_model[1] * (isy - semiTemplateSize) - (imx - semisearchSize);
	lsm_model[5] += lsm_model[3] * (isx - semiTemplateSize) + lsm_model[4] * (isy - semiTemplateSize) - (imy - semisearchSize);

	if (bDebug)
	{
		cv::Mat img_outImage;
		drawTogether(sTemplateMat, mSearchMat, img_outImage);
		cv::imwrite("lsm0.png", img_outImage);
	}
	leastSquareMatching(sTemplateMat, mSearchMat, lsm_model);
	if (bDebug)
	{
		cv::Mat img_outImage;
		drawTogether(sTemplateMat, mSearchMat, img_outImage);
		cv::imwrite("lsm1.png", img_outImage);
	}
	int ix = semiTemplateSize;
	int iy = semiTemplateSize;
	double mx = lsm_model[0] * ix + lsm_model[1] * iy + lsm_model[2];
	double my = lsm_model[3] * ix + lsm_model[4] * iy + lsm_model[5];
	//cout << mx << ", " << my << endl;
	ossimDpt mc = ossimDpt(imx - semisearchSize + mx, imy - semisearchSize + my);
	ossimDpt sc = ossimDpt(isx, isy);

	//ossimDpt mc = ossimDpt(mkeypoints[good_matches[maxId].trainIdx].pt.x, mkeypoints[good_matches[maxId].trainIdx].pt.y);
	//ossimDpt sc = ossimDpt(skeypoints[good_matches[maxId].queryIdx].pt.x, skeypoints[good_matches[maxId].queryIdx].pt.y);
	//////int isx = (int)(skeypoints[good_matches[maxId].queryIdx].pt.x + 0.5);
	//////int isy = (int)(skeypoints[good_matches[maxId].queryIdx].pt.y + 0.5);
	//////ossimDpt sc = ossimDpt(isx, isy);
	double c = mkeypoints[good_matches[maxId].trainIdx].response + skeypoints[good_matches[maxId].queryIdx].response;
	tDpt = ossimTDpt(mc, sc, mkeypoints[good_matches[maxId].trainIdx].response);


	if (bDebug)
	{
		std::vector<double> scale_ratioList;
		std::vector<double> scale_diffList;
		std::vector<double> angle_diffList;
		for (size_t i = 0; i < inliers.size(); i++)
		{
			float s1 = mkeypoints[good_matches[inliers[i]].trainIdx].size;
			float s2 = skeypoints[good_matches[inliers[i]].queryIdx].size;
			scale_ratioList.push_back(s1 / (s2 + FLT_EPSILON));
			scale_diffList.push_back((s1 - s2) / sqrt(s1*s1 + s2*s2 + FLT_EPSILON));

			float a1 = mkeypoints[good_matches[inliers[i]].trainIdx].angle;
			float a2 = skeypoints[good_matches[inliers[i]].queryIdx].angle;
			angle_diffList.push_back(a1 - a2);
		}
		//cout<<models[0]<<endl;
		vector< DMatch > final_matches;
		for (int i = 0; i < (int)inliers.size(); ++i)
		{
			final_matches.push_back(good_matches[inliers[i]]);
		}
		// Draw matches
		cv::Mat imgMatch;
		drawMatches(slaveMat, skeypoints, masterMat, mkeypoints, final_matches, imgMatch);
		cv::imwrite("result.png", imgMatch);

		cv::Mat outMat;
		_prepareImgAndDrawKeylines(slaveMat,
			masterMat,
			tDpt,
			outMat,
			cv::Scalar(0, 0, 255));
		cv::imwrite("matched.png", outMat);


		ossimFilename matchedFolder = "matched";
		if (!matchedFolder.exists())
		{
			_mkdir(matchedFolder.c_str());
		}
		char buf[256];
		sprintf_s(buf, "%s\\matched%04d.png\0", matchedFolder.c_str(), ++matched_counter);
		cv::imwrite(buf, outMat);
		sprintf_s(buf, "%s\\result%04d.png\0", matchedFolder.c_str(), matched_counter);
		cv::imwrite(buf, imgMatch);
	}
	return match_state::success;
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

string getRequestUrl(ossimProjection* sProjection, double res,
	ossimIrect sRect, double slaveAccuracy)
{
	ossimIpt sul = sRect.ul() + ossimIpt(-slaveAccuracy, -slaveAccuracy);
	ossimIpt sur = sRect.ur() + ossimIpt(slaveAccuracy, -slaveAccuracy);
	ossimIpt slr = sRect.lr() + ossimIpt(slaveAccuracy, slaveAccuracy);
	ossimIpt sll = sRect.ll() + ossimIpt(-slaveAccuracy, slaveAccuracy);
	ossimDpt mul, mlr;
	ossimGpt gpt;

	ossimDpt p[4];
	sProjection->lineSampleToWorld(sul, gpt);
	p[0].x = gpt.lon; p[0].y = gpt.lat;
	sProjection->lineSampleToWorld(sur, gpt);
	p[1].x = gpt.lon; p[1].y = gpt.lat;
	sProjection->lineSampleToWorld(slr, gpt);
	p[2].x = gpt.lon; p[2].y = gpt.lat;
	sProjection->lineSampleToWorld(sll, gpt);
	p[3].x = gpt.lon; p[3].y = gpt.lat;

	double xmin = p[0].x;
	double xmax = p[0].x;
	double ymin = p[0].y;
	double ymax = p[0].y;

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

	double clon = (xmin + xmax)*0.5;
	double clat = (ymin + ymax)*0.5;

	int nZoomLevel = getZoomLevel(res, clat);
	ossimDpt res_deg = gsd_degree(nZoomLevel);
	int nWidth = (xmax - xmin) / res_deg.x;
	int nHeight = (ymax - ymin) / res_deg.y;

	char buf[2048];
	sprintf_s(buf, "https://maps.googleapis.com/maps/api/staticmap?center=%lf,%lf&zoom=%d&size=%dx%d&maptype=satellite", 
		clon, clat, nZoomLevel, nWidth, nHeight);
	return string(buf);
}

ossimDrect getLonLatRect(ossimProjection* sProjection,
	ossimIrect sRect, double slaveAccuracy)
{
	ossimIpt sul = sRect.ul() + ossimIpt(-slaveAccuracy, -slaveAccuracy);
	ossimIpt sur = sRect.ur() + ossimIpt(slaveAccuracy, -slaveAccuracy);
	ossimIpt slr = sRect.lr() + ossimIpt(slaveAccuracy, slaveAccuracy);
	ossimIpt sll = sRect.ll() + ossimIpt(-slaveAccuracy, slaveAccuracy);
	ossimDpt mul, mlr;
	ossimGpt gpt;

	ossimDpt p[4];
	sProjection->lineSampleToWorld(sul, gpt);
	p[0].x = gpt.lon; p[0].y = gpt.lat;
	sProjection->lineSampleToWorld(sur, gpt);
	p[1].x = gpt.lon; p[1].y = gpt.lat;
	sProjection->lineSampleToWorld(slr, gpt);
	p[2].x = gpt.lon; p[2].y = gpt.lat;
	sProjection->lineSampleToWorld(sll, gpt);
	p[3].x = gpt.lon; p[3].y = gpt.lat;

	double xmin = p[0].x;
	double xmax = p[0].x;
	double ymin = p[0].y;
	double ymax = p[0].y;

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

	return ossimDrect(xmin, ymax, xmax, ymin);
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

struct MemoryStruct {
	char *memory;
	size_t size;
};

ossimGpt TileXy2LonLat(double x, double y, int nZoomLevel)
{
	double scale = 1.0/(double)(1 << (8 + nZoomLevel));
	double lon = x * 360.0 * scale - 180.0;
	double t = exp((0.5 - y * scale) * 4.0 * PI);
	double lat = asin((t - 1) / (t + 1))*180.0 / PI;

	return ossimGpt(lat, lon);
}

ossimDpt LonLat2TileXy(ossimGpt gpt, int nZoomLevel)
{
	double scale = (double)(1 << (8 + nZoomLevel));
	double sinLat = sin(gpt.lat * PI / 180.0);
	double x = (gpt.lon + 180.0) * scale / 360.0;
	double y = (0.5 - log((1.0 + sinLat) / (1.0 - sinLat)) / (4.0*PI))*scale;
	return ossimDpt(x, y);
}

ossimGpt masterLineSample2World(ossimDpt linesample, 
	ossimGpt ll_center, ossimDpt image_center, 
	int nZoomLevel, double scale = 1.0)
{
	ossimDpt ctileXy = LonLat2TileXy(ll_center, nZoomLevel);
	ossimDpt tileXy = ctileXy + (linesample - image_center) / scale;
	ossimGpt world = TileXy2LonLat(tileXy.x, tileXy.y, nZoomLevel);
	return world;
}


static size_t
WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
	size_t realsize = size * nmemb;
	struct MemoryStruct *mem = (struct MemoryStruct *)userp;

	mem->memory = (char*)realloc(mem->memory, mem->size + realsize + 1);
	if (mem->memory == NULL) {
		/* out of memory! */
		printf("not enough memory (realloc returned NULL)\n");
		return 0;
	}

	memcpy(&(mem->memory[mem->size]), contents, realsize);
	mem->size += realsize;
	mem->memory[mem->size] = 0;

	return realsize;
}


struct MemoryStruct getDataFromUrl(const char* url, const char* proxy="")
{
	CURL *curl_handle;
	CURLcode res;
	int timeout = 5;
	struct MemoryStruct chunk;

	chunk.memory = (char*)malloc(1);  /* will be grown as needed by the realloc above */
	chunk.size = 0;    /* no data at this point */

	curl_global_init(CURL_GLOBAL_ALL);

	/* init the curl session */
	curl_handle = curl_easy_init();

	/* specify URL to get */
	curl_easy_setopt(curl_handle, CURLOPT_URL, url);

	/* send all data to this function  */
	curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);

	/* we pass our 'chunk' struct to the callback function */
	curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)&chunk);

	/* some servers don't like requests that are made without a user-agent
	field, so we provide one */
	//curl_easy_setopt(curl_handle, CURLOPT_USERAGENT, "libcurl-agent/1.0");

	curl_easy_setopt(curl_handle, CURLOPT_CONNECTTIMEOUT, timeout);
	curl_easy_setopt(curl_handle, CURLOPT_PROXYAUTH, CURLAUTH_BASIC); //代理认证模式

	if (0 != strcmp(proxy, ""))
	{
		curl_easy_setopt(curl_handle, CURLOPT_PROXY, proxy);
		curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYPEER, 0L);
		curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYHOST, 0L);
	}

	curl_version_info_data *info = curl_version_info(CURLVERSION_NOW);

	/* get it! */
	res = curl_easy_perform(curl_handle);

	/* check for errors */
	if (res != CURLE_OK) {
		fprintf(stderr, "curl_easy_perform() failed: %s\n",
			curl_easy_strerror(res));
	}
	else {
		/*
		* Now, our chunk.memory points to a memory block that is chunk.size
		* bytes big and contains the remote file.
		*
		* Do something nice with it!
		*/

		//printf("%lu bytes retrieved\n", (long)chunk.size);
	}

	/* cleanup curl stuff */
	curl_easy_cleanup(curl_handle);

	/* we're done with libcurl, so clean it up */
	curl_global_cleanup();

	return chunk;
}

bool radiImageRegistration::getGridFeaturesParallel(const ossimIrect& rect, radiBlockTieGptSet& tSet)
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

	GdalRasterApp slaveApp;
	if (!slaveApp.open(getSlave().c_str()))
	{
		cerr<<"radiImageRegistration"<<"::execute can't open slave image  "<< getSlave().c_str() <<endl;
		return false;
	}
	// select only one band (if multiple)
	ossim_uint32 sbc = slaveApp.nBand();
	//add a band selector
	ossim_uint32 sb = theSlaveBand;
	if (sb>=sbc) 
	{
		cerr<<"radiImageRegistration"<<"::execute Warning not enough bands in slave, only "<< sbc <<endl;
		sb=0;
	}
		
	double sGSD = 1.0;// = min(theSlaveProjection->getMetersPerPixel().x, theSlaveProjection->getMetersPerPixel().y);
	if (NULL != theSlaveProjection)
	{
		sGSD = min(theSlaveProjection->getMetersPerPixel().x, theSlaveProjection->getMetersPerPixel().y);
	}
	//int ncore = omp_get_num_procs();//获取执行核的总数；  目前机器CPU的数量
	row_col_step = 5;
	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);	

	bool found = false;
	int search_size = (int)row_col_List.size();
	//search_size = min(search_size, 1);
	for (int i=0;i < search_size && !found;++i)
	{
		//int icol = int(row_col_List[i+j].col_idx+center_col);
		//int irow = int(row_col_List[i + j].row_idx + center_row);
		int icol = int(row_col_List[i].col_idx);
		int irow = int(row_col_List[i].row_idx);
		ossim_int32 samp = START_SAMP+icol*TILE_WIDTH;
		ossim_int32 line = START_LINE+irow*TILE_HEIGHT;

		ossim_int32 BufHeight = min(TILE_HEIGHT, nHeight);
		ossim_int32 BufWidth = min(TILE_WIDTH, nWidth);
		//行末尾小块处理
		if (irow == tilerows-1)
		{
			BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;
			BufHeight = min(BufHeight, TILE_HEIGHT);
		}
		//列末尾小块处理
		if (icol == tilecols-1)
		{
			BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;
			BufWidth = min(BufWidth, TILE_WIDTH);
		}

		// slave
		cv::Mat slaveMat;
		ossimDpt scenter((samp + samp + BufWidth - 1)*0.5, (line + line + BufHeight - 1)*0.5);


		// get slave rect
		ossimGpt ll_center;
		theSlaveProjection->lineSampleToWorld(scenter, ll_center);
		int nearestZoomLevel = getZoomLevel(sGSD, ll_center.lat);

		vector<ossimTDpt> tp(1);
		ossimDpt res = gsd_meter(ll_center.lat, nearestZoomLevel);
		double mGSD = res.x;
		double scale = mGSD / sGSD;
		nWidth = TILE_WIDTH / scale;
		nHeight = TILE_HEIGHT / scale;
		nWidth = min(nWidth, 600);
		nHeight = min(nHeight, 600);

		ossimDpt mimage_center = ossimDpt(nWidth*0.5, nHeight*0.5);
		ossimGpt mul_ll = masterLineSample2World(ossimDpt(0.0,0.0),
			ll_center, mimage_center, nearestZoomLevel, 1.0);
		ossimGpt mlr_ll = masterLineSample2World(ossimDpt(nWidth, nHeight),
			ll_center, mimage_center, nearestZoomLevel, 1.0);
		ossimGrect mllRect(mul_ll, mlr_ll);
		//double ulLon = mul_ll.lon;
		//double ulLat = mul_ll.lat;
		//double lrLon = mlr_ll.lon;
		//double lrLat = mlr_ll.lat;
		//ossimGrect mllRect(ossimGpt(ulLat, ulLon), ossimGpt(lrLat, lrLon));

		// get slave rect
		ossimDpt p[4];
		theSlaveProjection->worldToLineSample(mllRect.ul(), p[0]);
		theSlaveProjection->worldToLineSample(mllRect.ur(), p[1]);
		theSlaveProjection->worldToLineSample(mllRect.lr(), p[2]);
		theSlaveProjection->worldToLineSample(mllRect.ll(), p[3]);

		double xmin = p[0].x;
		double xmax = p[0].x;
		double ymin = p[0].y;
		double ymax = p[0].y;

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
		ossimIpt delta_lr((ossim_int32)(ceil(theSlaveAccuracy)), (ossim_int32)(ceil(theSlaveAccuracy)));
		ossimDrect srect(ossimDpt(xmin, ymin) - delta_lr, ossimDpt(xmax, ymax) + delta_lr);

		//srect = srect.clipToRect(handlerS->getBoundingRect(0));
		srect = worldRectIntersection(srect, handlerS->getBoundingRect(0));

		ossimDpt sul = srect.ul();
		ossimDpt slr = srect.lr();

		if (srect.hasNans() || !slaveApp.getCombinedRect2CvMatByte(srect, slaveMat, { 2, 3, 1 }, { .2126, .7152, .0722 }, ossimDpt(1.0, 1.0), 0.015))
		//if (srect.hasNans() || !slaveApp.getRect2CvMatByte(srect, slaveMat, theSlaveBand, 1.0, 0.015))
		{
			continue;
		}
		if (countNonZero(slaveMat) < 1)
		{
			continue;
		}

		char url_buf[2048];
		//string type = "Google";
		string type = "Bing";
		struct MemoryStruct chunk;
		if (0 == stricmp(type.c_str(), "Google"))
		{
			sprintf_s(url_buf, "https://maps.googleapis.com/maps/api/staticmap?center=%lf,%lf&zoom=%d&size=%dx%d&maptype=satellite",
				ll_center.lat, ll_center.lon, nearestZoomLevel, nWidth, nWidth);
			chunk = getDataFromUrl(url_buf, "http://127.0.0.1:8087");
		}
		else if (0 == stricmp(type.c_str(), "Bing"))
		{
			sprintf_s(url_buf, "http://dev.virtualearth.net/REST/v1/Imagery/Map/Aerial/%lf,%lf/%d?mapSize=%d,%d&key=AopmRHqd_NETfVAUdDQL9KTpveOTs5zX8T2nVS23AOVjbLnYRCoCZLb14jEhWbgj",
				ll_center.lat, ll_center.lon, nearestZoomLevel, nWidth, nHeight);
			chunk = getDataFromUrl(url_buf);
		}

		cv::Mat imgbuf = cv::Mat(nWidth, nHeight, CV_8UC3, chunk.memory);
		//masterMat = cv::imdecode(imgbuf, CV_LOAD_IMAGE_COLOR);

		cv::Mat masterMat;
		found = false;
		masterMat = cv::imdecode(imgbuf, CV_LOAD_IMAGE_GRAYSCALE);
		if (countNonZero(masterMat) < 1)
		{
			continue;
		}

		double nNewWidth = masterMat.cols * scale;
		double nNewHeight = masterMat.rows * scale;

		cv::resize(masterMat, masterMat, cv::Size(nNewWidth, nNewHeight));//resize image

		if (chunk.memory)
			free(chunk.memory);


		//cv::Mat temp;
		//cv::GaussianBlur(slaveMat, temp, cv::Size(0, 0), 3);
		//cv::addWeighted(slaveMat, 1.5, temp, -0.5, 0, slaveMat);

		//cv::GaussianBlur(masterMat, temp, cv::Size(0, 0), 3);
		//cv::addWeighted(masterMat, 1.5, temp, -0.5, 0, masterMat);

		//masterMat = cv::imread(buf, 0);
		double cx = nNewWidth * 0.5;
		double cy = nNewHeight * 0.5;
		if (match_state::success == runMatchParallel(slaveMat, masterMat, tp[0], this, theDebug))
		{
			ossimDpt ctileXy = LonLat2TileXy(ll_center, nearestZoomLevel);
			ossimDpt tileXy = ctileXy + (tp[0].getMasterPoint() - ossimDpt(cx, cy))/scale;
			ossimGpt ll = TileXy2LonLat(tileXy.x, tileXy.y, nearestZoomLevel);
			//double mLon = ll_center.lon + (tp[0].getMasterPoint().x - cx) * res_deg.x;
			//double mLat = ll_center.lat + (tp[0].getMasterPoint().y - cy) * res_deg.y;
			//ossimGpt ll(mLat, mLon);

			tp[0].setSlavePoint(tp[0].getSlavePoint() + sul);
			// convert "Image to Image" to "Ground to Image" tie points    //TBC : use more generic tie points
			ossimRefPtr<radiBlockTieGpt> tgi(new radiBlockTieGpt);
			ll.hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(ll);
			tgi->setGroundPoint(ll);
			//set slave image position
			tgi->refImagePoint() = tp[0].getSlavePoint();
			tgi->setSlaveId(theSlaveId);
			tgi->setMasterId(theMasterId);
			tgi->setPointType(radiBlockTieGpt::known_ground_control_points);
			//set score
			tgi->setScore(tp[0].score);
			//add to list
			tSet.addTiePoint(tgi);

			found = true;
			break;
		}

		//// search in the reference library
		//vector<ossimFilename> masterFileList;
		//masterFileList.clear();
		//masterFileList.push_back(theMaster);

		//cv::Mat masterMat;
		//found = false;
		//for (int iFile = 0;iFile < (int)masterFileList.size();++iFile)
		//{
		//	ossimFilename lastMaster = masterFileList[iFile];
		//	theLastMaster = lastMaster;


		//	//string url = getRequestUrl(theSlaveProjection, sGSD,
		//	//	srect, theSlaveAccuracy);
		//	ossimDrect mllRect = getLonLatRect(theSlaveProjection, srect, theSlaveAccuracy);

		//	double clon = mllRect.midPoint().x;
		//	double clat = mllRect.midPoint().y;

		//	int nZoomLevel = getZoomLevel(sGSD, clat);
		//	ossimDpt res_deg = gsd_degree(nZoomLevel);
		//	double aa = mllRect.lr().x - mllRect.ul().x;
		//	double bb = mllRect.ul().y - mllRect.lr().y;
		//	int nWidth = aa / res_deg.x;
		//	int nHeight = bb / res_deg.y;

		//	char url_buf[2048];
		//	//sprintf_s(url_buf, "http://maps.googleapis.com/maps/api/staticmap?center=%lf,%lf&zoom=%d&size=%dx%d&maptype=satellite",
		//	//	clat, clon, nZoomLevel, nWidth, nHeight);
		//	sprintf_s(url_buf, "http://dev.virtualearth.net/REST/v1/Imagery/Map/Aerial?mapArea=%lf,%lf,%lf,%lf&zoomLevel=%d&key=AopmRHqd_NETfVAUdDQL9KTpveOTs5zX8T2nVS23AOVjbLnYRCoCZLb14jEhWbgj",
		//		mllRect.lr().y, mllRect.ul().x, mllRect.ul().y, mllRect.lr().x, nZoomLevel);

		//

		//	struct MemoryStruct chunk = getDataFromUrl(url_buf);
		//	cv::Mat imgbuf = cv::Mat(nWidth, nHeight, CV_8UC3, chunk.memory);
		//	//masterMat = cv::imdecode(imgbuf, CV_LOAD_IMAGE_COLOR);
		//	masterMat = cv::imdecode(imgbuf, CV_LOAD_IMAGE_GRAYSCALE);

		//	if (chunk.memory)
		//		free(chunk.memory);

		//	//masterMat = cv::imread(buf, 0);
		//	double cx = nWidth * 0.5;
		//	double cy = nHeight * 0.5;
		//	if (match_state::success == runMatchParallel(slaveMat, masterMat, tp[0], this, theDebug))
		//	{
		//		double mLon = clon + (tp[0].getMasterPoint().x - cx) * res_deg.x;
		//		double mLat = clat + (tp[0].getMasterPoint().y - cy) * res_deg.y;
		//		ossimGpt ll(mLon, mLat);

		//		tp[0].setSlavePoint(tp[0].getSlavePoint() + sul);
		//		// convert "Image to Image" to "Ground to Image" tie points    //TBC : use more generic tie points
		//		ossimRefPtr<radiBlockTieGpt> tgi(new radiBlockTieGpt);
		//		ll.hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(ll);
		//		tgi->setGroundPoint(ll);
		//		//set slave image position
		//		tgi->refImagePoint() = tp[0].getSlavePoint();
		//		tgi->setSlaveId(theSlaveId);
		//		tgi->setMasterId(theMasterId);
		//		tgi->setPointType(radiBlockTieGpt::known_ground_control_points);
		//		//set score
		//		tgi->setScore(tp[0].score);
		//		//add to list
		//		tSet.addTiePoint(tgi);

		//		found = true;
		//		break;
		//	}
		//}
	}
	slaveApp.close();
	return found;
}

//bool radiImageRegistration::getGridFeaturesParallel(const ossimIrect& rect, radiBlockTieGpt& tSet)
//{
//	if (!theSlaveBandSelector)
//	{
//		ossimNotify(ossimNotifyLevel_WARN)
//			<< "WARN ossimTieGenerator::scanForEdges():"
//			<< "\nInput source is not a ossimImageChip.  Returning..." << std::endl;
//		return false;
//	}
//
//	const ossim_int32 TILE_HEIGHT	= theTileSize;
//	const ossim_int32 TILE_WIDTH	= theTileSize;
//	const ossim_int32 START_LINE = rect.ul().y;
//	const ossim_int32 STOP_LINE  = rect.lr().y;
//	const ossim_int32 START_SAMP = rect.ul().x;
//	const ossim_int32 STOP_SAMP  = rect.lr().x;
//	int nWidth = STOP_SAMP-START_SAMP;
//	int nHeight = STOP_LINE-START_LINE;
//
//	// For percent complete status.
//	ossim_int32 tilerows=(STOP_LINE-START_LINE+TILE_HEIGHT) / TILE_HEIGHT; //ceil : (stop-start+1+size-1)/size
//	ossim_int32 tilecols=(STOP_SAMP-START_SAMP+TILE_WIDTH) / TILE_WIDTH;
//	double total_tiles = ((double)tilerows)*tilecols;
//	double tiles_processed = 0.0;
//
//	// loop through all tiles
//	// need to use a sequencer for parallelism in the future TBD
//	ossim_int32 line=START_LINE;
//	ossim_int32 i,j;
//
//	vector<row_col> row_col_List;
//	double center_row = tilerows * 0.5;
//	double center_col = tilecols * 0.5;
//	for (i=0;(i<tilerows);++i)
//	{
//		for (j=0;(j<tilecols);++j )
//		{
//			row_col_List.push_back(row_col((center_row-i-1),(center_col-j-1)));
//		}
//	}
//
//	GdalRasterApp slaveApp;
//	if (!slaveApp.open(getSlave().c_str()))
//	{
//		cerr<<"radiImageRegistration"<<"::execute can't open slave image  "<< getSlave().c_str() <<endl;
//		return false;
//	}
//	// select only one band (if multiple)
//	ossim_uint32 sbc = slaveApp.nBand();
//
//	//add a band selector
//	ossim_uint32 sb = theSlaveBand;
//	if (sb>=sbc) 
//	{
//		cerr<<"radiImageRegistration"<<"::execute Warning not enough bands in slave, only "<< sbc <<endl;
//		sb=0;
//	}
//
//	double sGSD = min(theSlaveProjection->getMetersPerPixel().x, theSlaveProjection->getMetersPerPixel().y);
//	//int ncore = omp_get_num_procs();//获取执行核的总数；  目前机器CPU的数量
//	row_col_step = 5;
//	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);	
//
//	//GdalRasterApp slaveApp;
//	//slaveApp.open(theSlave.c_str());
//	bool bDebug = true;
//	//if (theTset.size() == 20)
//	//{
//	//	bDebug = true;
//	//}
//	bool found = false;
//	int N_PARALLEL = 1;//ncore*2;
//	int search_size = (int)row_col_List.size();
//	//search_size = min(search_size, 1);
//	for (int i=0;i < search_size && !found;)
//		//for (int i=0;i < (int)row_col_List.size() && !found;)
//	{
//		int nParallel = min(N_PARALLEL, (int)row_col_List.size()-i);	// 保证末尾不越界
//
//		std::vector<cv::Mat> slaveMatList;
//		std::vector<ossimIrect> srectList;
//		std::vector<bool>  slaveValidList(nParallel, false);
//		for (int j = 0; j < nParallel;j++)
//		{
//			int icol = floor(row_col_List[i+j].col_idx+center_col+0.5);
//			int irow = floor(row_col_List[i+j].row_idx+center_row+0.5);
//			//int ii = (5 * i) % (int)row_col_List.size();
//			//int icol = floor(row_col_List[ii+j].col_idx+center_col+0.5);
//			//int irow = floor(row_col_List[ii+j].row_idx+center_row+0.5);
//			ossim_int32 samp=START_SAMP+icol*TILE_WIDTH;
//			ossim_int32 line=START_LINE+irow*TILE_HEIGHT;
//
//
//			ossim_int32 BufHeight = TILE_HEIGHT;
//			ossim_int32 BufWidth = TILE_WIDTH;
//			//行末尾小块处理
//			if (irow == tilerows-1)
//			{
//				BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;//得出当前块的宽度Bufsizex，高度Bufsizey
//				BufHeight = min(BufHeight, TILE_HEIGHT);
//			}
//			//列末尾小块处理
//			if (icol == tilecols-1)
//			{
//				BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;//得出当前块的宽度Bufsizex，高度Bufsizey
//				BufWidth = min(BufWidth, TILE_WIDTH);
//			}
//
//			// slave
//			cv::Mat slaveMat;
//			ossimIrect srect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1));
//			//if (slaveApp.getPrincipalRect2CvMatByte( srect, slaveMat))
//			if (slaveApp.getRect2CvMatByte( srect, slaveMat, theSlaveBand))
//			{
//				slaveValidList[j] = true;
//			}
//			{
//				//if (createTileMat(caster[0], srect, slaveMat, 0))
//				//{
//				//	slaveValidList[j] = true;
//				//}
//			}
//			slaveMatList.push_back(slaveMat);
//			srectList.push_back(srect);
//		}
//		for (int j = 0; j < nParallel;j++)
//		{
//			if (found || !slaveValidList[j])
//			{
//				continue;
//			}
//			ossimIpt sul = srectList[j].ul();
//			ossimIpt slr = srectList[j].lr();
//			ossimIpt delta_lr((ossim_int32)(ceil(theSlaveAccuracy)), (ossim_int32)(ceil(theSlaveAccuracy)) );
//
//			vector<ossimTDpt> tp(1);
//			// search in the reference library
//			ossimGpt ul_latlon, lr_latlon;
//			theSlaveProjection->lineSampleToWorld(sul-delta_lr, ul_latlon);
//			theSlaveProjection->lineSampleToWorld(slr+delta_lr, lr_latlon);
//			//ossimDpt tempDpt;
//			//slaveApp.linesample2lonlat(sul-delta_lr, tempDpt);
//			//ul_latlon = ossimGpt(tempDpt.y, tempDpt.x);
//			//slaveApp.linesample2lonlat(slr+delta_lr, tempDpt);
//			//lr_latlon = ossimGpt(tempDpt.y, tempDpt.x);
//			vector<ossimFilename> masterFileList;
//			if (theMaster.ext().upcase() == "SHP")
//			{
//				getMasterList(theMaster, masterFileList, ul_latlon, lr_latlon);
//			}
//			else{
//				masterFileList.clear();
//				masterFileList.push_back(theMaster);
//			}
//
//			//ossimRefPtr<ossimImageHandler> master_handler = NULL;
//			//ossimProjection* master_projection = NULL;
//			//ossimRefPtr<ossimBandSelector> master_bandselector = NULL;
//			cv::Mat masterMat;
//			//ossimRefPtr<ossimCastTileSourceFilter> master_caster = new ossimCastTileSourceFilter();
//			//master_caster->setOutputScalarType(OSSIM_FLOAT64);
//			//master_caster->setOutputScalarType(OSSIM_UCHAR);
//
//			int match_result = match_state::success;
//			for (int iFile = 0;iFile < (int)masterFileList.size();++iFile)
//			{
//				if (found || match_result == match_state::slave_faild)
//				{
//					continue;
//				}
//
//				ossimFilename lastMaster = masterFileList[iFile];
//				theLastMaster = lastMaster;
//				GdalRasterApp masterApp;
//				if (!masterApp.open(lastMaster.c_str()))
//				{
//					cerr<<"radiImageRegistration"<<"::execute can't open master image  "<< lastMaster <<endl;
//					continue;
//				}
//				// select only one band (if multiple)
//				ossim_uint32 mbc = masterApp.nBand();
//
//				//add a band selector
//				ossim_uint32 mb = theMasterBand;
//				if (mb>=mbc) 
//				{
//					//cerr<<"radiImageRegistration"<<"::execute Warning not enough bands in master, only "<< mbc <<endl;
//					mb=0;
//				}
//
//				double gsd_scale;
//				// master
//				ossimIrect mrect;
//				if (0 == strcmp(masterApp.getGetProjectionRef(), ""))
//				{
//					// 无投影
//					mrect = ossimIrect(sul - delta_lr, slr + delta_lr);
//					gsd_scale = 1.0;
//				}
//				else
//				{
//					double mGSD = masterApp.getGeoTransform()[1];
//					gsd_scale = mGSD / sGSD;
//					//gsd_scale = 1.0;
//					// master
//					mrect = getMasterRect(theSlaveProjection, masterApp,
//						srectList[j], theSlaveAccuracy);
//				}
//				//ossimIpt mul, mlr;
//				//ossimGpt tempGpt;
//				//theSlaveProjection->lineSampleToWorld(sul - delta_lr, tempGpt);
//				////theSlaveProjection->lineSampleToWorld(sul, tempGpt);
//				//masterApp.lonlat2linesample(ossimDpt(tempGpt.lon, tempGpt.lat), mul);
//				//theSlaveProjection->lineSampleToWorld(slr + delta_lr, tempGpt);
//				////theSlaveProjection->lineSampleToWorld(slr, tempGpt);
//				//masterApp.lonlat2linesample(ossimDpt(tempGpt.lon, tempGpt.lat), mlr);
//
//				//delta_lr = delta_lr / gsd_scale;
//				//mul = mul - delta_lr;
//				//mlr = mlr + delta_lr;
//
//				//ossimDpt tempDpt;
//				//slaveApp.linesample2lonlat(sul - delta_lr, tempDpt);
//				//masterApp.lonlat2linesample(tempDpt, mul);
//				//slaveApp.linesample2lonlat(slr + delta_lr, tempDpt);
//				//masterApp.lonlat2linesample(tempDpt, mlr);
//				//ossimIrect mrect = ossimIrect(mul, mlr);
//				//if(!masterApp.getBoundary().pointWithin(mul) 
//				//	|| !masterApp.getBoundary().pointWithin(mlr))
//				//{
//				//	masterApp.close();
//				//	continue;
//				//}
//				//if (!masterApp.getPrincipalRect2CvMatByte( mrect, masterMat, gsd_scale))
//				if (!masterApp.getRect2CvMatByte( mrect, masterMat, mb, gsd_scale))
//				{
//					masterApp.close();
//					continue;
//				}
//				match_result = runMatchParallel(slaveMatList[j], masterMat, tp[0], this, bDebug);
//				if (match_result == match_state::slave_faild)
//				{
//					masterApp.close();
//					continue;
//				}
//				if (theFilename != ossimFilename::NIL && match_state::success == match_result)
//				{
//					//ossimTieGpt tiePt;
//					//tiePt.setImagePoint(ossimDpt(tp[0].getSlavePoint() + sul));
//					//ossimGpt gpt;
//					//masterApp.linesample2lonlat(tp[0].getMasterPoint() + mul, gpt);
//					//tiePt.setGroundPoint(gpt);
//					//theTset.addTiePoint(&tiePt);
//
//					if (0 == strcmp(masterApp.getGetProjectionRef(), "")
//						|| thePointType == point_type::tie)
//					{
//						tp[0].setMasterPoint(tp[0].getMasterPoint() + mrect.ul());
//						tp[0].setSlavePoint(tp[0].getSlavePoint() + sul);
//						// 无投影
//						ossimRefPtr<ossimTieGpt> tgi(new ossimTieGpt);
//						ossimDpt dpt = tp[0].getMasterPoint();
//						tgi->setGroundPoint(ossimGpt(dpt.y, dpt.x, 0.0));
//						//set slave image position
//						tgi->refImagePoint() = tp[0].getSlavePoint();
//						//set score
//						tgi->setScore(tp[0].score);
//
//						//add to list
//						tSet.addTiePoint(tgi);
//					}
//					else
//					{
//						tp[0].setMasterPoint(tp[0].getMasterPoint()/gsd_scale + mrect.ul());
//						tp[0].setSlavePoint(tp[0].getSlavePoint() + sul);
//						// convert "Image to Image" to "Ground to Image" tie points    //TBC : use more generic tie points
//						ossimRefPtr<ossimTieGpt> tgi(new ossimTieGpt);
//						//set master ground pos
//						ossimDpt lonlat;
//						masterApp.linesample2lonlat(tp[0].getMasterPoint(), lonlat);
//						double hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(ossimGpt(lonlat.y, lonlat.x));
//						//ossimDpt dpt = outProj->forward(*tgi);
//						tgi->setGroundPoint(ossimGpt(lonlat.y, lonlat.x, hgt));
//						//set slave image position
//						tgi->refImagePoint() = tp[0].getSlavePoint();
//						//set score
//						tgi->setScore(tp[0].score);
//						//add to list
//						tSet.addTiePoint(tgi);
//					}
//				}
//				if (match_state::success == match_result)
//				{
//					masterApp.close();
//					found = true;
//					continue;
//				}
//				masterApp.close();
//			}
//		}
//		slaveMatList.clear();
//		slaveValidList.clear();
//		srectList.clear();
//		i += N_PARALLEL;
//	}
//	slaveApp.close();
//	return found;
//}

//bool radiImageRegistration::getGridFeaturesParallel(const ossimIrect& rect)
//{
//	if (!theSlaveBandSelector)
//	{
//		ossimNotify(ossimNotifyLevel_WARN)
//			<< "WARN ossimTieGenerator::scanForEdges():"
//			<< "\nInput source is not a ossimImageChip.  Returning..." << std::endl;
//		return false;
//	}
//
//	const ossim_int32 TILE_HEIGHT	= theTileSize;
//	const ossim_int32 TILE_WIDTH	= theTileSize;
//	const ossim_int32 START_LINE = rect.ul().y;
//	const ossim_int32 STOP_LINE  = rect.lr().y;
//	const ossim_int32 START_SAMP = rect.ul().x;
//	const ossim_int32 STOP_SAMP  = rect.lr().x;
//	int nWidth = STOP_SAMP-START_SAMP;
//	int nHeight = STOP_LINE-START_LINE;
//
//	// For percent complete status.
//	ossim_int32 tilerows=(STOP_LINE-START_LINE+TILE_HEIGHT) / TILE_HEIGHT; //ceil : (stop-start+1+size-1)/size
//	ossim_int32 tilecols=(STOP_SAMP-START_SAMP+TILE_WIDTH) / TILE_WIDTH;
//	double total_tiles = ((double)tilerows)*tilecols;
//	double tiles_processed = 0.0;
//
//	// loop through all tiles
//	// need to use a sequencer for parallelism in the future TBD
//	ossim_int32 line=START_LINE;
//	ossim_int32 i,j;
//
//	vector<row_col> row_col_List;
//	double center_row = tilerows * 0.5;
//	double center_col = tilecols * 0.5;
//	for (i=0;(i<tilerows);++i)
//	{
//		for (j=0;(j<tilecols);++j )
//		{
//			row_col_List.push_back(row_col((center_row-i-1),(center_col-j-1)));
//		}
//	}
//
//	int ncore = omp_get_num_procs();//获取执行核的总数；  目前机器CPU的数量
//	row_col_step = 5;
//	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);	
//
//	//GdalRasterApp slaveApp;
//	//slaveApp.open(theSlave.c_str());
//	bool bDebug = true;
//	//if (theTset.size() == 20)
//	//{
//	//	bDebug = true;
//	//}
//	bool found = false;
//	int N_PARALLEL = 2;//ncore*2;
//	for (int i=0;i < (int)row_col_List.size() && !found;)
//	{
//		int nParallel = min(N_PARALLEL, (int)row_col_List.size()-1-i);	// 保证末尾不越界
//
//		std::vector<cv::Mat> slaveMatList;
//		std::vector<ossimIrect> srectList;
//		std::vector<bool>  slaveValidList(nParallel, false);
//		for (int j = 0; j < nParallel;j++)
//		{
//			int icol = floor(row_col_List[i+j].col_idx+center_col+0.5);
//			int irow = floor(row_col_List[i+j].row_idx+center_row+0.5);
//			ossim_int32 samp=START_SAMP+icol*TILE_WIDTH;
//			ossim_int32 line=START_LINE+irow*TILE_HEIGHT;
//
//
//			ossim_int32 BufHeight = TILE_HEIGHT;
//			ossim_int32 BufWidth = TILE_WIDTH;
//			//行末尾小块处理
//			if (irow == tilerows-1)
//			{
//				BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;//得出当前块的宽度Bufsizex，高度Bufsizey
//				BufHeight = min(BufHeight, TILE_HEIGHT);
//			}
//			//列末尾小块处理
//			if (icol == tilecols-1)
//			{
//				BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;//得出当前块的宽度Bufsizex，高度Bufsizey
//				BufWidth = min(BufWidth, TILE_WIDTH);
//			}
//
//			// slave
//			cv::Mat slaveMat;
//			ossimIrect srect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1));
//			//if (slaveApp.getRect2CvMat( srect, slaveMat, theSlaveBand))
//			//{
//			//	slaveValidList[j] = true;
//			//}
//			if (createTileMat(caster[0], srect, slaveMat, 0))
//			{
//				slaveValidList[j] = true;
//			}
//			slaveMatList.push_back(slaveMat);
//			srectList.push_back(srect);
//		}
//
//		ossimFilename lastMaster = theMaster;
//		theLastMaster = lastMaster;
//		GdalRasterApp masterApp;
//		if (!masterApp.open(lastMaster.c_str()))
//		{
//			cerr<<"radiImageRegistration"<<"::execute can't open master image  "<< lastMaster <<endl;
//			continue;
//		}
//		// select only one band (if multiple)
//		ossim_uint32 mbc = masterApp.nBand();
//
//		//add a band selector
//		ossim_uint32 mb = theMasterBand;
//		if (mb>=mbc) 
//		{
//			cerr<<"radiImageRegistration"<<"::execute Warning not enough bands in master, only "<< mbc <<endl;
//			mb=0;
//		}
//
//		int match_result = match_state::match_failed;
//#pragma omp parallel for num_threads(nParallel)
//		for (int j = 0; j < nParallel;j++)
//		{
//			if (found || !slaveValidList[j])
//			{
//				continue;
//			}
//			ossimIpt sul = srectList[j].ul();
//			ossimIpt slr = srectList[j].lr();
//			ossimIpt delta_lr((ossim_int32)(ceil(theSlaveAccuracy)), (ossim_int32)(ceil(theSlaveAccuracy)) );
//
//			vector<ossimTDpt> tp(1);
//			// search in the reference library
//			ossimGpt ul_latlon, lr_latlon;
//			theSlaveProjection->lineSampleToWorld(sul-delta_lr, ul_latlon);
//			theSlaveProjection->lineSampleToWorld(slr+delta_lr, lr_latlon);
//			
//
//			cv::Mat masterMat;
//			
//
//			// master
//			ossimIpt mul, mlr;
//			ossimGpt tempGpt;
//			theSlaveProjection->lineSampleToWorld(sul - delta_lr, tempGpt);
//			masterApp.lonlat2linesample(ossimDpt(tempGpt.lon, tempGpt.lat), mul);
//			theSlaveProjection->lineSampleToWorld(slr + delta_lr, tempGpt);
//			masterApp.lonlat2linesample(ossimDpt(tempGpt.lon, tempGpt.lat), mlr);
//			//ossimDpt tempDpt;
//			//slaveApp.linesample2lonlat(sul - delta_lr, tempDpt);
//			//masterApp.lonlat2linesample(tempDpt, mul);
//			//slaveApp.linesample2lonlat(slr + delta_lr, tempDpt);
//			//masterApp.lonlat2linesample(tempDpt, mlr);
//			ossimIrect mrect = ossimIrect(mul, mlr);
//			if(!masterApp.getBoundary().pointWithin(mul) 
//				|| !masterApp.getBoundary().pointWithin(mlr))
//			{
//				continue;
//			}
//			if (!masterApp.getRect2CvMat( mrect, masterMat, mb))
//			{
//				continue;
//			}
//			match_result = runMatchParallel(slaveMatList[j], masterMat, tp[0], this, bDebug);
//			if (match_result == match_state::slave_faild)
//			{
//				continue;
//			}
//			if (thedTieptFilename != ossimFilename::NIL && match_state::success == match_result)
//			{
//				tp[0].setMasterPoint(tp[0].getMasterPoint() + mul);
//				tp[0].setSlavePoint(tp[0].getSlavePoint() + sul);
//				//write on stream
//				writeTiePoints(tp);
//				if (getStoreFlag())
//				{
//					theTiePoints.push_back(tp[0]);
//					// convert "Image to Image" to "Ground to Image" tie points    //TBC : use more generic tie points
//
//					ossimRefPtr<ossimTieGpt> tgi(new ossimTieGpt);
//					//set master ground pos
//					ossimDpt lonlat;
//					masterApp.linesample2lonlat(tp[0].getMasterPoint(), lonlat);
//					double hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(ossimGpt(lonlat.y, lonlat.x));
//					//ossimDpt dpt = outProj->forward(*tgi);
//					tgi->setGroundPoint(ossimGpt(lonlat.y, lonlat.x, hgt));
//					//set slave image position
//					tgi->refImagePoint() = tp[0].getSlavePoint();
//					//set score
//					tgi->setScore(tp[0].score);
//
//					//add to list
//					theTset.addTiePoint(tgi);
//				}
//			}
//			if (match_state::success == match_result)
//			{
//				found = true;
//				continue;
//			}
//			
//		}
//		slaveMatList.clear();
//		slaveValidList.clear();
//		srectList.clear();
//		i += N_PARALLEL;
//	}
//	//slaveApp.close();
//	return found;
//}

//bool radiImageRegistration::getGridFeaturesParallel(const ossimIrect& rect, void *pData)
//{
//	radiImageRegistration* pThis = (radiImageRegistration*)pData;
//	if (!pThis->theSlaveBandSelector)
//	{
//		ossimNotify(ossimNotifyLevel_WARN)
//			<< "WARN ossimTieGenerator::scanForEdges():"
//			<< "\nInput source is not a ossimImageChip.  Returning..." << std::endl;
//		return false;
//	}
//
//	const ossim_int32 TILE_HEIGHT	= pThis->theTileSize;
//	const ossim_int32 TILE_WIDTH	= pThis->theTileSize;
//	const ossim_int32 START_LINE = rect.ul().y;
//	const ossim_int32 STOP_LINE  = rect.lr().y;
//	const ossim_int32 START_SAMP = rect.ul().x;
//	const ossim_int32 STOP_SAMP  = rect.lr().x;
//	int nWidth = STOP_SAMP-START_SAMP;
//	int nHeight = STOP_LINE-START_LINE;
//
//	// For percent complete status.
//	ossim_int32 tilerows=(STOP_LINE-START_LINE+TILE_HEIGHT) / TILE_HEIGHT; //ceil : (stop-start+1+size-1)/size
//	ossim_int32 tilecols=(STOP_SAMP-START_SAMP+TILE_WIDTH) / TILE_WIDTH;
//	double total_tiles = ((double)tilerows)*tilecols;
//	double tiles_processed = 0.0;
//
//	// loop through all tiles
//	// need to use a sequencer for parallelism in the future TBD
//	ossim_int32 line=START_LINE;
//	ossim_int32 i,j;
//
//	vector<row_col> row_col_List;
//	double center_row = tilerows * 0.5;
//	double center_col = tilecols * 0.5;
//	for (i=0;(i<tilerows);++i)
//	{
//		for (j=0;(j<tilecols);++j )
//		{
//			row_col_List.push_back(row_col((center_row-i-1),(center_col-j-1)));
//		}
//	}
//
//	row_col_step = 5;
//	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);	
//
//	bool bDebug = false;
//	//if (theTset.size() == 20)
//	//{
//	//	bDebug = true;
//	//}
//	bool found = false;
//	int ncore = omp_get_num_procs();//获取执行核的总数；  目前机器CPU的数量
//	//#pragma omp parallel for num_threads(ncore*2) shared(found)
////#pragma omp single
////{
//	for (int i=0;i < (int)row_col_List.size();i+=1)
//	{
//		if (found)
//		{
//			continue;
//		}
//		int icol = floor(row_col_List[i].col_idx+center_col+0.5);
//		int irow = floor(row_col_List[i].row_idx+center_row+0.5);
//		ossim_int32 samp=START_SAMP+icol*TILE_WIDTH;
//		ossim_int32 line=START_LINE+irow*TILE_HEIGHT;
//
//
//		ossim_int32 BufHeight = TILE_HEIGHT;
//		ossim_int32 BufWidth = TILE_WIDTH;
//		//行末尾小块处理
//		if (irow == tilerows-1)
//		{
//			BufHeight = nHeight - (tilerows-1) * TILE_HEIGHT;//得出当前块的宽度Bufsizex，高度Bufsizey
//			BufHeight = min(BufHeight, TILE_HEIGHT);
//		}
//		//列末尾小块处理
//		if (icol == tilecols-1)
//		{
//			BufWidth = nWidth - (tilecols-1) * TILE_WIDTH;//得出当前块的宽度Bufsizex，高度Bufsizey
//			BufWidth = min(BufWidth, TILE_WIDTH);
//		}
//		
//		// slave
//		cv::Mat slaveMat;
//		ossimIrect srect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1));
//		ossimIpt sul = srect.ul();
//		ossimIpt slr = srect.lr();
//		ossimIpt delta_lr((ossim_int32)(ceil(pThis->theSlaveAccuracy)), (ossim_int32)(ceil(pThis->theSlaveAccuracy)) );
//		if (!pThis->createTileMat(pThis->caster[0], srect, slaveMat, 0))
//		{
//			continue;
//		}
//
//		vector<ossimTDpt> tp(1);
//		// search in the reference library
//		ossimGpt ul_latlon, lr_latlon;
//		pThis->theSlaveProjection->lineSampleToWorld(sul-delta_lr, ul_latlon);
//		pThis->theSlaveProjection->lineSampleToWorld(slr+delta_lr, lr_latlon);
//		vector<ossimFilename> masterFileList;
//		if (pThis->theMaster.ext().upcase() == "SHP")
//		{
//			pThis->getMasterList(pThis->theMaster, masterFileList, ul_latlon, lr_latlon);
//		}
//		else{
//			masterFileList.clear();
//			masterFileList.push_back(pThis->theMaster);
//		}
//
//		ossimRefPtr<ossimImageHandler> master_handler = NULL;
//		ossimProjection* master_projection = NULL;
//		ossimRefPtr<ossimBandSelector> master_bandselector = NULL;
//		cv::Mat masterMat;
//		ossimRefPtr<ossimCastTileSourceFilter> master_caster = new ossimCastTileSourceFilter();
//		//master_caster->setOutputScalarType(OSSIM_FLOAT64);
//		master_caster->setOutputScalarType(OSSIM_UCHAR);
//		int match_result = match_state::success;
//		int nThreads = min(ncore*2, (int)masterFileList.size());
//		//int match_result = match_state::success;
////#pragma omp critical
////#pragma omp parallel for num_threads(nThreads) firstprivate(slaveMat) \
//		shared(found, match_result)\
//		private(master_handler, master_projection, master_bandselector, masterMat, master_caster)
//#pragma omp parallel for num_threads(nThreads)
//		for (int iFile = 0;iFile < (int)masterFileList.size();++iFile)
//		{
//			if (found || match_result == match_state::slave_faild)
//			{
//				continue;
//			}
//			//if (!master_caster)
//			//{
//			//	master_caster = new ossimCastTileSourceFilter();
//			//	master_caster->setOutputScalarType(OSSIM_FLOAT64);
//			//}
//			ossimFilename theLastMaster = masterFileList[iFile];
//			master_handler = ossimImageHandlerRegistry::instance()->open(theLastMaster);
//			pThis->theLastMaster = theLastMaster;
//			if (!master_handler)
//			{
//				cerr<<"radiImageRegistration"<<"::execute can't create handler for master image  "<< theLastMaster <<endl;
//				continue;
//			}
//			// select only one band (if multiple)
//			ossim_uint32 mbc = master_handler->getNumberOfOutputBands();
//
//			//add a band selector
//			ossim_uint32 mb = pThis->theMasterBand;
//			if (mb>=mbc) 
//			{
//				cerr<<"radiImageRegistration"<<"::execute Warning not enough bands in master, only "<< mbc <<endl;
//				mb=0;
//			}
//			//cout<<"Using band "<<mb<<" for master"<<endl; //TBR
//			master_bandselector = new ossimBandSelector;
//			master_bandselector->connectMyInputTo(0, master_handler.get());
//			master_bandselector->setOutputBandList(vector<ossim_uint32>(1,mb));
//			master_caster->connectMyInputTo(0, master_bandselector.get());
//			master_projection = master_handler->getImageGeometry()->getProjection();
//
//			// master
//			ossimIpt mul = pThis->slave2master(pThis->theSlaveProjection, master_projection, sul - delta_lr);
//			ossimIpt mlr = pThis->slave2master(pThis->theSlaveProjection, master_projection, slr + delta_lr);
//			ossimIrect mrect = ossimIrect(mul, mlr);
//			if(!master_handler->getBoundingRect(0).pointWithin(mul) 
//				|| !master_handler->getBoundingRect(0).pointWithin(mlr))
//			{
//				continue;
//			}
//			if (!createTileMat(master_caster, mrect, masterMat, 0))
//			{
//				//master_caster->disconnect();
//				//if (master_handler != NULL)
//				//{
//				//	master_handler->disconnect();
//				//	master_handler = NULL;
//				//}
//				continue;
//			}
//			match_result = pThis->runMatchParallel(slaveMat, masterMat, tp[0], pData, bDebug);
//			if (match_result == match_state::slave_faild)
//			{
//				//master_caster->disconnect();
//				//if (master_handler != NULL)
//				//{
//				//	master_handler->disconnect();
//				//	master_handler = NULL;
//				//}
//				continue;
//			}
//			if (pThis->thedTieptFilename != ossimFilename::NIL && match_state::success == match_result)
//			{
//				tp[0].setMasterPoint(tp[0].getMasterPoint() + mul);
//				tp[0].setSlavePoint(tp[0].getSlavePoint() + sul);
//				//write on stream
//				pThis->writeTiePoints(tp);
//				if (pThis->getStoreFlag())
//				{
//					pThis->theTiePoints.push_back(tp[0]);
//					// convert "Image to Image" to "Ground to Image" tie points    //TBC : use more generic tie points
//
//					ossimRefPtr<ossimTieGpt> tgi(new ossimTieGpt);
//					//set master ground pos
//					master_projection->lineSampleToWorld( tp[0].getMasterPoint() , *tgi ); //TBC : is it always lon/lat WGS84?
//					double hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(*tgi);
//					//ossimDpt dpt = outProj->forward(*tgi);
//					tgi->setGroundPoint(ossimGpt(tgi->lat, tgi->lon, hgt));
//					//set slave image position
//					tgi->refImagePoint() = tp[0].getSlavePoint();
//					//set score
//					tgi->setScore(tp[0].score);
//
//					//add to list
//					pThis->theTset.addTiePoint(tgi);
//				}
//			}
//			if (match_state::success == match_result)
//			{
//				master_caster->disconnect();
//				if (master_handler != NULL)
//				{
//					master_handler->disconnect();
//					master_handler = NULL;
//				}
//
//				found = true;
//				//return true;
//				continue;
//				//exit(0);
//				//found = true;
//				//break;
//				//return true;
//			}
//
//			//master_caster->disconnect();
//			//if (master_handler != NULL)
//			//{
//			//	master_handler->disconnect();
//			//	master_handler = NULL;
//			//}
//		}
//
//		master_caster->disconnect();
//		if (master_handler != NULL)
//		{
//			master_handler->disconnect();
//			master_handler = NULL;
//		}
//	}
////}
//	//setPercentComplete(100.0);
//	return found;
//}

bool radiImageRegistration::getAllFeatures()
{
	int  nTasks, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
	if (0 == rank)
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
	ossimRefPtr<ossimCastTileSourceFilter> slaveCaster = caster[0];
	totalBlocks = (int)row_col_List.size();
	if (theThreadNum == 0)
	{
		theThreadNum = OpenThreads::GetNumberOfProcessors() * 2;
	}
	//if(num_threads > totalBlocks) num_threads = totalBlocks;
	//num_threads = 1;

	int numForTask = 0;
	int startNum = 0;
	if (rank == 0)
	{
		// assign tasks
		int portion = 0;
		portion = totalBlocks / nTasks;
		startNum = 0;
		numForTask = portion - startNum;

		for (int i = 1; i < nTasks; i++)
		{
			// calculate the data for each thread.
			int curStartNum = i * portion;
			int curEndNum = (i + 1) * portion;
			if (i == nTasks - 1) { curEndNum = totalBlocks - 1; }
			if (curStartNum < 0) { curStartNum = 0; }

			// we need to send a thread the number of characters it will be receiving.
			int curLength = curEndNum - curStartNum;
			MPI_Send(&curLength, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&curStartNum, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
		}
	}
	else
	{
		// We are not the thread that read the file.
		// We need to receive data from whichever thread 
		MPI_Status status;
		MPI_Recv(&numForTask, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&startNum, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
		// Thread 0 is responsible for making sure the startNum and endNum calculated here are valid.
		// This is because thread 0 tells us exactly how many characters we were send.
	}

	theThreadNum = min(theThreadNum, numForTask);
	if (rank == 0)
	{
		cout << "using " << theThreadNum << " threads..." << endl;
	}
	//cout << "using " << theThreadNum << " threads..." << endl;
	GLOBAL_NUM_THREADS = theThreadNum + 1;
	std::vector<radiMatchRectThread *> threads(theThreadNum);
	OpenThreads::Thread::SetConcurrency(theThreadNum);
	OpenThreads::Thread::Init();


	for (int i = 0; i < theThreadNum; ++i) {
		threads[i] = new radiMatchRectThread(this);
		vector<ossimIrect> rectList;
		for (int j = i + startNum; j < i + startNum + numForTask; j += theThreadNum)
		{
			int irow = (int)row_col_List[j].row_idx;

			ossim_int32 line = START_LINE + irow*TILE_HEIGHT;
			ossim_int32 BufHeight = TILE_HEIGHT;
			ossim_int32 BufWidth = TILE_WIDTH;

			//ossim_int32 BufHeight = min(TILE_HEIGHT, nHeight);
			//ossim_int32 BufWidth = min(TILE_WIDTH, nWidth);
			//列末尾小块处理
			if (irow == tilerows - 1)
			{
				BufHeight = nHeight - (tilerows - 1) * TILE_HEIGHT;
				BufHeight = min(BufHeight, TILE_HEIGHT);
			}
			//cout<<BufWidth<<"\t"<<BufHeight<<endl;

			int icol = (int)row_col_List[j].col_idx;
			ossim_int32 samp = START_SAMP + icol*TILE_WIDTH;
			//列末尾小块处理
			if (icol == tilecols - 1)
			{
				BufWidth = nWidth - (tilecols - 1) * TILE_WIDTH;
				BufWidth = min(BufWidth, TILE_WIDTH);
			}
			rectList.push_back(ossimIrect(ossimIpt(samp, line), ossimIpt(samp + BufWidth - 1, line + BufHeight - 1)));
		}
		threads[i]->m_totalBlocks = numForTask;
		threads[i]->setRect(rectList);
		threads[i]->start();
	}

	bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished

		
	// mergeTiePoints
	theTset.clearTiePoints();
	for (int i = 0;i < theThreadNum;++i)
	{
		//cout<<"thread "<<i+1<<": "<<threads[i]->m_theTset.getTiePoints().size()<<endl;
		threads[i]->mergeTiePoints(theTset);
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("\r%d%%\n", 100);
		//setPercentComplete(100.0);
	}

#define MAX_PATH 260
	//reduce
	if (rank == 0)
	{
		// The master thread will need to receive all computations from all other threads.
		MPI_Status status;

		// MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
		// We need to go and receive the data from all other threads.
		// The arbitrary tag we choose is 1, for now.

		for (int i = 1; i < nTasks; i++)
		{
			char gcpfilename[MAX_PATH];
			MPI_Recv(gcpfilename, MAX_PATH, MPI_CHAR, i, 3, MPI_COMM_WORLD, &status);
			//printf("%s\n", gcpfilename);
			radiBlockTieGptSet* blockTieGptSet = new radiBlockTieGptSet;
			ossimFilename blockTieGptSetFile(gcpfilename);
			ossimXmlDocument gmlDoc(blockTieGptSetFile);
			ossimRefPtr<ossimXmlNode> aGmlNode = gmlDoc.getRoot();
			if (aGmlNode.valid())
			{
				blockTieGptSet->importFromGmlNode(aGmlNode);
				mergeTiePoints(theTset, blockTieGptSet);
			}
			if (blockTieGptSetFile.exists())
			{
				boost::filesystem::remove(blockTieGptSetFile.c_str());
			}
		}
	}
	else
	{
		// We are finished with the results in this thread, and need to send the data to thread 1.
		// MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
		// The destination is thread 0, and the arbitrary tag we choose for now is 1.
		char gcpfilename[MAX_PATH];
		sprintf_s(gcpfilename, "%s\\%s_%03d.xml", theFilename.path().c_str(), theFilename.fileNoExtension().c_str(), rank);
		//printf("%s\n", gcpfilename);
		//printf("%d\n", theTset.size());
		writePoints(gcpfilename);
		MPI_Send(gcpfilename, MAX_PATH, MPI_CHAR, 0, 3, MPI_COMM_WORLD);
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
