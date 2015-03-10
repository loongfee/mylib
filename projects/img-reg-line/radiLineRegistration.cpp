// radiLineRegistration.cpp
#include <iostream>

#include "radiLineRegistration.h"
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

#include "LineExtract.h"
#include "ecmlrAffine.h"

#include "lineFeature.h"
#include "LineMatchFLANN.h"

//#include "lbd/EDLineDetector.hh"
//#include "lbd/LineDescriptor.hh"
//#include "lbd/PairwiseLineMatching.hh"
#include "mySpectralMatching.h"

#include <ctime>
#include <iostream>


extern "C"{
#include <vl/generic.h>
#include <vl/stringop.h>
#include <vl/pgm.h>
#include <vl/sift.h>
#include <vl/getopt_long.h>
};

#include <QString>
#include <QDir>
#include <QXmlStreamWriter>
#include <QMapIterator>
using namespace std;
#pragma comment(lib, "vl.lib")

RTTI_DEF2(radiLineRegistration, "radiLineRegistration", ossimOutputSource, ossimProcessInterface);

radiLineRegistration::radiLineRegistration()
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

radiLineRegistration::radiLineRegistration(const radiLineRegistration& a)
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
	theUseGeographic(a.theUseGeographic),
	theDebug(a.theDebug),
	theStoreFlag(a.theStoreFlag)
{
}

radiLineRegistration::~radiLineRegistration()
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

ossimString radiLineRegistration::getRole() const
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
      cerr<<"radiLineRegistration::getRole unknown output projection, need to supply it"<<endl;
   }

   return role;
}

ossimImageHandler*  radiLineRegistration::getProjectionHandler()
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
			   cerr<<"radiLineRegistration"<<"::execute can't create handler for slave image  "<< theSlave <<endl;
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
      cerr<<"radiLineRegistration::getProjectionHandler cannot get handler for " << getRole() <<endl;
   }
   return projHandler;
}

ossimRefPtr<ossimImageGeometry> radiLineRegistration::getOutputImageGeometry()
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
ossimMapProjection* radiLineRegistration::getOutputProjection()
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
         cerr << "radiLineRegistration::getOutputProjection cannot create projection from " << getRole() <<" geometry." <<endl;
      }
   }
   else
   {
      cerr << "radiLineRegistration::getOutputProjection cannot get "
           <<getRole() << " geometry." << endl;
   }

   return mop;
}

// buildRenerer() - builds renderer for an imageSource
// accounts for :
//  -scale factor
//  -required projection
bool radiLineRegistration::buildRenderer(
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
         cerr<<"radiLineRegistration"<<"::buildRenderer cannot get projection from master/slave"<<endl;
         return false;
      }      
   }
   else
   {
      cerr<<"radiLineRegistration"<<"::buildRenderer NULL source"<<endl;
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
ossimIpt radiLineRegistration::slave2master(ossimProjection* slaveProjection,
									ossimProjection* masterProjection,
									ossimIpt slaveDpt)
{
	ossimDpt masterDpt;
	ossimGpt gpt;
	slaveProjection->lineSampleToWorld(slaveDpt, gpt);
	masterProjection->worldToLineSample(gpt, masterDpt);
	return ossimIpt(masterDpt.x, masterDpt.y);
}

ossimDpt radiLineRegistration::slave2master(ossimProjection* slaveProjection,
									ossimProjection* masterProjection,
									ossimDpt slaveDpt)
{
	ossimDpt masterDpt;
	ossimGpt gpt;
	slaveProjection->lineSampleToWorld(slaveDpt, gpt);
	masterProjection->worldToLineSample(gpt, masterDpt);
	return masterDpt;
}

ossimDpt radiLineRegistration::slave2master(ossimDpt slaveDpt)
{
	ossimDpt masterDpt;
	ossimGpt gpt;
	theSlaveProjection->lineSampleToWorld(slaveDpt, gpt);
	theMasterProjection->worldToLineSample(gpt, masterDpt);
	return masterDpt;
}

ossimIpt radiLineRegistration::slave2master(ossimIpt slaveDpt)
{
	ossimDpt masterDpt;
	ossimGpt gpt;
	theSlaveProjection->lineSampleToWorld(slaveDpt, gpt);
	theMasterProjection->worldToLineSample(gpt, masterDpt);
	return ossimIpt(masterDpt.x, masterDpt.y);
}

ossimDrect radiLineRegistration::imageRect2World(const ossimIrect& imageRect, ossimProjection* projection)
{
	ossimGpt p[4];
	if (point_type::control == thePointType)
	{
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

ossimIrect radiLineRegistration::worldRect2Image(const ossimDrect& wroldRect, ossimProjection* projection)
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

ossimDrect radiLineRegistration::worldRectIntersection(const ossimDrect& r1, const ossimDrect& r2)
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

bool radiLineRegistration::createTileMat(const ossimRefPtr<ossimCastTileSourceFilter>& cast, const ossimIrect& rect, cv::Mat& outMat, ossim_uint32 resLevel)
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

bool radiLineRegistration::getMasterList(ossimFilename spatial_index_file, vector<ossimFilename>& masterList,
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

bool radiLineRegistration::execute()
{
	bool result=true;
	// check point type
	ossimRefPtr<ossimImageHandler> mHandler;
	ossimRefPtr<ossimProjection> mProjection;
	if (theMaster.ext().upcase() != "SHP")
	{ 
		mHandler = ossimImageHandlerRegistry::instance()->open(theMaster);
		if (NULL == mHandler.get())
		{
			cerr<<"radiLineRegistration"<<"::execute can't open master image  "<< theMaster <<endl;
			return false;
		}
		if(thePointType == point_type::control &&
			(NULL == mHandler->getImageGeometry().get() || NULL == (mProjection = mHandler->getImageGeometry()->getProjection()).get()))
		{
			// 无投影
			cerr<<"radiLineRegistration::execute can't get control points as"
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
		cerr<<"radiLineRegistration"<<"::execute can't create handler for slave image  "<< theSlave <<endl;
		return false;
	}
	theSChain->add(handlerS.get());

	ossim_uint32 sbc = handlerS->getNumberOfOutputBands();
	//add a band selector
	ossim_uint32 sb = getSlaveBand();
	if (sb>=sbc) 
	{
		cerr<<"radiLineRegistration"<<"::execute Warning not enough bands in slave, only "<< sbc <<endl;
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
	if (theMaster.ext().upcase() != "SHP" &&
		!(NULL == mHandler->getImageGeometry().get() || NULL == (mProjection = mHandler->getImageGeometry()->getProjection()).get()))
	{
		ossimIrect mBoundary = mHandler->getBoundingRect(0);
		ossimDrect mWorldRect = imageRect2World(mBoundary, mProjection.get());
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
				mProjection->lineSampleToWorld(mBoundary.ul(), gpt);
				theSlaveProjection->worldToLineSample(gpt, mul);
				mProjection->lineSampleToWorld(mBoundary.lr(), gpt);
				theSlaveProjection->worldToLineSample(gpt, mlr);
			}
			else
			{
				mProjection->lineSampleHeightToWorld(mBoundary.ul(), 0.0, gpt);
				theSlaveProjection->worldToLineSample(gpt, mul);
				mProjection->lineSampleHeightToWorld(mBoundary.lr(), 0.0, gpt);
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

	m_TieLines.clear();
	// -- 4 -- run
	result = getAllFeatures();

	cout << m_TieLines.size() << " tie points are found." << endl;

	theHasRun = true;
	return true;
}

void radiLineRegistration::appendTiePoints(const ossimFilename& filename)
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

void radiLineRegistration::appendControlPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
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


void radiLineRegistration::writeTiePoints(const ossimFilename& filename)
{
	fstream ofs;
	ofs.open(theFilename.c_str(), ios_base::out);
	ofs.close();
	appendTiePoints(filename);
}
void radiLineRegistration::writeControlPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
{
	fstream ofs;
	ofs.open(theFilename.c_str(), ios_base::out);
	ofs.close();
	appendControlPoints(filename, pMapProjection);
}

void radiLineRegistration::appendPoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
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

void radiLineRegistration::writePoints(const ossimFilename& filename, ossimMapProjection* pMapProjection/* = NULL*/)
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


static void _prepareImgAndDrawKeylines(const cv::Mat& img1,
	const cv::Mat& img2,
	const vector<ossimTDpt> tptList,
	cv::Mat& outImg,
	const cv::Scalar& singlePointColor)
{
	Size size(img1.cols + img2.cols, MAX(img1.rows, img2.rows));

	outImg.create(size, CV_MAKETYPE(img1.depth(), 3));
	cv::Mat outImg1 = outImg(Rect(0, 0, img1.cols, img1.rows));
	cv::Mat outImg2 = outImg(Rect(img1.cols, 0, img2.cols, img2.rows));
	if (img1.type() == CV_8U)
		cvtColor(img1, outImg1, CV_GRAY2BGR);
	else
		img1.copyTo(outImg1);

	if (img2.type() == CV_8U)
		cvtColor(img2, outImg2, CV_GRAY2BGR);
	else
		img2.copyTo(outImg2);
	// draw keypoints
	for (size_t i = 0; i < tptList.size(); i++)
	{
		ossimDpt slavePt = tptList[i].getSlavePoint();
		ossimDpt masterPt = tptList[i].getMasterPoint();
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
}

//double point2polarline(const fPoint& p, const cvline_polar& l)
//{
//	double dist = -p.x * sin(l.theta) + p.y * cos(l.theta) - l.rho;
//	return dist;
//}
//
//
//double distance1(const cvline_polar& l1, const cvline_polar& l2)
//{
//	//return (point2polarline(l1.pt1, l2) + point2polarline(l1.pt2, l2))*0.5;
//	double d1 = point2polarline(l1.pt1, l2);
//	double d2 = point2polarline(l1.pt2, l2);
//	return sqrt(d1*d1 + d2*d2);
//}
//
//fPoint distance2(const cvline_polar& l1, const cvline_polar& l2)
//{
//
//	double centerx1 = (l1.pt1.x + l1.pt2.x)*0.5;
//	double centery1 = (l1.pt1.y + l1.pt2.y)*0.5;
//	double centerx2 = (l2.pt1.x + l2.pt2.x)*0.5;
//	double centery2 = (l2.pt1.y + l2.pt2.y)*0.5;
//	//fPoint shift = fPoint(centerx1-centerx2, centery1-centery2);
//	double shift = sqrt((centerx1 - centerx2)*(centerx1 - centerx2) + (centery1 - centery2)*(centery1 - centery2));
//	//double shift_weight = 1.0E-2;
//	//double shift_weight = 0.08;
//
//	fPoint dist;
//	// dist1
//	dist.x = distance1(l1, l2);
//	// dist2
//	dist.y = distance1(l2, l1);
//
//	//shift = max(dist.x,dist.y)*exp(-1.0/(shift*shift+DBL_EPSILON));
//	//shift = 3.0*exp(-m_delta_2/(shift*shift+DBL_EPSILON));
//	shift = 1.0*exp(-1.0 / (shift + DBL_EPSILON));
//	//if((dist.x*dist.x+dist.y*dist.y) < 4.0)
//	//{
//	dist.x += shift;
//	dist.y += shift;
//	//}
//
//	//dist.x = fabs(l1.rho - l2.rho);
//	//dist.y = scale_angle * fabs(l1.theta - l2.theta);
//	////dist.y = dist.y > lammda_rho*dist.x ? lammda_rho*dist.x : dist.y;
//	////dist.y = fabs(tan(l1.theta) - tan(l2.theta));
//	//return sqrt(dist.x*dist.x + dist.y*dist.y);
//	return dist;
//}
//
//cvline_polar forward(const cvline_polar& l, const vector<double>& parameters)
//{
//	cvline_polar outLine;
//	fPoint line[2];
//	line[0].x = parameters[0] + parameters[1] * l.pt1.x + parameters[2] * l.pt1.y;
//	line[0].y = parameters[3] + parameters[4] * l.pt1.x + parameters[5] * l.pt1.y;
//	line[1].x = parameters[0] + parameters[1] * l.pt2.x + parameters[2] * l.pt2.y;
//	line[1].y = parameters[3] + parameters[4] * l.pt2.x + parameters[5] * l.pt2.y;
//	outLine = line2polar(line);
//	return outLine;
//}
//
//struct dataStruct
//{
//	vector<cvline_polar> slave_vec_lines;
//	vector<cvline_polar> master_vec_lines;
//	vector<int> slaveAssign;
//	std::vector<double> parameters;
//};
//
//void levmar_function_fvec(double *param, double *hx, int nparameter, int nequation, void *adata)
//{
//	dataStruct *pThis = (dataStruct*)adata;
//	int nX = (int)pThis->slave_vec_lines.size();
//	int nY = (int)pThis->master_vec_lines.size();
//
//	int pos = 0;
//	int i;
//	for (i = 0; i<nparameter; ++i)
//	{
//		pThis->parameters[i] = param[i];
//	}
//
//	int c = 0;
//
//	for (int i = 0; i < nX; ++i)
//	{
//		cvline_polar outPt = forward(pThis->slave_vec_lines[i], pThis->parameters);
//
//		fPoint dist = distance2(outPt, pThis->master_vec_lines[pThis->slaveAssign[i]]);
//
//		hx[c++] = dist.x;
//		hx[c++] = dist.y;
//	}
//}
//
//void levmar_function_jac(double *param, double *jac, int nparameter, int nequation, void *adata)
//{
//	dataStruct *pThis = (dataStruct*)adata;
//	int nX = (int)pThis->slave_vec_lines.size();
//	int nY = (int)pThis->master_vec_lines.size();
//
//	int pos = 0;
//	int i;
//	for (i = 0; i<nparameter; ++i)
//	{
//		pThis->parameters[i] = param[i];
//	}
//
//	double pstep_scale = 1e-4;
//	static double den = 0.5 / pstep_scale;
//	int c = 0;
//
//	for (int i = 0; i < nX; ++i)
//	{
//		cvline_polar outPt = forward(pThis->slave_vec_lines[i], pThis->parameters);
//
//		for (int p = 0; p<nparameter; ++p)
//		{
//			double middle = pThis->parameters[p];
//			pThis->parameters[p] = middle + pstep_scale;
//			cvline_polar outLine1 = forward(pThis->slave_vec_lines[i], pThis->parameters);
//			//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
//			fPoint dist1 = distance2(outLine1, pThis->master_vec_lines[pThis->slaveAssign[i]]);
//
//			pThis->parameters[p] = middle - pstep_scale;
//			cvline_polar outLine2 = forward(pThis->slave_vec_lines[i], pThis->parameters);
//			//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
//			fPoint dist2 = distance2(outLine2, pThis->master_vec_lines[pThis->slaveAssign[i]]);
//
//			pThis->parameters[p] = middle;
//
//			double derivative_x = (dist1.x - dist2.x) * den;
//			double derivative_y = (dist1.y - dist2.y) * den;
//
//			jac[c*nparameter + p] = derivative_x;
//			jac[(c + 1)*nparameter + p] = derivative_y;
//		}
//		c += 2;
//	}
//}
//
//
//bool parameters_optimization(const vector<cvline_polar>& slave_vec_lines, const vector<cvline_polar>& master_vec_lines, 
//	const vector<int>& slaveAssign, std::vector<double>& parameters)
//{
//	int nX = (int)slave_vec_lines.size();
//	int nY = (int)master_vec_lines.size();
//	
//	VectorXd b(6);
//	MatrixXd A(6, 6);
//	b.fill(0.0);
//	A.fill(0.0);
//	
//	int nparam = 6;
//	std::vector<double> cparm(nparam);
//	
//	for (int i = 0; i<nparam; ++i) cparm[i] = parameters[i];
//	
//	// lm
//	double *p = &cparm[0];
//	
//	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
//	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
//	opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 
//	
//	dataStruct ds;
//	ds.slave_vec_lines = slave_vec_lines;
//	ds.master_vec_lines = master_vec_lines;
//	ds.parameters = parameters;
//	ds.slaveAssign = slaveAssign;
//	//int ret = dlevmar_dif(levmar_function_fvec, &cparm[0], &x[0], nparam, 2*nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
//	int ret = dlevmar_der(levmar_function_fvec, levmar_function_jac, p, NULL, nparam, 2*nX, 1000, opts, info, NULL, NULL, &ds); // with analytic Jacobian
//	//delete []x;
//	
//	for (int i = 0; i<nparam; ++i) parameters[i] = p[i];
//	return true;
//}
//
//
//bool parameters_optimization(const vector<cvline_polar>& slave_vec_lines, const vector<cvline_polar>& master_vec_lines, std::vector<double>& parameters)
//{
//	int nX = (int)slave_vec_lines.size();
//	vector<int> slaveAssign;
//	for (int i = 0; i < nX; i++)
//	{
//		slaveAssign.push_back(i);
//	}
//
//	return parameters_optimization(slave_vec_lines, master_vec_lines, slaveAssign, parameters);
//}


//vector<double> parameters_optimization(const vector<cvline_polar>& slave_vec_lines, const vector<cvline_polar>& master_vec_lines, const vector<int>& slaveAssign)
//{
//	int nX = (int)slave_vec_lines.size();
//	int nY = (int)master_vec_lines.size();
//
//	VectorXd b(6);
//	MatrixXd A(6, 6);
//	b.fill(0.0);
//	A.fill(0.0);
//
//	int c = 0;
//	for (int i = 0; i < (int)nX; ++i)
//	{
//		double x1 = slave_vec_lines[i].pt1.x;
//		double y1 = slave_vec_lines[i].pt1.y;
//		double x2 = slave_vec_lines[i].pt2.x;
//		double y2 = slave_vec_lines[i].pt2.y;
//
//		double sin_theta = sin(master_vec_lines[slaveAssign[i]].theta);
//		double cos_theta = cos(master_vec_lines[slaveAssign[i]].theta);
//		double rho = master_vec_lines[slaveAssign[i]].rho;
//		double len = segment_length(master_vec_lines[slaveAssign[i]]);
//
//		vector<double> coeff(6);
//		coeff[0] = -sin_theta;
//		coeff[1] = -sin_theta * x1;
//		coeff[2] = -sin_theta * y1;
//		coeff[3] = cos_theta;
//		coeff[4] = cos_theta * x1;
//		coeff[5] = cos_theta * y1;
//
//		for (int p1 = 0; p1<6; ++p1)
//		{
//			b[p1] += coeff[p1] * rho;
//			for (int p2 = 0; p2<6; ++p2)
//			{
//				A(p1, p2) += coeff[p1] * coeff[p2];
//			}
//		}
//
//		coeff[0] = -sin_theta;
//		coeff[1] = -sin_theta * x2;
//		coeff[2] = -sin_theta * y2;
//		coeff[3] = cos_theta;
//		coeff[4] = cos_theta * x2;
//		coeff[5] = cos_theta * y2;
//		for (int p1 = 0; p1<6; ++p1)
//		{
//			b[p1] += coeff[p1] * rho;
//			for (int p2 = 0; p2<6; ++p2)
//			{
//				A(p1, p2) += coeff[p1] * coeff[p2];
//			}
//		}
//	}
//
//	Eigen::MatrixXd damper = MatrixXd::Identity(6, 6);
//	double tmp = A.determinant();
//	if (fabs(tmp) > DBL_EPSILON)
//	{
//		Eigen::VectorXd parameter = A.inverse() * b;
//		vector<double> affineParameters(6);
//		for (int i = 0; i < 6; ++i) affineParameters[i] = parameter[i];
//		return affineParameters;
//	}
//	else
//	{
//		cout << "ERROR : cannot find a solution" << endl;
//		return vector<double>(6);
//	}
//	return vector<double>(6);
//}

void icl(const vector<cvline_polar>& slave_vec_lines, const vector<cvline_polar>& master_vec_lines,
	vector<cvline_polar>& matched_slave_lines, vector<cvline_polar>& matched_master_lines,
	vector<double>& parameters, double threshould = 5.0)
{
	int MAX_ITERATION = 200;
	int iTime = 0;
	double epsilon = 1e-6;
	double delta = 1e10;
	int nSlave = (int)slave_vec_lines.size();
	int nMaster = (int)master_vec_lines.size();
	if (nSlave < 1 || nMaster < 1)
	{
		return;
	}

	vector<int> slaveAssign(nSlave, 0);
	vector<double> distanceAssign(nSlave, 0);
	double rms = 1e10;
	do
	{
		vector<cvline_polar> trans_lines(nSlave);
		// transform
		for (int i = 0; i < nSlave; i++)
		{
			trans_lines[i] = forward(slave_vec_lines[i], parameters);
		}

		// assign
		for (int i = 0; i < nSlave; i++)
		{
			double minDis = 1e20;
			int minIdx = 0;
			for (int j = 0; j < nMaster; j++)
			{
				fPoint dis2 = distance2(trans_lines[i], master_vec_lines[j]);
				double dis = sqrt(dis2.x*dis2.x + dis2.y*dis2.y);
				if (minDis > dis)
				{
					minDis = dis;
					minIdx = j;
				}
			}
			slaveAssign[i] = minIdx;
			distanceAssign[i] = minDis;
		}

		// check errors
		double rms_new = 0.0;
		for (int i = 0; i < nSlave; i++)
		{
			rms_new += distanceAssign[i] * distanceAssign[i];
		}
		rms_new = sqrt(rms_new / nSlave);
		delta = fabs(rms - rms_new);
		rms = rms_new;

		// estimation
		//parameters = affine_estimation2(slave_vec_lines, master_vec_lines, slaveAssign);
		parameters_optimization(slave_vec_lines, master_vec_lines, slaveAssign, parameters);
		printf("%lf ", rms);
	} while (iTime < MAX_ITERATION && delta > epsilon);

	if (iTime == MAX_ITERATION)
	{
		printf("The maximum iteration time is reached, and the result may be not accurate.\n");
	}

	// nearest neighbour assign
	for (int i = 0; i < nSlave; i++)
	{
		if (distanceAssign[i] < threshould)
		{
			matched_slave_lines.push_back(slave_vec_lines[i]);
			matched_master_lines.push_back(master_vec_lines[slaveAssign[i]]);
		}
	}

	printf("rms: %lf\n", rms);
	printf("%d matched.\n\n", matched_slave_lines.size());
}

//void lbd_match(cv::Mat slaveMat, cv::Mat masterMat)
//{
//	//extract lines, compute their descriptors and match lines
//	LineDescriptor lineDesc;
//	PairwiseLineMatching lineMatch;
//	
//	ScaleLines   linesInReference;
//	ScaleLines   linesInSource;
//	std::vector<unsigned int> matchResult;
//	
//	cvtColor(masterMat, masterMat, CV_RGB2GRAY);
//	cvtColor(slaveMat, slaveMat, CV_RGB2GRAY);
//	
//	//cv::Mat debugImage;
//	//cv::Sobel(srcImage, debugImage, CV_8U, 1, 0, 3, 1.0, 0.0, cv::BORDER_CONSTANT);
//	//cv::imwrite("E:\\debug.png", debugImage);
//	
//	lineDesc.GetLineDescriptor(masterMat, linesInReference);
//	lineDesc.GetLineDescriptor(slaveMat, linesInSource);
//	
//	lineMatch.LineMatching(linesInReference,linesInSource,matchResult);
//		
//	vector<cvline_polar> matched_ref_lines_list;
//	vector<cvline_polar> matched_src_lines_list;
//	
//	
//	int nMatched = (int)matchResult.size()/2;
//	for(int i = 0; i < nMatched; i++)
//	{
//		cvline_polar refLine;
//		cvline_polar srcLine;
//	
//		fPoint line[2];
//		//line[0].x = linesInReference[matchRefIndex[i]][0].startPointX;
//		//line[0].y = linesInReference[matchRefIndex[i]][0].startPointY;
//		//line[1].x = linesInReference[matchRefIndex[i]][0].endPointX;
//		//line[1].y = linesInReference[matchRefIndex[i]][0].endPointY;
//		line[0].x = linesInReference[matchResult[2*i]][0].sPointInOctaveX;
//		line[0].y = linesInReference[matchResult[2*i]][0].sPointInOctaveY;
//		line[1].x = linesInReference[matchResult[2*i]][0].ePointInOctaveX;
//		line[1].y = linesInReference[matchResult[2*i]][0].ePointInOctaveY;
//		refLine = line2polar(line);
//	
//		//line[0].x = linesInSource[matchSrcIndex[i]][0].startPointX;
//		//line[0].y = linesInSource[matchSrcIndex[i]][0].startPointY;
//		//line[1].x = linesInSource[matchSrcIndex[i]][0].endPointX;
//		//line[1].y = linesInSource[matchSrcIndex[i]][0].endPointY;
//		line[0].x = linesInSource[matchResult[2*i+1]][0].sPointInOctaveX;
//		line[0].y = linesInSource[matchResult[2*i+1]][0].sPointInOctaveY;
//		line[1].x = linesInSource[matchResult[2*i+1]][0].ePointInOctaveX;
//		line[1].y = linesInSource[matchResult[2*i+1]][0].ePointInOctaveY;
//		srcLine = line2polar(line);
//	
//		matched_ref_lines_list.push_back(refLine);
//		matched_src_lines_list.push_back(srcLine);
//	}
//	
//	//std::vector<short> matchRefIndex;
//	//std::vector<short> matchSrcIndex;
//	//lineDesc.MatchLineByDescriptor(linesInReference, 	linesInSource,
//	//	matchRefIndex, matchSrcIndex,
//	//	LineDescriptor::NearestNeighbor);
//	
//	//int nMatched = (int)matchRefIndex.size();
//	//for(int i = 0; i < nMatched; i++)
//	//{
//	//	cvline_polar refLine;
//	//	cvline_polar srcLine;
//	
//	//	fPoint line[2];
//	//	//line[0].x = linesInReference[matchRefIndex[i]][0].startPointX;
//	//	//line[0].y = linesInReference[matchRefIndex[i]][0].startPointY;
//	//	//line[1].x = linesInReference[matchRefIndex[i]][0].endPointX;
//	//	//line[1].y = linesInReference[matchRefIndex[i]][0].endPointY;
//	//	line[0].x = linesInReference[matchRefIndex[i]][0].sPointInOctaveX;
//	//	line[0].y = linesInReference[matchRefIndex[i]][0].sPointInOctaveY;
//	//	line[1].x = linesInReference[matchRefIndex[i]][0].ePointInOctaveX;
//	//	line[1].y = linesInReference[matchRefIndex[i]][0].ePointInOctaveY;
//	//	refLine = line2polar(line);
//	
//	//	//line[0].x = linesInSource[matchSrcIndex[i]][0].startPointX;
//	//	//line[0].y = linesInSource[matchSrcIndex[i]][0].startPointY;
//	//	//line[1].x = linesInSource[matchSrcIndex[i]][0].endPointX;
//	//	//line[1].y = linesInSource[matchSrcIndex[i]][0].endPointY;
//	//	line[0].x = linesInSource[matchSrcIndex[i]][0].sPointInOctaveX;
//	//	line[0].y = linesInSource[matchSrcIndex[i]][0].sPointInOctaveY;
//	//	line[1].x = linesInSource[matchSrcIndex[i]][0].ePointInOctaveX;
//	//	line[1].y = linesInSource[matchSrcIndex[i]][0].ePointInOctaveY;
//	//	srcLine = line2polar(line);
//	
//	//	matched_ref_lines_list.push_back(refLine);
//	//	matched_src_lines_list.push_back(srcLine);
//	//}
//	
//	//LineMatch(linesInReference, linesInSource, matched_ref_lines_list, matched_src_lines_list);
//}


void radiLineRegistration::VLFeatSift(const cv::Mat& inMat, vector<KeyPoint>& kpts, cv::Mat& descriptors)
{
	int noctaves = 2, nlevels = 4, o_min = 0;
	// noctaves=(int)(log(min)/log(2));
	vl_sift_pix *ImageData = new vl_sift_pix[inMat.rows * inMat.cols];
	unsigned char *Pixel;
	for (int i = 0; i<inMat.rows; i++)
	{
		for (int j = 0; j<inMat.cols; j++)
		{
			Pixel = (unsigned char*)(inMat.data + i*inMat.cols + j);
			ImageData[i*inMat.cols + j] = *(Pixel);
		}
	}
	VlSiftFilt *SiftFilt = NULL;
	SiftFilt = vl_sift_new(inMat.cols, inMat.rows, noctaves, nlevels, o_min);
	double edge_thresh = 30;  //-1 will use the default (as in matlab)
	double peak_thresh = 0.001;
	double norm_thresh = -1;
	double magnif = -1;
	double window_size = -1;
	if (peak_thresh >= 0) vl_sift_set_peak_thresh(SiftFilt, peak_thresh);
	if (edge_thresh >= 0) vl_sift_set_edge_thresh(SiftFilt, edge_thresh);
	if (norm_thresh >= 0) vl_sift_set_norm_thresh(SiftFilt, norm_thresh);
	if (magnif >= 0) vl_sift_set_magnif(SiftFilt, magnif);
	if (window_size >= 0) vl_sift_set_window_size(SiftFilt, window_size);
	int nKeyPoint = 0;
	int idx = 0;

	//descriptor = cv::Mat(cv::Size(128, 1), CV_32FC1);
	if (vl_sift_process_first_octave(SiftFilt, ImageData) != VL_ERR_EOF)
	{
		while (TRUE)
		{
			//计算每组中的关键点
			vl_sift_detect(SiftFilt);
			//遍历并绘制每个点			
			nKeyPoint += SiftFilt->nkeys;
			VlSiftKeypoint *pKeyPoint = SiftFilt->keys;
			for (int i = 0; i<SiftFilt->nkeys; i++)
			{
				VlSiftKeypoint TemptKeyPoint = *pKeyPoint;
				pKeyPoint++;
				//cv::KeyPoint kpt(float x, float y, float _size, float _angle=-1,
				//	float _response=0, int _octave=0, int _class_id=-1);
				cv::KeyPoint kpt(TemptKeyPoint.x, TemptKeyPoint.y, TemptKeyPoint.sigma / 2);
				//kpts.push_back(kpt);

				//cvDrawCircle(Image, cvPoint(TemptKeyPoint.x,TemptKeyPoint.y),TemptKeyPoint.sigma/2,CV_RGB(255,0,0));
				idx++;
				//计算并遍历每个点的方向
				double angles[4];
				int angleCount = vl_sift_calc_keypoint_orientations(SiftFilt, angles, &TemptKeyPoint);
				for (int j = 0; j<angleCount; j++)
				{
					double TemptAngle = angles[j];
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
			if (vl_sift_process_next_octave(SiftFilt) == VL_ERR_EOF)
			{
				break;
			}
			//free(pKeyPoint);
			nKeyPoint = NULL;
		}
	}
	vl_sift_delete(SiftFilt);
	delete[]ImageData;
	ImageData = NULL;
}

double affine_invariant_sim(const cvline_polar& slave_line, const cvline_polar& master_line, const vector<ossimTDpt>& tptList)
{
	double minSim = 1e30;
	int np = min(tptList.size(), 10);
	for (size_t i = 0; i < np - 1; i++)
	{
		for (size_t j = i + 1; j < np; j++)
		{
			ossimDpt slave_dpt1 = tptList[i].getSlavePoint();
			ossimDpt slave_dpt2 = tptList[j].getSlavePoint();
			ossimDpt master_dpt1 = tptList[i].getMasterPoint();
			ossimDpt master_dpt2 = tptList[j].getMasterPoint();
			double slave_dist1 = fabs(-slave_dpt1.x * sin(slave_line.theta) + slave_dpt1.y * cos(slave_line.theta) - slave_line.rho);
			double slave_dist2 = fabs(-slave_dpt2.x * sin(slave_line.theta) + slave_dpt2.y * cos(slave_line.theta) - slave_line.rho);
			double master_dist1 = fabs(-master_dpt1.x * sin(master_line.theta) + master_dpt1.y * cos(master_line.theta) - master_line.rho);
			double master_dist2 = fabs(-master_dpt2.x * sin(master_line.theta) + master_dpt2.y * cos(master_line.theta) - master_line.rho);
			double sim = exp(-fabs(slave_dist1 / master_dist1 - slave_dist2 / master_dist2));
			//double dis = fabs(log(slave_dist1 * master_dist2) - log(slave_dist2 * master_dist1));

			// position
			ossimDpt slavePt1 = (slave_dpt1 + slave_dpt2) * 0.5;
			ossimDpt slavePt2 = ossimDpt((slave_line.pt1.x + slave_line.pt2.x)*0.5, (slave_line.pt1.y + slave_line.pt2.y)*0.5);
			ossimDpt masterPt1 = (master_dpt1 + master_dpt2) * 0.5;
			ossimDpt masterPt2 = ossimDpt((master_line.pt1.x + master_line.pt2.x)*0.5, (master_line.pt1.y + master_line.pt2.y)*0.5);
			double slave_a = sqrt((slavePt1.x - slavePt2.x)*(slavePt1.x - slavePt2.x) + (slavePt1.y - slavePt2.y)*(slavePt1.y - slavePt2.y));
			double master_a = sqrt((masterPt1.x - masterPt2.x)*(masterPt1.x - masterPt2.x) + (masterPt1.y - masterPt2.y)*(masterPt1.y - masterPt2.y));
			sim *= exp(-fabs(fabs(slave_a / master_a)-1.0)*1.0);
			if (minSim > sim)
			{
				minSim = sim;
			}
		}
	}
	return minSim;
}

double affine_invariant(const cvline_polar& line_polar, const ossimDpt& dpt1, const ossimDpt& dpt2)
{
	double dist1 = fabs(-dpt1.x * sin(line_polar.theta) + dpt1.y * cos(line_polar.theta) - line_polar.rho);
	double dist2 = fabs(-dpt2.x * sin(line_polar.theta) + dpt2.y * cos(line_polar.theta) - line_polar.rho);
	return log(dist1 / (dist2 + 1e-30));
}

int radiLineRegistration::runMatchParallel(const cv::Mat& slaveMat, const cv::Mat& masterMat, vector<ossimTDpt>& theTies, void* pData, bool bDebug)
{
	radiLineRegistration* pThis = (radiLineRegistration*)pData;
	// detect corners
	//cv::initModule_nonfree();
	cv::initModule_features2d();
	std::vector<KeyPoint> skeypoints, mkeypoints;
	cv::Mat sdescriptors, mdescriptors;

	if (bDebug)
	{
		//cv::imwrite("slave.png", slaveMat);
		//cv::imwrite("master.png", masterMat);

		int thickness = 2;
		cv::Mat img_outImage;
		drawLineMatches(slaveMat, vector<cvline_polar>(), masterMat, vector<cvline_polar>(),
			img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
			vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
		cv::imwrite("image.png", img_outImage);
	}
	//vector<cvline_polar> slave_vec_lines;
	//vector<cvline_polar> master_vec_lines;
	//vector<Vec4i> lines;
	//cv::Mat dst, cdst;
	//cv::blur(slaveMat, slaveMat, cv::Size(3, 3));
	//cv::Canny(slaveMat, dst, 50, 200, 3);
	//cvtColor(dst, cdst, CV_GRAY2BGR);
	//cv::HoughLinesP(dst, lines, 1, CV_PI / 180, 50, 50, 10);
	//for (size_t i = 0; i < lines.size(); i++)
	//{
	//	Vec4i l = lines[i];

	//	cvline_polar outLine;
	//	fPoint line[2];
	//	line[0].x = l[0];
	//	line[0].y = l[1];
	//	line[1].x = l[2];
	//	line[1].y = l[3];
	//	outLine = line2polar(line);
	//	slave_vec_lines.push_back(outLine);

	//	cv::line(cdst, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0, 0, 255), 3, CV_AA);
	//}
	//cv::imwrite("slave_lines.png", cdst);


	//lines.clear();
	//cv::blur(masterMat, masterMat, cv::Size(3, 3));
	//cv::Canny(masterMat, dst, 50, 200, 3);
	//cvtColor(dst, cdst, CV_GRAY2BGR);
	//cv::HoughLinesP(dst, lines, 1, CV_PI / 180, 50, 10, 2);
	//for (size_t i = 0; i < lines.size(); i++)
	//{
	//	Vec4i l = lines[i];

	//	cvline_polar outLine;
	//	fPoint line[2];
	//	line[0].x = l[0];
	//	line[0].y = l[1];
	//	line[1].x = l[2];
	//	line[1].y = l[3];
	//	outLine = line2polar(line);
	//	master_vec_lines.push_back(outLine);

	//	cv::line(cdst, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0, 0, 255), 3, CV_AA);
	//}
	//cv::imwrite("master_lines.png", cdst);

	IplImage slaveIplImageMat = slaveMat;
	IplImage masterIplImageMat = masterMat;
	vector<cvline_polar> slave_vec_lines;
	vector<cvline_polar> master_vec_lines;
	double scale = 0.9;
	int selection_num = 500;
	IplImage* slave_line = get_lines(&slaveIplImageMat, slave_vec_lines, scale, selection_num);
	IplImage* master_line = get_lines(&masterIplImageMat, master_vec_lines, scale, selection_num);
	if (!slave_line || !master_line)
	{
		return match_state::match_failed;
	}
	// detect straight line
	if (slave_vec_lines.size() < 10)
	{
		return match_state::slave_faild;
	}
	if (master_vec_lines.size() < 10)
	{
		return match_state::master_faild;
	}

	if (bDebug)
	{
		//cvSaveImage("slave_lines.png", slave_line);
		//cvSaveImage("master_lines.png", master_line);

		int thickness = 2;
		cv::Mat img_outImage;
		drawLineMatches(cv::cvarrToMat(slave_line), vector<cvline_polar>(), cv::cvarrToMat(master_line), vector<cvline_polar>(),
			img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
			vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
		cv::imwrite("lines.png", img_outImage);
	}

	cvReleaseImage(&slave_line);
	cvReleaseImage(&master_line);

	//// icl
	//vector<double> parameters(6);
	//parameters[0] = 25.0;
	//parameters[1] = 1.0;
	//parameters[2] = 0.0;
	//parameters[3] = 25.0;
	//parameters[4] = 0.0;
	//parameters[5] = 1.0;

	//vector<cvline_polar> matched_slave_lines;
	//vector<cvline_polar> matched_master_lines;
	//icl(slave_vec_lines, master_vec_lines,
	//	matched_slave_lines, matched_master_lines,
	//	parameters, 5.0);

	//int nMatched = (int)matched_slave_lines.size();

	//if (nMatched < 4)
	//{
	//	return match_state::match_failed;
	//}

	//vector<double> modelParameters = affine_estimation2(matched_slave_lines, matched_master_lines);

	//double delta_energy = 1.0 - fabs(modelParameters[1] * modelParameters[5] - modelParameters[2] * modelParameters[4]);
	//if (fabs(delta_energy) > 0.5)
	//{
	//	return match_state::match_failed;
	//}

	//if (bDebug)
	//{
	//	int thickness = 2;
	//	cv::Mat img_outImage;
	//	drawLineMatches(slaveMat, matched_slave_lines, masterMat, matched_master_lines,
	//		img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
	//		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//	cv::imwrite("matched.png", img_outImage);
	//}

	//// lbd_geo
	//vector<vector<double>> master_LineFeature;
	//vector<vector<double>> slave_LineFeature;
	//vector<int> master_LineNum;
	//vector<int> slave_LineNum;
	//const int R = 23;
	//LineFeature(&masterIplImageMat, master_vec_lines, master_LineFeature, master_LineNum, R);
	//LineFeature(&slaveIplImageMat, slave_vec_lines, slave_LineFeature, slave_LineNum, R);

	//vector<cvline_polar> matched_slave_lines;
	//vector<cvline_polar> matched_master_lines;
	//LineMatch(masterMat, slaveMat, master_LineFeature, slave_LineFeature,
	//	master_LineNum, slave_LineNum, master_vec_lines, slave_vec_lines,
	//	matched_master_lines, matched_slave_lines);

	//int nMatched = (int)matched_slave_lines.size();

	//if (nMatched < 4)
	//{
	//	return match_state::match_failed;
	//}
	//
	//vector<double> modelParameters = affine_estimation2(matched_slave_lines, matched_master_lines);
	////vector<double> modelParameters(6);
	////modelParameters[0] = 25.0;
	////modelParameters[1] = 1.0;
	////modelParameters[2] = 0.0;
	////modelParameters[3] = 25.0;
	////modelParameters[4] = 0.0;
	////modelParameters[5] = 1.0;
	////parameters_optimization(slave_vec_lines, master_vec_lines, modelParameters);

	//double delta_energy = 1.0 - fabs(modelParameters[1] * modelParameters[5] - modelParameters[2] * modelParameters[4]);
	//if (fabs(delta_energy) > 0.5)
	//{
	//	return match_state::match_failed;
	//}

	//if (bDebug)
	//{
	//	int thickness = 2;
	//	cv::Mat img_outImage;
	//	drawLineMatches(slaveMat, matched_slave_lines, masterMat, matched_master_lines,
	//		img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
	//		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//	cv::imwrite("matched.png", img_outImage);
	//}


	//return match_state::success;

	//// emlr
	//vector<cvline_polar> matched_master_lines;
	//vector<cvline_polar> matched_slave_lines;
	//vector<double> parameters(6);
	//parameters[0] = 25.0;
	//parameters[1] = 1.0;
	//parameters[2] = 0.0;
	//parameters[3] = 25.0;
	//parameters[4] = 0.0;
	//parameters[5] = 1.0;

	//ecmlrAffine* ecmlr = new ecmlrAffine(slave_vec_lines, master_vec_lines);
	//ecmlr->setParameters(parameters);
	//vector<double> modelParameters;
	//if(!ecmlr->solve(modelParameters))
	//{
	//	return match_state::match_failed;
	//}

	//int nMatched = (int)ecmlr->m_matched_src_ls.size();

	//if (nMatched < 6)
	//{
	//	return match_state::match_failed;
	//}

	//double delta_energy = 1.0 - fabs(modelParameters[1] * modelParameters[5] - modelParameters[2] * modelParameters[4]);
	//if (fabs(delta_energy) > 0.5)
	//{
	//	return match_state::match_failed;
	//}

	//double maxAngle = 0.0;

	//if (bDebug)
	//{
	//	int thickness = 2;
	//	cv::Mat img_outImage;
	//	drawLineMatches(slaveMat, ecmlr->m_matched_src_ls, masterMat, ecmlr->m_matched_ref_ls,
	//		img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
	//		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//	cv::imwrite("matched.png", img_outImage);
	//}


	//return match_state::success;




	//// lbd
	////extract lines, compute their descriptors and match lines
	//LineDescriptor lineDesc;
	//PairwiseLineMatching lineMatch;

	//ScaleLines   linesInReference;
	//ScaleLines   linesInSource;
	//std::vector<unsigned int> matchResult;

	////cvtColor(masterMat, masterMat, CV_RGB2GRAY);
	////cvtColor(slaveMat, slaveMat, CV_RGB2GRAY);

	////cv::Mat debugImage;
	////cv::Sobel(srcImage, debugImage, CV_8U, 1, 0, 3, 1.0, 0.0, cv::BORDER_CONSTANT);
	////cv::imwrite("E:\\debug.png", debugImage);

	//lineDesc.GetLineDescriptor(masterMat, linesInReference);
	//lineDesc.GetLineDescriptor(slaveMat, linesInSource);

	//lineMatch.LineMatching(linesInReference, linesInSource, matchResult);

	//vector<cvline_polar> matched_master_lines;
	//vector<cvline_polar> matched_slave_lines;


	//int nMatched = (int)matchResult.size() / 2;
	//for (int i = 0; i < nMatched; i++)
	//{
	//	cvline_polar refLine;
	//	cvline_polar srcLine;

	//	fPoint line[2];
	//	//line[0].x = linesInReference[matchRefIndex[i]][0].startPointX;
	//	//line[0].y = linesInReference[matchRefIndex[i]][0].startPointY;
	//	//line[1].x = linesInReference[matchRefIndex[i]][0].endPointX;
	//	//line[1].y = linesInReference[matchRefIndex[i]][0].endPointY;
	//	line[0].x = linesInReference[matchResult[2 * i]][0].sPointInOctaveX;
	//	line[0].y = linesInReference[matchResult[2 * i]][0].sPointInOctaveY;
	//	line[1].x = linesInReference[matchResult[2 * i]][0].ePointInOctaveX;
	//	line[1].y = linesInReference[matchResult[2 * i]][0].ePointInOctaveY;
	//	refLine = line2polar(line);

	//	//line[0].x = linesInSource[matchSrcIndex[i]][0].startPointX;
	//	//line[0].y = linesInSource[matchSrcIndex[i]][0].startPointY;
	//	//line[1].x = linesInSource[matchSrcIndex[i]][0].endPointX;
	//	//line[1].y = linesInSource[matchSrcIndex[i]][0].endPointY;
	//	line[0].x = linesInSource[matchResult[2 * i + 1]][0].sPointInOctaveX;
	//	line[0].y = linesInSource[matchResult[2 * i + 1]][0].sPointInOctaveY;
	//	line[1].x = linesInSource[matchResult[2 * i + 1]][0].ePointInOctaveX;
	//	line[1].y = linesInSource[matchResult[2 * i + 1]][0].ePointInOctaveY;
	//	srcLine = line2polar(line);

	//	matched_master_lines.push_back(refLine);
	//	matched_slave_lines.push_back(srcLine);
	//}

	//nMatched = (int)matched_slave_lines.size();

	//if (nMatched < 4)
	//{
	//	return match_state::match_failed;
	//}
	//
	//vector<double> modelParameters(6);
	//modelParameters[0] = 25.0;
	//modelParameters[1] = 1.0;
	//modelParameters[2] = 0.0;
	//modelParameters[3] = 25.0;
	//modelParameters[4] = 0.0;
	//modelParameters[5] = 1.0;
	////affine_estimation2(matched_slave_lines, matched_master_lines);
	//parameters_optimization(slave_vec_lines, master_vec_lines, modelParameters);

	//double delta_energy = 1.0 - fabs(modelParameters[1] * modelParameters[5] - modelParameters[2] * modelParameters[4]);
	//if (fabs(delta_energy) > 0.5)
	//{
	//	return match_state::match_failed;
	//}

	//if (bDebug)
	//{
	//	int thickness = 2;
	//	cv::Mat img_outImage;
	//	drawLineMatches(slaveMat, matched_slave_lines, masterMat, matched_master_lines,
	//		img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
	//		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//	cv::imwrite("matched.png", img_outImage);
	//}

	//return match_state::success;

	//// spectral matching
	//mySpectralMatching spectralMatching;
	//vector<unsigned int> matches;
	//spectralMatching.Matching(slave_vec_lines, master_vec_lines, matches);
	//vector<cvline_polar> matched_slave_lines;
	//vector<cvline_polar> matched_master_lines;

	//int nMatched = (int)matches.size() / 2;
	//for (int i = 0; i < nMatched; i++)
	//{
	//	matched_slave_lines.push_back(slave_vec_lines[matches[2 * i]]);
	//	matched_master_lines.push_back(master_vec_lines[matches[2 * i+1]]);
	//}

	//if (nMatched < 4)
	//{
	//	return match_state::match_failed;
	//}

	//vector<double> modelParameters = affine_estimation2(matched_slave_lines, matched_master_lines);

	//double delta_energy = 1.0 - fabs(modelParameters[1] * modelParameters[5] - modelParameters[2] * modelParameters[4]);
	//if (fabs(delta_energy) > 0.5)
	//{
	//	return match_state::match_failed;
	//}

	//if (bDebug)
	//{
	//	int thickness = 2;
	//	cv::Mat img_outImage;
	//	drawLineMatches(slaveMat, matched_slave_lines, masterMat, matched_master_lines,
	//		img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
	//		vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
	//	cv::imwrite("matched.png", img_outImage);
	//}


	//return match_state::success;

	// affine invariant
	VLFeatSift(slaveMat, skeypoints, sdescriptors);
	if (skeypoints.size() < 10)
	{
		return match_state::slave_faild;

	}
	VLFeatSift(masterMat, mkeypoints, mdescriptors);
	if (mkeypoints.size() < 10)
	{
		//handlerM->close();
		return match_state::master_faild;
	}

	BFMatcher matcher(NORM_L1, false);
	vector< vector< DMatch >  > matches;
	matcher.knnMatch(sdescriptors, mdescriptors, matches, 2);

	vector< DMatch > good_matches;
	findGoodMatches(matches, good_matches, 0.70f);

	if (good_matches.size() < 4)
	{
		//handlerM->close();
		return match_state::match_failed;
	}

	//-- Create input data
	Eigen::MatrixXd dataPoints((int)good_matches.size(), 4);
	for (unsigned int i = 0; i < good_matches.size(); ++i)
	{
		dataPoints(i, 0) = skeypoints[good_matches[i].queryIdx].pt.x;
		dataPoints(i, 1) = skeypoints[good_matches[i].queryIdx].pt.y;
		dataPoints(i, 2) = mkeypoints[good_matches[i].trainIdx].pt.x;
		dataPoints(i, 3) = mkeypoints[good_matches[i].trainIdx].pt.y;
	}

	// RANSAC detect outliers
	auto_ptr< estimators::Solver<Eigen::MatrixXd, Eigen::VectorXd> > ptrSolver(
		new estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>);
	vector<int> inliers;
	//for (int i = 0; i < (int)good_matches.size(); i++) inliers.push_back(i);
	vector<Eigen::VectorXd> models;

	ransac::Ransac_Handler ransac_fun_Handler;
	bool result = ransac::Ransac_RobustEstimator
		(
		dataPoints, // the input data
		estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::extractor, // How select sampled point from indices
		dataPoints.rows(),  // the number of putatives data
		*(ptrSolver.get()),  // compute the underlying model given a sample set
		estimators::affineSolver<Eigen::MatrixXd, Eigen::VectorXd>::defaultEvaluator,  // the function to evaluate a given model
		//Ransac Object that contain function:
		// CandidatesSelector, Sampler and TerminationFunction
		ransac_fun_Handler, // the basic ransac object
		1000,  // the maximum rounds for RANSAC routine
		inliers, // inliers to the final solution
		models, // models array that fit input data
		0.95 // the confidence want to achieve at the end
		);
	if (inliers.size() < 6)
	{
		//handlerM->close();
		return match_state::match_failed;
	}

	//if (fabs(models[0][1] * models[0][5] - models[0][2] * models[0][4]) < 0.5)
	//{
	//	//handlerM->close();
	//	return match_state::match_failed;
	//}
	double delta_energy = 1.0 - fabs(models[0][1] * models[0][5] - models[0][2] * models[0][4]);
	if (fabs(delta_energy) > 0.5)
	{
		//handlerM->close();
		return match_state::match_failed;
	}

	//// maybe a lot of correspondences are found, but we need only one correspondence for a tile
	//float minId = inliers[0];
	//float minDist = good_matches[inliers[0]].distance;
	//for (int i = 1; i < (int)inliers.size(); ++i)
	//{
	//	if (minDist > good_matches[inliers[i]].distance)
	//	{
	//		minDist = good_matches[inliers[i]].distance;
	//		minId = inliers[i];
	//	}
	//}
	vector<ossimTDpt> tpts;
	for (size_t i = 0; i < inliers.size(); i++)
	{
		ossimDpt sc = ossimDpt(dataPoints(inliers[i], 0), dataPoints(inliers[i], 1));
		ossimDpt mc = ossimDpt(dataPoints(inliers[i], 2), dataPoints(inliers[i], 3));
		tpts.push_back(ossimTDpt(mc, sc, good_matches[inliers[i]].distance));
	}


	if (bDebug)
	{
		cv::Mat outMat;
		_prepareImgAndDrawKeylines(slaveMat,
			masterMat,
			tpts,
			outMat,
			cv::Scalar(0, 0, 255));
		cv::imwrite("matched_points.png", outMat);
	}

	// match lines
	// compute affine invariant
	int nSlave = (int)slave_vec_lines.size();
	int nMaster = (int)master_vec_lines.size();
	cv::Mat slave_affine_invariant(nSlave, 1, CV_32FC1);
	cv::Mat master_affine_invariant(nMaster, 1, CV_32FC1);
	for (int i = 0; i < nSlave; i++)
	{
		double v = affine_invariant(slave_vec_lines[i],
			tpts[0].getSlavePoint(), tpts[1].getSlavePoint());
		slave_affine_invariant.at<float>(i) = (float)v;
	}
	for (int i = 0; i < nMaster; i++)
	{ 
		double v = affine_invariant(master_vec_lines[i],
			tpts[0].getMasterPoint(), tpts[1].getMasterPoint());
		master_affine_invariant.at<float>(i) = (float)v;
	}

	//FlannBasedMatcher flannMatcher;
	////matches.clear();
	//std::vector<DMatch> allMatches;
	//good_matches.clear();
	//flannMatcher.match(slave_affine_invariant, master_affine_invariant, allMatches);
	//FILE *pf = fopen("ratio.txt", "w+");
	//for (size_t i = 0; i < allMatches.size(); i++)
	//{
	//	fprintf(pf, "%lf\n", allMatches[i].distance);
	//	if (allMatches[i].distance < 0.001)
	//	{
	//		good_matches.push_back(allMatches[i]);
	//	}
	//}
	//fclose(pf);
	////findGoodMatches(matches, good_matches, 0.70f);

	vector<cvline_polar> matched_slave_lines;
	vector<cvline_polar> matched_master_lines;

	//for (unsigned int i = 0; i < good_matches.size(); ++i)
	//{
	//	matched_slave_lines.push_back(slave_vec_lines[good_matches[i].queryIdx]);
	//	matched_master_lines.push_back(master_vec_lines[good_matches[i].trainIdx]);
	//}

	FILE *pf = fopen("ratio.txt", "w+");
	for (size_t i = 0; i < slave_vec_lines.size(); i++)
	{
		double maxSim = 0;
		int iMaster = 0;
		for (size_t j = 0; j < master_vec_lines.size(); j++)
		{
			double sim = affine_invariant_sim(slave_vec_lines[i], master_vec_lines[j], tpts);
			fprintf(pf, "%lf\n", sim);
			if (maxSim < sim)
			{
				maxSim = sim;
				iMaster = j;
			}
		}
		if (maxSim > 0.95)
		{
			matched_slave_lines.push_back(slave_vec_lines[i]);
			matched_master_lines.push_back(master_vec_lines[iMaster]);
		}
	}
	fclose(pf);

	int nMatched = (int)matched_slave_lines.size();
	if (nMatched < 4)
	{
		return match_state::match_failed;
	}


	vector<double> modelParameters = affine_estimation2(matched_slave_lines, matched_master_lines);

	//delta_energy = 1.0 - fabs(modelParameters[1] * modelParameters[5] - modelParameters[2] * modelParameters[4]);
	//if (fabs(delta_energy) > 0.5)
	//{
	//	return match_state::match_failed;
	//}

	// choose two lines
	int idx1 = 0;
	int idx2 = 1;
	double maxSine = 0.0;
	for (int i = 0; i < nMatched-1; i++)
	{
		for (int j = i + 1; j < nMatched; j++)
		{
			double sine = fabs(sin(matched_master_lines[i].theta - matched_master_lines[j].theta));
			if (maxSine < sine)
			{
				idx1 = i;
				idx2 = j;
				maxSine = sine;
			}
		}
	}

	theTies.clear();
	ossimDpt sc = ossimDpt(matched_slave_lines[idx1].pt1.x, matched_slave_lines[idx1].pt1.y);
	ossimDpt mc = ossimDpt(matched_master_lines[idx1].pt1.x, matched_master_lines[idx1].pt1.y);
	theTies.push_back(ossimTDpt(mc, sc, 0.0));

	sc = ossimDpt(matched_slave_lines[idx1].pt2.x, matched_slave_lines[idx1].pt2.y);
	mc = ossimDpt(matched_master_lines[idx1].pt2.x, matched_master_lines[idx1].pt2.y);
	theTies.push_back(ossimTDpt(mc, sc, 0.0));

	sc = ossimDpt(matched_slave_lines[idx2].pt1.x, matched_slave_lines[idx2].pt1.y);
	mc = ossimDpt(matched_master_lines[idx2].pt1.x, matched_master_lines[idx2].pt1.y);
	theTies.push_back(ossimTDpt(mc, sc, 0.0));

	sc = ossimDpt(matched_slave_lines[idx2].pt2.x, matched_slave_lines[idx2].pt2.y);
	mc = ossimDpt(matched_master_lines[idx2].pt2.x, matched_master_lines[idx2].pt2.y);
	theTies.push_back(ossimTDpt(mc, sc, 0.0));


	vector<cvline_polar> matched_slave_lines1;
	vector<cvline_polar> matched_master_lines1;
	matched_slave_lines1.push_back(matched_slave_lines[idx1]);
	matched_slave_lines1.push_back(matched_slave_lines[idx2]);
	matched_master_lines1.push_back(matched_master_lines[idx1]);
	matched_master_lines1.push_back(matched_master_lines[idx2]);

	if (bDebug)
	{
		cv::Mat s = cv::imread("D:\\workspace\\Landsat\\xinjiang\\s.tif");
		cv::Mat m = cv::imread("D:\\workspace\\Landsat\\xinjiang\\m.tif");
		int thickness = 2;
		cv::Mat img_outImage;
		//drawLineMatches(slaveMat, matched_slave_lines, masterMat, matched_master_lines,
		//	img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
		//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
		drawLineMatches(s, matched_slave_lines, m, matched_master_lines,
			img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
			vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
		cv::imwrite("matched_lines.png", img_outImage);



		cv::Mat img_outImage1;
		//drawLineMatches(slaveMat, matched_slave_lines, masterMat, matched_master_lines,
		//	img_outImage, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
		//	vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
		drawLineMatches(s, matched_slave_lines1, m, matched_master_lines1,
			img_outImage1, cv::Scalar(0, 0, 255), cv::Scalar::all(-1),
			vector<vector<char> >(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, thickness);
		cv::imwrite("matched_lines1.png", img_outImage1);
	}

	return match_state::success;
}

ossimIrect radiLineRegistration::getMasterRect(ossimProjection* sProjection, const GdalRasterApp& mGdalApp,
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

ossimIrect radiLineRegistration::getMasterRect(ossimProjection* sProjection, ossimProjection* mProjection,
												ossimIrect sRect, double slaveAccuracy)
{
	//ossimIpt sul = sRect.ul() + ossimIpt(-slaveAccuracy, -slaveAccuracy);
	//ossimIpt sll = sRect.ll() + ossimIpt(-slaveAccuracy, slaveAccuracy);

	ossimIpt sul = sRect.ul() + ossimIpt(-slaveAccuracy,-slaveAccuracy); 
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

bool radiLineRegistration::getGridFeaturesParallel(const ossimIrect& rect, radiBlockTieGptSet& tSet)
{

	//cout<<rect.ul()<<"\t"<<rect.lr()<<endl;
	if (!theSlaveBandSelector)
	{
		ossimNotify(ossimNotifyLevel_WARN)
			<< "WARN ossimTieGenerator::scanForEdges():"
			<< "\nInput source is not a ossimImageChip.  Returning..." << std::endl;
		return false;
	}

	const ossim_int32 TILE_HEIGHT	= theTileSize;
	const ossim_int32 TILE_WIDTH	= theTileSize;
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
	center_row = 0;
	center_col = 0;
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
		cerr<<"radiLineRegistration"<<"::execute can't open slave image  "<< getSlave().c_str() <<endl;
		return false;
	}
	// select only one band (if multiple)
	ossim_uint32 sbc = slaveApp.nBand();
	//add a band selector
	ossim_uint32 sb = theSlaveBand;
	if (sb>=sbc) 
	{
		cerr<<"radiLineRegistration"<<"::execute Warning not enough bands in slave, only "<< sbc <<endl;
		sb=0;
	}
		
	double sGSD = min(theSlaveProjection->getMetersPerPixel().x, theSlaveProjection->getMetersPerPixel().y);
	//int ncore = omp_get_num_procs();//获取执行核的总数；  目前机器CPU的数量
	row_col_step = 5;
	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);	

	bool found = false;
	int search_size = (int)row_col_List.size();
	double image_scale = 1.0/1.0;
	//search_size = min(search_size, 1);
	for (int i=0;i < search_size && !found;++i)
	{
		int icol = int(row_col_List[i+j].col_idx+center_col);
		int irow = int(row_col_List[i+j].row_idx+center_row);
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
		ossimIrect srect(ossimIpt(samp, line),ossimIpt(samp+BufWidth-1,line+BufHeight-1));
		//cout<<srect.ul()<<"\t"<<srect.lr()<<endl;
		//cout<<BufWidth<<"\t"<<BufHeight<<endl;
		if (!slaveApp.getRect2CvMatByte(srect, slaveMat, theSlaveBand, image_scale))
		{
			//cerr<<"radiLineRegistration"<<"::execute slave tile is not valid.  "<<endl;
			//cout<<srect.ul()<<"\t"<<srect.lr()<<endl;
			continue;
		}

		ossimIpt sul = srect.ul();
		ossimIpt slr = srect.lr();
		ossimIpt delta_lr((ossim_int32)(ceil(theSlaveAccuracy)), (ossim_int32)(ceil(theSlaveAccuracy)) );

		vector<ossimTDpt> tps;
		// search in the reference library
		ossimGpt ul_latlon, lr_latlon;
		theSlaveProjection->lineSampleToWorld(sul-delta_lr, ul_latlon);
		theSlaveProjection->lineSampleToWorld(slr+delta_lr, lr_latlon);
		vector<ossimFilename> masterFileList;
		if (theMaster.ext().upcase() == "SHP")
		{
			getMasterList(theMaster, masterFileList, ul_latlon, lr_latlon);
		}
		else{
			masterFileList.clear();
			masterFileList.push_back(theMaster);
		}

		cv::Mat masterMat;
		found = false;
		for (int iFile = 0;iFile < (int)masterFileList.size();++iFile)
		{
			ossimFilename lastMaster = masterFileList[iFile];
			theLastMaster = lastMaster;
			GdalRasterApp masterApp;
			ossimRefPtr<ossimImageHandler> mHandler = ossimImageHandlerRegistry::instance()->open(lastMaster);
			if (!masterApp.open(lastMaster.c_str()) || NULL == mHandler.get())
			{
				cerr<<"radiLineRegistration"<<"::execute can't open master image  "<< lastMaster <<endl;
				continue;
			}
			ossimRefPtr<ossimProjection> mProjection;
			// select only one band (if multiple)
			ossim_uint32 mbc = masterApp.nBand();
			//add a band selector
			ossim_uint32 mb = theMasterBand;
			if (mb>=mbc) 
			{
				//cerr<<"radiLineRegistration"<<"::execute Warning not enough bands in master, only "<< mbc <<endl;
				mb=0;
			}

			double gsd_scale;
			// master
			ossimIrect mrect;
			//if (0 == strcmp(masterApp.getGetProjectionRef(), ""))
			if(NULL == mHandler->getImageGeometry().get() || NULL == (mProjection = mHandler->getImageGeometry()->getProjection()).get())
			{
				// 无投影
				mrect = ossimIrect(sul - delta_lr, slr + delta_lr);
				gsd_scale = 1.0;
			}
			else
			{
				double mGSD = mProjection->getMetersPerPixel().x;
				//double mGSD = masterApp.getGeoTransform()[1];
				gsd_scale = mGSD / sGSD;
				//gsd_scale = 1.0;
				// master
				mrect = getMasterRect(theSlaveProjection, mProjection.get(),
					srect, theSlaveAccuracy);
			}
			mrect = mrect.clipToRect(masterApp.getBoundary());
			if (!masterApp.getRect2CvMatByte(mrect, masterMat, mb, gsd_scale * image_scale))
			{
				//cerr<<"radiLineRegistration"<<"::execute master tile is not valid.  "<<endl;
				masterApp.close();
				continue;
			}
			if (match_state::success == runMatchParallel(slaveMat, masterMat, tps, this, theDebug) && tps.size() == 4)
			{
				for (size_t i = 0; i < 4; i++)
				{
					tps[i].setSlavePoint(tps[i].getSlavePoint() + sul);
					tps[i].setMasterPoint(tps[i].getMasterPoint() / gsd_scale + mrect.ul());

				}
				if(NULL == mHandler->getImageGeometry().get()
					|| NULL == (mProjection = mHandler->getImageGeometry()->getProjection()).get()
					|| thePointType == point_type::tie)
				{
					LinePolarPair lp1, lp2;
					lp1.line = line2polar(tps[0].getSlavePoint(), tps[1].getSlavePoint());
					lp1.line_prime = line2polar(tps[0].getMasterPoint(), tps[1].getMasterPoint());
					lp2.line = line2polar(tps[2].getSlavePoint(), tps[3].getSlavePoint());
					lp2.line_prime = line2polar(tps[2].getMasterPoint(), tps[3].getMasterPoint());

					m_TieLines.push_back(lp1);
					m_TieLines.push_back(lp2);
				}
				else
				{
					for (size_t i = 0; i < 4; i++)
					{
						ossimGpt ll;
						mProjection->lineSampleToWorld(tps[i].getMasterPoint(), ll);
						tps[i].setMasterPoint(ossimDpt(ll.lon, ll.lat));
					}
					LinePolarPair lp1, lp2;
					lp1.line = line2polar(tps[0].getSlavePoint(), tps[1].getSlavePoint());
					lp1.line_prime = line2polar(tps[0].getMasterPoint(), tps[1].getMasterPoint());
					lp2.line = line2polar(tps[2].getSlavePoint(), tps[3].getSlavePoint());
					lp2.line_prime = line2polar(tps[2].getMasterPoint(), tps[3].getMasterPoint());
					m_TieLines.push_back(lp1);
					m_TieLines.push_back(lp2);
				}
				masterApp.close();
				found = true;
				break;
			}
			else
			{
				//cerr<<"this tile does not matched."<<endl;
				masterApp.close();
				continue;
			}
		}
	}
	slaveApp.close();
	return found;
}

//bool radiLineRegistration::getGridFeaturesParallel(const ossimIrect& rect, radiBlockTieGpt& tSet)
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
//		cerr<<"radiLineRegistration"<<"::execute can't open slave image  "<< getSlave().c_str() <<endl;
//		return false;
//	}
//	// select only one band (if multiple)
//	ossim_uint32 sbc = slaveApp.nBand();
//
//	//add a band selector
//	ossim_uint32 sb = theSlaveBand;
//	if (sb>=sbc) 
//	{
//		cerr<<"radiLineRegistration"<<"::execute Warning not enough bands in slave, only "<< sbc <<endl;
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
//					cerr<<"radiLineRegistration"<<"::execute can't open master image  "<< lastMaster <<endl;
//					continue;
//				}
//				// select only one band (if multiple)
//				ossim_uint32 mbc = masterApp.nBand();
//
//				//add a band selector
//				ossim_uint32 mb = theMasterBand;
//				if (mb>=mbc) 
//				{
//					//cerr<<"radiLineRegistration"<<"::execute Warning not enough bands in master, only "<< mbc <<endl;
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

//bool radiLineRegistration::getGridFeaturesParallel(const ossimIrect& rect)
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
//			cerr<<"radiLineRegistration"<<"::execute can't open master image  "<< lastMaster <<endl;
//			continue;
//		}
//		// select only one band (if multiple)
//		ossim_uint32 mbc = masterApp.nBand();
//
//		//add a band selector
//		ossim_uint32 mb = theMasterBand;
//		if (mb>=mbc) 
//		{
//			cerr<<"radiLineRegistration"<<"::execute Warning not enough bands in master, only "<< mbc <<endl;
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

//bool radiLineRegistration::getGridFeaturesParallel(const ossimIrect& rect, void *pData)
//{
//	radiLineRegistration* pThis = (radiLineRegistration*)pData;
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
//				cerr<<"radiLineRegistration"<<"::execute can't create handler for master image  "<< theLastMaster <<endl;
//				continue;
//			}
//			// select only one band (if multiple)
//			ossim_uint32 mbc = master_handler->getNumberOfOutputBands();
//
//			//add a band selector
//			ossim_uint32 mb = pThis->theMasterBand;
//			if (mb>=mbc) 
//			{
//				cerr<<"radiLineRegistration"<<"::execute Warning not enough bands in master, only "<< mbc <<endl;
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

bool radiLineRegistration::getAllFeatures()
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
	const ossim_int32 TILE_HEIGHT    = ceil(nHeight / (double)(tilerows-0.5));
	const ossim_int32 TILE_WIDTH     = ceil(nWidth / (double)(tilecols-0.5));
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

void radiLineRegistration::writeTiePoints(const radiBlockTieGptSet& tp)
{
	ossimMapProjection* outProjection = radiLineRegistration::getOutputProjection();
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

void radiLineRegistration::writeTieLines(ossimFilename outPath, ossimMapProjection* pMapProjection/* = NULL*/)
{
	ossimMapProjection* outProjection = radiLineRegistration::getOutputProjection();
	ossimFilename slaveLineFile = outPath + "\\source.txt";
	ossimFilename masterLineFile = outPath + "\\reference.txt";
	FILE *pf_slave = fopen(slaveLineFile.c_str(), "w+");
	FILE *pf_master = fopen(masterLineFile.c_str(), "w+");
	for (size_t i = 0; i < m_TieLines.size(); i++)
	{
		if (i != 0)
		{
			fprintf(pf_slave, "\n");
			fprintf(pf_master, "\n");
		}

		fprintf(pf_slave, "2%03d\t%s\t%d", i + 1, "StraightLine", 2);
		fprintf(pf_master, "2%03d\t%s\t%d", i + 1, "StraightLine", 2);
		
		ossimDpt dpt = ossimDpt(m_TieLines[i].line.pt1.x, m_TieLines[i].line.pt1.y);
		ossimGpt gpt = ossimGpt(m_TieLines[i].line_prime.pt1.y, m_TieLines[i].line_prime.pt1.x);
		if (NULL!=pMapProjection && !pMapProjection->isGeographic())
		{
			ossimDpt d1 = outProjection->forward(gpt);
			gpt = ossimGpt(d1.x, d1.y);
		}
		else
		{
			gpt = ossimGpt(gpt.lon, gpt.lat);
		}
		fprintf(pf_slave, "\n%20.8lf\t%20.8lf", dpt.x, dpt.y);
		fprintf(pf_master, "\n%20.8lf\t%20.8lf", gpt.lat, gpt.lon);

		dpt = ossimDpt(m_TieLines[i].line.pt2.x, m_TieLines[i].line.pt2.y);
		gpt = ossimGpt(m_TieLines[i].line_prime.pt2.y, m_TieLines[i].line_prime.pt2.x);
		if (NULL != pMapProjection && !pMapProjection->isGeographic())
		{
			ossimDpt d1 = outProjection->forward(gpt);
			gpt = ossimGpt(d1.x, d1.y);
		}
		else
		{
			gpt = ossimGpt(gpt.lon, gpt.lat);
		}
		fprintf(pf_slave, "\n%20.8lf\t%20.8lf", dpt.x, dpt.y);
		fprintf(pf_master, "\n%20.8lf\t%20.8lf", gpt.lat, gpt.lon);

	}
	fclose(pf_slave);
	fclose(pf_master);
}

void radiLineRegistration::setOutputName(const ossimString& filename)
{
	ossimOutputSource::setOutputName(filename);

	if (isOpen()) close();

	if (filename != "")
	{
		theFilename = filename;
	}
}

void radiLineRegistration::setAreaOfInterest(const ossimIrect& rect)
{
	theAreaOfInterest = rect;
}

bool radiLineRegistration::isOpen()const
{
#if OSSIM_HAS_MPI
	return (theFileStream != NULL);
#else
	return const_cast<fstream*>(&theFileStream)->is_open();
#endif
}

bool radiLineRegistration::open()
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

void radiLineRegistration::close()
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
