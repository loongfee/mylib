#pragma warning(push)
#pragma warning(disable : 4482)

#include "RpcModel.h"

#define OVERLAP 0

RpcModel::RpcModel(void)
{
	m_FileOutType = "tiff";
	//m_FileOutType = "pix";
	m_SampleType = "BICUBIC";
	m_ctrSet = new ossimTieGptSet;
	m_chkSet = new ossimTieGptSet;

	m_bUseGcps = true;
	m_bReplace = true;
	m_bInitState = false;

	m_progress = new OrthProcessListener;
	m_TileCount = 1;
	m_CurrentTile = 0;
}

bool RpcModel::setProjection()
{
	m_MapProjection.clear();
	ossimString tmp;
	//double a,b;
	char* prefix="";
	//中央纬度
	double centerLat = 0.0;	
	//半球
	char theHemisphere = m_UTMZone[m_UTMZone.size() - 1];
	//投影带
	ossim_uint32 theZone = m_UTMZone.beforePos(m_UTMZone.size() - 1).toInt();
	//中央经线
	double centerLon = (double)((theZone - 30) * 6 - 3);

	m_MapProjection.add(prefix,ossimKeywordNames::ORIGIN_LATITUDE_KW, centerLat,true);
	m_MapProjection.add(prefix,ossimKeywordNames::CENTRAL_MERIDIAN_KW,centerLon,true);

	//使用椭球模型
	m_EllipsoidModel = "WE";
	const 	ossimEllipsoid * theelli=ossimEllipsoidFactory::instance()->create(m_EllipsoidModel.c_str());
	theelli->saveState(m_MapProjection, prefix);

	ossimDpt	theMetersPerPixel;
	//X方向分辨率
	theMetersPerPixel.lon = m_PixelSize;
	//Y方向分辨率
	theMetersPerPixel.lat = m_PixelSize;

	//设置分辨率 单位为米
	m_MapProjection.add(prefix,	ossimKeywordNames::PIXEL_SCALE_XY_KW,theMetersPerPixel.toString().c_str(),true);
	m_MapProjection.add(prefix,	ossimKeywordNames::PIXEL_SCALE_UNITS_KW,	ossimUnitTypeLut::instance()->getEntryString(OSSIM_METERS),true); 

	ossimDpt	theFalseEastingNorthing;
	//东向偏移值
	theFalseEastingNorthing.lon = 500000.0;
	//北向偏移值
	theFalseEastingNorthing.lat = 0.0;

	m_MapProjection.add(prefix,ossimKeywordNames::FALSE_EASTING_NORTHING_KW,theFalseEastingNorthing.toString().c_str(), true);
	m_MapProjection.add(prefix,	ossimKeywordNames::FALSE_EASTING_NORTHING_UNITS_KW,	ossimUnitTypeLut::instance()->getEntryString(OSSIM_METERS), true);
	bool theElevationLookupFlag=true;
	m_MapProjection.add(prefix,	ossimKeywordNames::ELEVATION_LOOKUP_FLAG_KW,	ossimString::toString(theElevationLookupFlag), true);

	//设置UTM投影
	m_MapProjection.add(prefix,ossimKeywordNames::ZONE_KW,theZone,true);
	m_MapProjection.add(prefix,ossimKeywordNames::HEMISPHERE_KW,theHemisphere,true);
	m_MapProjection.add(prefix,	"type","ossimUtmProjection",true);



	m_MapPar=PTR_CAST(ossimMapProjection,
		ossimMapProjectionFactory::instance()->createProjection(m_MapProjection));
	if(m_MapPar==NULL)
	{
		cout<<"您的投影信息设置有误！"<<endl;
		return false;
	}
	return true;
}

void RpcModel::ReadGcpAndProjection(ossimFilename strFilename,ossimTieGptSet* &m_gptset,ossimTieGptSet* &mcheck_gptset)
{
	ossimGpt tieg;
	ossimDpt tied;
	ossimTieGpt *aTiePt;
	ossimString index;
	m_gptset->clearTiePoints();
	mcheck_gptset->clearTiePoints();


	ossimTempFilename pp;
	bool flag=true;
	ossimString str;

	char strtmp[255];
	std::vector<ossimString> gcpcon;
	const char* ellipse_code;
	pp.generateRandomFile(); 

	fstream os;
	os.open(strFilename.c_str(), ios_base::in);
	fstream oss;
	oss.open(pp.c_str(),ios::out|ios_base::app);
	os>>str;
	while(os.getline(strtmp,255) ) {
		str=strtmp;
		oss<< str.c_str();oss<<endl;
		if(str.contains("MAP_PROJECTION_END")) {flag=false;break;}

	}
	oss.close();

	if (!flag) {////////判断控制点文件中是否有投影信息
		m_MapProjection.clear();
		m_MapPar=NULL;
		m_MapProjection.addFile(pp.c_str());
		m_MapPar=PTR_CAST(ossimMapProjection,
			ossimMapProjectionFactory::instance()->createProjection(m_MapProjection));
		ellipse_code = m_MapProjection.find(ossimKeywordNames::DATUM_KW);

	}
	else
	{
		if(m_MapPar ==NULL)  {
			return;
		}

		ellipse_code=m_MapProjection.find(ossimKeywordNames::DATUM_KW);
		os.clear();
		//os.seekg(0,ios::beg);
		os.close();
		os.open(strFilename.c_str());

	}
	mcheck_gptset->clearTiePoints();
	m_gptset->clearTiePoints();
	while(os.getline(strtmp,255))
	{
		int star;
		str=strtmp;
		if (str.empty()) continue;
		gcpcon=str.split("	", true);
		if (gcpcon.size()==1) {
			str=gcpcon[0];
			gcpcon=str.split(" ", true);
		}
		if (gcpcon.size()<5) continue;
		aTiePt=new ossimTieGpt(tieg,tied,0);
		if(gcpcon[0].empty()) {
			aTiePt->GcpNumberID = gcpcon[1];
			star=1;
		}
		else
		{
			aTiePt->GcpNumberID = gcpcon[0];
			star=0;
		}

		tied.x=gcpcon[star+1].toDouble();tied.y=gcpcon[star+2].toDouble();
		tieg.lat=gcpcon[star+3].toDouble();tieg.lon=gcpcon[star+4].toDouble();
		tieg.hgt=gcpcon[star+5].toDouble();
		tieg.datum(ossimDatumFactory::instance()->create(ossimString(ellipse_code))); 

		aTiePt->setGroundPoint(tieg);
		aTiePt->setImagePoint(tied);
		index=aTiePt->GcpNumberID.c_str();
		if (index.substr(0,1).compare("-")) m_gptset->addTiePoint(aTiePt);
		else
		{
			aTiePt->GcpNumberID=index.after("-").c_str();
			mcheck_gptset->addTiePoint(aTiePt);
		}

	}
	os.close();
}


bool RpcModel::GetElevations(ossimTieGptSet* &ctrlSet, double default_hgt/* = 0.0*/)
{
	vector<ossimRefPtr<ossimTieGpt> >&    theGcp =  ctrlSet->refTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::iterator iter,tit;

	ossimMapProjection* theMapProjection = m_MapPar;

	if(m_theElevManager == NULL)
	{
		cout<<"高程源还未选取！"<<endl;
		return false;
	}
	if(theMapProjection==NULL)
	{
		cout<<"您投影设置错误或没有设置投影信息！"<<endl;
		return false;
	}


	ossimGpt tGpt;
	for (iter = theGcp.begin() ; iter != theGcp.end() ; ++iter)
	{
		tGpt = theMapProjection->inverse(ossimDpt((*iter)->lat,(*iter)->lon));
		(*iter)->hgt = m_theElevManager->getHeightAboveEllipsoid(tGpt);
		if(ossim::isnan((*iter)->hgt)) (*iter)->hgt = 0.0;
		if((*iter)->hgt == 0.0)
		{
			theGcp.erase(iter);
			iter--;
		}
	}

	return true;
}


void RpcModel::UpdateSensorModel(ossimTieGptSet tieGptSet,
								 ossimSensorModel* &sensorModel,
								 ossimKeywordlist& geom)
{
	ossimDpt imagepoint,cimagepoint;
	ossimGpt goundpoint,tGpt;
	int i;

	int num = static_cast<int>(tieGptSet.size());

	vector<ossimRefPtr<ossimTieGpt> >& theTPV = tieGptSet.refTiePoints();

	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;

	for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
	{
		imagepoint=(*tit)->getImagePoint();
		goundpoint=(*tit)->getGroundPoint();
		tGpt = m_MapPar->inverse(ossimDpt(goundpoint.lat, goundpoint.lon));
		if(!(*tit)->isHgtNan()) tGpt.hgt = (*tit)->hgt;
		(*tit)->setGroundPoint(tGpt);
	}

	sensorModel->optimizeFit(tieGptSet);
	sensorModel->updateModel();
	sensorModel->saveState(geom);

	for(i = 0;i < static_cast<int>(tieGptSet.getTiePoints().size());i++)
	{
		ossimDpt dpt = sensorModel->m_proj->forward(*tieGptSet.getTiePoints()[i]);
		ossimGpt gpt(dpt.x,dpt.y);
		tieGptSet.refTiePoints()[i]->setGroundPoint(ossimGpt(dpt.x,dpt.y,tieGptSet.getTiePoints()[i]->hgt));
	}
}

bool RpcModel::Orthograph(ossimFilename outfile)
{
	ossimInit::instance()->initialize();
	m_TileCount = 1;
	m_CurrentTile = 0;
	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(m_imageFileName);
	if(!handler) return false;   //应该弹出警告对话框
	ossimImageGeometry geom;
	geom.loadState(m_geom);
	handler->setImageGeometry(&geom);
	ossimKeywordlist tt_geom;

	ossimImageFileWriter* writer;

	if (m_FileOutType.contains("tiff")) 
		writer= ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
	if (m_FileOutType.contains("pix")) 
		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_pcidsk"));
	if (m_FileOutType.contains("img")) 
		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_imagine_hfa"));

	if(writer==NULL) return false;

	ossimDpt imagesize(handler->getImageRectangle().width(), handler->getImageRectangle().height());
	ossimImageRenderer* renderer = new ossimImageRenderer;
	ossimPolyCutter* theCutter;
	ossimBandSelector* theBandSelector;
	vector<ossimDpt> polygon;
	ossimIrect bound,boundw;
	theCutter = new ossimPolyCutter;
	///////////////////以下四行应该从界面取得
	int starline,starpixel,endpixel,endline;
	starline	=	0;
	starpixel	=	0;
	endpixel	=	0;
	endline		=	0;
	if (0 == endpixel)
	{
		endpixel = imagesize.x - 1;
	}
	if (0 == endline)
	{
		endline = imagesize.y - 1;
	}/////////////////////////////////////////////////////////////////////

	theBandSelector = new ossimBandSelector;
	theBandSelector->connectMyInputTo(0, handler);
	theBandSelector->setOutputBandList(m_OutBandList);


	ossimDpt ps(starpixel,starline),p2(endpixel,starline),p3(endpixel,endline),p4(starpixel,endline),p5(starpixel,starline);

	polygon.push_back(ps);
	polygon.push_back(p2);
	polygon.push_back(p3);
	polygon.push_back(p4);
	polygon.push_back(p5);
	theCutter->connectMyInputTo(theBandSelector);
	theCutter->setPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);


	tt_geom.clear();
	writer->saveState(tt_geom);

	tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA",true);

	if(m_SampleType.contains("NEAREST_NEIGHBOR")) renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_NEAREST_NEIGHBOR);
	if(m_SampleType.contains("BILINEAR"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	if(m_SampleType.contains("BICUBIC"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);


	writer->loadState(tt_geom);

	renderer->connectMyInputTo(theCutter);
	renderer->setView(m_MapPar);

	writer->setFilename(outfile);
	writer->addListener((ossimProcessListener*)m_progress);
	writer->connectMyInputTo(0,renderer);

	//bool bResult = true;
	bool bResult = writer->execute();
	writer->disconnectAllInputs();
	renderer->disconnectAllInputs();
	theCutter->disconnectAllInputs();

	//handler->close();
	return bResult;
}

bool RpcModel::Orthograph_Tile(ossimFilename outfile, int tileSize)
{
	ossimInit::instance()->initialize();
	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(m_imageFileName);
	if(!handler) return false;   //应该弹出警告对话框
	ossimImageGeometry geom;
	geom.loadState(m_geom);
	handler->setImageGeometry(&geom);
	ossimKeywordlist tt_geom;

	ossimImageFileWriter* writer;

	if (m_FileOutType.contains("tiff")) 
		writer= ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
	if (m_FileOutType.contains("pix")) 
		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_pcidsk"));
	if (m_FileOutType.contains("img")) 
		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_imagine_hfa"));

	if(writer==NULL) return false;

	ossimDpt imagesize(handler->getImageRectangle().width(), handler->getImageRectangle().height());
	ossimImageRenderer* renderer = new ossimImageRenderer;
	ossimPolyCutter* theCutter;
	ossimBandSelector* theBandSelector;
	vector<ossimDpt> polygon;
	theCutter = new ossimPolyCutter;
	int starline,starpixel,endpixel,endline;
	starline	=	0;
	starpixel	=	0;
	endpixel	=	0;
	endline		=	0;
	if (0 == endpixel)
	{
		endpixel = imagesize.x - 1;
	}
	if (0 == endline)
	{
		endline = imagesize.y - 1;
	}/////////////////////////////////////////////////////////////////////

	theBandSelector = new ossimBandSelector;
	theBandSelector->connectMyInputTo(0, handler);
	theBandSelector->setOutputBandList(m_OutBandList);

	ossimDpt ps(starpixel,starline),p2(endpixel,starline),p3(endpixel,endline),p4(starpixel,endline),p5(starpixel,starline);
	polygon.push_back(ps);
	polygon.push_back(p2);
	polygon.push_back(p3);
	polygon.push_back(p4);
	polygon.push_back(p5);
	theCutter->connectMyInputTo(theBandSelector);
	theCutter->setPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	ossimIrect bound = theCutter->getBoundingRect();

	renderer->connectMyInputTo(theCutter);
	renderer->setView(m_MapPar);


	tt_geom.clear();
	writer->saveState(tt_geom);
	tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA",true);

	if(m_SampleType.contains("NEAREST_NEIGHBOR")) renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_NEAREST_NEIGHBOR);
	if(m_SampleType.contains("BILINEAR"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	if(m_SampleType.contains("BICUBIC"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);


	ossimStdOutProgress progress(0,true);
	writer->addListener((ossimProcessListener*)m_progress);
	writer->connectMyInputTo(0,renderer);


	ossimImageViewTransform*  IVT = renderer->getImageViewTransform();
	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, IVT);
	ossimDrect tempRect = IVT->getImageToViewBounds(bound);

	ossimIrect outrect(tempRect),blockrect;
	int row,col;
	ossimFilename outname;
	row=(int)(static_cast<int>(outrect.height()) / tileSize);
	col=(int)(static_cast<int>(outrect.width()) / tileSize);
	if(row * tileSize < static_cast<int>(outrect.height()) ) row =row +1;
	if(col * tileSize < static_cast<int>(outrect.width()) ) col =col +1;

	ossimString m_tilfile;
	m_tilfile = outfile.setExtension("TIL");
	fstream fout;
	fout.open(m_tilfile.c_str(),ios::out);
	ossimFilename file;
	/////////////////////////////////////////////////////////////////////////////
	fout<<"bandId = \"RGB\""<<endl;
	fout<<"numTiles = "<<col*row<<";"<<endl;
	fout<<"tileSizeX = "<<col<<";"<<endl;
	fout<<"tileSizeY = "<<row<<";"<<endl;
	fout<<"tileUnits = \"Pixels\";"<<endl;
	fout<<"tileOverlap = 0;"<<endl;

	int starx = outrect.ul().x - OVERLAP;
	int stary = outrect.ul().y - OVERLAP;
	m_TileCount = row * col;
	m_TileCount = (m_TileCount < 1)?1:m_TileCount;
	m_CurrentTile = 0;
	for(int i=1;i<=row;i++)
	{
		for(int j=1;j<=col;j++)
		{
			//m_CurrentTile = (i - 1) * col + j - 1;
			blockrect.set_ulx(outrect.ul().x+(j-1) * tileSize-OVERLAP);blockrect.set_uly(outrect.ul().y+(i-1)*tileSize-OVERLAP);
			blockrect.set_lrx(blockrect.ul().x + tileSize+OVERLAP-1);blockrect.set_lry(blockrect.ul().y+tileSize+OVERLAP-1);

			if(blockrect.lr().x > outrect.lr().x) blockrect.set_lrx(outrect.lr().x);
			if(blockrect.lr().y > outrect.lr().y) blockrect.set_lry(outrect.lr().y);

			//tmpname = outfile.noExtension() + "_R" + ossimString::toString(i) + "C" + ossimString::toString(j) + "_tmp.TIF";
			outname = outfile.noExtension() + "_R" + ossimString::toString(i) + "C" + ossimString::toString(j) + ".TIF";

			//		str="BEGIN_GROUP = TILE_"+str_filenum+"\n"+"filename = "+str_line+"_"+str_pixel+".tif"+"\n"+"ULColOffset = "+str_ULCO+";\n"+"ULRowOffset = "+str_ULRO+";\n"+"LRColOffset = "+str_LRCO+";\n"+"LRRowOffset = "+str_LRRO+";\n"+"END_GROUP = TILE_"+str_filenum+"\n";
			fout<<"BEGIN_GROUP = TILE_"<<col*(i-1)+j<<endl;
			fout<<"filename = "<<outname.file().c_str()<<endl;
			fout<<"ULColOffset = "<<blockrect.ul().x-starx<<";"<<endl;
			fout<<"ULRowOffset = "<<blockrect.ul().y-stary<<";"<<endl;
			fout<<"LRColOffset = "<<blockrect.lr().x-starx<<";"<<endl;
			fout<<"LRRowOffset = "<<blockrect.lr().y-stary<<";"<<endl;
			fout<<"END_GROUP = TILE_"<<col*(i-1)+j<<endl;

			if (outname.exists())
			{
				m_CurrentTile++;
				continue;
			}
			writer->setAreaOfInterest(blockrect);
			writer->setFilename(outname);
			if(!writer->execute())
			{
				m_CurrentTile++;
				continue;
			}
			m_CurrentTile++;
		}
	}
	fout.close();
	handler->close();

	theCutter->clear();
	theCutter->disableSource();
	writer->disableListener();
	writer->removeListener(m_progress);

	return true;
}

bool RpcModel::Clone_Tile(ossimFilename outfile, int tileSize)
{
	ossimInit::instance()->initialize();
	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(m_imageFileName);
	if(!handler) return false;   //应该弹出警告对话框
	ossimImageGeometry geom;
	ossimKeywordlist tt_geom;

	ossimImageFileWriter* writer;

	if (m_FileOutType.contains("tiff")) 
		writer= ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
	if (m_FileOutType.contains("pix")) 
		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_pcidsk"));
	if (m_FileOutType.contains("img")) 
		writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("gdal_imagine_hfa"));

	if(writer==NULL) return false;

	ossimDpt imagesize(handler->getImageRectangle().width(), handler->getImageRectangle().height());
	ossimImageRenderer* renderer = new ossimImageRenderer;
	ossimPolyCutter* theCutter;
	ossimBandSelector* theBandSelector;
	vector<ossimDpt> polygon;
	theCutter = new ossimPolyCutter;
	///////////////////以下四行应该从界面取得
	int starline,starpixel,endpixel,endline;
	starline	=	0;
	starpixel	=	0;
	endpixel	=	0;
	endline		=	0;
	if (0 == endpixel)
	{
		endpixel = imagesize.x - 1;
	}
	if (0 == endline)
	{
		endline = imagesize.y - 1;
	}/////////////////////////////////////////////////////////////////////

	theBandSelector = new ossimBandSelector;
	theBandSelector->connectMyInputTo(0, handler);
	m_OutBandList.clear();
	m_OutBandList.push_back(0);
	theBandSelector->setOutputBandList(m_OutBandList);

	ossimDpt ps(starpixel,starline),p2(endpixel,starline),p3(endpixel,endline),p4(starpixel,endline),p5(starpixel,starline);
	polygon.push_back(ps);
	polygon.push_back(p2);
	polygon.push_back(p3);
	polygon.push_back(p4);
	polygon.push_back(p5);
	theCutter->connectMyInputTo(theBandSelector);
	theCutter->setPolygon(polygon);
	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
	theCutter->setNumberOfPolygons(1);

	ossimIrect bound = theCutter->getBoundingRect();

	renderer->connectMyInputTo(theCutter);////////    handler  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	renderer->setView(m_MapPar);


	tt_geom.clear();
	writer->saveState(tt_geom);
	tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA",true);

	if(m_SampleType.contains("NEAREST_NEIGHBOR")) renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_NEAREST_NEIGHBOR);
	if(m_SampleType.contains("BILINEAR"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
	if(m_SampleType.contains("BICUBIC"))  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_CUBIC);


	ossimStdOutProgress progress(0,true);
	writer->addListener((ossimProcessListener*)m_progress);
	writer->connectMyInputTo(0,renderer);


	ossimImageViewTransform*  IVT = renderer->getImageViewTransform();
	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, IVT);
	ossimDrect tempRect = IVT->getImageToViewBounds(bound);

	ossimIrect outrect(tempRect),blockrect;
	int row,col;
	ossimFilename outname;
	row=(int)(outrect.height() / tileSize);
	col=(int)(outrect.width() / tileSize);
	if(row * tileSize < static_cast<int>(outrect.height()) ) row =row +1;
	if(col * tileSize < static_cast<int>(outrect.width()) ) col =col +1;

	ossimString m_tilfile;
	m_tilfile = outfile.setExtension("TIL");
	fstream fout;
	fout.open(m_tilfile.c_str(),ios::out);
	ossimFilename file;
	/////////////////////////////////////////////////////////////////////////////
	fout<<"bandId = \"RGB\""<<endl;
	fout<<"numTiles = "<<col*row<<";"<<endl;
	fout<<"tileSizeX = "<<col<<";"<<endl;
	fout<<"tileSizeY = "<<row<<";"<<endl;
	fout<<"tileUnits = \"Pixels\";"<<endl;
	fout<<"tileOverlap = 0;"<<endl;

	int starx = outrect.ul().x - OVERLAP;
	int stary = outrect.ul().y - OVERLAP;
	m_TileCount = row * col;
	m_TileCount = (m_TileCount < 1)?1:m_TileCount;
	for(int i=1;i<=row;i++)
	{
		for(int j=1;j<=col;j++)
		{
			m_CurrentTile = (i - 1) * col + j - 1;
			blockrect.set_ulx(outrect.ul().x+(j-1) * tileSize-OVERLAP);blockrect.set_uly(outrect.ul().y+(i-1)*tileSize-OVERLAP);
			blockrect.set_lrx(blockrect.ul().x + tileSize+OVERLAP-1);blockrect.set_lry(blockrect.ul().y+tileSize+OVERLAP-1);

			if(blockrect.lr().x > outrect.lr().x) blockrect.set_lrx(outrect.lr().x);
			if(blockrect.lr().y > outrect.lr().y) blockrect.set_lry(outrect.lr().y);

			//tmpname = outfile.noExtension() + "_R" + ossimString::toString(i) + "C" + ossimString::toString(j) + "_tmp.TIF";
			outname = outfile.noExtension() + "_R" + ossimString::toString(i) + "C" + ossimString::toString(j) + ".TIF";

			//		str="BEGIN_GROUP = TILE_"+str_filenum+"\n"+"filename = "+str_line+"_"+str_pixel+".tif"+"\n"+"ULColOffset = "+str_ULCO+";\n"+"ULRowOffset = "+str_ULRO+";\n"+"LRColOffset = "+str_LRCO+";\n"+"LRRowOffset = "+str_LRRO+";\n"+"END_GROUP = TILE_"+str_filenum+"\n";
			fout<<"BEGIN_GROUP = TILE_"<<col*(i-1)+j<<endl;
			fout<<"filename = "<<outname.file().c_str()<<endl;
			fout<<"ULColOffset = "<<blockrect.ul().x-starx<<";"<<endl;
			fout<<"ULRowOffset = "<<blockrect.ul().y-stary<<";"<<endl;
			fout<<"LRColOffset = "<<blockrect.lr().x-starx<<";"<<endl;
			fout<<"LRRowOffset = "<<blockrect.lr().y-stary<<";"<<endl;
			fout<<"END_GROUP = TILE_"<<col*(i-1)+j<<endl;

			if (outname.exists()) continue;
			writer->setAreaOfInterest(blockrect);
			writer->setFilename(outname);
			if(!writer->execute()) return false;
		}
	}
	fout.close();
	handler->close();

	theCutter->clear();
	theCutter->disableSource();
	writer->disableListener();
	writer->removeListener(&progress);

	return true;
}

int RpcModel::GetProgress()
{
	return floor((m_CurrentTile * 100 + m_progress->getPercentComplete()) / (double)m_TileCount + 0.5);
}

bool RpcModel::OutputReport(ossimFilename reportfile, ossimSensorModel* sensorModel, ossimTieGptSet* ctrlSet, ossimTieGptSet* chkSet)
{
	int i;
	fstream fs;
	fs.open(reportfile.c_str(),ios_base::out);
	fs.setf(ios::fixed, ios::floatfield);
	fs.precision(2);

	fs<<"残差报告\n\n"<<"残差单位：米\n";

	NEWMAT::ColumnVector residue1 = CalcResidue(sensorModel, *ctrlSet);
	NEWMAT::ColumnVector residue2 = CalcResidue(sensorModel, *chkSet);

	fs<<"控制点点数："<<ctrlSet->size();
	if(static_cast<int>(ctrlSet->size()) > 0)
	{
		double ResidueX = 0.0;
		double ResidueY = 0.0;

		for(i = 0;i < static_cast<int>(ctrlSet->size());i++)
		{
			ResidueX += residue1.element(i * 3 + 0) * residue1.element(i * 3 + 0);
			ResidueY += residue1.element(i * 3 + 1) * residue1.element(i * 3 + 1);
		}
		ResidueX = sqrt(ResidueX / static_cast<int>(ctrlSet->size()));
		ResidueY = sqrt(ResidueY / static_cast<int>(ctrlSet->size()));
		fs<<"\t\t"<<"均方差：X "<<ResidueX<<"\tY "<<ResidueY<<endl;
	}
	else
		fs<<endl;

	fs<<"检查点点数："<<chkSet->size();
	if(static_cast<int>(chkSet->size()) > 0)
	{
		double ResidueX = 0.0;
		double ResidueY = 0.0;

		for(i = 0;i < static_cast<int>(chkSet->size());i++)
		{
			ResidueX += residue2.element(i * 3 + 0) * residue2.element(i * 3 + 0);
			ResidueY += residue2.element(i * 3 + 1) * residue2.element(i * 3 + 1);
		}
		ResidueX = sqrt(ResidueX / static_cast<int>(chkSet->size()));
		ResidueY = sqrt(ResidueY / static_cast<int>(chkSet->size()));
		fs<<"\t\t"<<"均方差：X "<<ResidueX<<"\tY "<<ResidueY<<endl;
	}
	else
		fs<<endl;

	fs<<"清单："<<endl;

	fs<<"标示\t"<<"误差\t"<<"残差X\t"<<"残差Y\t"<<"类型\t"<<"幅号\t"<<"实际X\t"<<"实际Y\t"<<"对照X\t"<<"对照Y\t"<<endl;

	fs<<"用于计算模型的控制点点数："<<ctrlSet->size()<<endl;

	//residue = CalcResidue(*m_CtrTieGptSet);
	for(i = 0;i < static_cast<int>(ctrlSet->size());i++)
	{
		fs<<ctrlSet->getTiePoints()[i]->GcpNumberID<<"\t"
			<<sqrt(residue1.element(i * 3 + 0)*residue1.element(i * 3 + 0) + residue1.element(i * 3 + 1)*residue1.element(i * 3 + 1))<<"\t"
			<<residue1.element(i * 3 + 0)<<"\t"
			<<residue1.element(i * 3 + 1)<<"\t"
			<<"GCP"<<"\t"
			<<""<<"\t"
			<<ctrlSet->getTiePoints()[i]->getGroundPoint().lat<<"\t"
			<<ctrlSet->getTiePoints()[i]->getGroundPoint().lon<<"\t"
			<<ctrlSet->getTiePoints()[i]->getGroundPoint().lat + residue1.element(i * 3 + 0)<<"\t"
			<<ctrlSet->getTiePoints()[i]->getGroundPoint().lon + residue1.element(i * 3 + 1)<<"\n";
	}

	//fs<<endl;

	//residue = CalcResidue(*m_ChkTieGptSet);
	for(i = 0;i < static_cast<int>(chkSet->size());i++)
	{
		fs<<chkSet->getTiePoints()[i]->GcpNumberID<<"\t"
			<<sqrt(residue2.element(i * 3 + 0)*residue2.element(i * 3 + 0) + residue2.element(i * 3 + 1)*residue2.element(i * 3 + 1))<<"\t"
			<<residue2.element(i * 3 + 0)<<"\t"
			<<residue2.element(i * 3 + 1)<<"\t"
			<<"CHECK"<<"\t"
			<<""<<"\t"
			<<chkSet->getTiePoints()[i]->getGroundPoint().lat<<"\t"
			<<chkSet->getTiePoints()[i]->getGroundPoint().lon<<"\t"
			<<chkSet->getTiePoints()[i]->getGroundPoint().lat + residue2.element(i * 3 + 0)<<"\t"
			<<chkSet->getTiePoints()[i]->getGroundPoint().lon + residue2.element(i * 3 + 1)<<"\n";
	}

	fs.close();

	return true;
}

NEWMAT::ColumnVector RpcModel::CalcResidue(ossimSensorModel* sensorModel,ossimTieGptSet gptSet)
{
	ossimDpt imagepoint,cimagepoint;
	ossimGpt goundpoint,tGpt;
	ossimDpt residue1,residue2;
	int i;

	int num = static_cast<int>(gptSet.size());
	NEWMAT::ColumnVector residue;
	residue.ReSize(num*3);

	vector<ossimRefPtr<ossimTieGpt> >& theTPV = gptSet.refTiePoints();

	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;

	for(i = 0;i < num;++i)
	{
		imagepoint = gptSet.getTiePoints()[i]->getImagePoint();
		tGpt.hgt = gptSet.getTiePoints()[i]->hgt;
		sensorModel->lineSampleToWorld(imagepoint,tGpt);
		residue1 = sensorModel->m_proj->forward(tGpt);
		ossimGpt gpt = gptSet.getTiePoints()[i]->getGroundPoint();
		residue2 = ossimDpt(gpt.lat,gpt.lon);

		residue.element(i * 3 + 0) = residue2.x - residue1.x;
		residue.element(i * 3 + 1) = residue2.y - residue1.y;
		residue.element(i * 3 + 2) = 0;
	}

	return residue;
}


void RpcModel::SavePointToFile(ossimFilename filenametosave,ossimTieGptSet* m_gptset,ossimTieGptSet* mcheck_gptset)
{
	//////////////////////////////////////////////////输出txt
	if (filenametosave.ext().contains("txt")) {

		std::fstream     theTieFileStream;
		theTieFileStream.open(filenametosave.c_str(), ios_base::out);

		theTieFileStream<<"MAP_PROJECTION_BEGIN"<<endl;
		m_MapProjection.print(theTieFileStream);
		theTieFileStream<<"MAP_PROJECTION_END"<<endl;


		for (int i = 0; i<(int)m_gptset->refTiePoints().size(); i++)
		{
			theTieFileStream<<m_gptset->refTiePoints()[i]->GcpNumberID.c_str()<<"	"<< setiosflags(ios::fixed)<<
				setprecision(9)<<m_gptset->refTiePoints()[i]->refImagePoint().x<<"	"<<
				setprecision(9)<<m_gptset->refTiePoints()[i]->refImagePoint().y<<"	"<<
				setprecision(9)<<m_gptset->refTiePoints()[i]->refGroundPoint().lat<<"	"<<
				setprecision(9)<<m_gptset->refTiePoints()[i]->refGroundPoint().lon<<"	"<<
				setprecision(9)<<m_gptset->refTiePoints()[i]->refGroundPoint().hgt<<endl;

		}
		for (int i = 0; i<(int)mcheck_gptset->refTiePoints().size(); i++)
		{
			theTieFileStream<<"-"<<mcheck_gptset->refTiePoints()[i]->GcpNumberID.c_str()<<"	"<< setiosflags(ios::fixed)<<
				setprecision(9)<<mcheck_gptset->refTiePoints()[i]->refImagePoint().x<<"	"<<
				setprecision(9)<<mcheck_gptset->refTiePoints()[i]->refImagePoint().y<<"	"<<
				setprecision(9)<<mcheck_gptset->refTiePoints()[i]->refGroundPoint().lat<<"	"<<
				setprecision(9)<<mcheck_gptset->refTiePoints()[i]->refGroundPoint().lon<<"	"<<
				setprecision(9)<<mcheck_gptset->refTiePoints()[i]->refGroundPoint().hgt<<endl;

		}
		theTieFileStream.close();

	}
}

#pragma warning(pop)
//void outputTile(ossimFilename shpfilename,ossimFilename imgfilename,ossimFilename outfilename)
//{
//	vector<ossimDpt>  polygon;
//	ossimKeywordlist geo_geom, tt_geom;
//	ossimImageHandler *handler   = ossimImageHandlerRegistry::instance()->open(imgfilename);
//
//
//	ossimRefPtr<ossimImageGeometry>  dd = handler->getImageGeometry();
//	dd.get()->getProjection()->saveState(geo_geom);
//
//	//ossimImageGeometry *pl=handler->getImageGeometry();
//	//pl->saveState(geo_geom);
//	//ossimProjection *proj=pl->getProjection();
//
//	ossimMapProjection*	proj = NULL;
//	//handler->getImageGeometry(geo_geom);
//
//	proj=PTR_CAST(ossimMapProjection,ossimMapProjectionFactory::instance()->createProjection(geo_geom));
//	if (proj==NULL) {handler->close();return;}
//
//	ossimMapProjection*	MapPar=NULL;
//	tt_geom.clear();
//
//	tt_geom.addFile(PROVLambertCC);
//	MapPar=PTR_CAST(ossimMapProjection,ossimMapProjectionFactory::instance()->createProjection(tt_geom));
//	if (MapPar==NULL) {handler->close();return;}
//
//	tt_geom.clear();
//
//	const char* lookup;
//	ossimDpt tie;
//	//  lookup = geo_geom.find("projection.", ossimKeywordNames::TIE_POINT_XY_KW);
//	lookup = geo_geom.find("", ossimKeywordNames::TIE_POINT_XY_KW);
//	if (lookup)
//	{
//
//		tie.toPoint(std::string(lookup));
//		//theTranslateX  = tie.x;
//		//theTranslateY  = tie.y;
//	}
//	else {
//		handler->close();return;
//	}
//
//
//	//std::vector<ossim_uint32> m_OutBandList;
//	//m_OutBandList.clear();
//	//m_OutBandList.push_back(4);
//	//m_OutBandList.push_back(3);
//	//m_OutBandList.push_back(2);
//	//m_OutBandList.push_back(0);
//	//m_OutBandList.push_back(1);
//	//m_OutBandList.push_back(2);
//
//
//	OGRDataSource       *poDS;
//	int pos;
//	ossimString m_IdString;
//	poDS = OGRSFDriverRegistrar::Open( shpfilename.c_str(), FALSE );
//	if( poDS == NULL )
//	{
//		printf( "Open failed.\n" );
//		exit( 1 );
//	}
//	OGRLayer  *poLayer;
//	poLayer = poDS->GetLayer(0);//GetLayerByName( "point" );
//	OGRFeature *poFeature;
//	poLayer->ResetReading();
//	ossimGpt point_ll;
//	ossimDpt pt;
//	double MINX,MAXX,MINY,MAXY;		
//	MINX=999999999.0;MAXX=-999999999.0;
//	MINY=999999999.0;MAXY=-999999999.0;
//	while( (poFeature = poLayer->GetNextFeature()) != NULL )
//	{
//
//		OGRGeometry *poGeometry;
//		poGeometry = poFeature->GetGeometryRef();
//		pos=poFeature->GetFieldAsInteger(0);
//
//		if( poGeometry != NULL 	&& wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon && (pos==PROVINCE_ID))
//		{
//
//
//			OGRPolygon* pol=(OGRPolygon *) poGeometry;
//			OGRLinearRing* ring=pol->getExteriorRing();
//			polygon.clear();
//			for (int i=0;i<ring->getNumPoints();i++) {
//				OGRPoint poPoint ;
//
//				ring->getPoint(i,&poPoint);
//				//point_ll.lat=poPoint.getY();point_ll.lon=poPoint.getX();
//				//pt=proj->forward(point_ll);
//				pt.x=poPoint.getX();pt.y=poPoint.getY();
//
//				if (pt.x < MINX ) MINX=pt.x;
//				if (pt.x > MAXX ) MAXX=pt.x;
//				if (pt.y < MINY ) MINY=pt.y;
//				if (pt.y > MAXY ) MAXY=pt.y;
//			}
//		}
//		else
//		{
//			//  printf( "no point geometry\n" );
//		}       
//		OGRFeature::DestroyFeature( poFeature );
//	}
//
//	OGRDataSource::DestroyDataSource( poDS );
//
//	polygon.clear();
//	pt.x=(MINX-tie.x)/RESOLUTION;pt.y=(tie.y-MAXY)/RESOLUTION;
//	polygon.push_back(pt);//左上
//	pt.x=(MAXX-tie.x)/RESOLUTION;pt.y=(tie.y-MAXY)/RESOLUTION;
//	polygon.push_back(pt);//右上
//	pt.x=(MAXX-tie.x)/RESOLUTION;pt.y=(tie.y-MINY)/RESOLUTION;
//	polygon.push_back(pt);//右下
//	pt.x=(MINX-tie.x)/RESOLUTION;pt.y=(tie.y-MINY)/RESOLUTION;
//	polygon.push_back(pt);//左下
//	pt.x=(MINX-tie.x)/RESOLUTION;pt.y=(tie.y-MAXY)/RESOLUTION;
//	polygon.push_back(pt);
//
//
//	//ossimBandSelector* theBandSelector;
//	//theBandSelector = new ossimBandSelector;
//	//theBandSelector->connectMyInputTo(0, handler);
//	//theBandSelector->setOutputBandList(m_OutBandList);
//
//
//	ossimPolyCutter* theCutter;
//
//	ossimIrect bound,boundw;
//	theCutter = new ossimPolyCutter;
//	theCutter->setNumberOfPolygons(0);
//
//
//	theCutter->connectMyInputTo(handler);
//	theCutter->addPolygon(polygon);
//	theCutter->setCutType(ossimPolyCutter::ossimPolyCutterCutType::OSSIM_POLY_NULL_OUTSIDE);
//	theCutter->setNumberOfPolygons(1);
//
//	bound = theCutter->getBoundingRect();
//
//
//	ossimImageFileWriter* writer = ossimImageWriterFactoryRegistry::instance()->createWriter(ossimString("tiff_strip"));
//	ossimImageRenderer* renderer = new ossimImageRenderer;
//  renderer->getResampler()->setFilterType(ossimFilterResampler::ossimFilterResamplerType::ossimFilterResampler_BILINEAR);
//	renderer->connectMyInputTo(handler);////////    handler  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//	renderer->setView(MapPar);
//
//	tt_geom.clear();
//	writer->saveState(tt_geom);
//	tt_geom.add("",ossimKeywordNames::PIXEL_TYPE_KW,"PIXEL_IS_AREA",true);
//
//
//
//	ossimStdOutProgress progress(0,true);
//	writer->addListener(&progress);
//	writer->connectMyInputTo(0,renderer);	
//
//
//	ossimImageViewTransform*  IVT=renderer->getImageViewTransform();
//
//	ossimImageViewProjectionTransform* IVP = PTR_CAST(ossimImageViewProjectionTransform, IVT);
//
//	ossimDrect tempRect = IVT->getImageToViewBounds(bound);
//	ossimIrect outrect(tempRect),blockrect;
//	int row,col;
//	ossimFilename tmpname,outname;
//	row=(int)(outrect.height()/blocksize);
//	col=(int)(outrect.width()/blocksize);
//	if(row*blocksize < outrect.height() ) row =row +1;
//	if(col*blocksize < outrect.width() ) col =col +1;
//	ossimString m_tilfile;
//	m_tilfile=PROVLambertTIL;
//	ofstream fout(m_tilfile.c_str(),ios::out);
//	ossimFilename file;   
//
//	//////////////////////////////////////////////////////////////////////////////
//
//
//
//
//	/////////////////////////////////////////////////////////////////////////////
//
//	fout<<"bandId = \"RGB\""<<endl;
//	fout<<"numTiles = "<<col*row<<";"<<endl;
//	fout<<"tileSizeX = "<<col<<";"<<endl;
//	fout<<"tileSizeY = "<<row<<";"<<endl;
//	fout<<"tileUnits = \"Pixels\";"<<endl;
//	fout<<"tileOverlap = 0;"<<endl;
//
//
//	int starx=outrect.ul().x-OVERLAP;
//	int stary=outrect.ul().y-OVERLAP;
//	std::vector<ossimFilename> tmpnameList;
//	for(int i=1;i<=row;i++)
//		for(int j=1;j<=col;j++)
//		{
//
//			blockrect.set_ulx(outrect.ul().x+(j-1)*blocksize-OVERLAP);blockrect.set_uly(outrect.ul().y+(i-1)*blocksize-OVERLAP);
//			blockrect.set_lrx(blockrect.ul().x+blocksize+OVERLAP-1);blockrect.set_lry(blockrect.ul().y+blocksize+OVERLAP-1);
//
//			tmpname=outfilename+"_R"+ossimString::toString(i)+"C"+ossimString::toString(j)+"_tmp.TIF";
//			outname=outfilename+"_R"+ossimString::toString(i)+"C"+ossimString::toString(j)+".TIF";
//
//
//			tmpnameList.push_back(tmpname);
//
//			//		str="BEGIN_GROUP = TILE_"+str_filenum+"\n"+"filename = "+str_line+"_"+str_pixel+".tif"+"\n"+"ULColOffset = "+str_ULCO+";\n"+"ULRowOffset = "+str_ULRO+";\n"+"LRColOffset = "+str_LRCO+";\n"+"LRRowOffset = "+str_LRRO+";\n"+"END_GROUP = TILE_"+str_filenum+"\n";
//			fout<<"BEGIN_GROUP = TILE_"<<col*(i-1)+j<<endl;
//			fout<<"filename = "<<tmpname.file().c_str()<<endl;
//			fout<<"ULColOffset = "<<blockrect.ul().x-starx<<";"<<endl;
//			fout<<"ULRowOffset = "<<blockrect.ul().y-stary<<";"<<endl;
//			fout<<"LRColOffset = "<<blockrect.lr().x-starx<<";"<<endl;
//			fout<<"LRRowOffset = "<<blockrect.lr().y-stary<<";"<<endl;
//			fout<<"END_GROUP = TILE_"<<col*(i-1)+j<<endl;
//
//
//			if (tmpname.exists()) continue;
//		/*	if (tmpname.exists()) cutbyOUTshpfile(tmpname,outname); 
//			else
//			{*/
//				writer->setAreaOfInterest(blockrect);
//				writer->setFilename(tmpname);
//				writer->execute();
//			//}
//			////////////////////////////////////////////////////////////////////////////
//			//	cutbyOUTshpfile(tmpname,outname);
//			/////////////////////////////////////////////////////////////////////////////////////
//		}
//		fout.close();
//		handler->close();
//
//				theCutter->clear();
//		theCutter->disableSource();
//		writer->disableListener();
//	//	writer->disconnectAllInputs();
//	//	writer->disconnectAllOutputs();
//		writer->removeListener(&progress);
//	//	writer->r
//		//	delete theBandSelector;
//		//	delete theCutter;
//		//	delete renderer;
//		//delete writer;
//////////////////////////////////////////////////////////////////////////////////
//		ossimFilename shp;
//		shp=PROVLambertTILSHP;
//		if(shp.exists()) return;
//
//
//			OGRDataSource *poWriteDS;
//			const char *pszDriverName = "ESRI Shapefile";
//			OGRSFDriver *poDriver;
//			poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName( pszDriverName );
//		
//			if( poDriver == NULL )
//			{
//				printf( "%s driver not available.\n", pszDriverName );
//				exit( 1 );
//			}
//			poWriteDS = poDriver->CreateDataSource( PROVLambertTILSHP, NULL );
//			if( poWriteDS == NULL )
//			{
//				printf( "Creation of output file failed.\n" );
//				exit( 1 );
//			}
//		
//			OGRLayer *poLayerWrite;
//		
//			poLayerWrite = poWriteDS->CreateLayer( "border", NULL, wkbPolygon, NULL );
//			if( poLayerWrite == NULL )
//			{
//				printf( "Layer creation failed.\n" );
//				exit( 1 );
//			}
//		
//		
//			OGRFieldDefn oField2( "NAME", OFTString );
//		
//			oField2.SetWidth(100);
//		
//			if( poLayerWrite->CreateField( &oField2 ) != OGRERR_NONE )
//			{
//				printf( "Creating Name field failed.\n" );
//				exit( 1 );
//			}
//				
//
//			const char* lookuptmp;
//			ossimDpt tietmp;
//			ossimKeywordlist geo_tmp;
//			ossimImageHandler *handlertmp;
//			ossimIrect recttmp;
//		for (int num=0;num<tmpnameList.size();num++) {
//
//			OGRFeature *WritepoFeature;
//
//			WritepoFeature = OGRFeature::CreateFeature( poLayerWrite->GetLayerDefn() );
//
//			WritepoFeature->SetField(0, tmpnameList[num].c_str() );
//
//			OGRLinearRing* pRing=NULL;
//			OGRPolygon* pPolygon=new OGRPolygon();
//			pRing=new OGRLinearRing();
//			
//			handlertmp   = ossimImageHandlerRegistry::instance()->open(tmpnameList[num]);
//			if(!handlertmp) return;
//			recttmp=handlertmp->getImageRectangle();
//			geo_tmp.clear();
//
//			ossimRefPtr<ossimImageGeometry>  ddd=handlertmp->getImageGeometry();
//			ddd.get()->getProjection()->saveState(geo_tmp);
//
//			//handlertmp->getImageGeometry(geo_tmp);
//
//			lookuptmp = geo_tmp.find("", ossimKeywordNames::TIE_POINT_XY_KW);
//			if (lookuptmp)
//			{
//
//				tietmp.toPoint(std::string(lookuptmp));
//
//				//blocksize
//			}
//			else {
//				handlertmp->close();return;
//			}
//
//					pRing->addPoint(tietmp.x,tietmp.y);
//					pRing->addPoint(tietmp.x+recttmp.width()*RESOLUTION90,tietmp.y);
//					pRing->addPoint(tietmp.x+recttmp.width()*RESOLUTION90,tietmp.y-recttmp.height()*RESOLUTION90);
//					pRing->addPoint(tietmp.x,tietmp.y-recttmp.height()*RESOLUTION90);
//					pRing->addPoint(tietmp.x,tietmp.y);
//
//	/*			tietmp.x=tietmp.x+OVERLAP*RESOLUTION90;
//				tietmp.y=tietmp.y-OVERLAP*RESOLUTION90;
//					pRing->addPoint(tietmp.x,tietmp.y);
//					pRing->addPoint(tietmp.x+blocksize*RESOLUTION90,tietmp.y);
//					pRing->addPoint(tietmp.x+blocksize*RESOLUTION90,tietmp.y-blocksize*RESOLUTION90);
//					pRing->addPoint(tietmp.x,tietmp.y-blocksize*RESOLUTION90);
//					pRing->addPoint(tietmp.x,tietmp.y);*/
//
//
//			pPolygon->addRing(pRing);
//			WritepoFeature->SetGeometry( pPolygon ); 
//
//			if( poLayerWrite->CreateFeature( WritepoFeature ) != OGRERR_NONE )
//				{
//					printf( "Failed to create feature in shapefile.\n" );
//					exit( 1 );
//				}
//
//			OGRFeature::DestroyFeature( WritepoFeature );	
//
//		}
//
//	OGRDataSource::DestroyDataSource( poWriteDS );
//
//
//
//
//}