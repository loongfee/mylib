#include "gcpUtil.h"
#include <ossim/projection/ossimEquDistCylProjection.h>

namespace mylib{

void readTextGcpFile(ossimFilename strFilename,
							ossimTieGptSet* &gcpSet,
							ossimTieGptSet* &chkSet,
							ossimKeywordlist* prjKwl/* = NULL*/)
{
	// 首先读取投影文件
	ossimFilename strProjectionFile = strFilename;
	strProjectionFile.setExtension("geom");
	const char* ellipse_code = "";
	ossimRefPtr<ossimMapProjection> MapProjection;
	if (prjKwl)
	{
		// 如果需要读取投影
		if (strProjectionFile.exists())
		{
			prjKwl->addFile(strProjectionFile);
			ellipse_code = prjKwl->find(ossimKeywordNames::DATUM_KW);
			// 测试投影有效性
			if(!(MapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(*prjKwl))))
			{
				// 如果无效，则清空
				cout<<"Warning: invalid projection, use geographic coordinates as default."<<endl;
				MapProjection = new ossimEquDistCylProjection();
				MapProjection->saveState(*prjKwl);
				//prjKwl->clear();

			}
		}
		else{
			// 无投影文件
			//cout<<"Warning: 未找到相应的投影文件投影文件("<<strProjectionFile.file()<<")"<<endl;
			prjKwl->clear();
		}
	}

	if (prjKwl && !(MapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(*prjKwl))))
	{
		// 兼容旧控制点文件投影格式
		ossimTempFilename pp;
		char strtmp[255];
		std::vector<ossimString> gcpcon;
		const char* ellipse_code;
		pp.generateRandomFile();
		fstream os;
		os.open(strFilename.c_str(), ios_base::in);
		fstream oss;
		oss.open(pp.c_str(),ios::out|ios_base::app);
		ossimString str;
		os>>str;
		while(os.getline(strtmp,255) ) {
			str=strtmp;
			oss<< str.c_str();oss<<endl;
			if(str.contains("MAP_PROJECTION_END"))
			{
				break;
			}
		}
		oss.close();

		prjKwl->addFile(pp.c_str());
		ellipse_code = prjKwl->find(ossimKeywordNames::DATUM_KW);
		// 测试投影有效性
		if(!(MapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(*prjKwl))))
		{
			// 如果无效，则清空
			cout << "Warning: invalid projection, use geographic coordinates as default." << endl;
			MapProjection = new ossimEquDistCylProjection();
			MapProjection->saveState(*prjKwl); 
		}
	}

	ossimGpt tieg;
	ossimDpt tied;
	ossimTieGpt *aTiePt;
	gcpSet->clearTiePoints();
	chkSet->clearTiePoints();

	fstream os;
	os.open(strFilename.c_str(), ios_base::in);
	char strtmp[256];
	vector<string> delimiterList;
	delimiterList.push_back(" ");
	delimiterList.push_back("\t");
	while(os.getline(strtmp,255))
	{
		vector<string> strList;
		string str(strtmp);
		int commentPos = str.find_first_of('#');
		if (commentPos >= 0)
		{
			str = str.substr(0, commentPos);
		}

		splitString(str, delimiterList, strList);

		if (strList.size() < 5) continue;

		aTiePt=new ossimTieGpt(tieg,tied,0);
		if(strList[0].empty()) {
			aTiePt->GcpNumberID=strList[1];
		}
		else
		{
			aTiePt->GcpNumberID=strList[0];
		}
		try
		{
			aTiePt->GcpNumberID=strList[0];
			tied.x = ossimString(strList[1]).toDouble();
			tied.y = ossimString(strList[2]).toDouble();
			if (!MapProjection || MapProjection->isGeographic())
			{
				tieg.lat = ossimString(strList[4]).toDouble();
				tieg.lon = ossimString(strList[3]).toDouble();
			}
			else
			{
				tieg.lat = ossimString(strList[3]).toDouble();
				tieg.lon = ossimString(strList[4]).toDouble();
			}
			tieg.hgt = ossimString(strList[5]).toDouble();
		}
		catch (...)
		{
			continue;
		}
		tieg.datum(ossimDatumFactory::instance()->create(ossimString(ellipse_code))); 

		aTiePt->setGroundPoint(tieg);
		aTiePt->setImagePoint(tied);
		if (0 != strcmp(aTiePt->GcpNumberID.substr(0, 1).c_str(), "-"))
		//if (0 != strId.left(1).compare("-"))
		{
			if (gcpSet)
			{
				gcpSet->addTiePoint(aTiePt);
			}
		}
		else
		{
			if (chkSet)
			{
				aTiePt->GcpNumberID = aTiePt->GcpNumberID.substr(1, aTiePt->GcpNumberID.size()-1);
				chkSet->addTiePoint(aTiePt);
			}
		}

	}
	os.close();
}

void readXmlGcpFile(ossimFilename strFilename,
						   ossimTieGptSet* &gcpSet,
						   ossimTieGptSet* &chkSet,
						   ossimKeywordlist* prjKwl/* = NULL*/)
{
	gcpSet->clearTiePoints();
	chkSet->clearTiePoints();
	//ossimTieGptSet* theTieSet;
	vector<ossimRefPtr<ossimTieGpt> >::iterator tit;

	ossimDpt imagepoint;
	ossimGpt tGpt;
	//theTieSet = new ossimTieGptSet;
	ossimXmlDocument gmlDoc;
	gmlDoc.openFile(strFilename);
	std::vector< ossimRefPtr< ossimXmlNode > > tieSetList;
	gmlDoc.findNodes(ossimString("/") + ossimTieGptSet::TIEPTSET_TAG, tieSetList);

	if (tieSetList.size() != 1)
	{
		ossimNotify(ossimNotifyLevel_WARN) << 
			"WARNING: ossimModelOptimizer::loadGMLTieSet need exactly one element of type "<<
			ossimTieGptSet::TIEPTSET_TAG<<", found "<<tieSetList.size()<<"\n";
		return;
	}
	//theTieSet->clearTiePoints();
	gcpSet->importFromGmlNode(tieSetList[0]);
	vector<ossimRefPtr<ossimTieGpt> >& theGcp = gcpSet->refTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::iterator iter;
	for (iter = theGcp.begin() ; iter != theGcp.end() && iter->valid() ;)
	{
		//if ((*iter)->getScore() < 0.95)
		//{
		//	theGcp.erase(iter);
		//	continue;
		//}
		tGpt = (*iter)->getGroundPoint();
		int index = iter - theGcp.begin();
		char strId[10];
		_itoa(index+1, strId, 10);
		(*iter)->GcpNumberID = ossimString(strId);
		++iter;
	}
}


void readGcpFile(ossimFilename strFilename,
						ossimTieGptSet* &gcpSet,
						ossimTieGptSet* &chkSet,
						ossimKeywordlist* prjKwl/* = NULL*/)
{
	if(strFilename.ext().upcase().contains("XML"))
	{
		readXmlGcpFile(strFilename, gcpSet, chkSet, prjKwl);
	}
	else
	{
		readTextGcpFile(strFilename, gcpSet, chkSet, prjKwl);
	}
}

void saveGcpFile(ossimFilename filenametosave,
						ossimTieGptSet* gcpSet, 
						ossimTieGptSet* chkSet/* = NULL*/, 
						ossimKeywordlist* prjKwl/* = NULL*/,
						bool extern_file/* = true*/)
{

	////////////////////////////////////////////////
	if (filenametosave.ext().contains("txt"))
	{
		std::fstream     theTieFileStream;
		theTieFileStream.open(filenametosave.c_str(), ios_base::out);

		if (prjKwl)
		{
			if (extern_file)
			{
				// 保存投影文件
				ossimFilename strProjectionFile = filenametosave.path() + "//" + filenametosave.fileNoExtension() + ".geom";
				fstream os;
				os.open(strProjectionFile.c_str(), ios_base::out);
				os<<*prjKwl;
				os.close();
			}
			else{

				theTieFileStream<<"MAP_PROJECTION_BEGIN"<<endl;
				prjKwl->print(theTieFileStream);
				theTieFileStream<<"MAP_PROJECTION_END"<<endl;
			}
		}

		// 输出txt
		theTieFileStream.setf(ios::fixed, ios::floatfield);
		if (gcpSet)
		{
			for (int i = 0; i<(int)gcpSet->refTiePoints().size(); i++)
			{
				if (i != 0)
				{
					theTieFileStream<<endl;
				}
				theTieFileStream<<setiosflags(ios::left)<<
					setprecision(0)<<setw(10)<<gcpSet->refTiePoints()[i]->GcpNumberID.c_str()<<
					setprecision(10)<<setw(25)<<gcpSet->refTiePoints()[i]->refImagePoint().x<<
					setprecision(10)<<setw(25)<<gcpSet->refTiePoints()[i]->refImagePoint().y<<
					setprecision(10)<<setw(25)<<gcpSet->refTiePoints()[i]->refGroundPoint().lon<<
					setprecision(10)<<setw(25)<<gcpSet->refTiePoints()[i]->refGroundPoint().lat<<
					setprecision(10)<<setw(25)<<gcpSet->refTiePoints()[i]->refGroundPoint().hgt;

			}
		}
		if (chkSet)
		{
			for (int i = 0; i<(int)chkSet->refTiePoints().size(); i++)
			{
				//if (i != 0)
				//{
				theTieFileStream<<endl;
				//}
				theTieFileStream<<setiosflags(ios::left)<<
					setprecision(0)<<setw(10)<<"-"+chkSet->refTiePoints()[i]->GcpNumberID<<
					setprecision(10)<<setw(25)<<chkSet->refTiePoints()[i]->refImagePoint().x<<
					setprecision(10)<<setw(25)<<chkSet->refTiePoints()[i]->refImagePoint().y<<
					setprecision(10)<<setw(25)<<chkSet->refTiePoints()[i]->refGroundPoint().lon<<
					setprecision(10)<<setw(25)<<chkSet->refTiePoints()[i]->refGroundPoint().lat<<
					setprecision(10)<<setw(25)<<chkSet->refTiePoints()[i]->refGroundPoint().hgt;
			}
		}
		theTieFileStream.close();
	}
	////////////////////////////////////////////////////输出shpfile
	//if (filenametosave.ext().contains("shp")) {
	//	bool flag=false;
	//	flag=WriteShpPoint(filenametosave,m_gptset,mcheck_gptset);

	//}
}

bool get_elevation(ossimTieGptSet* &ctrlSet, ossimKeywordlist prjKwl, double defaultElev/* = 0.0*/)
{
	vector<ossimRefPtr<ossimTieGpt> >& theGcp = ctrlSet->refTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::iterator iter,tit;

	ossimMapProjection* theMapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(prjKwl));
	if (!ossimElevManager::instance())
	{

		for (iter = theGcp.begin() ; iter != theGcp.end() ; ++iter)
		{
			(*iter)->hgt = defaultElev;
		}
		if (!ossimElevManager::instance())
		{
			cout<<"高程源还未选取！"<<endl;
			return false;
		}
	}

	ossimGpt tGpt;
	for (iter = theGcp.begin() ; iter != theGcp.end() && iter->valid() ;)
	{
		if (NULL == theMapProjection || theMapProjection->isGeographic())
		{
			tGpt = ossimGpt((*iter)->lat, (*iter)->lon);
		}
		else
		{
			tGpt = theMapProjection->inverse(ossimDpt((*iter)->lat, (*iter)->lon));
		}

		(*iter)->hgt=ossimElevManager::instance()->getHeightAboveEllipsoid(tGpt);
		if(ossim::isnan((*iter)->hgt))
		//if(ossim::isnan((*iter)->hgt) || (*iter)->hgt == 0.0)
			//if(ossim::isnan((*iter)->hgt) || (*iter)->hgt == 0.0)
		{
			//cout<<tGpt.lat<<" "<<tGpt.lon<<" "<<(*iter)->hgt<<endl;
			// 如果为空，则去掉该点
			//theGcp.erase(iter);
			// 如果为空，则赋为平均值
			(*iter)->hgt = defaultElev;
			++iter;
		}
		else{
			++iter;
		}
	}
	return true;
}

bool get_elevation(ossimTieGptSet* &ctrlSet, double defaultElev/* = 0.0*/)
{
	vector<ossimRefPtr<ossimTieGpt> >& theGcp = ctrlSet->refTiePoints();
	vector<ossimRefPtr<ossimTieGpt> >::iterator iter, tit;

	if (!ossimElevManager::instance())
	{

		for (iter = theGcp.begin(); iter != theGcp.end(); ++iter)
		{
			(*iter)->hgt = defaultElev;
		}
		if (!ossimElevManager::instance())
		{
			cout << "高程源还未选取！" << endl;
			return false;
		}
	}

	ossimGpt tGpt;
	for (iter = theGcp.begin(); iter != theGcp.end() && iter->valid();)
	{
		tGpt = ossimGpt((*iter)->lat, (*iter)->lon);

		(*iter)->hgt = ossimElevManager::instance()->getHeightAboveEllipsoid(tGpt);
		if (ossim::isnan((*iter)->hgt))
			//if(ossim::isnan((*iter)->hgt) || (*iter)->hgt == 0.0)
			//if(ossim::isnan((*iter)->hgt) || (*iter)->hgt == 0.0)
		{
			//cout<<tGpt.lat<<" "<<tGpt.lon<<" "<<(*iter)->hgt<<endl;
			// 如果为空，则去掉该点
			//theGcp.erase(iter);
			// 如果为空，则赋为平均值
			(*iter)->hgt = defaultElev;
			++iter;
		}
		else{
			++iter;
		}
	}
	return true;
}

void get_elevation(ossimFilename infile, ossimFilename outfile, ossimFilename elevationPath/* = ""*/, double defaultElev/* = 0.0*/)
{
	if(!ossimElevManager::instance()->loadElevationPath(elevationPath))
	{
		cout<<"warning: 加载DEM失败！"<<endl;
		return;
	}
	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	readGcpFile(infile, gcpSet, chkSet, &prjKwl);
	ossimElevManager::instance()->loadElevationPath(elevationPath);
	get_elevation(gcpSet, prjKwl, defaultElev);
	get_elevation(chkSet, prjKwl, defaultElev);
	saveGcpFile(outfile, gcpSet, chkSet, &prjKwl);
}

void batch_get_elevation(ossimFilename inDir, ossimFilename outDir, 
								ossimFilename filter/*="*.txt"*/, 
								ossimFilename elevationPath/*=""*/,
								double defaultElev/* = 0.0*/)
{
	if (!outDir.isDir())
	{
		_mkdir(outDir.c_str());
	}

	std::vector<ossimFilename> allFiles;
	//ossimString strReg = "/\.txt$/i ";
	ossimDirectory(inDir).findAllFilesThatMatch(allFiles, filter);

	for (int i = 0;i < (int)allFiles.size();++i)
	{
		ossimFilename infile = allFiles[i];
		ossimFilename outfile = outDir + "\\" + infile.file();
		get_elevation(infile, outfile, elevationPath, defaultElev);
	}
}


ossimTieGptSet* projection2ll(ossimTieGptSet* inGcpSet, ossimKeywordlist prjKwl)
{
	ossimTieGptSet* outGcpSet = new ossimTieGptSet;
	*outGcpSet = *inGcpSet;

	ossimMapProjection* theMapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(prjKwl));

	if (!theMapProjection)
	{
		cout<<"投影错误或没有投影信息！"<<endl;
		return NULL;
	}

	if (theMapProjection->isGeographic())
	{
		return inGcpSet;
	}

	for (vector<ossimRefPtr<ossimTieGpt> >::iterator iter = outGcpSet->refTiePoints().begin() ; iter != outGcpSet->refTiePoints().end() && iter->valid();)
	{
		ossimGpt oldGpt = (*iter)->getGroundPoint();
		ossimGpt gpt = theMapProjection->inverse(ossimDpt((*iter)->getGroundPoint().lat, (*iter)->getGroundPoint().lon));
		(*iter)->setGroundPoint(ossimGpt(gpt.lat, gpt.lon, oldGpt.hgt));
		iter++;
	}
	return outGcpSet;
}

ossimTieGptSet* datum_shift(ossimTieGptSet* inGcpSet, const char* outDatumString)
{
	const ossimDatum* outDatum =
		ossimDatumFactory::instance()->create(ossimString(outDatumString));

	ossimTieGptSet* outGcpSet = new ossimTieGptSet;
	*outGcpSet = *inGcpSet;

	
	if (!outDatum)
	{
		cout << "椭球体错误！" << endl;
		return NULL;
	}

	for (vector<ossimRefPtr<ossimTieGpt> >::iterator iter = outGcpSet->refTiePoints().begin(); iter != outGcpSet->refTiePoints().end() && iter->valid();)
	{
		(*iter)->refGroundPoint().changeDatum(outDatum);
		iter++;
	}
	return outGcpSet;
}

bool projection2ll(ossimFilename inFile, ossimFilename outFile, bool bGetElevation/* = false*/)
{
	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist prjKwl;
	readGcpFile(inFile, gcpSet, chkSet, &prjKwl);
	if (bGetElevation)
	{
		get_elevation(gcpSet, prjKwl, 0.0);
		get_elevation(chkSet, prjKwl, 0.0);
	}
	gcpSet = projection2ll(gcpSet, prjKwl);
	chkSet = projection2ll(chkSet, prjKwl);
	//saveGcpFile(outFile, gcpSet, chkSet);

	// 输出txt

	std::fstream     theTieFileStream;
	theTieFileStream.open(outFile.c_str(), ios_base::out);

	theTieFileStream.setf(ios::fixed, ios::floatfield);
	if (gcpSet)
	{
		for (int i = 0; i<(int)gcpSet->refTiePoints().size(); i++)
		{
			theTieFileStream<<setiosflags(ios::left)<<
				setprecision(0)<<setw(10)<<gcpSet->refTiePoints()[i]->GcpNumberID.c_str()<<
				setprecision(3)<<setw(15)<<gcpSet->refTiePoints()[i]->refImagePoint().x<<
				setprecision(3)<<setw(15)<<gcpSet->refTiePoints()[i]->refImagePoint().y<<
				setprecision(10)<<setw(25)<<gcpSet->refTiePoints()[i]->refGroundPoint().lon<<
				setprecision(10)<<setw(25)<<gcpSet->refTiePoints()[i]->refGroundPoint().lat<<
				setprecision(3)<<setw(15)<<gcpSet->refTiePoints()[i]->refGroundPoint().hgt<<endl;

		}
	}
	if (chkSet)
	{
		for (int i = 0; i<(int)chkSet->refTiePoints().size(); i++)
		{
			theTieFileStream<<setiosflags(ios::left)<<
				setprecision(0)<<setw(10)<<"-"+chkSet->refTiePoints()[i]->GcpNumberID<<
				setprecision(3)<<setw(15)<<chkSet->refTiePoints()[i]->refImagePoint().x<<
				setprecision(3)<<setw(15)<<chkSet->refTiePoints()[i]->refImagePoint().y<<
				setprecision(10)<<setw(25)<<chkSet->refTiePoints()[i]->refGroundPoint().lon<<
				setprecision(10)<<setw(25)<<chkSet->refTiePoints()[i]->refGroundPoint().lat<<
				setprecision(3)<<setw(15)<<chkSet->refTiePoints()[i]->refGroundPoint().hgt<<endl;
		}
	}
	theTieFileStream.close();

	//if (prjKwl)
	//{
	//	// 保存投影文件
	//	ossimFilename strProjectionFile = filenametosave + ".geom";
	//	fstream os;
	//	os.open(strProjectionFile.c_str(), ios_base::out);
	//	os<<*prjKwl;
	//	os.close();			
	//}

	return true;
}

ossimTieGptSet* reprojectionPoints(ossimTieGptSet* inGcpSet, ossimKeywordlist inPrjKwl, ossimKeywordlist outPrjKwl)
{
	ossimTieGptSet* outGcpSet = new ossimTieGptSet;
	*outGcpSet = *inGcpSet;

	ossimMapProjection* inMapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(inPrjKwl));
	ossimMapProjection* outMapProjection = PTR_CAST(ossimMapProjection, ossimMapProjectionFactory::instance()->createProjection(outPrjKwl));

	if (!inMapProjection || !outMapProjection)
	{
		cout<<"投影错误或没有投影信息！"<<endl;
		return NULL;
	}

	const char* in_ellipse_code;
	const char* out_ellipse_code;
	in_ellipse_code = inPrjKwl.find(ossimKeywordNames::DATUM_KW);
	out_ellipse_code = outPrjKwl.find(ossimKeywordNames::DATUM_KW);

	const ossimDatum* inDatum = NULL;
	const ossimDatum* outDatum = NULL;
	bool bDatumShift = false;
	if (0 != _strcmpi(in_ellipse_code, out_ellipse_code))
	{
		bDatumShift = true;
		inDatum =
			ossimDatumFactory::instance()->create(ossimString(in_ellipse_code));
		outDatum =
			ossimDatumFactory::instance()->create(ossimString(out_ellipse_code));
	}

	for (vector<ossimRefPtr<ossimTieGpt> >::iterator iter = outGcpSet->refTiePoints().begin() ; iter != outGcpSet->refTiePoints().end() && iter->valid() ;)
	{
		ossimGpt oldGpt = (*iter)->getGroundPoint();
		ossimGpt ll = inMapProjection->inverse(ossimDpt(oldGpt.lat, oldGpt.lon));
		ll.hgt = oldGpt.hgt;
		if (bDatumShift)
		{
			ll = ossimGpt(ll.lat, ll.lon, ll.hgt, inDatum);
			ll.changeDatum(outDatum);
		}
		ossimDpt newGpt = outMapProjection->forward(ll);
		(*iter)->setGroundPoint(ossimGpt(newGpt.x, newGpt.y, ll.hgt));
		iter++;
	}
	return outGcpSet;
}

bool reprojectionPoints(ossimFilename inFile, ossimKeywordlist outPrjKwl, ossimFilename outFile)
{
	ossimTieGptSet* gcpSet = new ossimTieGptSet;
	ossimTieGptSet* chkSet = new ossimTieGptSet;
	ossimKeywordlist inPrjKwl;
	readGcpFile(inFile, gcpSet, chkSet, &inPrjKwl);

	gcpSet = reprojectionPoints(gcpSet, inPrjKwl, outPrjKwl);
	chkSet = reprojectionPoints(chkSet, inPrjKwl, outPrjKwl);
	saveGcpFile(outFile, gcpSet, chkSet, &outPrjKwl);
	return true;
}

//bool drawPointsInThumbnail(ossimFilename originImageFile,
//						   ossimTieGptSet* ctrSet,
//						   ossimTieGptSet* chkSet,
//						   int xSize = 1024,
//						   int ySize = 0,
//						   vector<int> bandList = vector<int>(),
//						   )
//{
//
//}

} // end of namespace mylib