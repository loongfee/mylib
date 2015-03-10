#include "exteriorOptimize.h"

namespace Cbers04Optimize{
	namespace exteriorOptimize{		

		NEWMAT::ColumnVector
			getResidue(const vector< Cbers04OptStruct >& optStructList)
		{
			NEWMAT::ColumnVector residue;
			int nImages = (int)optStructList.size();
			int nTotalEquation = 0;
			int dimObs = 3;
			for(int i = 0;i < nImages;++i)
			{
				nTotalEquation += dimObs * (int)optStructList[i].gcpSet->getTiePoints().size();
			}
			if (nImages < 1)
			{
				return residue;
			}
			residue.ReSize(nTotalEquation);
			unsigned long c=1;
			for (int image_index = 0;image_index < nImages;++image_index)
			{
				const vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->getTiePoints();
				vector<ossimRefPtr<ossimTieGpt> >::const_iterator tit;
				//if(!optStructList[image_index].cbers04Model->m_proj) return residue;

				//// forward
				//ossimDpt resIm;
				//for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
				//{
				//	resIm = (*tit)->tie - optStructList[image_index].cbers04Model->forward(**tit);
				//	residue(c++) = resIm.x;
				//	residue(c++) = resIm.y;
				//}

				// inverse
				// ground observations
				ossimGpt gd;
				// loop on tie points
				for (tit = theTPV.begin(); tit != theTPV.end(); ++tit)
				{
					//compute residue
					gd = optStructList[image_index].cbers04Model->inverse((*tit)->tie);
					residue(c++) = ((*tit)->lon - gd.lon) * 100000.0; //approx meters //TBC TBD
					residue(c++) = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
					residue(c++) = (*tit)->hgt - gd.hgt; //meters
				}
			}
			return residue;
		}

		ossimDpt getForwardDeriv(radiCbers04Model* cbers04Model, const ossimGpt& gpos, int paramIdx,
			int ephIdx/* = 0*/, double hdelta/* = 1e-6*/)
		{
			double den = 0.5/hdelta;
			ossimDpt res;

			double middle = getParameter(cbers04Model, paramIdx, ephIdx);
			setParameter(middle + hdelta, cbers04Model, paramIdx, ephIdx);
			//res = inverse(gpos);
			res = cbers04Model->forward(gpos);
			setParameter(middle - hdelta, cbers04Model, paramIdx, ephIdx);
			//res -= inverse(gpos);
			res -= cbers04Model->forward(gpos);
			res = res*den;
			setParameter(middle, cbers04Model, paramIdx, ephIdx);

			return res;
		}

		ossimGpt getInverseDeriv(radiCbers04Model* cbers04Model, const ossimDpt& dpos, int paramIdx,
			int ephIdx/* = 0*/, double hdelta/* = 1e-6*/)
		{
			double den = 0.5/hdelta;
			ossimGpt res, gd;

			double middle = getParameter(cbers04Model, paramIdx, ephIdx);
			setParameter(middle + hdelta, cbers04Model, paramIdx, ephIdx);
			//res = inverse(gpos);
			res = cbers04Model->inverse(dpos);
			setParameter(middle - hdelta, cbers04Model, paramIdx, ephIdx);
			//res -= inverse(gpos);
			gd = cbers04Model->inverse(dpos);
			setParameter(middle, cbers04Model, paramIdx, ephIdx);

			res.lon = den*(res.lon - gd.lon) * 100000.0; //TBC : approx meters
			res.lat = den*(res.lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
			res.hgt = den*(res.hgt - gd.hgt);

			return res;
		}

		//ossimDpt getForwardDeriv(radiCbers04Model* cbers04Model,
		//	int iAtt, int iComponent, const ossimGpt& gpos, double hdelta)
		//{   
		//	double den = 0.5/hdelta;
		//	ossimDpt res;
		//	int idx = 0;
		//	if (idx++ == iComponent)
		//	{
		//		double middle = cbers04Model->theSupportData->theAttitudeSamples[iAtt].x;
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].x = middle + hdelta;
		//		//res = inverse(gpos);
		//		res = cbers04Model->forward(gpos);
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].x = middle - hdelta;
		//		//res -= inverse(gpos);
		//		res -= cbers04Model->forward(gpos);
		//		res = res*den;
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].x = middle;
		//	}
		//	else if(idx++ == iComponent)
		//	{
		//		double middle = cbers04Model->theSupportData->theAttitudeSamples[iAtt].y;
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].y = middle + hdelta;
		//		//res = inverse(gpos);
		//		res = cbers04Model->forward(gpos);
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].y = middle - hdelta;
		//		//res -= inverse(gpos);
		//		res -= cbers04Model->forward(gpos);
		//		res = res*den;
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].y = middle;
		//	}
		//	else if(idx++ == iComponent)
		//	{
		//		double middle = cbers04Model->theSupportData->theAttitudeSamples[iAtt].z;
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].z = middle + hdelta;
		//		//res = inverse(gpos);
		//		res = cbers04Model->forward(gpos);
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].z = middle - hdelta;
		//		//res -= inverse(gpos);
		//		res -= cbers04Model->forward(gpos);
		//		res = res*den;
		//		cbers04Model->theSupportData->theAttitudeSamples[iAtt].z = middle;
		//	}
		//	else if(idx++ == iComponent)
		//	{
		//		double middle = cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].x;
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].x = middle + hdelta;
		//		//res = inverse(gpos);
		//		res = cbers04Model->forward(gpos);
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].x = middle - hdelta;
		//		//res -= inverse(gpos);
		//		res -= cbers04Model->forward(gpos);
		//		res = res*den;
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].x = middle;
		//	}
		//	else if(idx++ == iComponent)
		//	{
		//		double middle = cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].y;
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].y = middle + hdelta;
		//		//res = inverse(gpos);
		//		res = cbers04Model->forward(gpos);
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].y = middle - hdelta;
		//		//res -= inverse(gpos);
		//		res -= cbers04Model->forward(gpos);
		//		res = res*den;
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].y = middle;
		//	}
		//	else if(idx++ == iComponent)
		//	{
		//		double middle = cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].z;
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].z = middle + hdelta;
		//		//res = inverse(gpos);
		//		res = cbers04Model->forward(gpos);
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].z = middle - hdelta;
		//		//res -= inverse(gpos);
		//		res -= cbers04Model->forward(gpos);
		//		res = res*den;
		//		cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].z = middle;
		//	}
		//	else if(idx++ == iComponent)
		//	{
		//		double middle = cbers04Model->theSupportData->theSampTimesBias[iAtt];
		//		cbers04Model->theSupportData->theSampTimesBias[iAtt] = middle + hdelta;
		//		//res = inverse(gpos);
		//		res = cbers04Model->forward(gpos);
		//		cbers04Model->theSupportData->theSampTimesBias[iAtt] = middle - hdelta;
		//		//res -= inverse(gpos);
		//		res -= cbers04Model->forward(gpos);
		//		res = res*den;
		//		cbers04Model->theSupportData->theSampTimesBias[iAtt] = middle;
		//	}
		//	return res;
		//}

		//void buildNormalEquation(const vector< Cbers04OptStruct >& optStructList, 
		//	const vector<int>& paramList,
		//	arma::mat& A,
		//	arma::vec& residue,
		//	double pstep_scale)
		//{
		//	//init
		//	int nImages = (int)optStructList.size();
		//	int np = getParameterNum(optStructList, paramList);

		//	int nTotalEquation = 0;
		//	for(int i = 0;i < nImages;++i)
		//	{
		//		nTotalEquation += 2 * (int)optStructList[i].gcpSet->getTiePoints().size();
		//	}

		//	//Zeroify matrices that will be accumulated
		//	A = arma::zeros(nTotalEquation, np);
		//	residue = arma::zeros(nTotalEquation);

		//	unsigned long c=1;
		//	//image observations 
		//	std::vector<ossimDpt> imDerp(np, ossimDpt(0.0, 0.0));
		//	ossimDpt resIm;
		//	for(int image_index = 0;image_index < nImages;++image_index)
		//	{
		//		const vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->getTiePoints();
		//		vector<ossimRefPtr<ossimTieGpt> >::const_iterator tit;
		//		// loop on tie points
		//		for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
		//		{
		//			//compute residue
		//			resIm = (*tit)->tie - optStructList[image_index].cbers04Model->forward(*(*tit));
		//			residue(c++) = resIm.x;
		//			residue(c++) = resIm.y;

		//			//compute all image derivatives regarding parametres for the tie point position
		//			int pos = 0;
		//			for (int j = 0;j < image_index;++j)
		//			{
		//				pos += getImageParameterNum(optStructList[image_index].cbers04Model, paramList);
		//			}

		//			int nEhp = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
		//			for (int ii = 0;ii < (int)paramList.size();++ii)
		//			{
		//				int iParam = paramList[ii];
		//				if (iParam < EPH_SPLIT)
		//				{
		//					// parameters for a scene
		//					imDerp[pos++] = getForwardDeriv(optStructList[image_index].cbers04Model, *(*tit), iParam, 0, pstep_scale);
		//				}
		//				else
		//				{
		//					// parameters for an ephemeris
		//					for (int iEhp = 0;iEhp < nEhp;++iEhp)
		//					{
		//						imDerp[pos++] = getForwardDeriv(optStructList[image_index].cbers04Model, *(*tit), iParam, iEhp, pstep_scale);
		//					}
		//				}
		//			}

		//			//int nAtt = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
		//			//for (int iAtt = 0;iAtt < nAtt;++iAtt)
		//			//{
		//			//	for (int ii = 0;ii < np_single;++ii)
		//			//	{
		//			//		imDerp[pos+iAtt*np_single+ii] = getForwardDeriv( optStructList[image_index].cbers04Model, iAtt , ii, *(*tit) , pstep_scale);
		//			//	}
		//			//}

		//			//compute influence of tie point on all sytem elements
		//			for(int p1=0;p1<np;++p1)
		//			{        
		//				//proj residue: J * residue
		//				projResidue.element(p1) += imDerp[p1].x * resIm.x + imDerp[p1].y * resIm.y;

		//				//normal matrix A = transpose(J)*J
		//				for(int p2=p1;p2<np;++p2)
		//				{
		//					A.element(p1,p2) += imDerp[p1].x * imDerp[p2].x + imDerp[p1].y * imDerp[p2].y;
		//				}
		//			}
		//		}
		//	}
		//}

		void buildNormalEquation(const vector< Cbers04OptStruct >& optStructList, 
			const vector<int>& paramList,
			NEWMAT::SymmetricMatrix& A,
			NEWMAT::ColumnVector& residue,
			NEWMAT::ColumnVector& projResidue,
			double pstep_scale)
		{
			//init
			int nImages = (int)optStructList.size();
			int np = getParameterNum(optStructList, paramList);

			int dimObs = 3;

			int nTotalEquation = 0;
			for(int i = 0;i < nImages;++i)
			{
				nTotalEquation += dimObs * (int)optStructList[i].gcpSet->getTiePoints().size();
			}
			A.ReSize(np);
			residue.ReSize(nTotalEquation);
			projResidue.ReSize(np);
			//Zeroify matrices that will be accumulated
			A = 0.0;
			projResidue = 0.0;


			unsigned long c=1;
			//image observations 
			std::vector<ossimDpt> imDerp(np, ossimDpt(0.0, 0.0));
			// ground observations
			std::vector<ossimGpt>  gdDerp(np);
			ossimDpt resIm;
			for(int image_index = 0;image_index < nImages;++image_index)
			{
				const vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[image_index].gcpSet->getTiePoints();
				vector<ossimRefPtr<ossimTieGpt> >::const_iterator tit;
				ossimGpt gd, resGd;
				// loop on tie points
				for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
				{
					// forward
					////compute residue
					//resIm = (*tit)->tie - optStructList[image_index].cbers04Model->forward(*(*tit));
					//residue(c++) = resIm.x;
					//residue(c++) = resIm.y;

					////compute all image derivatives regarding parametres for the tie point position
					//int pos = 0;
					//for (int j = 0;j < image_index;++j)
					//{
					//	pos += getImageParameterNum(optStructList[image_index].cbers04Model, paramList);
					//}

					//int nEhp = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
					//for (int ii = 0;ii < (int)paramList.size();++ii)
					//{
					//	int iParam = paramList[ii];
					//	if (iParam < EPH_SPLIT)
					//	{
					//		// parameters for a scene
					//		imDerp[pos++] = getForwardDeriv(optStructList[image_index].cbers04Model, *(*tit), iParam, 0, pstep_scale);
					//	}
					//	else
					//	{
					//		// parameters for an ephemeris
					//		for (int iEhp = 0;iEhp < nEhp;++iEhp)
					//		{
					//			imDerp[pos++] = getForwardDeriv(optStructList[image_index].cbers04Model, *(*tit), iParam, iEhp, pstep_scale);
					//		}
					//	}
					//}
					////compute influence of tie point on all sytem elements
					//for (int p1 = 0; p1<np; ++p1)
					//{
					//	//proj residue: J * residue
					//	projResidue.element(p1) += imDerp[p1].x * resIm.x + imDerp[p1].y * resIm.y;

					//	//normal matrix A = transpose(J)*J
					//	for (int p2 = p1; p2<np; ++p2)
					//	{
					//		A.element(p1, p2) += imDerp[p1].x * imDerp[p2].x + imDerp[p1].y * imDerp[p2].y;
					//	}
					//}

					// inverse
					//compute residue
					gd = optStructList[image_index].cbers04Model->inverse((*tit)->tie);
					residue(c++) = resGd.lon = ((*tit)->lon - gd.lon) * 100000.0;
					residue(c++) = resGd.lat = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
					residue(c++) = resGd.hgt = (*tit)->hgt - gd.hgt; //TBD : normalize to meters?

					//compute all image derivatives regarding parametres for the tie point position
					int pos = 0;
					for (int j = 0; j < image_index; ++j)
					{
						pos += getImageParameterNum(optStructList[image_index].cbers04Model, paramList);
					}

					int nEhp = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
					for (int ii = 0; ii < (int)paramList.size(); ++ii)
					{
						int iParam = paramList[ii];
						if (iParam < EPH_SPLIT)
						{
							// parameters for a scene
							gdDerp[pos++] = getInverseDeriv(optStructList[image_index].cbers04Model, (*tit)->tie, iParam, 0, pstep_scale);
						}
						else
						{
							// parameters for an ephemeris
							for (int iEhp = 0; iEhp < nEhp; ++iEhp)
							{
								gdDerp[pos++] = getInverseDeriv(optStructList[image_index].cbers04Model, (*tit)->tie, iParam, 0, pstep_scale);
							}
						}
					}

					//int nAtt = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
					//for (int iAtt = 0;iAtt < nAtt;++iAtt)
					//{
					//	for (int ii = 0;ii < np_single;++ii)
					//	{
					//		imDerp[pos+iAtt*np_single+ii] = getForwardDeriv( optStructList[image_index].cbers04Model, iAtt , ii, *(*tit) , pstep_scale);
					//	}
					//}

					//compute influence of tie point on all sytem elements
					for(int p1=0;p1<np;++p1)
					{        
						//proj residue: J * residue
						projResidue.element(p1) += gdDerp[p1].lon * resGd.lon + gdDerp[p1].lat * resGd.lat + gdDerp[p1].hgt * resGd.hgt;

						//normal matrix A = transpose(J)*J
						for(int p2=p1;p2<np;++p2)
						{
							A.element(p1, p2) += gdDerp[p1].lon * gdDerp[p2].lon + gdDerp[p1].lat * gdDerp[p2].lat + gdDerp[p1].hgt * gdDerp[p2].hgt;
						}
					}
				}
			}
		}


		void adjustment(const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList)
		{
			int nImages = (int)optStructList.size();
			int np = getParameterNum(optStructList, paramList);
			int dimObs = 3;

			int nTotalEquation = 0;
			for(int i = 0;i < nImages;++i)
			{
				nTotalEquation += dimObs * (int)optStructList[i].gcpSet->getTiePoints().size();
			}
			if (nImages < 1)
			{
				return;
			}
			int nobs = nTotalEquation / dimObs;

			//setup initail values
			int iter=0;
			int iter_max = 200;
			double minResidue = 1e-10; //TBC
			double minDelta = 1e-10; //TBC

			//build Least Squares initial normal equation
			// don't waste memory, add samples one at a time
			NEWMAT::SymmetricMatrix A;
			NEWMAT::ColumnVector residue;
			NEWMAT::ColumnVector projResidue;
			double deltap_scale = 1e-4; //step_Scale is 1.0 because we expect parameters to be between -1 and 1
			buildNormalEquation(optStructList, paramList,
				A, residue, projResidue, deltap_scale);
			double ki2=residue.SumSquare();

			NEWMAT::ColumnVector cparm(np), nparm(np);
			cparm = getParameters(optStructList, paramList);
			//int n = 1;
			//for (int image_index = 0;image_index < nImages;++image_index)
			//{
			//	int nAtt = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
			//	for (int iParam = 0;iParam < (int)paramList.size();++iParam)
			//	{
			//		if (iParam < EPH_SPLIT)
			//		{
			//			// parameters for a scene
			//		}
			//		else
			//		{
			//			// parameters for a ephemeris
			//		}
			//	}
			//	for (int iAtt = 0;iAtt < nAtt;++iAtt)
			//	{
			//		cparm(n++) = optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].x;
			//		cparm(n++) = optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].y;
			//		cparm(n++) = optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].z;
			//		cparm(n++) = optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].x;
			//		cparm(n++) = optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].y;
			//		cparm(n++) = optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].z;
			//		//cparm(n++) = optStructList[image_index].cbers04Model->theSupportData->theSampTimesBias[iAtt];
			//	}
			//}
			
			double damping_speed = 2.0;
			//find max diag element for A
			double maxdiag=0.0;
			for(int d=1;d<=np;++d) {
				if (maxdiag < A(d,d)) maxdiag=A(d,d);
			}
			double damping = 1e-3 * maxdiag;
			double olddamping = 0.0;
			bool found = false;

			while ( (!found) && (iter < iter_max) ) //non linear optimization loop
			{
				bool decrease = false;

				do
				{
					//add damping update to normal matrix
					for(int d=1;d<=np;++d) A(d,d) += damping - olddamping;
					olddamping = damping;

					NEWMAT::ColumnVector deltap = optStructList[0].cbers04Model->solveLeastSquares(A, projResidue);

					if (deltap.NormFrobenius() <= minDelta) 
					{
						found = true;
					} else {
						//update adjustment
						nparm = cparm + deltap;
						setParameters(nparm, optStructList, paramList);
						for (int i = 0; i < (int)optStructList.size(); i++)
						{
							optStructList[i].cbers04Model->updateModel();
						}
						//int n = 1;
						//for (int image_index = 0;image_index < nImages;++image_index)
						//{
						//	int nAtt = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
						//	for (int iAtt = 0;iAtt < nAtt;++iAtt)
						//	{
						//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].x = nparm(n++);
						//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].y = nparm(n++);
						//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].z = nparm(n++);
						//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].x = nparm(n++);
						//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].y = nparm(n++);
						//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].z = nparm(n++);
						//		//optStructList[image_index].cbers04Model->theSupportData->theSampTimesBias[iAtt] = nparm(n++);
						//	}
						//	optStructList[image_index].cbers04Model->updateModel();
						//}

						//check residue is reduced
						NEWMAT::ColumnVector newresidue = getResidue(optStructList);
						double newki2=newresidue.SumSquare();
						double res_reduction = (ki2 - newki2) / (deltap.t()*(deltap*damping + projResidue)).AsScalar();
						//DEBUG TBR
						cout<<sqrt(newki2/nobs)<<" ";
						cout.flush();

						if (res_reduction > 0)
						{
							//accept new parms
							cparm = nparm;
							ki2=newki2;

							deltap_scale = max(1e-15, deltap.NormInfinity()*1e-4);


							buildNormalEquation(optStructList, paramList,
								A, residue, projResidue, deltap_scale);
							olddamping = 0.0;

							found = ( projResidue.NormInfinity() <= minResidue );
							//update damping factor
							damping *= max( 1.0/3.0, 1.0-std::pow((2.0*res_reduction-1.0),3));
							damping_speed = 2.0;
							decrease = true;
						} else {
							//cancel parameter update
							setParameters(nparm, optStructList, paramList);
							for (int i = 0; i < (int)optStructList.size(); i++)
							{
								optStructList[i].cbers04Model->updateModel();
							}
							//int n = 1;
							//for (int image_index = 0;image_index < nImages;++image_index)
							//{
							//	int nAtt = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
							//	for (int iAtt = 0;iAtt < nAtt;++iAtt)
							//	{
							//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].x = nparm(n++);
							//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].y = nparm(n++);
							//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].z = nparm(n++);
							//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].x = nparm(n++);
							//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].y = nparm(n++);
							//		optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].z = nparm(n++);
							//		//optStructList[image_index].cbers04Model->theSupportData->theSampTimesBias[iAtt] = nparm(n++);
							//	}
							//	optStructList[image_index].cbers04Model->updateModel();
							//}

							damping *= damping_speed;
							damping_speed *= 2.0;
						}
					}
				} while (!decrease && !found);
				++iter;
			}
			
			cout<<endl;
			cout<<getParameters(optStructList, paramList).t();
			//cout<<getParameters(optStructList, paramList);
			//for (int image_index = 0;image_index < nImages;++image_index)
			//{
			//	int nAtt = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
			//	for (int iAtt = 0;iAtt < nAtt;++iAtt)
			//	{
			//		cout<<optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].x
			//		<<" "<<optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].y
			//		<<" "<<optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples[iAtt].z
			//		<<" "<<optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].x
			//		<<" "<<optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].y
			//		<<" "<<optStructList[image_index].cbers04Model->theSupportData->theAttitudeVelSamples[iAtt].z
			//		//<<" "<<optStructList[image_index].cbers04Model->theSupportData->theSampTimesBias[iAtt]
			//		<<endl;
			//	}
			//}
			//DEBUG TBR
			cout<<endl;
		}

		
		double getParameter(radiCbers04Model* cbers04Model, int paramIdx, int ephIdx/* = 0*/)
		{
			if (ROLL_OFFSET == paramIdx)
			{
				return cbers04Model->theRollOffset;
			}
			else if (PITCH_OFFSET == paramIdx)
			{
				return cbers04Model->thePitchOffset;
			}
			else if (YAW_OFFSET == paramIdx)
			{
				return cbers04Model->theYawOffset;
			}
			else if (ATT_PITCH_OFFSET == paramIdx)
			{
				return cbers04Model->theSupportData->theAttitudeSamples[ephIdx].x;
			}
			else if (ATT_ROLL_OFFSET == paramIdx)
			{
				return cbers04Model->theSupportData->theAttitudeSamples[ephIdx].y;
			}
			else if (ATT_YAW_OFFSET == paramIdx)
			{
				return cbers04Model->theSupportData->theAttitudeSamples[ephIdx].z;
			}
			else if (ATT_VEL_PITCH_OFFSET == paramIdx)
			{
				return cbers04Model->theSupportData->theAttitudeVelSamples[ephIdx].x;
			}
			else if (ATT_VEL_ROLL_OFFSET == paramIdx)
			{
				return cbers04Model->theSupportData->theAttitudeVelSamples[ephIdx].y;
			}
			else if (ATT_VEL_YAW_OFFSET == paramIdx)
			{
				return cbers04Model->theSupportData->theAttitudeVelSamples[ephIdx].z;
			}

			return 0;
		}

		void setParameter(double v, radiCbers04Model* cbers04Model, int paramIdx, int ephIdx/* = 0*/)
		{
			if (ROLL_OFFSET == paramIdx)
			{
				cbers04Model->theRollOffset = v;
			}
			else if (PITCH_OFFSET == paramIdx)
			{
				cbers04Model->thePitchOffset = v;
			}
			else if (YAW_OFFSET == paramIdx)
			{
				cbers04Model->theYawOffset = v;
			}
			else if (ATT_PITCH_OFFSET == paramIdx)
			{
				cbers04Model->theSupportData->theAttitudeSamples[ephIdx].x = v;
			}
			else if (ATT_ROLL_OFFSET == paramIdx)
			{
				cbers04Model->theSupportData->theAttitudeSamples[ephIdx].y = v;
			}
			else if (ATT_YAW_OFFSET == paramIdx)
			{
				cbers04Model->theSupportData->theAttitudeSamples[ephIdx].z = v;
			}
			else if (ATT_VEL_PITCH_OFFSET == paramIdx)
			{
				cbers04Model->theSupportData->theAttitudeVelSamples[ephIdx].x = v;
			}
			else if (ATT_VEL_ROLL_OFFSET == paramIdx)
			{
				cbers04Model->theSupportData->theAttitudeVelSamples[ephIdx].y = v;
			}
			else if (ATT_VEL_YAW_OFFSET == paramIdx)
			{
				cbers04Model->theSupportData->theAttitudeVelSamples[ephIdx].z = v;
			}
		}

		int getImageParameterNum(radiCbers04Model* cbers04Model, const vector<int>& paramList)
		{
			int np = 0;
			for (int ii = 0;ii < (int)paramList.size();++ii)
			{
				int iParam = paramList[ii];
				if (iParam < EPH_SPLIT)
				{
					// parameters for a scene
					np++;
				}
				else
				{
					// parameters for an ephemeris
					np += (int)cbers04Model->theSupportData->theAttitudeSamples.size();
				}
			}
			return np;
		}

		int getParameterNum(const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList)
		{
			int nImages = (int)optStructList.size();
			int np = 0;
			for (int ii = 0;ii < (int)paramList.size();++ii)
			{
				int iParam = paramList[ii];
				if (iParam < EPH_SPLIT)
				{
					// parameters for a scene
					np += nImages;
				}
				else
				{
					for (int i = 0;i < nImages;++i)
					{
						// parameters for an ephemeris
						np += (int)optStructList[i].cbers04Model->theSupportData->theAttitudeSamples.size();
					}
				}
			}
			return np;
		}


		NEWMAT::ColumnVector getParameters(const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList)
		{
			int nImages = (int)optStructList.size();

			int np = getParameterNum(optStructList, paramList);

			NEWMAT::ColumnVector cparm(np);
			int n = 1;
			for (int image_index = 0;image_index < nImages;++image_index)
			{
				int nEhp = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
				for (int ii = 0;ii < (int)paramList.size();++ii)
				{
					int iParam = paramList[ii];
					if (iParam < EPH_SPLIT)
					{
						// parameters for a scene
						cparm(n++) = getParameter(optStructList[image_index].cbers04Model, iParam);
					}
					else
					{
						// parameters for an ephemeris
						for (int iEhp = 0;iEhp < nEhp;++iEhp)
						{
							cparm(n++) = getParameter(optStructList[image_index].cbers04Model, iParam, iEhp);
						}
					}
				}
			}
			return cparm;
		}

		void setParameters(NEWMAT::ColumnVector nparm, const vector< Cbers04OptStruct >& optStructList, const vector<int>& paramList)
		{
			int nImages = (int)optStructList.size();

			int np = getParameterNum(optStructList, paramList);

			int n = 1;
			for (int image_index = 0;image_index < nImages;++image_index)
			{
				int nEhp = (int)optStructList[image_index].cbers04Model->theSupportData->theAttitudeSamples.size();
				for (int ii = 0;ii < (int)paramList.size();++ii)
				{
					int iParam = paramList[ii];
					if (iParam < EPH_SPLIT)
					{
						// parameters for a scene
						setParameter(nparm(n++), optStructList[image_index].cbers04Model, iParam);
					}
					else
					{
						// parameters for an ephemeris
						for (int iEhp = 0;iEhp < nEhp;++iEhp)
						{
							setParameter(nparm(n++), optStructList[image_index].cbers04Model, iParam, iEhp);
						}
					}
				}
			}
		}
	}
}