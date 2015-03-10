#include "innerOptimize.h"

namespace Cbers04Optimize{
	namespace innerOptimize
	{
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

void buildNormalEquation(const vector< Cbers04OptStruct >& optStructList, 
			const vector<int>& parameterList,
	NEWMAT::SymmetricMatrix& A,
	NEWMAT::ColumnVector& residue,
	NEWMAT::ColumnVector& projResidue,
	double pstep_scale)
{
	//init
	int nImages = (int)optStructList.size();
	int np = (int)parameterList.size();
	int nTotalEquation = 0;
	int dimObs = 3;
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
	std::vector<ossimDpt> imDerp(np);
	// ground observations
	std::vector<ossimGpt>  gdDerp(np);
	ossimDpt resIm;
	ossimGpt gd, resGd;
	for(int i = 0;i < nImages;++i)
	{
		const vector<ossimRefPtr<ossimTieGpt> >& theTPV = optStructList[i].gcpSet->getTiePoints();
		vector<ossimRefPtr<ossimTieGpt> >::const_iterator tit;
		// loop on tie points
		for (tit = theTPV.begin() ; tit != theTPV.end() ; ++tit)
		{
			//// forward
			//// compute residue
			//resIm = (*tit)->tie - optStructList[i].cbers04Model->forward(*(*tit));
			////optStructList[i].cbers04Model->worldToLineSample(*(*tit), resIm);
			////resIm = (*tit)->tie - resIm;
			////optStructList[i].cbers04Model->worldToLineSample(*(*tit), resIm);
			////resIm = (*tit)->tie - resIm;
			//residue(c++) = resIm.x;
			//residue(c++) = resIm.y;

			////compute all image derivatives regarding parametres for the tie point position
			//for(int p=0;p<np;p++)
			//{
			//	imDerp[p] = optStructList[i].cbers04Model->getForwardDeriv( parameterList[p] , *(*tit) , pstep_scale);
			//}

			////compute influence of tie point on all sytem elements
			//for(int p1=0;p1<np;++p1)
			//{        
			//	//proj residue: J * residue
			//	projResidue.element(p1) += imDerp[p1].x * resIm.x + imDerp[p1].y * resIm.y;

			//	//normal matrix A = transpose(J)*J
			//	for(int p2=p1;p2<np;++p2)
			//	{
			//		A.element(p1,p2) += imDerp[p1].x * imDerp[p2].x + imDerp[p1].y * imDerp[p2].y;
			//	}
			//}
			

			// inverse
			//compute residue
			gd = optStructList[i].cbers04Model->inverse((*tit)->tie);
			residue(c++) = resGd.lon = ((*tit)->lon - gd.lon) * 100000.0;
			residue(c++) = resGd.lat = ((*tit)->lat - gd.lat) * 100000.0 * cos(gd.lat / 180.0 * M_PI);
			residue(c++) = resGd.hgt = (*tit)->hgt - gd.hgt; //TBD : normalize to meters?

			//compute all image derivatives regarding parametres for the tie point position
			for (int p = 0; p<np; p++)
			{
				gdDerp[p] = optStructList[i].cbers04Model->getInverseDeriv(parameterList[p], (*tit)->tie, pstep_scale);
			}

			//compute influence of tie point on all sytem elements
			for (int p1 = 0; p1<np; ++p1)
			{
				//proj residue: J * residue
				projResidue.element(p1) += gdDerp[p1].lon * resGd.lon + gdDerp[p1].lat * resGd.lat + gdDerp[p1].hgt * resGd.hgt;

				//normal matrix A = transpose(J)*J
				for (int p2 = p1; p2<np; ++p2)
				{
					A.element(p1, p2) += gdDerp[p1].lon * gdDerp[p2].lon + gdDerp[p1].lat * gdDerp[p2].lat + gdDerp[p1].hgt * gdDerp[p2].hgt;
				}
			}
		}
	}	
}

void adjustment(const vector< Cbers04OptStruct >& optStructList, const vector<int>& parameterList)
{
	int nImages = (int)optStructList.size();
	int nTotalParam = (int)parameterList.size();
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

	int np = (int)parameterList.size();
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
	buildNormalEquation(optStructList, parameterList,
		A, residue, projResidue, deltap_scale);
	double ki2=residue.SumSquare();

	//get current adjustment (between -1 and 1 normally) and convert to ColumnVector
	ossimAdjustmentInfo cadj;
	optStructList[0].cbers04Model->getAdjustment(cadj);
	std::vector< ossimAdjustableParameterInfo >& parmlist = cadj.getParameterList();
	NEWMAT::ColumnVector cparm(np), nparm(np);
	for(int n=0;n<np;++n)
	{
		cparm(n+1) = parmlist[parameterList[n]].getParameter();
	}

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
				for (int image_index = 0;image_index < nImages;++image_index)
				{
					for(int n=0;n<np;++n)
					{
						optStructList[image_index].cbers04Model->setAdjustableParameter(
							parameterList[n], nparm(n+1), false); //do not update now, wait
					}
					optStructList[image_index].cbers04Model->updateModel();
				}

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


					buildNormalEquation(optStructList, parameterList,
						A, residue, projResidue, deltap_scale);
					olddamping = 0.0;

					found = ( projResidue.NormInfinity() <= minResidue );
					//update damping factor
					damping *= max( 1.0/3.0, 1.0-std::pow((2.0*res_reduction-1.0),3));
					damping_speed = 2.0;
					decrease = true;
				} else {
					//cancel parameter update
					for (int image_index = 0;image_index < nImages;++image_index)
					{
						for(int n=0;n<np;++n)
						{
							optStructList[image_index].cbers04Model->setAdjustableParameter(
								parameterList[n], nparm(n+1), false); //do not update now, wait
						}
						optStructList[image_index].cbers04Model->updateModel();
					}

					damping *= damping_speed;
					damping_speed *= 2.0;
				}
			}
		} while (!decrease && !found);
		++iter;
	}

	//DEBUG TBR
	cout<<endl;
}

	}
}