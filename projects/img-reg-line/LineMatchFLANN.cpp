#include "lineFeature.h"
#include "cv.h"
#include <stdio.h>
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <math.h>
#include "func.h"
#include "quick_selection.h"
#include "LineMatchFLANN.h"
using namespace cv;

cv::Mat affine_estimation(const std::vector<cvline_polar>& src_line_polar, const std::vector<cvline_polar>& ref_line_polar)
{
	int nLines = (int)src_line_polar.size();
	if (nLines != (int)ref_line_polar.size())
	{
		return cv::Mat(2, 3, CV_64F);
	}

	double bb[6] = {0.0};
	double AA[6][6] = {0.0};
	//b *= 0.0;
	//A *= 0.0;

	int c = 0;
	for(int i = 0;i < nLines;++i)
	{
		double x1 = src_line_polar[i].pt1.x;
		double y1 = src_line_polar[i].pt1.y;
		double x2 = src_line_polar[i].pt2.x;
		double y2 = src_line_polar[i].pt2.y;

		double sin_theta = sin(ref_line_polar[i].theta);
		double cos_theta = cos(ref_line_polar[i].theta);
		double rho = ref_line_polar[i].rho;
		double len = segment_length(ref_line_polar[i]);

		vector<double> coeff(6);
		coeff[0] = -sin_theta;
		coeff[1] = -sin_theta * x1;
		coeff[2] = -sin_theta * y1;
		coeff[3] = cos_theta;
		coeff[4] = cos_theta * x1;
		coeff[5] = cos_theta * y1;

		for(int p1=0;p1<6;++p1)
		{        
			bb[p1] += coeff[p1] * rho;
			for(int p2=0;p2<6;++p2)
			{
				AA[p1][p2] += coeff[p1] * coeff[p2];
			}
		}

		coeff[0] = -sin_theta;
		coeff[1] = -sin_theta * x2;
		coeff[2] = -sin_theta * y2;
		coeff[3] = cos_theta;
		coeff[4] = cos_theta * x2;
		coeff[5] = cos_theta * y2;
		for(int p1=0;p1<6;++p1)
		{        
			bb[p1] += coeff[p1] * rho;
			for(int p2=0;p2<6;++p2)
			{
				AA[p1][p2] += coeff[p1] * coeff[p2];
			}
		}
	}

	cv::Mat b(6, 1, CV_64F, bb);
	cv::Mat A(6, 6, CV_64F, AA);
	cv::Mat paramTemp = A.inv()*b;
	cv::Mat parameter(2, 3, CV_64F);
	double a0 = parameter.at<double>(0, 2) = paramTemp.at<double>(0, 0);
	double a1 = parameter.at<double>(0, 0) = paramTemp.at<double>(1, 0);
	double a2 = parameter.at<double>(0, 1) = paramTemp.at<double>(2, 0);
	double b0 = parameter.at<double>(1, 2) = paramTemp.at<double>(3, 0);
	double b1 = parameter.at<double>(1, 0) = paramTemp.at<double>(4, 0);
	double b2 = parameter.at<double>(1, 1) = paramTemp.at<double>(5, 0);
	return parameter;
}

vector<double> affine_estimation2(const std::vector<cvline_polar>& src_line_polar,
	const std::vector<cvline_polar>& ref_line_polar)
{
	int nLines = (int)src_line_polar.size();
	if (nLines != (int)ref_line_polar.size())
	{
		return cv::Mat(2, 3, CV_64F);
	}

	double bb[6] = { 0.0 };
	double AA[6][6] = { 0.0 };
	//b *= 0.0;
	//A *= 0.0;

	int c = 0;
	for (int i = 0; i < nLines; ++i)
	{
		double x1 = src_line_polar[i].pt1.x;
		double y1 = src_line_polar[i].pt1.y;
		double x2 = src_line_polar[i].pt2.x;
		double y2 = src_line_polar[i].pt2.y;

		double sin_theta = sin(ref_line_polar[i].theta);
		double cos_theta = cos(ref_line_polar[i].theta);
		double rho = ref_line_polar[i].rho;
		double len = segment_length(ref_line_polar[i]);

		vector<double> coeff(6);
		coeff[0] = -sin_theta;
		coeff[1] = -sin_theta * x1;
		coeff[2] = -sin_theta * y1;
		coeff[3] = cos_theta;
		coeff[4] = cos_theta * x1;
		coeff[5] = cos_theta * y1;

		for (int p1 = 0; p1<6; ++p1)
		{
			bb[p1] += coeff[p1] * rho;
			for (int p2 = 0; p2<6; ++p2)
			{
				AA[p1][p2] += coeff[p1] * coeff[p2];
			}
		}

		coeff[0] = -sin_theta;
		coeff[1] = -sin_theta * x2;
		coeff[2] = -sin_theta * y2;
		coeff[3] = cos_theta;
		coeff[4] = cos_theta * x2;
		coeff[5] = cos_theta * y2;
		for (int p1 = 0; p1<6; ++p1)
		{
			bb[p1] += coeff[p1] * rho;
			for (int p2 = 0; p2<6; ++p2)
			{
				AA[p1][p2] += coeff[p1] * coeff[p2];
			}
		}
	}

	cv::Mat b(6, 1, CV_64F, bb);
	cv::Mat A(6, 6, CV_64F, AA);
	cv::Mat paramTemp = A.inv()*b;
	//cv::Mat parameter(2, 3, CV_64F);
	vector<double> parameter(6);
	double a0 = parameter[0] = paramTemp.at<double>(0, 0);
	double a1 = parameter[1] = paramTemp.at<double>(1, 0);
	double a2 = parameter[2] = paramTemp.at<double>(2, 0);
	double b0 = parameter[3] = paramTemp.at<double>(3, 0);
	double b1 = parameter[4] = paramTemp.at<double>(4, 0);
	double b2 = parameter[5] = paramTemp.at<double>(5, 0);
	return parameter;
}


vector<double> affine_estimation2(const std::vector<cvline_polar>& src_line_polar,
	const std::vector<cvline_polar>& ref_line_polar, const vector<int>& slaveAssign)
{
	int nLines = (int)src_line_polar.size();
	if (nLines != (int)slaveAssign.size())
	{
		return cv::Mat(2, 3, CV_64F);
	}

	double bb[6] = { 0.0 };
	double AA[6][6] = { 0.0 };
	//b *= 0.0;
	//A *= 0.0;

	int c = 0;
	for (int i = 0; i < nLines; ++i)
	{
		double x1 = src_line_polar[i].pt1.x;
		double y1 = src_line_polar[i].pt1.y;
		double x2 = src_line_polar[i].pt2.x;
		double y2 = src_line_polar[i].pt2.y;

		double sin_theta = sin(ref_line_polar[slaveAssign[i]].theta);
		double cos_theta = cos(ref_line_polar[slaveAssign[i]].theta);
		double rho = ref_line_polar[slaveAssign[i]].rho;
		double len = segment_length(ref_line_polar[slaveAssign[i]]);

		vector<double> coeff(6);
		coeff[0] = -sin_theta;
		coeff[1] = -sin_theta * x1;
		coeff[2] = -sin_theta * y1;
		coeff[3] = cos_theta;
		coeff[4] = cos_theta * x1;
		coeff[5] = cos_theta * y1;

		for (int p1 = 0; p1<6; ++p1)
		{
			bb[p1] += coeff[p1] * rho;
			for (int p2 = 0; p2<6; ++p2)
			{
				AA[p1][p2] += coeff[p1] * coeff[p2];
			}
		}

		coeff[0] = -sin_theta;
		coeff[1] = -sin_theta * x2;
		coeff[2] = -sin_theta * y2;
		coeff[3] = cos_theta;
		coeff[4] = cos_theta * x2;
		coeff[5] = cos_theta * y2;
		for (int p1 = 0; p1<6; ++p1)
		{
			bb[p1] += coeff[p1] * rho;
			for (int p2 = 0; p2<6; ++p2)
			{
				AA[p1][p2] += coeff[p1] * coeff[p2];
			}
		}
	}

	cv::Mat b(6, 1, CV_64F, bb);
	cv::Mat A(6, 6, CV_64F, AA);
	cv::Mat paramTemp = A.inv()*b;
	//cv::Mat parameter(2, 3, CV_64F);
	vector<double> parameter(6);
	double a0 = parameter[0] = paramTemp.at<double>(0, 0);
	double a1 = parameter[1] = paramTemp.at<double>(1, 0);
	double a2 = parameter[2] = paramTemp.at<double>(2, 0);
	double b0 = parameter[3] = paramTemp.at<double>(3, 0);
	double b1 = parameter[4] = paramTemp.at<double>(4, 0);
	double b2 = parameter[5] = paramTemp.at<double>(5, 0);
	return parameter;
}

cvline_polar affineTransLine(cv::Mat affineTrans, const cvline_polar& inLine)
{
	double a0 = affineTrans.at<double>(0, 2);
	double a1 = affineTrans.at<double>(0, 0);
	double a2 = affineTrans.at<double>(0, 1);
	double b0 = affineTrans.at<double>(1, 2);
	double b1 = affineTrans.at<double>(1, 0);
	double b2 = affineTrans.at<double>(1, 1);
	cvline_polar outLine;
	fPoint line[2];
	line[0].x = a0 + a1 * inLine.pt1.x + a2 * inLine.pt1.y;
	line[0].y = b0 + b1 * inLine.pt1.x + b2 * inLine.pt1.y;
	line[1].x = a0 + a1 * inLine.pt2.x + a2 * inLine.pt2.y;
	line[1].y = b0 + b1 * inLine.pt2.x + b2 * inLine.pt2.y;
	outLine = line2polar(line);
	return outLine;
}

cvline_polar affineTransLine(vector<double> affineTrans, const cvline_polar& inLine)
{
	double a0 = affineTrans[0];
	double a1 = affineTrans[1];
	double a2 = affineTrans[2];
	double b0 = affineTrans[3];
	double b1 = affineTrans[4];
	double b2 = affineTrans[5];
	cvline_polar outLine;
	fPoint line[2];
	line[0].x = a0 + a1 * inLine.pt1.x + a2 * inLine.pt1.y;
	line[0].y = b0 + b1 * inLine.pt1.x + b2 * inLine.pt1.y;
	line[1].x = a0 + a1 * inLine.pt2.x + a2 * inLine.pt2.y;
	line[1].y = b0 + b1 * inLine.pt2.x + b2 * inLine.pt2.y;
	outLine = line2polar(line);
	return outLine;
}

fPoint dist2(const cvline_polar& l1, const cvline_polar& l2)
{
	fPoint dist;
	// dist1
	double d1 = -l1.pt1.x * sin(l2.theta) + l1.pt1.y * cos(l2.theta) - l2.rho;
	double d2 = -l1.pt2.x * sin(l2.theta) + l1.pt2.y * cos(l2.theta) - l2.rho;
	dist.x = sqrt(d1*d1+d2*d2);
	// dist2
	d1 = -l2.pt1.x * sin(l1.theta) + l2.pt1.y * cos(l1.theta) - l1.rho;
	d2 = -l2.pt2.x * sin(l1.theta) + l2.pt2.y * cos(l1.theta) - l1.rho;
	dist.y = sqrt(d1*d1+d2*d2);
	return dist;
}

double center_dist(const cvline_polar& l1, const cvline_polar& l2)
{
	fPoint ct1;
	ct1.x = (l1.pt1.x+l1.pt2.x)*0.5;
	ct1.y = (l1.pt1.y+l1.pt2.y)*0.5;
	fPoint ct2;
	ct2.x = (l2.pt1.x+l2.pt2.x)*0.5;
	ct2.y = (l2.pt1.y+l2.pt2.y)*0.5;
	return sqrt((ct1.x-ct2.x)*(ct1.x-ct2.x)+(ct1.y-ct2.y)*(ct1.y-ct2.y));
}



double point2polarline(const fPoint& p, const cvline_polar& l)
{
	double dist = -p.x * sin(l.theta) + p.y * cos(l.theta) - l.rho;
	return dist;
}


double distance1(const cvline_polar& l1, const cvline_polar& l2)
{
	//return (point2polarline(l1.pt1, l2) + point2polarline(l1.pt2, l2))*0.5;
	double d1 = point2polarline(l1.pt1, l2);
	double d2 = point2polarline(l1.pt2, l2);
	return sqrt(d1*d1 + d2*d2);
}

fPoint distance2(const cvline_polar& l1, const cvline_polar& l2)
{

	double centerx1 = (l1.pt1.x + l1.pt2.x)*0.5;
	double centery1 = (l1.pt1.y + l1.pt2.y)*0.5;
	double centerx2 = (l2.pt1.x + l2.pt2.x)*0.5;
	double centery2 = (l2.pt1.y + l2.pt2.y)*0.5;
	//fPoint shift = fPoint(centerx1-centerx2, centery1-centery2);
	double shift = sqrt((centerx1 - centerx2)*(centerx1 - centerx2) + (centery1 - centery2)*(centery1 - centery2));
	//double shift_weight = 1.0E-2;
	//double shift_weight = 0.08;

	fPoint dist;
	// dist1
	dist.x = distance1(l1, l2);
	// dist2
	dist.y = distance1(l2, l1);

	//shift = max(dist.x,dist.y)*exp(-1.0/(shift*shift+DBL_EPSILON));
	//shift = 3.0*exp(-m_delta_2/(shift*shift+DBL_EPSILON));
	shift = 0.0*exp(-1.0 / (shift + DBL_EPSILON));
	//if((dist.x*dist.x+dist.y*dist.y) < 4.0)
	//{
	dist.x += shift;
	dist.y += shift;
	//}

	//dist.x = fabs(l1.rho - l2.rho);
	//dist.y = scale_angle * fabs(l1.theta - l2.theta);
	////dist.y = dist.y > lammda_rho*dist.x ? lammda_rho*dist.x : dist.y;
	////dist.y = fabs(tan(l1.theta) - tan(l2.theta));
	//return sqrt(dist.x*dist.x + dist.y*dist.y);
	return dist;
}

cvline_polar forward(const cvline_polar& l, const vector<double>& parameters)
{
	cvline_polar outLine;
	fPoint line[2];
	line[0].x = parameters[0] + parameters[1] * l.pt1.x + parameters[2] * l.pt1.y;
	line[0].y = parameters[3] + parameters[4] * l.pt1.x + parameters[5] * l.pt1.y;
	line[1].x = parameters[0] + parameters[1] * l.pt2.x + parameters[2] * l.pt2.y;
	line[1].y = parameters[3] + parameters[4] * l.pt2.x + parameters[5] * l.pt2.y;
	outLine = line2polar(line);
	return outLine;
}


void levmar_function_fvec(double *param, double *hx, int nparameter, int nequation, void *adata)
{
	dataStruct *pThis = (dataStruct*)adata;
	int nX = (int)pThis->slave_vec_lines.size();
	int nY = (int)pThis->master_vec_lines.size();

	int pos = 0;
	int i;
	for (i = 0; i<nparameter; ++i)
	{
		pThis->parameters[i] = param[i];
	}

	int c = 0;

	for (int i = 0; i < nX; ++i)
	{
		cvline_polar outPt = forward(pThis->slave_vec_lines[i], pThis->parameters);

		fPoint dist = distance2(outPt, pThis->master_vec_lines[pThis->slaveAssign[i]]);

		hx[c++] = dist.x;
		hx[c++] = dist.y;
	}
}

void levmar_function_jac(double *param, double *jac, int nparameter, int nequation, void *adata)
{
	dataStruct *pThis = (dataStruct*)adata;
	int nX = (int)pThis->slave_vec_lines.size();
	int nY = (int)pThis->master_vec_lines.size();

	int pos = 0;
	int i;
	for (i = 0; i<nparameter; ++i)
	{
		pThis->parameters[i] = param[i];
	}

	double pstep_scale = 1e-4;
	static double den = 0.5 / pstep_scale;
	int c = 0;

	for (int i = 0; i < nX; ++i)
	{
		cvline_polar outPt = forward(pThis->slave_vec_lines[i], pThis->parameters);

		for (int p = 0; p<nparameter; ++p)
		{
			double middle = pThis->parameters[p];
			pThis->parameters[p] = middle + pstep_scale;
			cvline_polar outLine1 = forward(pThis->slave_vec_lines[i], pThis->parameters);
			//fPoint dist1 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine1);
			fPoint dist1 = distance2(outLine1, pThis->master_vec_lines[pThis->slaveAssign[i]]);

			pThis->parameters[p] = middle - pstep_scale;
			cvline_polar outLine2 = forward(pThis->slave_vec_lines[i], pThis->parameters);
			//fPoint dist2 = pThis->distance2(pThis->m_observe_lines_polar[j], outLine2);
			fPoint dist2 = distance2(outLine2, pThis->master_vec_lines[pThis->slaveAssign[i]]);

			pThis->parameters[p] = middle;

			double derivative_x = (dist1.x - dist2.x) * den;
			double derivative_y = (dist1.y - dist2.y) * den;

			jac[c*nparameter + p] = derivative_x;
			jac[(c + 1)*nparameter + p] = derivative_y;
		}
		c += 2;
	}
}


bool parameters_optimization(const vector<cvline_polar>& slave_vec_lines, const vector<cvline_polar>& master_vec_lines,
	const vector<int>& slaveAssign, std::vector<double>& parameters)
{
	int nX = (int)slave_vec_lines.size();
	int nY = (int)master_vec_lines.size();

	VectorXd b(6);
	MatrixXd A(6, 6);
	b.fill(0.0);
	A.fill(0.0);

	int nparam = 6;
	std::vector<double> cparm(nparam);

	for (int i = 0; i<nparam; ++i) cparm[i] = parameters[i];

	// lm
	double *p = &cparm[0];

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0] = LM_INIT_MU; opts[1] = 1E-15; opts[2] = 1E-15; opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	dataStruct ds;
	ds.slave_vec_lines = slave_vec_lines;
	ds.master_vec_lines = master_vec_lines;
	ds.parameters = parameters;
	ds.slaveAssign = slaveAssign;
	//int ret = dlevmar_dif(levmar_function_fvec, &cparm[0], &x[0], nparam, 2*nX, 1000, opts, info, NULL, NULL, this);  // no Jacobian
	int ret = dlevmar_der(levmar_function_fvec, levmar_function_jac, p, NULL, nparam, 2 * nX, 1000, opts, info, NULL, NULL, &ds); // with analytic Jacobian
	//delete []x;

	for (int i = 0; i<nparam; ++i) parameters[i] = p[i];
	return true;
}


bool parameters_optimization(const vector<cvline_polar>& slave_vec_lines, const vector<cvline_polar>& master_vec_lines, std::vector<double>& parameters)
{
	int nX = (int)slave_vec_lines.size();
	vector<int> slaveAssign;
	for (int i = 0; i < nX; i++)
	{
		slaveAssign.push_back(i);
	}

	return parameters_optimization(slave_vec_lines, master_vec_lines, slaveAssign, parameters);
}


void LineMatch(const cv::Mat& ref_seg, const cv::Mat& src_seg,
	vector< vector<double> > &ref_LineFeature, vector< vector<double> > &src_LineFeature,
	vector<int>&ref_LineNum, vector<int>&src_LineNum, 
	vector<cvline_polar>&ref_lines_list, vector<cvline_polar>&src_lines_list,
	vector<cvline_polar>& matched_ref_ls, vector<cvline_polar>& matched_src_ls)
{
	int N_ref=ref_LineFeature.size();//参考影像不靠近边缘的直线个数
	int N_src=src_LineFeature.size();//源影像不靠近边缘的直线个数
	if (N_ref < 1 || N_src < 1) return;
	int M=ref_LineFeature[0].size();//描述子维数
	//cout<<N_ref<<"\t"<<N_src<<endl;
	int i,j;
	FlannBasedMatcher matcher;
	std::vector< DMatch > matches;
	//CvMat* descriptors_ref=cvCreateMat(N_ref,M, CV_32FC1);
	//CvMat* descriptors_src=cvCreateMat(N_src,M, CV_32FC1);
	cv::Mat descriptors_ref(N_ref, M, CV_32FC1);
	cv::Mat descriptors_src(N_src, M, CV_32FC1);
	for(i=0;i!=N_ref;i++)
	{
		for(j=0;j!=M;j++)
		{
			//cvmSet(descriptors_ref,i,j, ref_LineFeature[i][j]);
			descriptors_ref.at<float>(i, j) = (float)ref_LineFeature[i][j];
		}
	}

	for(i=0;i!=N_src;i++)
	{
		for(j=0;j!=M;j++)
		{
			//cvmSet(descriptors_src, i, j, src_LineFeature[i][j]);
			descriptors_src.at<float>(i, j) = (float)src_LineFeature[i][j];
		}
	}
	matcher.match(descriptors_src,  descriptors_ref, matches );
	int nMatches = (int)matches.size();
	std::vector<double> distList(nMatches);
	for (int i = 0;i <nMatches;i++)
	{
		distList[i] = matches[i].distance;
	}
	int selection_num = 3;
	double min_dist = quick_select(&distList[0], 0, nMatches-1, selection_num);
	vector<int> indexList;
	for (int i = 0;i < nMatches;++i)
	{
		if (matches[i].distance <= min_dist)
		{
			indexList.push_back(i);
		}
	}
	std::vector<cvline_polar> srcGoodList;
	std::vector<cvline_polar> refGoodList;
	for (int i = 0;i<(int)indexList.size();++i)
	{
		int newRefIndex = ref_LineNum[matches[indexList[i]].trainIdx];
		refGoodList.push_back(ref_lines_list[newRefIndex]);
		int newSrcIndex = src_LineNum[matches[indexList[i]].queryIdx];
		srcGoodList.push_back(src_lines_list[newSrcIndex]);

		//matched_src_ls.push_back(src_lines_list[newSrcIndex]);
		//matched_ref_ls.push_back(ref_lines_list[newRefIndex]);
	}
	//return;
	cv::Mat affine_parameters = affine_estimation(srcGoodList, refGoodList);
	//std::vector<double> affine_parameters(6);
	//affine_parameters[0] = 25.0;
	//affine_parameters[1] = 1.0;
	//affine_parameters[2] = 0.0;
	//affine_parameters[3] = 25.0;
	//affine_parameters[4] = 0.0;
	//affine_parameters[5] = 1.0;
	//parameters_optimization(srcGoodList, refGoodList, affine_parameters);

	std::vector<cvline_polar> trans_line_list;
	for (int i  = 0;i <(int)src_lines_list.size();++i)
	{
		trans_line_list.push_back(affineTransLine(affine_parameters, src_lines_list[i]));
	}

	std::vector< DMatch > good_matches;
	double distance_threshold =10.0;
	double center_distance_threshold = 10;
	std::vector<int> transOpenList((int)src_lines_list.size());
	for (int i = 0;i < (int)src_lines_list.size();++i)
	{
		transOpenList[i] = i;
	}
	for ( i = 0; i < (int)ref_lines_list.size(); i++ )
	{
		int ref_i=i;

		int minIndex = 0;
		double minDistance = DBL_MAX;
		for (int j = 0;j < (int)transOpenList.size() && j >=0;++j)
		{
			int trans_j=transOpenList[j];
			fPoint dist_pt = dist2(ref_lines_list[ref_i], trans_line_list[trans_j]);
			double dist = sqrt(dist_pt.x*dist_pt.x+dist_pt.y*dist_pt.y);

			if (dist < minDistance)
			{
				if (center_dist(trans_line_list[trans_j], ref_lines_list[ref_i]) <= center_distance_threshold)
				{
					minDistance = dist;
					minIndex = j;
					continue;
				}
			}
		}
		if (minDistance < distance_threshold)
		{
			DMatch matc;
			matc.queryIdx = transOpenList[minIndex];
			matc.trainIdx = i;
			//if (center_dist(trans_line_list[matc.queryIdx], ref_lines_list[matc.trainIdx]) > center_distance_threshold)
			//{
			//	continue;
			//}
			good_matches.push_back(matc);
			transOpenList.erase(transOpenList.begin()+minIndex);
		}
	}
	cout<<"total matches:"<<good_matches.size()<<endl;

	matched_src_ls.clear();
	matched_ref_ls.clear();	
	for(i = 0; i < good_matches.size(); i++)
	{
		int ref_index = good_matches[i].trainIdx;
		int src_index = good_matches[i].queryIdx;

		matched_src_ls.push_back(src_lines_list[src_index]);
		matched_ref_ls.push_back(ref_lines_list[ref_index]);

		//char temp[64];
		//sprintf_s(temp, "%d", i+1);
		//string strName(temp);
		//cv::Scalar lineColor = cv::Scalar(0, 0, 255);
		//int nThickness = 1;
		//cv::line(ref_seg, cv::Point(ref_lines_list[ref_index].pt1.x,ref_lines_list[ref_index].pt1.y),
		//	cv::Point(ref_lines_list[ref_index].pt2.x,ref_lines_list[ref_index].pt2.y),lineColor,nThickness, CV_AA);
		//_drawLabel( ref_seg, ref_lines_list[ref_index], cv::Scalar::all(-1), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, strName);
		//cv::line(src_seg, cv::Point(src_lines_list[src_index].pt1.x,src_lines_list[src_index].pt1.y),
		//	cv::Point(src_lines_list[src_index].pt2.x,src_lines_list[src_index].pt2.y),lineColor,nThickness, CV_AA);
		//_drawLabel( src_seg, src_lines_list[src_index], cv::Scalar::all(-1), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS, strName);
	}
	////clock_t clockEnd;
	////clockEnd = clock();
	////cout<<"time consuming:"<<(clockEnd - clockBegin)*1e-3<<"second"<<endl;
	//clock_t clockEnd2;
	//clockEnd2 = clock();
	//cout<<"Line match time consuming:"<<(clockEnd2 - clockBegin2)*1e-3<<"second"<<endl;
	//cv::imshow( "ref", ref_seg);
	//cv::imshow( "src", src_seg);

	//string RefImageOut3 = "data\\out_3.png";//用于匹配后图像
	//string SrcImageOut4 = "data\\out_4.png";
	//cv::imwrite(RefImageOut3.c_str(), ref_seg);
	//cv::imwrite(SrcImageOut4.c_str(), src_seg);
	//cvWaitKey();
	//cvDestroyAllWindows(); 


	//int n = (int)m_inliers.size() / 2;
}