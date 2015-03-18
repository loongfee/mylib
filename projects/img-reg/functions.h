#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <ossim/base/ossimString.h>
#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimOutputSource.h>
#include <ossim/base/ossimProcessInterface.h>
#include <ossim/base/ossimProcessProgressEvent.h>
#include <ossim/imaging/ossimFilterResampler.h>
#include <ossim/imaging/ossimImageChain.h>
#include <ossim/imaging/ossimImageHandler.h>
#include <ossim/imaging/ossimBandSelector.h>
#include <ossim/imaging/ossimImageRenderer.h>
#include <ossim/imaging/ossimCastTileSourceFilter.h>
#include <ossim/base/ossimTDpt.h>
#include <ossim_plugin/radi/radiBlockTieGpt.h>
#include <ossim_plugin/radi/radiBlockTieGptSet.h>
#include <vector>
//#include <mpi.h>
//#include <ossim/parallel/ossimMpi.h>

#include <ossim/elevation/ossimElevManager.h>

#include <ogrsf_frmts.h>
#include <gdal.h>
#include <ogr_api.h>
#include <ogr_geometry.h>
#include <ogr_feature.h>

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/xfeatures2d/nonfree.hpp"
#include <algorithm>
#include <GdalRasterApp.h>
#include <boost/filesystem.hpp>
#include <curl/curl.h>

#include <ez_ransac.h>

#include <levmar.h>

extern "C"{
#include <vl/generic.h>
#include <vl/stringop.h>
#include <vl/pgm.h>
#include <vl/sift.h>
#include <vl/getopt_long.h>
#include <vl/covdet.h>
};

using namespace mylib;
using namespace std;
using namespace cv;
using namespace ossimplugins;


struct row_col{
	double row_idx;
	double col_idx;
	row_col(double r, double c)
	{
		row_idx = r;
		col_idx = c;
	}
};

enum block_position{
	left = 0,
	right,
	top,
	bottom,
	left_top,
	left_bottom,
	right_top,
	right_bottom,
	center,
};

struct image_block
{
	ossimIrect rect;
	block_position position;
};


static int row_col_step = 7;
// multiple keywords
static bool RowColCompare(row_col rc1, row_col rc2)
{
	return (fabs(rc1.row_idx) + fabs(rc1.col_idx)) < (fabs(rc2.row_idx) + fabs(rc2.col_idx));
	//int fd1 = fabs(rc1.row_idx)+fabs(rc1.col_idx);
	//int fd2 = fabs(rc2.row_idx)+fabs(rc2.col_idx);
	//int d1 = (int)(fd1 + 0.5);
	//int d2 = (int)(fd2 + 0.5);
	//int m1 = d1%row_col_step;
	//int m2 = d2 %row_col_step;
	//if (m1 == m2)
	//{
	//	return fd1 < fd2;
	//}
	//return m1 < m2;
}

static vector<row_col> assign_row_col_list(int tilerows, int tilecols, block_position position = block_position::center)
{
	vector<row_col> row_col_List;
	double center_row = (tilerows - 1) * 0.5;
	double center_col = (tilecols - 1) * 0.5;
	double start_row = center_row;
	double start_col = center_col;

	//bool row_first = true;

	//if (position == block_position::left
	//	|| position == block_position::left_bottom
	//	|| position == block_position::left_top)
	//{
	//	start_col = 0.0;
	//	if (position == block_position::left)
	//	{
	//		row_first = false;
	//	}
	//}
	//else if (position == block_position::right
	//	|| position == block_position::right_bottom
	//	|| position == block_position::right_top)
	//{
	//	start_col = tilecols;
	//	if (position == block_position::right)
	//	{
	//		row_first = false;
	//	}
	//}
	//else
	//{
	//	start_col = center_col;
	//}

	//if (position == block_position::top
	//	|| position == block_position::left_top
	//	|| position == block_position::right_top)
	//{
	//	start_row = 0;
	//	if (position == block_position::top)
	//	{
	//		row_first = true;
	//	}
	//}
	//else if (position == block_position::bottom
	//	|| position == block_position::left_bottom
	//	|| position == block_position::right_bottom)
	//{
	//	start_row = tilerows;
	//	if (position == block_position::bottom)
	//	{
	//		row_first = true;
	//	}
	//}
	//else
	//{
	//	start_row = center_row;
	//}

	if (position == block_position::left
		|| position == block_position::left_bottom
		|| position == block_position::left_top)
	{
		start_col = center_col;
		start_col = 0.0;
	}
	else if (position == block_position::right
		|| position == block_position::right_bottom
		|| position == block_position::right_top)
	{
		start_col = tilecols;
	}
	else
	{
		start_col = center_col;
		start_col = 0.0;
	}

	if (position == block_position::top
		|| position == block_position::left_top
		|| position == block_position::right_top)
	{
		start_row = center_row;
		start_row = 0;
	}
	else if (position == block_position::bottom
		|| position == block_position::left_bottom
		|| position == block_position::right_bottom)
	{
		start_row = tilerows;
	}
	else
	{
		start_row = center_row;
		start_row = 0;
	}

	//double weight = 2.0; 
	//if (row_first)
	//{
	//	weight = tilecols*0.5;
	//}
	//else
	//{
	//	weight = tilerows*0.5;

	//}
	for (int i = 0; (i<tilerows); ++i)
	{
		for (int j = 0; (j<tilecols); ++j)
		{
			row_col_List.push_back(row_col(i - start_row, j - start_col));
			//if (row_first)
			//{
			//	row_col_List.push_back(row_col((i - start_row)*weight, j - start_col));
			//}
			//else
			//{
			//	row_col_List.push_back(row_col((i - start_row)*weight, (j - start_col)*weight));

			//}
		}
	}
	std::sort(row_col_List.begin(), row_col_List.end(), RowColCompare);
	for (size_t i = 0; i < row_col_List.size(); i++)
	{
		//if (row_first)
		//{
		//	row_col_List[i].row_idx /= weight;
		//}
		//else
		//{
		//	row_col_List[i].col_idx /= weight;
		//}
		row_col_List[i].col_idx = int(row_col_List[i].col_idx + start_col);
		row_col_List[i].row_idx = int(row_col_List[i].row_idx + start_row);
		if (row_col_List[i].col_idx < 0)
		{
			row_col_List[i].col_idx = 0;
		}
		else if (row_col_List[i].col_idx >= tilecols)
		{
			row_col_List[i].col_idx = tilecols - 1;
		}

		if (row_col_List[i].row_idx < 0)
		{
			row_col_List[i].row_idx = 0;
		}
		else if (row_col_List[i].row_idx >= tilerows)
		{
			row_col_List[i].row_idx = tilerows - 1;
		}
	}
	return row_col_List;
}

static ossimDpt gsd_degree(int nZoolLevel = 0)
{
	//Map resolution = 156543.04 meters / pixel * cos(latitude) / (2 ^ zoomlevel)
	double x = 1.40625 / (1 << nZoolLevel);
	double y = 0.703125 / (1 << nZoolLevel);
	return ossimDpt(x, y);
}

static ossimDpt gsd_meter(double lat = 0, int nZoolLevel = 0)
{
	//Map resolution = 156543.04 meters / pixel * cos(latitude) / (2 ^ zoomlevel)
	double d = 156543.033900000 * cos(lat*0.0174532925) / (1 << nZoolLevel);
	return ossimDpt(d, d);
}

static int getZoomLevel(double res, double lat = 0)
{
	//return (int)(log(156543.033900000 * cos(lat*0.0174532925) / res) / log(2) + 0.5);
	//return (int)(log(156543.033900000 * cos(lat*0.0174532925) / res) / log(2));
	//return ceil(log(156543.033900000 * cos(lat*0.0174532925) / res) / log(2.0));
	return int(log(156543.033900000 * cos(lat*0.0174532925) / res) / log(2.0) + 0.5);
}


static bool kpt_compare(const cv::KeyPoint& d1, const cv::KeyPoint& d2)
{
	return d1.response > d2.response;
}

static void normalization(const cv::Mat& inMat, cv::Mat& outMat)
{
	cv::Scalar mean_value;// = cv::mean(inMat);
	cv::Mat stdDevMat;
	cv::meanStdDev(inMat, mean_value, stdDevMat);
	cv::divide(inMat - mean_value, stdDevMat, outMat);
}

static void findGoodMatches(vector< vector< DMatch >  > all_matches_2, vector< DMatch >& good_matches, float nndrRatio = 0.80f)
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

		if (m1.distance <= nndrRatio * m2.distance)
			good_matches.push_back(m1);
	}
}

static void findGoogdMatches(vector< vector< DMatch >  > all_matches_1, vector< vector< DMatch >  > all_matches_2,
	vector< DMatch >& good_matches, float nndrRatio = 0.80f, bool bUseCrossMatch = true)
{
	// "cross-matching" and "first and second minimum distances ratio test"
	//vector< DMatch > good_matches;
	good_matches.clear();
	for (size_t i = 0; i < all_matches_1.size(); i++)
	{
		if (all_matches_1[i].size() != 2)
		{
			continue;
		}
		const DMatch &m1 = all_matches_1[i][0];
		const DMatch &m2 = all_matches_1[i][1];

		if (m1.distance <= nndrRatio * m2.distance)
			//	good_matches.push_back(m1);
			//if (matches[i][0].distance / (matches[i][1].distance + FLT_EPSILON) < 0.6)
		{
			good_matches.push_back(m1);
			continue;
		}
		if (bUseCrossMatch)
		{
			int queryIdx = all_matches_1[i][0].queryIdx;
			int trainIdx = all_matches_1[i][0].trainIdx;
			for (size_t j = 0; j < all_matches_2.size(); j++)
			{
				int queryIdx2 = all_matches_2[j][0].trainIdx;
				int trainIdx2 = all_matches_2[j][0].queryIdx;

				if (queryIdx == queryIdx2 && trainIdx == trainIdx2)
				{
					good_matches.push_back(all_matches_1[i][0]);
					break;
				}
			}
		}
	}
}

static void removeRepeated(const std::vector<KeyPoint>& skeypoints, const std::vector<KeyPoint>&  mkeypoints, vector< DMatch >& good_matches, double pos_threshold = 2.0)
{
	// eliminating repeated points
	vector< DMatch > existed_matches;
	for (size_t i = 0; i < good_matches.size();)
	{
		bool bExisted = false;
		for (size_t j = 0; j < existed_matches.size(); j++)
		{
			if (fabs(mkeypoints[existed_matches[j].trainIdx].pt.x - mkeypoints[good_matches[i].trainIdx].pt.x) < pos_threshold
				&& fabs(skeypoints[existed_matches[j].queryIdx].pt.x - skeypoints[good_matches[i].queryIdx].pt.x) < pos_threshold)
			{
				bExisted = true;
				break;
			}
		}

		if (!bExisted)
		{
			existed_matches.push_back(good_matches[i]);
		}
		else
		{
			good_matches.erase(good_matches.begin() + i);
			continue;
		}

		i++;
	}
}


static void _prepareImgAndDrawKeylines(const cv::Mat& img1,
	const cv::Mat& img2,
	ossimTDpt tpt,
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

	outImg = cv::Mat(size, CV_MAKETYPE(img1.depth(), 3), cv::Scalar(0, 0, 0));

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

static void drawTogether(const cv::Mat& img1,
	const cv::Mat& img2,
	const cv::Mat& img3,
	cv::Mat& outImg)
{
	int nBlank = 0;
	cv::Mat outImg1;
	cv::Mat outImg2;
	cv::Mat outImg3;
	Size size(img1.cols + img2.cols + img3.cols + nBlank * 2, MAX(MAX(img1.rows, img2.rows), img3.rows));

	outImg = cv::Mat(size, CV_MAKETYPE(img1.depth(), 3), cv::Scalar(0, 0, 0));

	outImg1 = outImg(Rect(0, 0, img1.cols, img1.rows));
	outImg2 = outImg(Rect(img1.cols + nBlank, 0, img2.cols, img2.rows));
	outImg3 = outImg(Rect(img1.cols + img2.cols + nBlank * 2, 0, img3.cols, img3.rows));

	if (img1.type() == CV_8U)
		cvtColor(img1, outImg1, CV_GRAY2BGR);
	else
		img1.copyTo(outImg1);

	if (img2.type() == CV_8U)
		cvtColor(img2, outImg2, CV_GRAY2BGR);
	else
		img2.copyTo(outImg2);

	if (img3.type() == CV_8U)
		cvtColor(img3, outImg3, CV_GRAY2BGR);
	else
		img3.copyTo(outImg3);
}

static void drawTogether(const cv::Mat& img1,
	const cv::Mat& img2,
	const cv::Mat& img3,
	const cv::Mat& img4,
	cv::Mat& outImg)
{
	int nBlank = 0;
	cv::Mat outImg1;
	cv::Mat outImg2;
	cv::Mat outImg3;
	cv::Mat outImg4;
	Size size(img1.cols + img2.cols + img3.cols + img4.cols + nBlank * 3, MAX(MAX(MAX(img1.rows, img2.rows), img3.rows), img4.rows));

	outImg = cv::Mat(size, CV_MAKETYPE(img1.depth(), 3), cv::Scalar(0, 0, 0));

	outImg1 = outImg(Rect(0, 0, img1.cols, img1.rows));
	outImg2 = outImg(Rect(img1.cols + nBlank, 0, img2.cols, img2.rows));
	outImg3 = outImg(Rect(img1.cols + img2.cols + nBlank * 2, 0, img3.cols, img3.rows));
	outImg4 = outImg(Rect(img1.cols + img2.cols + img3.cols + nBlank * 3, 0, img4.cols, img4.rows));

	if (img1.type() == CV_8U)
		cvtColor(img1, outImg1, CV_GRAY2BGR);
	else
		img1.copyTo(outImg1);

	if (img2.type() == CV_8U)
		cvtColor(img2, outImg2, CV_GRAY2BGR);
	else
		img2.copyTo(outImg2);

	if (img3.type() == CV_8U)
		cvtColor(img3, outImg3, CV_GRAY2BGR);
	else
		img3.copyTo(outImg3);
	if (img4.type() == CV_8U)
		cvtColor(img4, outImg4, CV_GRAY2BGR);
	else
		img4.copyTo(outImg4);
}


static float get_sub_pix(cv::Mat const &img, cv::Point2f const &pt)
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

static void funcErrorEquation(double *param, double *hx, int nparameter, int nequation, void *adata)
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

static void jacErrorEquation(double *param, double *j, int nparameter, int nequation, void *adata)
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

static bool leastSquareMatching(cv::Mat templateMat, cv::Mat searchMat, double lsm_model[8])
{
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

	if (ret >= 0 &&
		(1 == info[6] || 2 == info[6] || 6 == info[6]))
	{
		return true;
	}
	else
	{
		return false;
	}
}


struct MemoryStruct {
	char *memory;
	size_t size;
};

static ossimGpt TileXy2LonLat(double x, double y, int nZoomLevel)
{
	double scale = 1.0 / (double)(1 << (8 + nZoomLevel));	// 2^(-n-3)
	double lon = x * 360.0 * scale - 180.0;
	double t = exp((0.5 - y * scale) * 4.0 * PI);
	double lat = asin((t - 1) / (t + 1))*180.0 / PI;

	return ossimGpt(lat, lon);
}

static ossimDpt LonLat2TileXy(ossimGpt gpt, int nZoomLevel)
{
	double scale = (double)(1 << (8 + nZoomLevel));	// 2^(n+3)
	double sinLat = sin(gpt.lat * PI / 180.0);
	double x = (gpt.lon + 180.0) * scale / 360.0;
	double y = (0.5 - log((1.0 + sinLat) / (1.0 - sinLat)) / (4.0*PI))*scale;
	return ossimDpt(x, y);
}

static ossimGpt masterLineSample2World(ossimDpt linesample,
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


static struct MemoryStruct getDataFromUrl(const char* url, const char* proxy = "")
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


static ossimDrect getLonLatRect(ossimProjection* sProjection,
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


static string getRequestUrl(ossimProjection* sProjection, double res,
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


static VlSiftKeypoint cov2sift(const VlCovDetFeature &feature)
{
	VlSiftKeypoint siftKeyPoint;
	siftKeyPoint.x = feature.frame.x;
	siftKeyPoint.y = feature.frame.y;
	siftKeyPoint.sigma = feature.laplacianScaleScore;
	siftKeyPoint.contrast = feature.peakScore;
	siftKeyPoint.s = feature.edgeScore;
	return siftKeyPoint;
}

static void VLFeatCovdet(const cv::Mat& inMat, vector<KeyPoint>& kpts, cv::Mat& descriptors)
{
	//// create a detector object
	//VlCovDet * covdet = vl_covdet_new(VlCovDetMethod::VL_COVDET_METHOD_HESSIAN_LAPLACE);
	//// set various parameters (optional)
	//vl_covdet_set_first_octave(covdet, -1); // start by doubling the image resolution
	////vl_covdet_set_octave_resolution(covdet, octaveResolution);
	////vl_covdet_set_peak_threshold(covdet, peakThreshold);
	////vl_covdet_set_edge_threshold(covdet, edgeThreshold);


	//vl_sift_pix *ImageData = new vl_sift_pix[inMat.rows * inMat.cols];

	//int noctaves = 2, nlevels = 4, o_min = 0;
	//unsigned char *Pixel;
	//for (int i = 0; i<inMat.rows; i++)
	//{
	//	for (int j = 0; j<inMat.cols; j++)
	//	{
	//		Pixel = (unsigned char*)(inMat.data + i*inMat.cols + j);
	//		ImageData[i*inMat.cols + j] = *(Pixel);
	//	}
	//}

	//// process the image and run the detector
	//vl_covdet_put_image(covdet, ImageData, inMat.cols, inMat.rows);
	//vl_covdet_detect(covdet);
	//// drop features on the margin (optional)
	////vl_covdet_drop_features_outside(covdet, boundaryMargin);
	//// compute the affine shape of the features (optional)
	//vl_covdet_extract_affine_shape(covdet);
	//// compute the orientation of the features (optional)
	//vl_covdet_extract_orientations(covdet);
	//// get feature frames back
	//vl_size numFeatures = vl_covdet_get_num_features(covdet);
	//VlCovDetFeature const * feature = vl_covdet_get_features(covdet);
	//// get normalized feature appearance patches (optional)
	//vl_size w = 2 * patchResolution + 1;
	//for (i = 0; i < numFeatures; ++i) {
	//	float * patch = malloc(w*w*sizeof(*desc));
	//	vl_covdet_extract_patch_for_frame(covdet,
	//		patch,
	//		patchResolution,
	//		patchRelativeExtent,
	//		patchRelativeSmoothing,
	//		feature[i].frame);
	//	// do something with patch
	//}




	int noctaves = 2, octaveResolution = 4, o_min = 0;

	VlSiftFilt *SiftFilt = NULL;
	SiftFilt = vl_sift_new(inMat.cols, inMat.rows, noctaves, octaveResolution, o_min);
	//VlCovDet *vl_covDet = vl_covdet_new(VlCovDetMethod::VL_COVDET_METHOD_HESSIAN_LAPLACE);
	VlCovDet *vl_covDet = vl_covdet_new(VlCovDetMethod::VL_COVDET_METHOD_HESSIAN_LAPLACE);
	vl_sift_pix *ImageData = new vl_sift_pix[inMat.rows * inMat.cols];

	double peakThreshold = 400.0;
	double edgeThreshold = 10.0;
	unsigned char *Pixel;
	for (int i = 0; i<inMat.rows; i++)
	{
		for (int j = 0; j<inMat.cols; j++)
		{
			Pixel = (unsigned char*)(inMat.data + i*inMat.cols + j);
			ImageData[i*inMat.cols + j] = *(Pixel);
		}
	}

	vl_size patchResolution = 5;
	double patchRelativeExtent = 3;
	double patchRelativeSmoothing = 1.6;

	// set various parameters (optional)
	vl_covdet_set_first_octave(vl_covDet, -1); // start by doubling the image resolution
	vl_covdet_set_first_octave(vl_covDet, o_min); // start by doubling the image resolution
	//vl_covdet_set_octave_resolution(vl_covDet, octaveResolution);
	vl_covdet_set_peak_threshold(vl_covDet, peakThreshold);
	vl_covdet_set_edge_threshold(vl_covDet, edgeThreshold);

	int nKeyPoint = 0;
	int idx = 0;
	vl_covdet_put_image(vl_covDet, ImageData, inMat.cols, inMat.rows);


	vl_sift_process_first_octave(SiftFilt, ImageData);

	vl_covdet_detect(vl_covDet);
	int nkeys = vl_covdet_get_num_features(vl_covDet);
	VlCovDetFeature const * features = (VlCovDetFeature const *)vl_covdet_get_features(vl_covDet);
	for (int i = 0; i<nkeys; i++)
	{
		VlCovDetFeature feature = features[i];
		cv::KeyPoint kpt(feature.frame.x, feature.frame.y, feature.peakScore);
		idx++;
		//计算并遍历每个点的方向
		vl_size angleCount;
		VlCovDetFeatureOrientation * angles = vl_covdet_extract_orientations_for_frame(vl_covDet, &angleCount, feature.frame);

		for (int j = 0; j<angleCount; j++)
		{
			double TemptAngle = angles[j].angle;
			VlSiftKeypoint TemptKeyPoint = cov2sift(feature);

			//printf("%d: %f\n",j,TemptAngle);
			//计算每个方向的描述
			cv::Mat aDescriptor = cv::Mat(cv::Size(128, 1), CV_32FC1);
			//float *Descriptors=new float[128];
			vl_sift_calc_keypoint_descriptor(SiftFilt, (float*)aDescriptor.data, &TemptKeyPoint, TemptAngle);
			//vl_size w = 2 * patchResolution + 1;
			//cv::Mat aDescriptor = cv::Mat(cv::Size(w*w, 1), CV_32FC1);
			////cv::Mat aDescriptor = cv::Mat(cv::Size(128, 1), CV_32FC1);
			vl_covdet_extract_patch_for_frame(vl_covDet,
				(float*)aDescriptor.data,
				patchResolution,
				patchRelativeExtent,
				patchRelativeSmoothing,
				feature.frame);


			//float *Descriptors=new float[128];
			//vl_sift_calc_keypoint_descriptor(SiftFilt, (float*)aDescriptor.data, &TemptKeyPoint, TemptAngle);
			descriptors.push_back(aDescriptor);
			kpt.angle = TemptAngle;
			kpt.response = TemptKeyPoint.contrast;
			kpts.push_back(kpt);
		}
	}

	vl_sift_delete(SiftFilt);
	vl_covdet_delete(vl_covDet);
	delete[]ImageData;
	ImageData = NULL;
}

static void OrbDetector(const cv::Mat& inMat, vector<KeyPoint>& kpts, cv::Mat& descriptors)
{
	int numKeyPoints = 500;
	float distThreshold = 15.0;
	//"FAST" – FastFeatureDetector
	//"STAR" – StarFeatureDetector
	//"SIFT" – SiftFeatureDetector
	//"SURF" – SurfFeatureDetector
	//"ORB" – OrbFeatureDetector
	//"MSER" – MserFeatureDetector
	//"GFTT" – GoodFeaturesToTrackDetector
	//"HARRIS" – GoodFeaturesToTrackDetector with Harris detector enabled
	//"Dense" – DenseFeatureDetector
	//"SimpleBlob" – SimpleBlobDetector
	Ptr<FeatureDetector> detector = FeatureDetector::create("ORB");
	detector->detect(inMat, kpts);
	//cv::OrbFeatureDetector* detector = new cv::OrbFeatureDetector(numKeyPoints);
	//cv::OrbDescriptorExtractor* extractor = new cv::OrbDescriptorExtractor;
	//detector->detect(inMat, kpts);
	//extractor->compute(inMat, kpts, descriptors);

	//"SIFT" – SiftDescriptorExtractor
	//"SURF" – SurfDescriptorExtractor
	//"ORB" – OrbDescriptorExtractor
	//"BRIEF" – BriefDescriptorExtractor
	Ptr<DescriptorExtractor> extractor = DescriptorExtractor::create("SIFT");
	extractor->compute(inMat, kpts, descriptors);

	//delete detector;
	//delete extractor;
}

static void VLFeatSift(const cv::Mat& inMat, vector<KeyPoint>& kpts, cv::Mat& descriptors)
{
	//int noctaves = 2, nlevels = 4, o_min = 0;
	int noctaves = 1, nlevels = 5, o_min = 0;
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
	// default
	double edge_thresh = 5;  //-1 will use the default (as in matlab)
	double peak_thresh = 0.04;// 0.04;
	//double edge_thresh = -1;  //-1 will use the default (as in matlab)
	//double peak_thresh = -1;
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


static Eigen::VectorXd getAffineFromPyramid(const GdalRasterApp& slaveApp, const GdalRasterApp& masterApp, int sband = 0, int mband = 0, double scale = 16)
{
	Eigen::VectorXd affine_model;
	bool bStretch = true;
	cv::Mat slaveMat;
	ossimIrect srect(ossimIpt(0, 0), ossimIpt(slaveApp.width() - 1, slaveApp.height() - 1));
	if (!slaveApp.getRect2CvMatByte(srect, slaveMat, sband, ossimDpt(1.0 / scale, 1.0 / scale), 0.015, bStretch))
	{
		return affine_model;
	}

	cv::Mat masterMat;
	ossimIrect mrect(ossimIpt(0, 0), ossimIpt(masterApp.width() - 1, masterApp.height() - 1));
	if (!masterApp.getRect2CvMatByte(mrect, masterMat, sband, ossimDpt(1.0 / scale, 1.0 / scale), 0.015, bStretch))
	{
		return affine_model;
	}

	float nndrRatio = 0.75f;
	float angle_diff_threshold = 20.0f / 180.0f * VL_PI;
	double pos_threshold = 0.5;
	double affine_residual_threshold = 2.0; //pixel
	float std_scale_diff_threshold = 0.2f;
	float scale_ratio_threshold = 0.7f;	// 0.8f
	int count_pos = 0;
	int matched_counts[10];


	// detect corners
	//cv::initModule_features2d();
	std::vector<KeyPoint> skeypoints, mkeypoints;
	cv::Mat sdescriptors, mdescriptors;

	VLFeatSift(slaveMat, skeypoints, sdescriptors);
	if (skeypoints.size() < 10)
	{
		return affine_model;

	}
	VLFeatSift(masterMat, mkeypoints, mdescriptors);
	if (mkeypoints.size() < 10)
	{
		return affine_model;
	}
	
	BFMatcher matcher(NORM_L1, false);
	vector< vector< DMatch >  > matches;
	matcher.knnMatch(sdescriptors, mdescriptors, matches, 2);

	// inverse mathcing
	vector< vector< DMatch >  > matches_inverse;
	matcher.knnMatch(mdescriptors, sdescriptors, matches_inverse, 2);

	// "cross-matching" and "first and second minimum distances ratio test"
	vector< DMatch > good_matches;
	findGoogdMatches(matches, matches_inverse, good_matches, nndrRatio, true);
	vector< DMatch > good_matches0 = good_matches;

	matched_counts[count_pos++] = (int)good_matches.size();

	// find max scale
	const int scale_bins = 41;
	float scale_step = 1.2f;
	int scale_hist[scale_bins] = { 0 };
	float scale_bin_length = 1.0f;
	float scale_bin_start = (scale_bins / 2);
	std::vector<double> scaleList;
	for (size_t i = 0; i < good_matches.size(); i++)
	{
		float a1 = mkeypoints[good_matches[i].trainIdx].size;
		float a2 = skeypoints[good_matches[i].queryIdx].size;
		float scale_ratio = a1 / a2;

		float log_scale_ratio = log(scale_ratio) / log(scale_step);

		scale_hist[(int)(scale_bin_start + log_scale_ratio / scale_bin_length + 0.5)]++;
	}

	double scalePeak = 1.0;
	/* find the histogram maximum */
	int maxh = 0;
	for (int i = 0; i < scale_bins; ++i)
		maxh = max(maxh, scale_hist[i]);

	for (int i = 0; i < scale_bins; ++i) {
		double h0 = scale_hist[i];
		double hm = scale_hist[(i - 1 + scale_bins) % scale_bins];
		double hp = scale_hist[(i + 1 + scale_bins) % scale_bins];

		/* is this a peak? */
		if (h0 > 0.95*maxh && h0 > hm && h0 > hp) {

			/* quadratic interpolation */
			double di = -0.5 * (hp - hm) / (hp + hm - 2 * h0);
			double th = (i + di + 0.5) - scale_bin_start;
			scalePeak = pow(scale_step, th);
			break;
		}
	}

	// eliminating by scale
	for (size_t i = 0; i < good_matches.size();)
	{
		float s1 = mkeypoints[good_matches[i].trainIdx].size;
		float s2 = skeypoints[good_matches[i].queryIdx].size;
		float scale_ratio = fabs(s1 / s2);
		if (fabs(log(scale_ratio / scalePeak)) > fabs(log(scale_ratio_threshold)))
		//if (scale_ratio < scale_ratio_threshold || scale_ratio * scale_ratio_threshold > scalePeak)
		{
			good_matches.erase(good_matches.begin() + i);
			continue;
		}
		i++;
	}

	matched_counts[count_pos++] = (int)good_matches.size();
	if (good_matches.size() < 4)
	{
		return affine_model;
	}

	// find max rotation angle
	const int angle_bins = 36;
	int angle_hist[angle_bins] = { 0 };
	float angle_bin_length = 2.0f*VL_PI / (float)angle_bins;
	std::vector<double> angle_diffList;
	for (size_t i = 0; i < good_matches.size(); i++)
	{
		float a1 = mkeypoints[good_matches[i].trainIdx].angle;
		float a2 = skeypoints[good_matches[i].queryIdx].angle;
		float angle_diff = a1 - a2;
		if (angle_diff < 0.0)
		{
			angle_diff += 2.0*VL_PI;
		}

		angle_hist[(int)(angle_diff / angle_bin_length)]++;
	}

	double ratAngle = 0.0;
	/* find the histogram maximum */
	maxh = 0;
	for (int i = 0; i < angle_bins; ++i)
		maxh = max(maxh, angle_hist[i]);

	for (int i = 0; i < angle_bins; ++i) {
		double h0 = angle_hist[i];
		double hm = angle_hist[(i - 1 + angle_bins) % angle_bins];
		double hp = angle_hist[(i + 1 + angle_bins) % angle_bins];

		/* is this a peak? */
		if (h0 > 0.95*maxh && h0 > hm && h0 > hp) {

			/* quadratic interpolation */
			double di = -0.5 * (hp - hm) / (hp + hm - 2 * h0);
			double th = 2 * VL_PI * (i + di + 0.5) / angle_bins;
			ratAngle = th;
			break;
		}
	}

	// eliminating by rotation
	for (size_t i = 0; i < good_matches.size();)
	{
		float a1 = mkeypoints[good_matches[i].trainIdx].angle;
		float a2 = skeypoints[good_matches[i].queryIdx].angle;
		float angle_diff = a1 - a2;
		if (angle_diff < 0.0)
		{
			angle_diff += 2.0*VL_PI;
		}

		if (fabs(angle_diff - ratAngle) > angle_diff_threshold)
		{
			good_matches.erase(good_matches.begin() + i);
			continue;
		}
		i++;
	}
	matched_counts[count_pos++] = (int)good_matches.size();

	// eliminating repeated points
	removeRepeated(skeypoints, mkeypoints, good_matches, pos_threshold);


	if (good_matches.size() < 4)
	{
		return affine_model;
	}
	
	vector<int> inliers;
	vector<double> rigid_model = mylib::rigidRansac(skeypoints, mkeypoints, good_matches, inliers, 0.85);
	matched_counts[count_pos++] = (int)inliers.size();

	//double s = sqrt(rigid_model[0] * rigid_model[0] + rigid_model[1] * rigid_model[1]);
	////double angle = atan(rigid_model[1] / rigid_model[0]);
	//double arcsin = rigid_model[1] / s;
	//double arccos = rigid_model[0] / s;
	//double angle;
	//if (arccos >= 0)
	//{
	//	angle = asin(arcsin);
	//}
	//else
	//{
	//	if (arcsin >= 0)
	//	{
	//		angle = PI - asin(arcsin);
	//	}
	//	else if (arcsin < 0)
	//	{
	//		angle = -PI - asin(arcsin);
	//	}
	//}

	//double angle_threshold = 0.5; // rad -- 28.6478897565deg
	//if (fabs(angle) > PI*0.25)
	//{
	//	return affine_model;
	//}
	//for (size_t i = 0; i < inliers.size(); i++)
	//{
	//	float a1 = mkeypoints[good_matches[inliers[i]].trainIdx].angle * PI / 360.0;
	//	float a2 = skeypoints[good_matches[inliers[i]].queryIdx].angle * PI / 360.0;
	//	double angle_diff = a1 - a2;
	//	if (fabs(angle_diff + angle) > angle_threshold)
	//	{
	//		return affine_model;
	//		inliers.erase(inliers.begin() + i);
	//		i--;
	//	}
	//}

	//if (fabs(s - 1) > 0.5)
	//{
	//	return affine_model;
	//}

	// fit affine model
	//Eigen::VectorXd affine_model;
	while (inliers.size() >= 4)
	{
		// Build matrices to solve Ax = b problem:
		Eigen::VectorXd b(inliers.size() * 2);
		Eigen::MatrixXd A(inliers.size() * 2, 6);
		//b = candidates.col(6);
		for (int i = 0; i < (int)inliers.size(); ++i)
		{
			A(2 * i, 0) = 1.0;
			A(2 * i, 1) = skeypoints[good_matches[inliers[i]].queryIdx].pt.x;
			A(2 * i, 2) = skeypoints[good_matches[inliers[i]].queryIdx].pt.y;
			A(2 * i, 3) = 0.0;
			A(2 * i, 4) = 0.0;
			A(2 * i, 5) = 0.0;
			b(2 * i) = mkeypoints[good_matches[inliers[i]].trainIdx].pt.x;

			A(2 * i + 1, 0) = 0.0;
			A(2 * i + 1, 1) = 0.0;
			A(2 * i + 1, 2) = 0.0;
			A(2 * i + 1, 3) = 1.0;
			A(2 * i + 1, 4) = skeypoints[good_matches[inliers[i]].queryIdx].pt.x;
			A(2 * i + 1, 5) = skeypoints[good_matches[inliers[i]].queryIdx].pt.y;
			b(2 * i + 1) = mkeypoints[good_matches[inliers[i]].trainIdx].pt.y;
		}
		// Compute least-squares solution:
		affine_model = (A.transpose()*A).inverse()*A.transpose()*b;
		Eigen::VectorXd resMat = A*affine_model - b;
		double max_residual = 0.0;
		int max_pos = 0;
		for (size_t i = 0; i < inliers.size(); i++)
		{
			double residual = sqrt(resMat[2 * i] * resMat[2 * i] + resMat[2 * i + 1] * resMat[2 * i + 1]);
			if (residual > max_residual)
			{
				max_residual = residual;
				max_pos = i;
			}
		}

		if (max_residual > affine_residual_threshold)
		{
			inliers.erase(inliers.begin() + max_pos);
		}
		else
		{
			break;
		}
	}

	matched_counts[count_pos++] = (int)inliers.size();
	if (inliers.size() < 4)
	{
		return affine_model;
	}
	
	affine_model[0] *= scale;
	affine_model[3] *= scale;
	return affine_model;
}

#endif