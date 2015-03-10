#pragma once

#include "lineFeature.h"
#include "cv.h"
#include <stdio.h>
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <levmar.h>
using namespace cv;

vector<double> affine_estimation2(const std::vector<cvline_polar>& src_line_polar,
	const std::vector<cvline_polar>& ref_line_polar);
vector<double> affine_estimation2(const std::vector<cvline_polar>& src_line_polar,
	const std::vector<cvline_polar>& ref_line_polar, const vector<int>& slaveAssign);

void LineMatch(const cv::Mat& ref_seg, const cv::Mat& src_seg,vector< vector<double> > &ref_LineFeature,vector< vector<double> > &src_LineFeature,
	vector<int>&ref_LineNum, vector<int>&src_LineNum, vector<cvline_polar>&ref_lines_list, vector<cvline_polar>&src_lines_list,
	vector<cvline_polar>& matched_ref_ls, vector<cvline_polar>& matched_src_ls);




double point2polarline(const fPoint& p, const cvline_polar& l);

double distance1(const cvline_polar& l1, const cvline_polar& l2);

fPoint distance2(const cvline_polar& l1, const cvline_polar& l2);

cvline_polar forward(const cvline_polar& l, const vector<double>& parameters);



struct dataStruct
{
	vector<cvline_polar> slave_vec_lines;
	vector<cvline_polar> master_vec_lines;
	vector<int> slaveAssign;
	std::vector<double> parameters;
};

void levmar_function_fvec(double *param, double *hx, int nparameter, int nequation, void *adata);

void levmar_function_jac(double *param, double *jac, int nparameter, int nequation, void *adata);

bool parameters_optimization(const vector<cvline_polar>& slave_vec_lines, const vector<cvline_polar>& master_vec_lines,
	const vector<int>& slaveAssign, std::vector<double>& parameters);


bool parameters_optimization(const vector<cvline_polar>& slave_vec_lines, const vector<cvline_polar>& master_vec_lines, std::vector<double>& parameters);