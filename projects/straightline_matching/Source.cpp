#include <Windows.h>
#include <direct.h>
//// hough.cpp : Defines the entry point for the console application.
////
#include "LineExtract.h"
#include "lineConstant.h"
#include "func.h"

#include <iostream>
#include <Eigen/Eigen>
using namespace std;
//#pragma comment(lib,"LineMatcher.lib")

#include <GdalRasterApp.h>

#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "mylib.lib")
#pragma comment(lib, "ossim20.lib")
#pragma comment(lib, "opencv_ts300.lib")
#pragma comment(lib, "opencv_world300.lib")
#pragma comment(lib, "blas_win32_MT.lib")
#pragma comment(lib, "lapack_win32_MT.lib")
#pragma comment(lib, "mlpack.lib")
vector<double> parameters(6);

void gdal_getlines(const char* strImage, vector<cvline_polar>& vec_lines, double scale/* = 0.8*/, int seletion_num/* = 200*/)
{
	mylib::GdalRasterApp gdalApp;
	gdalApp.open(strImage);
	gdalApp.setTileWidth(2048);
	gdalApp.setTileHeight(2048);
	int nTileX = gdalApp.getTileCountX();
	int nTileY = gdalApp.getTileCountY();

	// Some constants needed throughout...
	const ossim_int32 START_LINE = theAreaOfInterest.ul().y;
	const ossim_int32 STOP_LINE = theAreaOfInterest.lr().y;
	const ossim_int32 START_SAMP = theAreaOfInterest.ul().x;
	const ossim_int32 STOP_SAMP = theAreaOfInterest.lr().x;

	// For percent complete status.
	ossim_int32 tilerows = ceil(sqrt(nPointRequired * nHeight / (double)nWidth));
	ossim_int32 tilecols = ceil(sqrt(nPointRequired * nWidth / (double)nHeight));
	const ossim_int32 TILE_HEIGHT = ceil(nHeight / (double)tilerows);
	const ossim_int32 TILE_WIDTH = ceil(nWidth / (double)tilecols);
	double total_tiles = ((double)tilerows)*tilecols;
	double total_tiles_processed;
	double tiles_processed = 0.0;


	// loop through all tiles
	// need to use a sequencer for parallelism in the future TBD
	theTiePoints.clear();

	vector<row_col> row_col_List;
	for (int i = 0; i<tilerows; ++i)
	{
		for (int j = 0; j<tilecols; ++j)
		{
			row_col_List.push_back(row_col(i, j));
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
	totalBlocks = (int)row_col_List.size();
	if (theThreadNum == 0)
	{
		theThreadNum = OpenThreads::GetNumberOfProcessors() * 2;
	}
	//if(num_threads > totalBlocks) num_threads = totalBlocks;
	//num_threads = 1;
	cout << "using " << theThreadNum << " threads..." << endl;
	GLOBAL_NUM_THREADS = theThreadNum + 1;
	std::vector<radiMatchRectThread *> threads(theThreadNum);
	OpenThreads::Thread::SetConcurrency(theThreadNum);
	OpenThreads::Thread::Init();
	//totalBlocks = 20; // for debug
	for (int i = 0; i<theThreadNum; ++i) {
		threads[i] = new radiMatchRectThread(this);
		vector<ossimIrect> rectList;
		for (int j = i; j < totalBlocks; j += theThreadNum)
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
		threads[i]->setRect(rectList);
		threads[i]->start();
	}
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till ready
	bar.block(GLOBAL_NUM_THREADS);  // Block 'till finished

	for (size_t i = 0; i < nTileY; i++)
	{
		for (size_t j = 0; j < nTileX; j++)
		{

		}
	}

	gdalApp.close();
}

void main()
{
	string sourceImage = "D:\\workspace\\Landsat\\xinjiang\\LS5_TM_20090721_050335_050401_146030_FASTB_L2\\image.tif";
	string sourceImage_line = "D:\\workspace\\Landsat\\xinjiang\\source_line.tif";
	string referImage = "D:\\workspace\\Landsat\\xinjiang\\reference.tif";
	string referImage_line = "D:\\workspace\\Landsat\\xinjiang\\reference_line.tif";
	vector<cvline_polar> src_lines_polar;
	vector<cvline_polar> ref_lines_polar;
	double sourceimage_p1 = 0.8;
	double sourceimage_p2 = 500;
	double referimage_p1 = 0.8;
	double referimage_p2 = 500;
	IplImage* src_line = get_lines(sourceImage.c_str(), src_lines_polar, sourceimage_p1, sourceimage_p2);
	IplImage* ref_line = get_lines(referImage.c_str(), ref_lines_polar, referimage_p1, referimage_p2);
	if (!src_line || !ref_line)
	{
		return;
	}
	//outLineAngle(line_angle, ref_lines_polar, src_lines_polar);
	cvSaveImage(sourceImage_line.c_str(), src_line, 0);
	cvSaveImage(referImage_line.c_str(), ref_line, 0);

	cvReleaseImage(&src_line);
	cvReleaseImage(&ref_line);
}