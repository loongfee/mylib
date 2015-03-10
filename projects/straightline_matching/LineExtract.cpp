#include "LineExtract.h"
#include "lsd.h"
#include "func.h"
//#include "LSWMS.h"

#include "quick_selection.h"
//
//IplImage* get_lines(IplImage* img,vector<cvline_polar>& vec_lines)
//{
//    //to grey
//    //IplImage* grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
//    //cvCvtColor(img, grey, CV_BGR2GRAY);
//
//    image_double image;
//    ntuple_list out;
//    unsigned int x,y,i,j;
//    image = new_image_double(img->width,img->height);
//    for(x=0;x</*grey*/img->width;x++)
//    for(y=0;y</*grey*/img->height;y++)
//    {
//      CvScalar s= cvGet2D(/*grey*/img,y,x);
//      double pix= s.val[0];
//      image->data[ x + y * image->xsize ]= pix;
//    }
//
//    /* call LSD */
//    out = lsd(image);
//    //out= lsd_scale(image,1);
//
//    /* print output */
//    //printf("%u line segments found:\n",out->size);
//    //vector<Line> vec;
//    for(i=0;i<out->size;i++)
//    {
//      //for(j=0;j<out->dim;j++)
//      {
//        //printf("%f ",out->values[ i * out->dim + j ]);
//
//		  cvline_polar outLine;
//		  fPoint line[2];
//		  line[0].x = out->values[ i * out->dim + 0];
//		  line[0].y = out->values[ i * out->dim + 1];
//		  line[1].x = out->values[ i * out->dim + 2];
//		  line[1].y = out->values[ i * out->dim + 3];
//
//		  outLine = line2polar(line);
//		  double len = segment_length(outLine);
//		  if(len > 5.0)
//          /*vec*/vec_lines.push_back(outLine);
//      }
//      //printf("\n");
//    }
//
//    IplImage* black= cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
//    cvZero(black);
//    for(int i=0;i<vec_lines.size();++i)
//    {
//        //if(vec[i].x1==vec[i].x2||vec[i].y1==vec[i].y2)
//        cvLine(black,cvPoint(vec_lines[i].pt1.x,vec_lines[i].pt1.y),cvPoint(vec_lines[i].pt2.x,vec_lines[i].pt2.y),CV_RGB(255,255,255),1, CV_AA);
//    }
//    /*cvNamedWindow("img", 0);
//    cvShowImage("img", img);*/
//    //cvSaveImage("lines_detect.png",black/*img*/);
//    /* free memory */
//    //cvReleaseImage(&grey);
//    free_image_double(image);
//    free_ntuple_list(out);
//
//    return black;
//}

bool isBoundary(IplImage *img, cvline_polar line, cv::Scalar bgColor/* = CV_RGB(0,0,0)*/)
{
	int X = img->width;
	int Y = img->height;
	int band = img->nChannels;
	int step = img->widthStep;
	uchar * data = (uchar *)img->imageData;

	fPoint pt1 = line.pt1;
	fPoint pt2 = line.pt2;

	int k = 3;
	// 外扩k个像素
	cvline_polar l1 = segment_move_vertical(line, (double)k);
	cvline_polar l2 = segment_move_vertical(line, -(double)k);

	double len = (int)(sqrt((l1.pt1.x-l1.pt2.x)*(l1.pt1.x-l1.pt2.x)+(l1.pt1.y-l1.pt2.y)*(l1.pt1.y-l1.pt2.y)) + 0.5);
	int i;
	for (i = 0;i <= len;++i)
	{
		double lammda = i / (double)len;
		int x = (int)(lammda*(l1.pt1.x)+(1.0-lammda)*l1.pt2.x+0.5);
		int y = (int)(lammda*(l1.pt1.y)+(1.0-lammda)*l1.pt2.y+0.5);
		if (x <0 || x>= X || y < 0 || y >=Y)
		{
			continue;
		}

		bool same_color = true;
		for (int k = 0;k < band;++k)
		{
			if (bgColor.val[k] != data[ y*step + x*band + k])
			{
				same_color = false;
				break;
			}
		}
		if (!same_color)
		{
			break;
		}
	}
	if (i > len)
	{
		return true;
	}

	len = (int)(sqrt((l2.pt1.x-l2.pt2.x)*(l2.pt1.x-l2.pt2.x)+(l2.pt1.y-l2.pt2.y)*(l2.pt1.y-l2.pt2.y)) + 0.5);
	for (i = 0;i <= len;++i)
	{
		double lammda = i / (double)len;
		int x = (int)(lammda*(l2.pt1.x)+(1.0-lammda)*l2.pt2.x+0.5);
		int y = (int)(lammda*(l2.pt1.y)+(1.0-lammda)*l2.pt2.y+0.5);
		if (x <0 || x>= X || y < 0 || y >=Y)
		{
			continue;
		}

		bool same_color = true;
		for (int k = 0;k < band;++k)
		{
			if (bgColor[k] != img->imageData[ y*step + x*band + k])
			{
				same_color = false;
				break;
			}
		}
		if (!same_color)
		{
			break;
		}
	}
	if (i > len)
	{
		return true;
	}

	return false;
}


bool isBoundary(cv::Mat img, cvline_polar line, cv::Scalar bgColor/* = CV_RGB(0,0,0)*/)
{
	int X = img.cols;
	int Y = img.rows;
	int band = img.channels();
	int step = img.step;
	uchar * data = (uchar *)img.data;

	fPoint pt1 = line.pt1;
	fPoint pt2 = line.pt2;

	int k = 5;
	// 外扩k个像素
	cvline_polar l1 = segment_move_vertical(line, (double)k);
	cvline_polar l2 = segment_move_vertical(line, -(double)k);

	double len = (int)(sqrt((l1.pt1.x-l1.pt2.x)*(l1.pt1.x-l1.pt2.x)+(l1.pt1.y-l1.pt2.y)*(l1.pt1.y-l1.pt2.y)) + 0.5);
	int i;
	for (i = 0;i <= len;++i)
	{
		double lammda = i / (double)len;
		int x = (int)(lammda*(l1.pt1.x)+(1.0-lammda)*l1.pt2.x+0.5);
		int y = (int)(lammda*(l1.pt1.y)+(1.0-lammda)*l1.pt2.y+0.5);
		if (x <0 || x>= X || y < 0 || y >=Y)
		{
			continue;
		}

		bool same_color = true;
		for (int k = 0;k < band;++k)
		{
			if (bgColor.val[k] != data[ y*step + x*band + k])
			{
				same_color = false;
				break;
			}
		}
		if (!same_color)
		{
			break;
		}
	}
	if (i > len)
	{
		return true;
	}

	len = (int)(sqrt((l2.pt1.x-l2.pt2.x)*(l2.pt1.x-l2.pt2.x)+(l2.pt1.y-l2.pt2.y)*(l2.pt1.y-l2.pt2.y)) + 0.5);
	for (i = 0;i <= len;++i)
	{
		double lammda = i / (double)len;
		int x = (int)(lammda*(l2.pt1.x)+(1.0-lammda)*l2.pt2.x+0.5);
		int y = (int)(lammda*(l2.pt1.y)+(1.0-lammda)*l2.pt2.y+0.5);
		if (x <0 || x>= X || y < 0 || y >=Y)
		{
			continue;
		}

		bool same_color = true;
		for (int k = 0;k < band;++k)
		{
			if (bgColor[k] != img.data[ y*step + x*band + k])
			{
				same_color = false;
				break;
			}
		}
		if (!same_color)
		{
			break;
		}
	}
	if (i > len)
	{
		return true;
	}

	return false;
}

IplImage* get_lines(IplImage* img, vector<cvline_polar>& vec_lines, double scale/* = 0.8*/, int seletion_num/* = 200*/)
{
    //to grey
    //IplImage* grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    //cvCvtColor(img, grey, CV_BGR2GRAY);

	double * image;
	double * out;

	int X = img->width;
	int Y = img->height;

	/* create a simple image: left half black, right half gray */
	image = (double *) malloc( X * Y * sizeof(double) );
	if( image == NULL )
	{
		fprintf(stderr,"error: not enough memory\n");
		exit(EXIT_FAILURE);
	}

    for(int i=0;i</*grey*/img->width;i++)
	{
		for(int j=0;j</*grey*/img->height;j++)
		{
		  CvScalar s= cvGet2D(/*grey*/img,j,i);
		  double pix= s.val[0];
		  image[ i + j * X ]= pix;
		}
	}

    /* call LSD */
	int n;
	//out = lsd(&n, image, X, Y);
	out = lsd_scale(&n, image, X, Y, scale);

	// x1,y1,x2,y2,width,p,-log10(NFA);
    //out= lsd_scale(image,1);

    /* print output */
    //printf("%u line segments found:\n",out->size);
    //vector<Line> vec;

	double minLen, minNFA;
	//double *nfas = new double[n];
	vector<double> nfas;
	minLen = 5.0;
	vector<int> notBoundary;
	for(int i=0;i<n;i++)
	{
		cvline_polar outLine;
		fPoint line[2];
		line[0].x = out[ i * 7 + 0];
		line[0].y = out[ i * 7 + 1];
		line[1].x = out[ i * 7 + 2];
		line[1].y = out[ i * 7 + 3];

		outLine = line2polar(line);
		if (isBoundary(img, outLine, cvScalar((0,0,0))))
		{
			//continue;
		}
		notBoundary.push_back(i);
		double len = segment_length(outLine);
		if(len >= minLen) nfas.push_back(-out[i*7+6]);
	}

	minNFA = -quick_select(&nfas[0], 0, (int)nfas.size()-1, seletion_num);
	//for(int i=0;i<(int)nfas.size();i++)
	for(int m=0;m<(int)notBoundary.size();m++)
    {
		int i = notBoundary[m];
      //for(j=0;j<out->dim;j++)
      {
        //printf("%f ",out->values[ i * out->dim + j ]);

		  cvline_polar outLine;
		  fPoint line[2];
		  line[0].x = out[ i * 7 + 0];
		  line[0].y = out[ i * 7 + 1];
		  line[1].x = out[ i * 7 + 2];
		  line[1].y = out[ i * 7 + 3];

		  //double w = out[i*7+4];
		  //double p = out[i*7+5];
		  //double NFA = out[i*7+6];

		  outLine = line2polar(line);

		  if(out[i*7+6] >= minNFA)
		  {
			  vec_lines.push_back(outLine);
		  }
      }
      //printf("\n");
    }

	IplImage* black = cvCloneImage(img);
    //IplImage* black= cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    //cvZero(black);
    for(int i=0;i<(int)vec_lines.size();++i)
    {
        //if(vec[i].x1==vec[i].x2||vec[i].y1==vec[i].y2)
		//cvLine(black, cvPoint(vec_lines[i].pt1.x, vec_lines[i].pt1.y), cvPoint(vec_lines[i].pt2.x, vec_lines[i].pt2.y), CV_RGB(255, 255, 255), 2, CV_AA);
		cvLine(black, cvPoint(vec_lines[i].pt1.x, vec_lines[i].pt1.y), cvPoint(vec_lines[i].pt2.x, vec_lines[i].pt2.y), CV_RGB(255, 0, 0), 1, CV_AA);
	}

	/* free memory */
	free( (void *) image );
	free( (void *) out );

    return black;
}

IplImage* get_lines(const char* strImage, vector<cvline_polar>& vec_lines, double scale/* = 0.8*/, int seletion_num/* = 200*/)
{
    //to grey
    //IplImage* grey = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    //cvCvtColor(img, grey, CV_BGR2GRAY);
	
	IplImage* img = cvLoadImage(strImage, 0 );
	if (!img)
	{
		cout<<"open image failed for \""<<strImage<<"\"!"<<endl;
		return NULL;
	}

	double * image;
	double * out;

	int X = img->width;
	int Y = img->height;


	/* create a simple image: left half black, right half gray */
	image = (double *) malloc( X * Y * sizeof(double) );
	if( image == NULL )
	{
		fprintf(stderr,"error: not enough memory\n");
		exit(EXIT_FAILURE);
	}

    for(int i=0;i</*grey*/img->width;i++)
	{
		for(int j=0;j</*grey*/img->height;j++)
		{
		  CvScalar s= cvGet2D(/*grey*/img,j,i);
		  double pix= s.val[0];
		  image[ i + j * X ]= pix;
		}
	}

    /* call LSD */
	int n;
	//out = lsd(&n, image, X, Y);
	out = lsd_scale(&n, image, X, Y, scale);

	cout<<"all detected lines:"<<n<<endl;

	// x1,y1,x2,y2,width,p,-log10(NFA);
    //out= lsd_scale(image,1);

    /* print output */
    //printf("%u line segments found:\n",out->size);
    //vector<Line> vec;

	double minLen, minNFA;
	//double *nfas = new double[n];
	vector<double> nfas;
	minLen = 5.0;
	vector<int> notBoundary;
	for(int i=0;i<n;i++)
	{
		cvline_polar outLine;
		fPoint line[2];
		line[0].x = out[ i * 7 + 0];
		line[0].y = out[ i * 7 + 1];
		line[1].x = out[ i * 7 + 2];
		line[1].y = out[ i * 7 + 3];

		outLine = line2polar(line);
		if (isBoundary(img, outLine, cvScalar((0,0,0))))
		{
			continue;
		}
		notBoundary.push_back(i);
		double len = segment_length(outLine);
		if(len >= minLen) nfas.push_back(-out[i*7+6]);
	}
	if (nfas.size() < 3)
	{
		return NULL;
	}
	minNFA = -quick_select(&nfas[0], 0, (int)nfas.size()-1, seletion_num);
	//for(int i=0;i<(int)nfas.size();i++)
	for(int m=0;m<(int)notBoundary.size();m++)
    {
		int i = notBoundary[m];
      //for(j=0;j<out->dim;j++)
      {
        //printf("%f ",out->values[ i * out->dim + j ]);

		  cvline_polar outLine;
		  fPoint line[2];
		  line[0].x = out[ i * 7 + 0];
		  line[0].y = out[ i * 7 + 1];
		  line[1].x = out[ i * 7 + 2];
		  line[1].y = out[ i * 7 + 3];

		  //double w = out[i*7+4];
		  //double p = out[i*7+5];
		  //double NFA = out[i*7+6];

		  outLine = line2polar(line);

		  if(out[i*7+6] >= minNFA)
		  {
			  vec_lines.push_back(outLine);
		  }
      }
      //printf("\n");
    }

    //IplImage* black= cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);
    //cvZero(black);
	IplImage* line_img = cvLoadImage(strImage, 1 );
    for(int i=0;i<(int)vec_lines.size();++i)
    {
        //if(vec[i].x1==vec[i].x2||vec[i].y1==vec[i].y2)
        cvLine(line_img,cvPoint(vec_lines[i].pt1.x,vec_lines[i].pt1.y),cvPoint(vec_lines[i].pt2.x,vec_lines[i].pt2.y),CV_RGB(255,0,0),2, CV_AA);
	}

	/* free memory */
	free( (void *) image );
	free( (void *) out );
	cvReleaseImage(&img);

    return line_img;
}

bool GetImGray(IplImage* im, double x, double y, double &gray)
{
	if (x < 1 || x > (im->width-2) || y < 1 || y > (im->height-2))
		return false;
	int x1 = (int)x;
	int y1 = (int)y;
	int x2 = x1 + 1;
	int y2 = y1 + 1;
	int step = im->widthStep;
	int _x1 = x1 < 0 ? 0 : x1;
	int _y1 = y1 < 0 ? 0 : y1;
	int _x2 = x2 > (im->width-1) ? (im->width-1) : x2;
	int _y2 = y2 > (im->height-1) ? (im->height-1) : y2;
	gray = (x2 - x) * (y2 - y) * ((BYTE*)im->imageData)[_y1*step+_x1]/1.0 + 
		(x - x1) * (y2 - y) * ((BYTE*)im->imageData)[_y1*step+_x2]/1.0 + 
		(x2 - x) * (y - y1) * ((BYTE*)im->imageData)[_y2*step+_x1]/1.0 + 
		(x - x1) * (y - y1) * ((BYTE*)im->imageData)[_y2*step+_x2]/1.0;
	return true;
}