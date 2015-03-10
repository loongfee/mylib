#include "lineFeature.h"
#include "math.h"
void LineFeature(IplImage* img,vector<cvline_polar>& lines_polar_list,vector<vector<double>> &AllLBD,vector<int>&LineNum,const int R)
{
	int i;//检测出的直线的下标
	int j;//检测出的直线上的点的下标
	int r;
	int m;//检测出的直线的平行邻域内直线的下标
	//const int R=23;//平行邻域的大小
	const int sz=5;//子区域的大小
	int M=2*R+1;
    int NumSubRegion=(2*R-1)/sz;
	int Mspt=NumSubRegion;
	int N;
	vector <vector<double>> Feature;//所有直线的灰度矩阵
	vector<vector<double>> LinesFeature;//所有直线的梯度描述矩阵
	vector<int>LineLength;
	for(i=0;i!=int(lines_polar_list.size());i++)
{       
	if( (lines_polar_list[i].pt1.x>=R)&&(lines_polar_list[i].pt1.x<=(img->width-R))&&
		(lines_polar_list[i].pt1.y>=R)&& (lines_polar_list[i].pt1.y<=(img->height-R))&&
		(lines_polar_list[i].pt2.x>=R)&&(lines_polar_list[i].pt2.x<=(img->width-R))&&
	   (lines_polar_list[i].pt2.y>=R)&& (lines_polar_list[i].pt2.y<=(img->height-R)))
 {
	    N=abs(floor(lines_polar_list[i].pt1.x)-floor(lines_polar_list[i].pt2.x) );//检测出的直线上点的个数
	    //N=sqrt((lines_polar_list[i].pt1.x-lines_polar_list[i].pt2.x)*(lines_polar_list[i].pt1.y-lines_polar_list[i].pt2.y));
		if(N>=6)
 {
		LineNum.push_back (i);
	    LineLength.push_back(N);
		vector<fPoint> PointInLine;//检测出的直线上的点
		vector<fPoint> PointsInRegion;//存放平行邻域内各点坐标的数组
		vector<double>FeatureValue;//平行域内各点灰度值
		double x0,y0,xj,yj,dx,dy;
	    if (lines_polar_list[i].theta<PI/2&&lines_polar_list[i].theta>=0)//theta在0~π/2之间时
    {  
		x0=(lines_polar_list[i].pt1.x<lines_polar_list[i].pt2.x)?lines_polar_list[i].pt1.x:lines_polar_list[i].pt2.x;//从具有较小x、y坐标的点开始
		y0=(lines_polar_list[i].pt1.y<lines_polar_list[i].pt2.y)?lines_polar_list[i].pt1.y:lines_polar_list[i].pt2.y;
		PointInLine.push_back(fPoint(x0,y0));
		for(j=1;j!=N;j++)
		{
			
			xj=PointInLine[0].x+j;
			yj=PointInLine[0].y+j*tan(lines_polar_list[i].theta);
			//xj=PointInLine[0].x+j*cos(lines_polar_list[i].theta);
			//yj=PointInLine[0].y+j*sin(lines_polar_list[i].theta);
			PointInLine.push_back(fPoint(xj,yj));
		}	
		 for(r=R;r>=-R;r--) //获得平行邻域内各点坐标值
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x-r*sin(lines_polar_list[i].theta);
			   dy=PointInLine[j].y+r*cos(lines_polar_list[i].theta);
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

	else if (lines_polar_list[i].theta>-PI/2&&lines_polar_list[i].theta<0)//theta在-π/2~0之间时
	{
		x0=(lines_polar_list[i].pt1.x>lines_polar_list[i].pt2.x)?lines_polar_list[i].pt1.x:lines_polar_list[i].pt2.x;
		y0=(lines_polar_list[i].pt1.y<lines_polar_list[i].pt2.y)?lines_polar_list[i].pt1.y:lines_polar_list[i].pt2.y;
		PointInLine.push_back(fPoint(x0,y0));
		for(j=1;j!=N;j++)
		{
			xj=PointInLine[0].x-j;
			yj=PointInLine[0].y-j*tan(lines_polar_list[i].theta);
			//xj=PointInLine[0].x+j*cos(lines_polar_list[i].theta);
			//yj=PointInLine[0].y+j*sin(lines_polar_list[i].theta);
			PointInLine.push_back(fPoint(xj,yj));
		}
		 for(r=R;r>=-R;r--) //获得平行邻域内各点坐标值
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x+r*sin(lines_polar_list[i].theta);
			   dy=PointInLine[j].y-r*cos(lines_polar_list[i].theta);
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

	else                                                                                             //theta等于π/2时
	{
        x0=lines_polar_list[i].pt1.x;
	    y0=(lines_polar_list[i].pt1.y<lines_polar_list[i].pt2.y)?lines_polar_list[i].pt1.y:lines_polar_list[i].pt2.y;//从直线中具有较小y坐标的点开始
		PointInLine.push_back(fPoint(x0,y0));
		for(j=1;j!=N;j++)
		{
           xj=PointInLine[0].x;
		   yj=PointInLine[0].y+j;
		   PointInLine.push_back(fPoint(xj,yj));
		}	
		 for(r=R;r>=-R;r--) //获得平行邻域内各点坐标值
		{
		   for(j=0;j!=N;j++)
		   {
			   dx=PointInLine[j].x-r;
			   dy=PointInLine[j].y;
			   PointsInRegion.push_back(fPoint(dx,dy));
		   }
		}
	}

		  for(m=0;m!=M;m++)//获得每个点的灰度值
		{
		  for(j=0;j!=N;j++)
		   {
		      CvScalar s= cvGet2D(img,PointsInRegion[j +m*N].y,PointsInRegion[j +m*N].x);
		      double pix= s.val[0];
		      FeatureValue.push_back(pix);
		   }
		}
		Feature.push_back(FeatureValue);
		PointInLine.swap(vector<fPoint>());
		PointsInRegion.swap(vector<fPoint>());
		FeatureValue.swap(vector<double>());
	}
  }
	else
	{;}
}

	int dg=0;
	double thetag=R;
	double wg=1/(sqrt(2*PI)*thetag)*exp(-dg*dg/(2*thetag*thetag));//全局高斯权重系数
    //构造MSLD描述子，对灰度矩阵中的(M-1)*(N-1)求梯度
	vector <vector<double>> All_Grad_dL;//dL方向的梯度
	for(i=0;i!=Feature.size();i++)
	{
		vector<double>Each_Grad_dL;
		N=LineLength[i];
		for(m=0;m!=M-1;m++)
		{   
			dg=abs(R-m);
			for(j=0;j!=N-1;j++)
			{
				Each_Grad_dL.push_back(wg* (Feature[i][j+1+m*N]-Feature[i][j+m*N]));
			}
		}
		All_Grad_dL.push_back (Each_Grad_dL);
		Each_Grad_dL.swap(vector<double>());
	}
	vector <vector<double>> All_Grad_dT;//dT方向的梯度
	for(i=0;i!=Feature.size();i++)
	{
		vector<double>Each_Grad_dT;
		N=LineLength[i];
		for(j=0;j!=N-1;j++)
		{
			for(m=0;m!=M-1;m++)
			{
				dg=abs(R-m);
				Each_Grad_dT.push_back( wg*(Feature[i][j+(m+1)*N]-Feature[i][j+m*N]));
			}
		}
		All_Grad_dT.push_back (Each_Grad_dT);
		Each_Grad_dT.swap(vector<double>());
	}
	vector<vector<double>> All_Grad_dTL;//沿dT方向计算的梯度按dL方向排列，与dL方向的梯度对应。
	for(i=0;i!=Feature.size();i++)
	{
	     vector<double> Each_Grad_dTL;
		 N=LineLength[i];
	     for(m=0;m!=M-1;m++)
		  {
			for(j=0;j!=N-1;j++)
				{
					Each_Grad_dTL.push_back (All_Grad_dT[i][m+j*(M-1)]);
			    } 
	      }
	     All_Grad_dTL.push_back(Each_Grad_dTL);
	     Each_Grad_dTL.swap(vector<double>());
	}
	    All_Grad_dT.swap(vector<vector<double>>()); //得到两个方向上的梯度All_Grad_dL、All_Grad_dTL
      
	int dl=0;
	double thetal=5;
    double wl=1/(sqrt(2*PI)*thetal)*exp(-dl*dl/(2*thetal*thetal));
	//vector<vector<double>>AllLBD;
   for(i=0;i!=Feature.size();i++)
	{   
		vector<double>LBD;
		N=LineLength[i];
		for(int k=0;k<9;k++)
		{   
			if(k==0)
			{
				vector <double>BD1;
				for(m=1;m!=11;m++)
				{
					dl=abs(m-3);
					double SumV1=0.0;
			        double SumV2=0.0;
			        double SumV3=0.0;
			       double SumV4=0.0;
					for(j=0;j!=N-1;j++)
				    {
					   All_Grad_dL[i][j+m*(N-1)]=wl*All_Grad_dL[i][j+m*(N-1)];
					   All_Grad_dTL[i][j+m*(N-1)]=wl* All_Grad_dTL[i][j+m*(N-1)];
					    if(All_Grad_dTL[i][j+m*(N-1)]>=0)
						  {SumV1+=All_Grad_dTL[i][j+m*(N-1)];  }
					    if(All_Grad_dTL[i][j+m*(N-1)]<0)
						  { SumV2+=abs(All_Grad_dTL[i][j+m*(N-1)]);}
						if(All_Grad_dL[i][j+m*(N-1)]>=0)
						  {  SumV3+=All_Grad_dL[i][j+m*(N-1)];}
					    if(All_Grad_dL[i][j+m*(N-1)]<0)
						  {SumV4+=abs(All_Grad_dL[i][j+m*(N-1)]);}
				    }
					BD1.push_back(SumV1);
			        BD1.push_back(SumV2);
			        BD1.push_back(SumV3);
			        BD1.push_back(SumV4);
				}
            for(int v=0;v!=4;v++)
			 {   
				 double sum=0.0;
				 for(int n=v;n<40;n+=4)
				 {
					 sum+=BD1[n];
				 }
				 double mean=sum/10;
				 LBD.push_back(mean);
				 double sumsqr=0.0;
				 for(int n=v;n<40;n+=4)
				 {
					 sumsqr+=(BD1[n]-mean)* (BD1[n]-mean);
				 }
				double std=sqrt( sumsqr);
				LBD.push_back(std);
			 }
			 BD1.swap(vector<double>());
		}
			else if(k==8)
			{
				vector<double>BD9;
				for(m=36;m!=46;m++)
				{
					dl=abs(m-43);
					double SumV1=0.0;
			        double SumV2=0.0;
			        double SumV3=0.0;
			        double SumV4=0.0;
					for(j=0;j!=N-1;j++)
				    {
					   All_Grad_dL[i][j+m*(N-1)]=wl*All_Grad_dL[i][j+m*(N-1)];
					   All_Grad_dTL[i][j+m*(N-1)]=wl* All_Grad_dTL[i][j+m*(N-1)];
					    if(All_Grad_dTL[i][j+m*(N-1)]>=0)
						  {SumV1+=All_Grad_dTL[i][j+m*(N-1)];  }
					    if(All_Grad_dTL[i][j+m*(N-1)]<0)
						  { SumV2+=abs(All_Grad_dTL[i][j+m*(N-1)]);}
						if(All_Grad_dL[i][j+m*(N-1)]>=0)
						  {  SumV3+=All_Grad_dL[i][j+m*(N-1)];}
					    if(All_Grad_dL[i][j+m*(N-1)]<0)
						  {SumV4+=abs(All_Grad_dL[i][j+m*(N-1)]);}
				    }
					BD9.push_back(SumV1);
			        BD9.push_back(SumV2);
			        BD9.push_back(SumV3);
			        BD9.push_back(SumV4);
				}
            for(int v=0;v!=4;v++)
			 {   
				 double sum=0.0;
				 for(int n=v;n<40;n+=4)
				 {
					 sum+=BD9[n];
				 }
				 double mean=sum/10;
				 LBD.push_back(mean);
				 double sumsqr=0.0;
				 for(int n=v;n<40;n+=4)
				 {
					 sumsqr+=(BD9[n]-mean)* (BD9[n]-mean);
				 }
				double std=sqrt( sumsqr);
				LBD.push_back(std);
			 }
			 BD9.swap(vector<double>());
			}
			else
			{
			 vector<double>BDM;
             for(m=5*k-4;m!=5*k+11;m++)
		       {
				  dl=abs(m-(5*k+3));
				  double SumV1=0.0;
			      double SumV2=0.0;
			      double SumV3=0.0;
			      double SumV4=0.0;
				   for(j=0;j!=N-1;j++)
				   {
					   All_Grad_dL[i][j+m*(N-1)]=wl*All_Grad_dL[i][j+m*(N-1)];
					   All_Grad_dTL[i][j+m*(N-1)]=wl* All_Grad_dTL[i][j+m*(N-1)];
					    if(All_Grad_dTL[i][j+m*(N-1)]>=0)
						  {SumV1+=All_Grad_dTL[i][j+m*(N-1)];  }
					    if(All_Grad_dTL[i][j+m*(N-1)]<0)
						  { SumV2+=abs(All_Grad_dTL[i][j+m*(N-1)]);}
						if(All_Grad_dL[i][j+m*(N-1)]>=0)
						  {  SumV3+=All_Grad_dL[i][j+m*(N-1)];}
					    if(All_Grad_dL[i][j+m*(N-1)]<0)
						  {SumV4+=abs(All_Grad_dL[i][j+m*(N-1)]);}
				   }
              BDM.push_back(SumV1);
			  BDM.push_back(SumV2);
			  BDM.push_back(SumV3);
			  BDM.push_back(SumV4);
		       }
			 for(int v=0;v!=4;v++)
			 {   
				 double sum=0.0;
				 for(int n=v;n<60;n+=4)
				 {
					 sum+=BDM[n];
				 }
				 double mean=sum/15;
				 LBD.push_back(mean);
				 double sumsqr=0.0;
				 for(int n=v;n<60;n+=4)
				 {
					 sumsqr+=(BDM[n]-mean)* (BDM[n]-mean);
				 }
				double std=sqrt( sumsqr);
				LBD.push_back(std);
			 }
			 BDM.swap(vector<double>());
			}
		}
		AllLBD.push_back(LBD);
		LBD.swap(vector<double>());
    }

	 for(i=0;i!=Feature.size();i++)
	 {
		 for(m=0;m!=2;m++)
	  {
		  double summs=0.0;
		 for(j=m;j<72;j+=2)
		 {summs+=AllLBD[i][j]*AllLBD[i][j];}
		 double len=sqrt(summs);
		 for(j=m;j<72;j+=2)
		 {AllLBD[i][j]=AllLBD[i][j]/len;}
	  }
	 }

	  for(i=0;i!=Feature.size();i++)
	  {
		  for(j=0;j!=72;j++)
		  if(AllLBD[i][j]>0.4)
		  {AllLBD[i][j]=0.4;}
	  }

	   for(i=0;i!=Feature.size();i++)
	   {
		   double renom=0.0;
		   for(j=0;j!=72;j++)
		   {
			   renom+=AllLBD[i][j]*AllLBD[i][j];
		   }
		   double relen=sqrt(renom);
		   for(j=0;j!=72;j++)
		   {
			   AllLBD[i][j]=AllLBD[i][j]/relen;
		   }
	   }
}

