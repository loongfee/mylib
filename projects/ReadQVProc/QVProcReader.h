#ifndef QVProcReader_H
#define QVProcReader_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "pugixml\pugiconfig.hpp"
#include "pugixml\pugixml.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
using namespace std;

static void posEci2Ecr(double& x, double &y, double &z, double t)
{
	//ω=7292115×10-11rads-1±0.150×10-11rads-1
	double omega = 7.292115e-5;
	double a = t * omega;
	double sa = sin(a);
	double ca = cos(a);

	double xx = ca *x + sa * y;
	double yy = -sa * x + ca * y;

	x = xx;
	y = yy;
}

static void velEci2Ecr(double x, double y, double z,
	double &vx, double &vy, double &vz, double t)
{
	//ω=7292115×10-11rads-1±0.150×10-11rads-1
	double omega = 7.292115e-5;
	double a = t * omega;
	double sa = sin(a);
	double ca = cos(a);

	double vxx = ca *vx + sa * vy;
	double vyy = -sa * vx + ca * vy;

	double xx = -omega*sa*x + omega*ca * y;
	double yy = -omega*ca * x - omega*sa * y;

	vx = vxx + xx;
	vy = vyy + yy;
}

static void ecef2lla(double x, double y, double z,
	double& lat, double& lon, double& alt)
{
	double b, ep, p, th, n;

	/* Constants (WGS ellipsoid) */
	const double a = 6378137;
	const double e = 8.1819190842622e-2;
	const double pi = 3.141592653589793;

	/* Calculation */
	b = sqrt(pow(a, 2)*(1 - pow(e, 2)));
	ep = sqrt((pow(a, 2) - pow(b, 2)) / pow(b, 2));
	p = sqrt(pow(x, 2) + pow(y, 2));
	th = atan2(a*z, b*p);
	lon = atan2(y, x);
	lat = atan2((z + ep*ep*b*pow(sin(th), 3)), (p - e*e*a*pow(cos(th), 3)));
	n = a / sqrt(1 - e*e*pow(sin(lat), 2));
	alt = p / cos(lat) - n;
	lat = (lat * 180) / pi;
	lon = (lon * 180) / pi;
}

class QVProcReader
{
public:
	QVProcReader();

	typedef struct _tagBL
	{
		float		longitude;
		float		latitude;
	}LONLAT_POS;

	typedef struct _tagTimeInfo
	{
		int			year;
		int			mon;
		int			day;
		int			hour;
		int			min;
		int			sec;		// 秒
		int			nsec;
		int			millsec;	// 毫秒
		int			microsec;	// 微秒
	}_TIMESTAMP;

	struct SAT_POS
	{
		float x;
		float y;
		float z;
		float vx;
		float vy;
		float vz;
	};

	struct SAT_ATT
	{
		float roll;		//q1
		float pitch;	//q2
		float yaw;		//q3
		float vroll;
		float vpitch;
		float vyaw;
	};

	typedef struct 
	{
		char  station_id[16];		//	接收该数据的地面站的标识
		char	 satellite_id[16];	    //	对应的卫星标识
		char   sensor_id[16];		//	对应的传感器标识
		char	 work_mode[16];		//	对应的工作模式信息
		char	 JobTaskID [20];	    //	任务标识
		int	 orbit_num;   		//  轨道号
		int	 channel_num;	 	//  通道号
		int	 totalfiles;	 	    //  本任务单传输快视文件总数
		int	 files_num;	 	    //  本次传输的文件在所属任务中的文件序号
		int	 desample_num;	    //	进行快数处理时的降采样数,上一步的降采样数
		int	 data_width;			//	上一步传进来的宽度,(如果我在做降采样的话),这里会进一步变小
		int	 sample_bit_count;	//	快视处理后的图像量化比特数，只有8Bit和16Bit两种情况
		int	 gray_image_flag;	//	当此值为真时，该图像为单通道灰度图像，否则为多通道伪彩色图像
		_TIMESTAMP		start_time;	//	数据接收的开始时间
		_TIMESTAMP		end_time;		//	数据接收的结束时间
		LONLAT_POS		station_pos;	//	对应的地面站在地球上的位置信息	地面站的位信息

		int  sample_num;
		int  line_num;
		int  band_num;
	} QUIVIMAGE_HEAD_INFO;

	typedef struct _tagQuivAuxInfo
	{
		int		   valid_flag;		//	该值为=1时有效，其他情况下辅助信息为无效值
		LONLAT_POS		nadir_pos;	//	对应的星下点的位置
		LONLAT_POS		line_left_pos;//	该距离线对应的左侧的位置（精度和纬度）
		LONLAT_POS		line_right_pos;	//	该距离线对应的右侧的位置（精度和纬度）
		SAT_POS			satpos;     //卫星的x,y,z坐标以及vx,vy,vz
		SAT_ATT			satatt;     //卫星的姿态参数,俯仰滚动偏航，r,p,y
		unsigned short line_count;
		_TIMESTAMP		line_time;			//	照射该脉冲的时间信息
		double      line_num;
		double		att_time;			//	姿态数据的时间信息
		double		gps_time;			//	GPS定位数据时间
		double      star_time;			//  整星星务对时数据
		double      gps_c_time;			//  GPS对时码
	} QUIVIMAGE_AUX_INFO;

	virtual bool read(const char* filename, const char* outPath);
	virtual bool read_band_seperate(const char* filename, const char* outPath);
	virtual bool read_by_scene(const char* filename, const char* outPath);
	virtual bool read_by_scene_band_seperate(const char* filename, const char* outPath);

	bool readHeadInfo(ifstream& fs, QUIVIMAGE_HEAD_INFO& header);
	bool readTimeInfo(ifstream& fs, _TIMESTAMP& t);
	bool readSAT_POS(ifstream& fs, SAT_POS& pos);
	bool readSAT_ATT(ifstream& fs, SAT_ATT& att);
	bool readQuivAuxInfo(ifstream& fs, QUIVIMAGE_AUX_INFO& aux);

	int timeCompare(const _TIMESTAMP& t1, const _TIMESTAMP& t2)const;

	virtual void writeDescXML(const QUIVIMAGE_HEAD_INFO& theHeaderInfo, const char* descXML)const;
	virtual void writeDescXML(const QUIVIMAGE_HEAD_INFO& theHeaderInfo, const vector<QUIVIMAGE_AUX_INFO>& auxList, const char* descXML, int line_offset = 0)const;
	
	bool time_ready;

	void BitInverse(char* chrs, int nBit);

	template <class T>
	void BitInverse(T &data)
	{
		int nBit = sizeof(T);
		BitInverse((char*)&data, nBit);
	};

	string time2String(const _TIMESTAMP& t)const;
	void convertTimeStamp(const string& time_stamp,
		double& ti) const;
	double convertTimeStamp(const string& time_stamp) const;

protected:
	//int theTotalLines;
	//int theSamples;
	//int theBands;
	double theReferenceTime;
	double theReferenceLine;
};

#endif