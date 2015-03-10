#ifndef HJ1QVProcReader_H
#define HJ1QVProcReader_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "pugixml\pugiconfig.hpp"
#include "pugixml\pugixml.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
using namespace std;

namespace QVProc{
class HJ1QVProcReader
{
public:
	HJ1QVProcReader();

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
		int			sec;
		int			nsec;
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
		float roll;
		float pitch;
		float yaw;
		float vroll;
		float vpitch;
		float vyaw;
	};

	typedef struct 
	{
		char  station_id[16];		//	���ո����ݵĵ���վ�ı�ʶ
		char	 satellite_id[16];	    //	��Ӧ�����Ǳ�ʶ
		char   sensor_id[16];		//	��Ӧ�Ĵ�������ʶ
		char	 work_mode[16];		//	��Ӧ�Ĺ���ģʽ��Ϣ
		char	 JobTaskID [20];	    //	�����ʶ
		int	 orbit_num;   		//  �����
		int	 channel_num;	 	//  ͨ����
		int	 totalfiles;	 	    //  �����񵥴�������ļ�����
		int	 files_num;	 	    //  ���δ�����ļ������������е��ļ����
		int	 desample_num;	    //	���п�������ʱ�Ľ�������,��һ���Ľ�������
		int	 data_width;			//	��һ���������Ŀ���,(����������������Ļ�),������һ����С
		int	 sample_bit_count;	//	���Ӵ������ͼ��������������ֻ��8Bit��16Bit�������
		int	 gray_image_flag;	//	����ֵΪ��ʱ����ͼ��Ϊ��ͨ���Ҷ�ͼ�񣬷���Ϊ��ͨ��α��ɫͼ��
		_TIMESTAMP		start_time;	//	���ݽ��յĿ�ʼʱ��
		_TIMESTAMP		end_time;		//	���ݽ��յĽ���ʱ��
		LONLAT_POS		station_pos;	//	��Ӧ�ĵ���վ�ڵ����ϵ�λ����Ϣ	����վ��λ��Ϣ

		int  sample_num;
		int  line_num;
		int  band_num;
	} QUIVIMAGE_HEAD_INFO;

	typedef struct _tagQuivAuxInfo
	{
		int		   valid_flag;		//	��ֵΪ=1ʱ��Ч����������¸�����ϢΪ��Чֵ
		LONLAT_POS		nadir_pos;	//	��Ӧ�����µ��λ��
		LONLAT_POS		line_left_pos;//	�þ����߶�Ӧ������λ�ã����Ⱥ�γ�ȣ�
		LONLAT_POS		line_right_pos;	//	�þ����߶�Ӧ���Ҳ��λ�ã����Ⱥ�γ�ȣ�
		SAT_POS			satpos;     //���ǵ�x,y,z�����Լ�vx,vy,vz
		SAT_ATT			satatt;     //���ǵ���̬����,��������ƫ����r,p,y
		unsigned short line_count;
		_TIMESTAMP		line_time;			//	����������ʱ����Ϣ
		double      line_num;
		double		att_time;			//	��̬���ݵ�ʱ����Ϣ
		double		gps_time;			//	GPS��λ����ʱ��
		double      star_time;			//  ���������ʱ����
		double      gps_c_time;			//  GPS��ʱ��
	}QUIVIMAGE_AUX_INFO;

	bool read(const char* filename, const char* outPath);
	bool read_band_seperate(const char* filename, const char* outPath);
	bool read_by_scene(const char* filename, const char* outPath);
	bool read_by_scene_band_seperate(const char* filename, const char* outPath);

	bool readHeadInfo(ifstream& fs, QUIVIMAGE_HEAD_INFO& header);
	bool readTimeInfo(ifstream& fs, _TIMESTAMP& t);
	bool readSAT_POS(ifstream& fs, SAT_POS& pos);
	bool readSAT_ATT(ifstream& fs, SAT_ATT& att);
	bool readQuivAuxInfo(ifstream& fs, QUIVIMAGE_AUX_INFO& aux);

	int timeCompare(const _TIMESTAMP& t1, const _TIMESTAMP& t2)const;

	void writeDescXML(const QUIVIMAGE_HEAD_INFO& theHeaderInfo, const vector<QUIVIMAGE_AUX_INFO>& auxList, const char* descXML, int line_offset = 0)const;
	string time2String(const _TIMESTAMP& t)const;

	bool time_ready;

	void BitInverse(char* chrs, int nBit);

	template <class T>
	void BitInverse(T &data)
	{
		int nBit = sizeof(T);
		BitInverse((char*)&data, nBit);
	};

protected:
	//int theTotalLines;
	//int theSamples;
	//int theBands;
	double theReferenceTime;
	double theReferenceLine;
};
} // end of namespace QVHJ

#endif