#ifndef ZY3_QVProcReader_H
#define ZY3_QVProcReader_H

#include "QVProcReader.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

class ZY3QVProcReader : public QVProcReader
{
public:
	ZY3QVProcReader();

	typedef struct _tagRpcStruct
	{
		double LINE_OFF;
		double SAMP_OFF;
		double LAT_OFF;
		double LONG_OFF;
		double HEIGHT_OFF;
		double LINE_SCALE;
		double SAMP_SCALE;
		double LAT_SCALE;
		double LONG_SCALE;
		double HEIGHT_SCALE;
		double LINE_NUM_COEFF[20];
		double LINE_DEN_COEFF[20];
		double SAMP_NUM_COEFF[20];
		double SAMP_DEN_COEFF[20];
	} RPC_STRUCT;

	typedef struct 
	{
		char format_id[64];			// 格式Id
		char version_id[16];		// 版本号
		char station_id[16];		//	接收该数据的地面站的标识
		char satellite_id[64];	    //	对应的卫星标识全称
		char sensor_id[512];		//	对应的传感器标识全称
		char jobTaskID [32];	    //	任务标识
		int	 orbit_num;   		//  轨道号
		int	 channel_num;	 	//  通道号
		int	 totalfiles;	 	    //  本任务单传输快视文件总数
		int	 files_num;	 	    //  本次传输的文件在所属任务中的文件序号
		int	 desample_num;	    //	进行快数处理时的降采样数,上一步的降采样数
		int image_mode;			// 图像模式。当此值=1时，该图像为灰度图像，否则为彩色合成图像，如：0
		_TIMESTAMP start_time;	// 数据接收的开始时间，为UTC时间
		_TIMESTAMP end_time;	// 数据接收的结束时间，为UTC时间
		LONLAT_POS station_pos;	//	对应的地面站在地球上的位置信息	地面站的位信息（经纬度）
	} QUIVIMAGE_HEAD_INFO;

	typedef struct _tagSceneGlobalInfo
	{// 景全局辅助信息
		char info_header[4];		// 景全局辅助信息帧头
		int	 scene_count;			// 景计数。该成像模式下景的计数值，从1开始，变换成像模式之后重新从1开始计数
		char satellite_id[64];		// 对应的卫星标识。卫星标识按“原始数据记录与交换格式”的规定
		char sensor_id[512];		// 对应的传感器标识。传感器标识按照“原始数据记录与交换格式”的规定
		int work_mode;				// 该景图像工作模式信息，针对同一传感器可变工作模式情况，以1开始的整数代替
		int scene_lines;			// 该景图像对应快视图像行数
		int scene_samples;			// 该景图像对应快视图像列数
		int qv_line_info_length;	// 快视图像行辅助信息总字节数 单位为字节
		int bands;					// 波段数
		int sample_bit_count;		// 量化位数
		char interleave[4];			// 多波段图像排列方式 针对不同工作模式，进行约定，分别为BSQ（Band Sequential format 波段顺序格式）、BIL（Band Interleaved by Line format，
									// 波段按像元交叉格式）、BIP（Band Interleaved by Pixel format，波段按行交叉格式），前填0x20
		_TIMESTAMP first_time;		// 景第一行图像数据起始时间 该时刻为成像时间 UTC
		_TIMESTAMP last_time;		// 景最后一行图像数据起始时间 该时刻为成像时间 UTC
		int image_format;			// 图像格式，采用整数表示，1：Raw，2：jpeg2000，3：tiff，4：gtif
		int compress_algorithm;		// 该景图像使用的压缩算法，采用数字表示，1：jpeg2000
		float compress_rate;		// 该景图像的压缩比率信息，如2，3，4
		unsigned long long scene_size;	// 景图像数据大小，单位为字节，如果图像格式为Raw则为压缩前数据大小，其他则为压缩后数据大小，小字节序
		LONLAT_POS left_upper;		// 左上角经纬度
		LONLAT_POS right_upper;		// 右上角经纬度
		LONLAT_POS left_lower;		// 左下角经纬度
		LONLAT_POS right_lower;		// 右下角经纬度
		LONLAT_POS scene_center;	// 景中心纬度
		RPC_STRUCT rpc_info;		// RPC系数相关信息
		char valid_flag[2];		// 原始数据行辅助信息标识，该值为=1时有效，代表有原始数据行辅助信息，可供存放在成像过程中产生的辅助信息
		int line_info_size;			// 原始数据行辅助信息大小
		int line_info_length;		// 原始数据行辅助信息总字节数
	}SCENE_GLOBAL_INFO;

	typedef struct _tagQuivAuxInfo
	{
		int line_index;		//	快视图像行辅助信息索引 快视图像行计数从1开始，该快视图像行信息索引对应于后续快视图像行数据
		LONLAT_POS		nadir_pos;	//	对应的星下点的位置
		LONLAT_POS		line_left_pos;//	该行图像左侧位置（经度和纬度）
		LONLAT_POS		line_right_pos;	//	该行图像右（经度和纬度）
		SAT_POS			satpos;     //卫星的x,y,z坐标以及vx,vy,vz
		int				att_model;	// 0: r,p,y   1: q1,q2,q3
		SAT_ATT			satatt;     //卫星的姿态参数,俯仰滚动偏航，r,p,y
		double		line_time;			//	该时刻为成像时间，UTC
		double		att_time;			//	，UTC
		double		star_time;			//	，UTC
		//_TIMESTAMP		line_time;			//	该时刻为成像时间，UTC
	}QUIVIMAGE_AUX_INFO;

	bool readHeader(const char* header);
	bool extractHeader(const char* filename);

	bool read(const char* filename, const char* outPath);
	bool read_by_scene(const char* filename, const char* outPath);

	bool readHeadInfo(FILE* pfIn, QUIVIMAGE_HEAD_INFO& header, FILE* pfOut = NULL);

	bool readHeadInfo(ifstream& fs, QUIVIMAGE_HEAD_INFO& header);
	bool readTimeInfo(FILE* pfIn, _TIMESTAMP& t, FILE* pfOut = NULL);
	bool readTimeInfo(ifstream& fs, _TIMESTAMP& t);
	bool readSAT_POS(FILE* pfIn, SAT_POS& pos, FILE* pfOut = NULL);
	bool readSAT_POS(ifstream& fs, SAT_POS& pos);
	bool readSAT_ATT(FILE* pfIn, SAT_ATT& att, FILE* pfOut = NULL);
	bool readSAT_ATT(ifstream& fs, SAT_ATT& att);
	bool readSceneGlobalInfo(FILE* pfIn, SCENE_GLOBAL_INFO& scene_info, FILE* pfOut = NULL);
	bool readSceneGlobalInfo(ifstream& fs, SCENE_GLOBAL_INFO& scene_info);
	bool readRPC_INFO(FILE* pfIn, RPC_STRUCT& rpc_info, FILE* pfOut = NULL);
	bool readRPC_INFO(ifstream& fs, RPC_STRUCT& rpc_info);
	bool readQuivAuxInfo(FILE* pfIn, QUIVIMAGE_AUX_INFO& aux, FILE* pfOut = NULL);
	bool readQuivAuxInfo(ifstream& fs, QUIVIMAGE_AUX_INFO& aux);

	bool writeRPC(const char* filename, const RPC_STRUCT& rpc_info);


	void writeDescXML(const QUIVIMAGE_HEAD_INFO& headerInfo, const char* descXML)const;
	void writeDescXML(const QUIVIMAGE_HEAD_INFO& headerInfo,
		const SCENE_GLOBAL_INFO& sceneInfo,
		const vector<QUIVIMAGE_AUX_INFO>& auxList,
		const char* descXML,
		int line_offset = 0)const;

	bool time_ready;

protected:
	QUIVIMAGE_HEAD_INFO theHeaderInfo;
	int theLines;
	int theSamples;
	int theBands;
};

#endif