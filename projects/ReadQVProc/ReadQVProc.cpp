#include <fstream>
#include <iostream>
#include <vector>
#include <string>

/////////////// gdal
#include "gdal_priv.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_multiproc.h"
#include "ogrsf_frmts.h"


#include "QVProcReader.h"
#include "HJ1QVProcReader.h"
#include "ZY3QVProcReader.h"
#include "CBERS04QVProcReader.h"

#include <QDir>

#include <argtable2/argtable2.h>
#include <assert.h>
#include <boost/filesystem.hpp>

/* for memory leak debugging */
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#include "strUtil.h"
#pragma comment(lib, "Qt5Core.lib")
#pragma comment(lib, "Qt5Gui.lib")
#pragma comment(lib, "gdal_i.lib")
#pragma comment(lib, "argtable2.lib")

using namespace std;
namespace fs = boost::filesystem;

enum sensorType{
	HJ1 = 0,
	ZY3,
	UnknownSensorType,
};
const char* strHJ1 = "HJ1";
const char* strZY3 = "ZY3";
const char* strCBERS04 = "CBERS04";
//const char* pszInputFile = "";
//const char* pszOutputDir = "";
string strInputFile = "";
string strOutputDir = "";
//const char* sensorType = "";


///*************************************************************************/
///*                               Usage()                                */
///************************************************************************/
//
//static void Usage(const char* pszErrorMsg = NULL)
//
//{
//	printf(
//		"Usage: ReadQVProc -i inputFile -t sensorType{HJ1|ZY3}\n");
//
//	if (pszErrorMsg != NULL)
//		fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);
//	exit(1);
//}
//
//#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
//	do { if (i + nExtraArg >= argc) \
//	Usage(CPLSPrintf("%s option requires %d argument(s)", argv[i], nExtraArg)); } while(0)


int mymain(const char** sensorType,
	const char **infile,
	const char **outpath)
{
	if (0 != _strcmpi(sensorType[0], strHJ1)
		&& 0 != _strcmpi(sensorType[0], strZY3)
		&& 0 != _strcmpi(sensorType[0], strCBERS04))
	{
		printf("Unknown Sensor Type.\n");
		return 0;
	}

	strInputFile = infile[0];
	if (0 == _strcmpi(strInputFile.c_str(), ""))
	{
		printf("inputFile can not be empty!\n");
		return 0;
	}

	if (!fs::exists(strInputFile))
	{
		printf("inputFile does not exist!\n");
	}
	else
	{
		fs::path inputPathname = fs::absolute(strInputFile);
		//char fullPath[_MAX_PATH];
		//_fullpath(fullPath, pszInputFile, _MAX_PATH);
		
		if (0 == _strcmpi(outpath[0], ""))
		{
			strOutputDir = inputPathname.parent_path().string();
			//pszOutputDir = SBeforeLast(string(fullPath), '\\').c_str();
		}
		else
		{
			strOutputDir = fs::absolute(outpath[0]).string();
		}

		if (!fs::exists(strOutputDir))
		{
			_mkdir(strOutputDir.c_str());
		}

		QVProcReader* qvreader = NULL;
		clock_t  clockBegin, clockEnd;
		clockBegin = clock();
		if (0 == _strcmpi(sensorType[0], strHJ1))
		{
			// HJ1
			qvreader = new HJ1QVProcReader;
			qvreader->read_by_scene(inputPathname.string().c_str(), strOutputDir.c_str());
		}
		else if (0 == _strcmpi(sensorType[0], strZY3))
		{
			// ZY3
			qvreader = new ZY3QVProcReader;
			qvreader->read_by_scene(inputPathname.string().c_str(), strOutputDir.c_str());
		}
		else if (0 == _strcmpi(sensorType[0], strCBERS04))
		{
			// CBERS04
			qvreader = new CBERS04QVProcReader;
			qvreader->read_by_scene(inputPathname.string().c_str(), strOutputDir.c_str());
		}
		else
		{
			// Unknown Type
			printf("Unknown Sensor Type.\n");
			return 0;
		}

		clockEnd = clock();
		printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
	}
	return 0;
}

int main(int argc, char** argv)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO"); // gdal 中文路径支持
	CPLSetConfigOption("SHAPE_ENCODING", "");	// shapefile 中文字段支持

	struct arg_str  *sensorType = arg_strn("t", "sensortype", NULL, 0, argc + 2, "sensor type {HJ1|ZY3|CBERS04} (default is \"HJ1\")");
	struct arg_lit  *help = arg_lit0("h", "help", "print this help and exit");
	struct arg_lit  *version = arg_lit0("v", "version", "print version information and exit");
	struct arg_file *infile = arg_file1(NULL, NULL, NULL, "input file");
	struct arg_file *outpath = arg_file0(NULL, NULL, "<output>", "output path");
	struct arg_end  *end = arg_end(20);


	//struct arg_lit  *list = arg_lit0("lL", NULL, "list files");
	//struct arg_lit  *recurse = arg_lit0("R", NULL, "recurse through subdirectories");
	//struct arg_int  *repeat = arg_int0("k", "scalar", NULL, "define scalar value k (default is 3)");
	//struct arg_str  *defines = arg_strn("D", "define", "MACRO", 0, argc + 2, "macro definitions");
	////struct arg_file *outfile = arg_file0("o", NULL, "<output>", "output file (default is \"-\")");
	//struct arg_lit  *verbose = arg_lit0("v", "verbose,debug", "verbose messages");
	//struct arg_lit  *help = arg_lit0(NULL, "help", "print this help and exit");
	//struct arg_lit  *version = arg_lit0(NULL, "version", "print version information and exit");
	//struct arg_file *infiles = arg_filen(NULL, NULL, "<input>", 1, argc + 2, "input file(s)");
	//struct arg_file *outfile = arg_file0(NULL, NULL, "<output>", "output file");
	//struct arg_end  *end = arg_end(20);
	////void* argtable[] = { list, recurse, repeat, defines, outfile, verbose, help, version, infiles, end };
	//void* argtable[] = { list, recurse, repeat, defines, verbose, help, version, infiles, outfile, end };
	void* argtable[] = { sensorType, help, version, infile, outpath, end };
	const char* progname = "ReadQVProc";
	int nerrors;
	int exitcode = 0;

	/* verify the argtable[] entries were allocated sucessfully */
	if (arg_nullcheck(argtable) != 0)
	{
		/* NULL entries were detected, some allocations must have failed */
		printf("%s: insufficient memory\n", progname);
		exitcode = 1;
		goto exit;
	}

	// default value
	sensorType->sval[0] = "HJ1";

	///* set any command line default values prior to parsing */
	//repeat->ival[0] = 3;
	//outfile->filename[0] = "-";

	/* Parse the command line as defined by argtable[] */
	nerrors = arg_parse(argc, argv, argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0)
	{
		printf("Usage: %s", progname);
		arg_print_syntax(stdout, argtable, "\n");
		printf("This program converts quicklook dat to tiff format with auxiliary data.\n");
		//printf("for parsing command line arguments. Argtable accepts integers\n");
		//printf("in decimal (123), hexadecimal (0xff), octal (0o123) and binary\n");
		//printf("(0b101101) formats. Suffixes KB, MB and GB are also accepted.\n");
		arg_print_glossary(stdout, argtable, "  %-25s %s\n");
		exitcode = 0;
		goto exit;
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0)
	{
		printf("'%s' Read Quicklook product.\n", progname);
		printf("December 2014, Long Tengfei.\n");
		exitcode = 0;
		goto exit;
	}

	/* If the parser returned any errors then display them and exit */
	if (nerrors > 0)
	{
		/* Display the error details contained in the arg_end struct.*/
		arg_print_errors(stdout, end, progname);
		printf("Try '%s -h or --help' for more information.\n", progname);
		exitcode = 1;
		goto exit;
	}

	/* special case: uname with no command line options induces brief help */
	if (argc == 1)
	{
		printf("Try '%s -h or --help' for more information.\n", progname);
		exitcode = 0;
		goto exit;
	}

	/* normal case: take the command line options at face value */
	exitcode = mymain(sensorType->sval, infile->filename, outpath->filename);

exit:
	/* deallocate each non-null entry in argtable[] */
	arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

	return exitcode;

	//if (argc > 1)
	//{
	//	/* -------------------------------------------------------------------- */
	//	/*      Parse arguments.                                                */
	//	/* -------------------------------------------------------------------- */
	//	for (int i = 1; i < argc; i++)
	//	{
	//		if (0 == _stricmp(argv[i], "-i"))
	//		{
	//			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
	//			pszInputFile = argv[++i];
	//		}
	//		else if (0 == _stricmp(argv[i], "-t"))
	//		{
	//			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
	//			sensorType = argv[++i];
	//			if (0 != _strcmpi(sensorType, strHJ1)
	//				&& 0 != _strcmpi(sensorType, strZY3))
	//			{
	//				printf("Unknown Sensor Type.\n");
	//				Usage();
	//			}
	//		}
	//		else if (0 == _stricmp(argv[i], "-op"))
	//		{
	//			CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(1);
	//			pszOutputDir = argv[++i];
	//		}
	//		else
	//		{
	//			Usage();
	//		}
	//	}

	//	if (0 == _strcmpi(pszInputFile, ""))
	//	{
	//		printf("inputFile can not be empty!\n");
	//		Usage();
	//	}
	//	else
	//	{
	//		char fullPath[_MAX_PATH];
	//		_fullpath(fullPath, pszInputFile, _MAX_PATH);
	//		if (0 == _strcmpi(pszOutputDir, ""))
	//		{
	//			pszOutputDir = SBeforeLast(string(fullPath), '\\').c_str();
	//		}
	//		if (!QDir(pszOutputDir).exists())
	//		{
	//			_mkdir(pszOutputDir);
	//		}

	//		clock_t  clockBegin, clockEnd;
	//		clockBegin = clock();
	//		if (0 == _strcmpi(sensorType, strHJ1))
	//		{
	//			// HJ1
	//			QVHJ::QVProcReader qv;
	//			qv.read_by_scene(fullPath, pszOutputDir);
	//		}
	//		else if (0 == _strcmpi(sensorType, strZY3))
	//		{
	//			// ZY3
	//			QVProc::ZY3QVProcReader zy3_qv;
	//			zy3_qv.read_by_scene(fullPath, pszOutputDir);
	//		}
	//		else
	//		{
	//			// Unknown Type
	//			printf("Unknown Sensor Type.\n");
	//			Usage();
	//		}

	//		clockEnd = clock();
	//		printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);
	//	}
	//}
	//else
	//{
	//	pszInputFile = "E:\\HJ1\\HJ-1A_CCD-1\\HJ-1A_CCD-1_SYC_201411250259_201411250305.dat";
	//	sensorType = strHJ1;
	//	char fullPath[_MAX_PATH];
	//	_fullpath(fullPath, pszInputFile, _MAX_PATH);
	//	if (0 == _strcmpi(pszOutputDir, ""))
	//	{
	//		pszOutputDir = SBeforeLast(string(fullPath), '\\').c_str();
	//	}
	//	if (!QDir(pszOutputDir).exists())
	//	{
	//		_mkdir(pszOutputDir);
	//	}

	//	clock_t  clockBegin, clockEnd;
	//	clockBegin = clock();
	//	if (0 == _strcmpi(sensorType, strHJ1))
	//	{
	//		// HJ1
	//		QVHJ::QVProcReader qv;
	//		qv.read_by_scene(fullPath, pszOutputDir);
	//	}
	//	else if (0 == _strcmpi(sensorType, strZY3))
	//	{
	//		// ZY3
	//		QVProc::ZY3QVProcReader zy3_qv;
	//		zy3_qv.read_by_scene(fullPath, pszOutputDir);
	//	}
	//	else
	//	{
	//		// Unknown Type
	//		printf("Unknown Sensor Type.\n");
	//		Usage();
	//	}

	//	clockEnd = clock();
	//	printf("Time consuming: %lf s\n", (clockEnd - clockBegin)*1e-3);

	//	Usage();
	//}
	//return 0;
}