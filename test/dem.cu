#include <iostream>
#include <fstream>
#include "dem.cuh"
#include <math.h>
#include <cstddef>
#include "structfile.cuh"
#include <sstream>
#include <string>
#include <iomanip>
#ifdef linux  
#include <unistd.h>  
#include <dirent.h>  
#endif  
#ifdef WIN32  
#include <direct.h>  //_mkdir fun
#include <io.h>  //_access fun
#endif 

#define NDEBUG
using namespace DEM;
using namespace STRUCTFILE;
using namespace std;
Dem::Dem(int rows, int cols, int mbc1, int mbc2, int mbc3, int mbc4)
{
#ifndef NDEBUG
	std::cerr << "DEM:begain parameters constructor" << std::endl;
	std::cerr << __func__ << std::endl;
#endif
	ynrows = rows;
	xncols = cols;
	xllcorner = 0;
	yllcorner = 0;
	Cellsize = 0;
	NODATA_VALUE = 0;
	mbc_1 = mbc1;
	mbc_2 = mbc2;
	mbc_3 = mbc3;
	mbc_4 = mbc4;
	DATA = new double *[cols];
	for (int i = 0; i < cols; i++)
	{
		DATA[i] = new double[rows];
	}
}

Dem::Dem(string file)
{
	#ifndef NDEBUG
		std::cerr << "Replace file:" << file << std::endl;
		std::cerr << "In file:" << file << std::endl;
		// std::cerr << __func__ << std::endl;
	#endif
		try
		{
			double title[6];
			char str[100];
			FILE *fin;
			char *fname = (char *)file.c_str();
			// int fileOK = 0;
			// double firstdata = 0.0;
	
			//打开dem读取文件头
			if ((fin= fopen(fname, "r"))== NULL)
			{
				// std::cout << fileOK << std::endl;
				std::cout << "Failed to open  file:" << fname << std::endl;
				// exit(1);
				return;
			}
			for (int i = 0; i < 6; i++)
			{
				fscanf(fin, "%s", &str);
				// printf("%s: ", &str);
				fscanf(fin, "%lf", &title[i]);
				// printf("%lf\n", title[i]);
			}
			set_xncols(int(title[0]));
			set_ynrows(int(title[1]));
			set_xllcorner(title[2]);
			set_yllcorner(title[3]);
			set_Cellsize(title[4]);
			set_NODATA_VALUE(title[5]);
			DATA = new double *[xncols];
			//DATA分配内存
			for (int i = 0; i < xncols; i++)
			{
				DATA[i] = new double[ynrows];
			}
			//读取网格，临时存在DATA
			for (int j = ynrows - 1; j >= 0; j--)
			{
				for (int i = 0; i <= xncols - 1; i++)
				{
					fscanf(fin, "%lf", &DATA[i][j]);
					if(DATA[i][j] == title[5])
					{
						DATA[i][j] = -9999.;
					}			
				}
			}
			// mbc_1 = 2;
			// mbc_2 = 2;
			// mbc_3 = 2;
			// mbc_4 = 2;
	#ifndef NDEBUG
			std::cerr << "Out test data:nx,ny,mbc" << std::endl;
			std::cerr << get_xncols() << std::endl;
			std::cerr << get_ynrows() << std::endl;
			std::cerr << get_xllcorner() << std::endl;
			std::cerr << get_yllcorner() << std::endl;
			std::cerr << get_Cellsize() << std::endl;
			std::cerr << get_NODATA_VALUE() << std::endl;
	#endif
			cout << "DEM constructor ok" << endl;
		}
		catch (...)
		{
			std::cout << "File error" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
Dem::Dem(FilePath &filepath)
{
#ifndef NDEBUG
	std::cerr << "File number construct DEM" << std::endl;
	std::cerr << __func__ << std::endl;
#endif
	try
	{
		double title[6];
		char str[100];
		FILE *demfile;
		char *fname = (char *)filepath.DEMPath.c_str();
		// int fileOK = 0;
		// double firstdata = 0.0;

		//打开dem读取文件头
		// fileOK = fopen_s(&demfile, fname, "r");
		if ((demfile= fopen(fname, "r"))== NULL)
		{
			// std::cout << fileOK << std::endl;
			std::cout << "Failed to open file:" << fname << std::endl;
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < 6; i++)
		{
			fscanf(demfile, "%s", &str);
			printf("%s: ", &str);
			fscanf(demfile, "%lf", &title[i]);
			printf("%lf\n", title[i]);
		}
		set_xncols(int(title[0]));
		set_ynrows(int(title[1]));
		set_xllcorner(title[2]);
		set_yllcorner(title[3]);
		set_Cellsize(title[4]);
		set_NODATA_VALUE(title[5]);
		DATA = new double *[xncols];
		//DATA分配内存
		for (int i = 0; i < xncols; i++)
		{
			DATA[i] = new double[ynrows];
		}
		//读取网格，临时存在DATA
		for (int j = ynrows - 1; j >= 0; j--)
		{
			for (int i = 0; i <= xncols - 1; i++)
			{
				fscanf(demfile, "%lf", &DATA[i][j]);
				if(DATA[i][j] == title[5])
				{
					DATA[i][j] = -9999.;
				}			
			}
		}
		// mbc_1 = 2;
		// mbc_2 = 2;
		// mbc_3 = 2;
		// mbc_4 = 2;
#ifndef NDEBUG
		std::cerr << "Out test data:nx,ny..." << std::endl;
		std::cerr << get_xncols() << std::endl;
		std::cerr << get_ynrows() << std::endl;
		std::cerr << get_xllcorner() << std::endl;
		std::cerr << get_yllcorner() << std::endl;
		std::cerr << get_Cellsize() << std::endl;
		std::cerr << get_NODATA_VALUE() << std::endl;
#endif
		cout << "DEM constructor ok" << endl;
	}
	catch (...)
	{
		std::cout << "File error" << std::endl;
		exit(EXIT_FAILURE);
	}
}

Dem::~Dem()
{
	// std::cout << "DEM in destructor" << std::endl;
	delete[] DATA;
}
void Dem::outDem()
{
#ifndef NDEBUG
	cout << "out DEM parameter free" << endl;
	cout << __func__ << endl;
#endif
}
void Dem::outDem(string filename, double **&x1, double cur_time)//二维输出
{
#ifndef NDEBUG
	cout << "out DEM with parameter" << endl;
	cout << __func__ << endl;
#endif
	stringstream ss;
	string cur_time_string;
	ss << cur_time;
	ss >> cur_time_string;
	// int fileOK;
	filename = filename + cur_time_string + ".txt";
	// filename = filename + to_string(cur_time) + ".txt";
	// if (!isFolderExist(filename))
	// {
	// 	// createDirectory(filename);
	// 	cout<<"no output file menu"<<endl;
	// 	exit(EXIT_FAILURE);
	// }
	const char *fname = filename.data();
	FILE *fp; //打开文件
	if (cur_time == 0.)
	{
		// fileOK = fopen_s(&fp, fname, "w");//windows系统
		fp = fopen(fname, "w");//linux系统
	}
	else if (cur_time > 0.)
	{
		// fileOK = fopen_s(&fp, fname, "w");//后续可选择直接写在文件后面
		fp = fopen(fname, "w");
	}
	if ((fp = fopen(fname, "w"))==NULL)
	{
		cout << "output " << fname << " failed" << endl;
		exit(EXIT_FAILURE);
	}
	fprintf(fp, "Time         ");
	fprintf(fp, "%-20.8lf\n", cur_time);
	fprintf(fp, "ncols        ");
	fprintf(fp, "%-8d\n", xncols);
	fprintf(fp, "nrows        ");
	fprintf(fp, "%-8d\n", ynrows);
	fprintf(fp, "xllcorner    ");
	fprintf(fp, "%-20.8lf\n", xllcorner);
	fprintf(fp, "yllcorner    ");
	fprintf(fp, "%-20.8lf\n", xllcorner);
	fprintf(fp, "cellsize     ");
	fprintf(fp, "%-20.5lf\n", Cellsize);
	fprintf(fp, "NODATA_value ");
	fprintf(fp, "%-20.5lf\n", NODATA_VALUE);
	for (int j = ynrows + mbc_3 - 1; j >= mbc_3; j--)
	{
		for (int i = mbc_1; i < xncols + mbc_1; i++)
		{
			if (x1[i][j] == NODATA_VALUE)
				fprintf(fp, "%22.10lf", NODATA_VALUE);
			else
				fprintf(fp, "%22.10lf", x1[i][j]);
		}
		fprintf(fp, " \n");
	}
	fclose(fp);
}
void Dem::outDem(string filename, double *&x1,double *&x2, double cur_time,int dem_type,int out_time_dem)//一维输出
{
#ifndef NDEBUG
	cout << "out DEM with parameter" << endl;
	cout << __func__ << endl;
	cout << filename << endl;
#endif
	int cur_time_int = floor(cur_time/out_time_dem) * floor(out_time_dem);

	if(dem_type == 0 )//系统平台输出
	{filename = filename + ".txt";}
	else if(dem_type == 1)//文件名输出带时间
	{filename = filename + std::to_string(cur_time_int) + ".txt";}
	else if(dem_type == 2)//文件名输出不带时间
	{filename = filename + ".txt";}
	
	ofstream file;
	file.open(filename,ios::out | ios::app);
	if(file.is_open())
	{
		if(dem_type == 0){file<<"Time         "<<fixed<<setw(-6)<<setprecision(1)<<cur_time/60.<<"\n";}
		file<<"ncols        "<<fixed<<setw(-8)<<xncols<<"\n";
		file<<"nrows        "<<fixed<<setw(-8)<<ynrows<<"\n";
		file<<"xllcorner    "<<fixed<<setw(-20)<<setprecision(8)<<xllcorner<<"\n";
		file<<"yllcorner    "<<fixed<<setw(-20)<<setprecision(8)<<yllcorner<<"\n";
		file<<"cellsize     "<<fixed<<setw(-20)<<setprecision(5)<<Cellsize<<"\n";
		if(dem_type < 3)
		{
			file<<"NODATA_value "<<fixed<<setw(-20)<<-9999.<<"\n";
			for (int j = ynrows + mbc_3 - 1;j >= mbc_3;  j--)
			{
				for (int i = mbc_1; i < xncols + mbc_1; i++)
				{
					if(x2[j + i * (ynrows + mbc_3 + mbc_4)] == -9999.)
					{
						file<<setw(8)<<setprecision(2)<<-9999.<<" ";
					}
					else
					{
						file<<setw(16)<<setprecision(10)<<x1[j + i * (ynrows + mbc_3 + mbc_4)]<<" ";
					}
				}
				file<<"\n";
			}
		}
		else
		{
			file<<"NODATA_value "<<fixed<<setw(-20)<<NODATA_VALUE<<"\n";
			for (int j = ynrows + mbc_3 - 1;j >= mbc_3;  j--)
			{
				for (int i = mbc_1; i < xncols + mbc_1; i++)
				{
					file<<setw(16)<<setprecision(10)<<x1[j + i * (ynrows + mbc_3 + mbc_4)]<<" ";
				}
				file<<"\n";
			}		
		}
		file.close();
	}
	else
	{
		cout << "output " << filename << " failed" << endl;
		exit(EXIT_FAILURE);
	}
}

void Dem::outDem(string filename, double *&x1,double *&x2, double cur_time,int dem_type)//一维输出
{
#ifndef NDEBUG
	cout << "out DEM with parameter" << endl;
	cout << __func__ << endl;
	cout << filename << endl;
#endif
	int cur_time_int = floor(cur_time/3600.) * 3600;
	
	// stringstream ss;
	// string cur_time_string;
	// ss << cur_time_int;
	// ss >> cur_time_string;

	if(dem_type == 0 )//系统平台输出
	{filename = filename + ".txt";}
	else if(dem_type == 1)//文件名输出带时间
	{filename = filename + std::to_string(cur_time_int) + ".txt";}
	else if(dem_type == 2)//文件名输出不带时间
	{filename = filename + ".txt";}
	
	ofstream file;
	file.open(filename,ios::out | ios::app);
	if(file.is_open())
	{
		if(dem_type == 0){file<<"Time         "<<fixed<<setw(-6)<<setprecision(1)<<cur_time/60.<<"\n";}
		file<<"ncols        "<<fixed<<setw(-8)<<xncols<<"\n";
		file<<"nrows        "<<fixed<<setw(-8)<<ynrows<<"\n";
		file<<"xllcorner    "<<fixed<<setw(-20)<<setprecision(8)<<xllcorner<<"\n";
		file<<"yllcorner    "<<fixed<<setw(-20)<<setprecision(8)<<yllcorner<<"\n";
		file<<"cellsize     "<<fixed<<setw(-20)<<setprecision(5)<<Cellsize<<"\n";
		if(dem_type < 3)
		{
			file<<"NODATA_value "<<fixed<<setw(-20)<<-9999.<<"\n";
			for (int j = ynrows + mbc_3 - 1;j >= mbc_3;  j--)
			{
				for (int i = mbc_1; i < xncols + mbc_1; i++)
				{
					if(x2[j + i * (ynrows + mbc_3 + mbc_4)] == -9999.)
					{
						file<<setw(8)<<setprecision(2)<<-9999.<<" ";
					}
					else
					{
						file<<setw(16)<<setprecision(10)<<x1[j + i * (ynrows + mbc_3 + mbc_4)]<<" ";
					}
				}
				file<<"\n";
			}
		}
		else
		{
			file<<"NODATA_value "<<fixed<<setw(-20)<<NODATA_VALUE<<"\n";
			for (int j = ynrows + mbc_3 - 1;j >= mbc_3;  j--)
			{
				for (int i = mbc_1; i < xncols + mbc_1; i++)
				{
					file<<setw(16)<<setprecision(10)<<x1[j + i * (ynrows + mbc_3 + mbc_4)]<<" ";
				}
				file<<"\n";
			}		}
		file.close();
	}
	else
	{
		cout << "output " << filename << " failed" << endl;
		exit(EXIT_FAILURE);
	}
}
// bool Dem::isFolderExist(string folder)
// {
// 	int ret = 0;

// 	ret = _access(folder.c_str(), 0);
// 	if (ret == 0)
// 		ret = true;
// 	else
// 		ret = false;

// 	return ret;
// }
// int Dem::createDirectory(string directoryPath)
// {
// 	uint32_t dirPathLen = 0;
// 	if (directoryPath != "") {
// 		dirPathLen = strlen(directoryPath.c_str());
// 	}
// 	if (dirPathLen > FILENAME_MAX)
// 	{
// 		return -1;
// 	}
// 	char tmpDirPath[FILENAME_MAX] = { 0 };
// 	for (uint32_t i = 0; i < dirPathLen; ++i)
// 	{
// 		tmpDirPath[i] = directoryPath[i];
// 		if (tmpDirPath[i] == '\\' || tmpDirPath[i] == '/')
// 		{
// 			if (!isFolderExist(tmpDirPath))
// 			{
// 				int ret = _mkdir(tmpDirPath);
// 				//BOOL ret = CreateDirectory(tmpDirPath, NULL);
// 				if (ret != 0)
// 					return -1;
// 			}
// 		}
// 	}
// 	return 0;
// }

void Dem::replace_DEM(double **data)
{
	for (int i = 0; i < xncols; i++)
	{
		for (int j = 0; j < ynrows; j++)
		{
			data[i][j] = 0.;
		}
	}
}
double *Dem::operator[](int i)
{
	return DATA[i];
}
bool Dem::equal(const Dem &obj)
{
	if (this->ynrows == obj.ynrows && this->xncols == obj.xncols && this->xllcorner == obj.xllcorner && this->yllcorner == obj.yllcorner)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Dem::readGrid(string filename, double **var)
{
	double **data_t;
	char *fname = (char *)filename.c_str();
	// char line[200];
	// int fileOK = 0;
	FILE *demfile;
	// fileOK = fopen_s(&demfile, fname, "r");
	if ((demfile = fopen(fname, "r"))==NULL)
	{
		std::cout << "Failed to open file:" << fname << std::endl;
	}

	for (int i = 0; i < 6; i++)
	{
		// getline(demfile,line);
	}
	for (int j = ynrows - 1; j >= 0; j--)
	{
		for (int k = xncols - 1; k >= 0; k--)
		{
			fscanf(demfile, "%lf", &data_t[k][j]);
			var[k][j] = data_t[k][j];
		}
	}
}
void Dem::set_mbc(int diff, int mbc_1, int mbc_2, int mbc_3, int mbc_4, int nx, int ny)
{
#ifndef NDEBUG
	std::cerr << "after set_mbc" << std::endl;
#endif
	int tBlock_x, tBlock_y;
	diff = 0;
	int m, n;
	if (diff == 0)
	{
		tBlock_x = 8;
		tBlock_y = 8;
	}
	else if (diff == 1)
	{
		tBlock_x = 8;
		tBlock_y = 16;
	}
	else if (diff == 2)
	{
		tBlock_x = 16;
		tBlock_y = 8;
	}
	else if (diff == 3)
	{
		tBlock_x = 16;
		tBlock_y = 16;
	}
	else if (diff == 4)
	{
		tBlock_x = 16;
		tBlock_y = 32;
	}
	else if (diff == 5)
	{
		tBlock_x = 16;
		tBlock_y = 32;
	}
	else if (diff == 6)
	{
		tBlock_x = 32;
		tBlock_y = 32;
	}
	else
	{
		tBlock_x = 8;
		tBlock_y = 8;
	}
	m = tBlock_x - nx % tBlock_x;
	n = tBlock_y - ny % tBlock_y;
	if (m < 4)
	{
		m = m + tBlock_x;
	}
	if (n < 4)
	{
		n = n + tBlock_y;
	}
	set_mbc_1(int(m / 2));
	set_mbc_2(m - get_mbc_1());
	set_mbc_3(int(n / 2));
	set_mbc_4(n - get_mbc_3());
#ifndef NDEBUG
	std::cerr << "IN DEM mbc1/2/3/4:" << get_mbc_1() << get_mbc_2() << get_mbc_3() << get_mbc_4() << endl;
#endif
}