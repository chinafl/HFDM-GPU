#include<iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include"wind_cuda.cuh"
#include"windmath_cuda.cuh"
#include"dem.cuh" 
#include <cuda_runtime.h>
#include "include/rapidjson/document.h"//读取截面点文件
#include "include/rapidjson/writer.h"
#include "include/rapidjson/stringbuffer.h"
#include "include/rapidjson/filereadstream.h"
#include "include/rapidjson/error/en.h"
#include <cstdio>
#include <vector>
#include <map>
#define CHECK(res) if(res!=cudaSuccess){exit(-1);}
#define NDEBUG
using namespace WINDMATH_CUDA;
using namespace WIND_CUDA;
using namespace DEM;
using namespace STRUCTFILE;
struct point
{
	int pointX;
	int pointY;
};
using namespace rapidjson;
// void diffusive_model();
clock_t FDM_start, FDM_end;//统计计算时间
void printDeviceProp(const cudaDeviceProp &prop);
int main()
{
	cudaError_t res;
	//计算时间及输出时间
	// double cur_time = 0.;
	// float total_time = 5.1;
	// float out_time_plt = 3600.;
	// float out_time_dem = 600.;
	// float out_time_point = 30.;
	// float temp_time_plt = out_time_plt;
	// float temp_time_dem = out_time_dem;
	// float temp_time_point = out_time_point;
	int json_file = 0;
	int rain_hours = 72;//降雨总时间
	double EPS = 1.e-15;
	FDM_start = clock();//计时开始
	

	FilePath filepath;
	//读取配置文件
	// setPath(filePath);
	//Wind实体化，整个计算区域
	Wind_Cuda wind(filepath, rain_hours);
	//Windmath计算模块1实体化，整个计算区域
	Windmath_Cuda diff_all(wind.get_nx(), wind.get_ny(), 
						   wind.get_mbc_1(), wind.get_mbc_2(), wind.get_mbc_3(), wind.get_mbc_4());
	//Windmath计算模块1实体化，分区
	// Windmath diff_regions(2,2,2,2);

	FILE* fp ;
	if ((fp = fopen("./input/point.json", "rb"))==NULL)
	{
		json_file = 1;
		cout<<"without point.json,no output point data"<<endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		char readBuffer[65536];
		FileReadStream is(fp, readBuffer, sizeof(readBuffer));

		Document d;
		d.ParseStream(is);
		if (d.HasParseError())
		{
			fprintf(stderr, "\nError(offset %u): %s\n",
				(unsigned)d.GetErrorOffset(),
				d.GetParseError());
		}
		std::map<int, int> mapgroup;
		int i = 1;
		std::string g = "group1";
		char *group = new char[6];// "group1";
		strcpy(group, g.c_str());
		//std::cout<<d.HasParseError()<<std::endl;
		while (d.HasMember(group))
		{
			std::cout << group << std::endl;
			//group = cc.c_str();
			//const char* group="group"+i.toString();
			//std::cout<<d.HasMember(group)<<std::endl;
			// if (d.HasMember(group))
			// {
			Value &jvobject = d[group];
			int j = 1;
			char *point = new char[6];// "point1";
			std::string p = "point1";
			strcpy(point, p.c_str());
			while (jvobject.HasMember(point))
			{

				//group=cc.c_str();
				//point = dd.c_str();
				//const char* point="point"+j.toString();
				//std::cout<<jvobject.HasMember(point)<<std::endl;
				// if (jvobject.HasMember(point))
				// {
				Value &jvobject2 = jvobject[point];
				if (jvobject2.HasMember("X"))
				{
					Value &jvobject3 = jvobject2["X"];
					Value &jvobject4 = jvobject2["Y"];
					std::cout << group << "  " << point << std::endl;
					std::cout << "x:" << jvobject3.GetInt();
					std::cout << "Y:" << jvobject4.GetInt() << std::endl;
					wind.get_point_x()[j-1 + (i-1)*1000] = jvobject3.GetInt();
					wind.get_point_y()[j-1 + (i-1)*1000] = jvobject4.GetInt();
				}
				j++;
				std::string const &dd = std::string("point") + std::to_string(j);
				int ddlen = dd.length();
				strcpy(point, dd.c_str());
				//point = (char *)malloc((ddlen + 1) * sizeof(char));
				//dd.copy(point, ddlen, 0);
				// }
			}
			//mapgroup[i] = j;
			//mapStudent.insert(pair<int, string>(1, "student_one"));

			//}
			i++;
			//std::cout << "openDOM" << std::endl;
			std::string const &cc = std::string("group") + std::to_string(i);
			int cclen = cc.length();
			strcpy(group, cc.c_str());
			//group = (char *)malloc((cclen + 1) * sizeof(char));
			//cc.copy(group, cclen, 0);
			mapgroup.insert(std::pair<int, int>(i - 1, j - 1));
		}

		std::map<int, int>::iterator iter;
		iter = mapgroup.begin();
		while (iter != mapgroup.end())
		{
			std::cout << iter->first << "-" << iter->second << std::endl;
			wind.set_point_num(iter->first-1 ,iter->second );
			iter++;
		}
	}
	//liWB AV nux
	wind.set_xy(wind.get_x_data(), wind.get_y_data());//计算x,y坐标并存储
	//计算预处理（初始化边界、通量）
	// diff_all.pre_cul(wind, rain_hours);//开边界和没有初始流速流深时不启用
	int nx=wind.get_nx();
	int ny=wind.get_ny();
	diff_all.set_wind(wind);
	wind.set_R_ALL(rain_hours,wind.get_area(),3);//
	cout<<"in main rain_hours:"<<rain_hours<<endl;
	cout<<"channel_num:"<<wind.get_channel_num()<<endl;

	int istat;
	int nDevices;

	istat = cudaGetDeviceCount(&nDevices);
	std::cout<<"devices: "<< nDevices <<endl;

	for (int i = 0; i < nDevices; ++i)
	{
		cudaDeviceProp prop;
		if (cudaGetDeviceProperties(&prop, i) == cudaSuccess) 
		{
			if (prop.major >= 1) 
			{
				printDeviceProp(prop);
				break;
			}
		}
	}
	//exit(EXIT_FAILURE);
	//填洼处理
	// diff_all.wind_fill();

	//处理基础流量
	diff_all.wind_baseflow(rain_hours,json_file);

	// FILE *discharge_file;
	// if ((discharge_file= fopen("./output/Q1.txt", "r"))== NULL)
    // {
    //     std::cout << "Failed to open file:" << "Q1" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
	// else
	// {
	// 	string H_file;
	// 	string U_file;
	// 	string V_file;
		
	// 	double discharge = 0. ;
	// 	double discharge_time = 0. ;
	// 	double max_time = 0.;
	// 	int max_time_int = 0;
	// 	double max_discharge = 0.;
	// 	//查找最大流量及其出现时间
	// 	while(!feof(discharge_file))
	// 	{
	// 		fscanf(discharge_file, "%lf", &discharge_time);
	// 		printf("%lf: ", discharge_time);
	// 		fscanf(discharge_file, "%lf", &discharge);
	// 		printf("%lf\n", discharge);
	// 		if(max_discharge < discharge)
	// 		{
	// 			max_discharge = discharge;
	// 			max_time = discharge_time;
	// 		}
	// 		// if(discharge_time < 1.){discharge = -10.;}
	// 	}
	// 	std::cout<<"max_time: "<<max_time<<endl;
	// 	std::cout<<"max_discharge: "<<max_discharge<<endl;
	// 	//根据最大流量出现时间点，找到指定文件，读取内容并替换h/u/v
	// 	max_time_int = floor(max_time / 3600) * 3600 ;
	// 	H_file = "./output/H/H" + std::to_string(max_time_int) + std::string(".txt");
	// 	U_file = "./output/U/U" + std::to_string(max_time_int) + std::string(".txt");
	// 	V_file = "./output/V/V" + std::to_string(max_time_int) + std::string(".txt");
	// 	wind.replace_H(H_file);
	// 	wind.replace_U(U_file);
	// 	wind.replace_V(V_file);
	// 	//将最大流量时h/U/V输出到当前工作目录
	// 	wind.outDem("H_C",wind.get_h_data(),wind.get_z_data(),cur_time,2);
	// 	wind.outDem("U_C",wind.get_u_data(),wind.get_z_data(),cur_time,2);
	// 	wind.outDem("V_C",wind.get_v_data(),wind.get_z_data(),cur_time,2);
	// }

	FDM_end = clock();
	printf("%lf seconds\n", (double)(FDM_end - FDM_start) / CLOCKS_PER_SEC);
	std::cout<<"Compution is done!"<<endl;
	// system("pause");//捕获停止
	return 0;
}
void printDeviceProp(const cudaDeviceProp &prop)
{
	printf("Device Name : %s.\n", prop.name);
	printf("totalGlobalMem : %d.\n", prop.totalGlobalMem);
	printf("sharedMemPerBlock : %d.\n", prop.sharedMemPerBlock);
	printf("regsPerBlock : %d.\n", prop.regsPerBlock);
	printf("warpSize : %d.\n", prop.warpSize);
	printf("memPitch : %d.\n", prop.memPitch);
	printf("maxThreadsPerBlock : %d.\n", prop.maxThreadsPerBlock);    printf("maxThreadsDim[0 - 2] : %d %d %d.\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
	printf("maxGridSize[0 - 2] : %d %d %d.\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
	printf("totalConstMem : %d.\n", prop.totalConstMem);
	printf("major.minor : %d.%d.\n", prop.major, prop.minor);
	printf("clockRate : %d.\n", prop.clockRate);
	printf("textureAlignment : %d.\n", prop.textureAlignment);
	printf("deviceOverlap : %d.\n", prop.deviceOverlap);
	printf("multiProcessorCount : %d.\n", prop.multiProcessorCount);
}
