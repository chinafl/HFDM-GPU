#ifndef DEM_H_
#define DEM_H_
#include <iostream>
#include"structfile.cuh"
using namespace STRUCTFILE;
namespace DEM
{
	class Dem
	{
	private:
		int mbc_1, mbc_2, mbc_3, mbc_4;
		int ynrows, xncols; //行列号
		double xllcorner, yllcorner;      //左下角点的坐标值
		double Cellsize;      //栅格单元大小
		double NODATA_VALUE;
		double **DATA;         //矩阵

	public:
		/**
		 * 传入文件名。通过读取DEM文件进行调用
		*/
		Dem(int x1, int x2, int = 2, int = 2, int = 2, int = 2);
		Dem(string file);
		Dem(FilePath &filepath);
		bool isFolderExist(string folder);
		int createDirectory(string directoryPath);
		virtual ~Dem();
		void outDem();//输出方法
		void outDem(string filename, double **&x1, double cur_time);
		void outDem(string filename, double *&x1,double *&x2, double cur_time, int dem_type, int out_time_dem);//一维输出
		void outDem(string filename, double *&x1,double *&x2, double cur_time,int dem_type);
		void readDem(string file);//更新方法
		void readGrid(string file, double **var);//读取网格
		void replace_DEM(double **data);
		/**
		 * 通过[]获取该DEM数据的DATA值
		*/
		double * operator[](int i);
		/**
		 * 函数比对-
		 * 只比对两个DEM的网格大小和左下角坐标
		 * 该注释勿删
		*/
		bool equal(const Dem &DEM);
		//以下是内联函数。不需要重写
		int get_mbc_1() { return mbc_1; };
		int get_mbc_2() { return mbc_2; };
		int get_mbc_3() { return mbc_3; };
		int get_mbc_4() { return mbc_4; };
		int get_xncols() { return xncols; };
		int get_ynrows() { return ynrows; };
		double get_xllcorner()const { return xllcorner; };
		double get_yllcorner()const { return yllcorner; };
		double get_Cellsize()const { return Cellsize; };
		double get_NODATA_VALUE()const { return NODATA_VALUE; };
		void set_mbc_1(int x1) { mbc_1 = x1; };
		void set_mbc_2(int x1) { mbc_2 = x1; };
		void set_mbc_3(int x1) { mbc_3 = x1; };
		void set_mbc_4(int x1) { mbc_4 = x1; };
		void set_xncols(int x1) { xncols = x1; };
		void set_ynrows(int x1) { ynrows = x1; };
		void set_xllcorner(double x1) { xllcorner = x1; };
		void set_yllcorner(double x1) { yllcorner = x1; };
		void set_Cellsize(double x1) { Cellsize = x1; };
		void set_NODATA_VALUE(double x1) { NODATA_VALUE = x1; };
		void set_DATA(double **data) { DATA = data; };//单独更新DATA数据
		void set_mbc(int diff, int mbc_1, int mbc_2, int mbc_3, int mbc_4, int nx, int ny);
		double **& get_Data() { return DATA; };
		double** get_Data_COPY() { return DATA; };
	};
}
#endif