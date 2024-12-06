#ifndef WINDMATH_CUDA_H_
#define WINDMATH_CUDA_H_
#include <iostream>
#include"wind_cuda.cuh"
#include <cuda_runtime.h>
using namespace std;
using namespace WIND_CUDA;
namespace WINDMATH_CUDA
{
	class Windmath_Cuda
	{
	private:
		int nx, ny;
		int mbc_1, mbc_2, mbc_3, mbc_4;
		// double cur_time;
		// double total_time;
		// double cur_time;
		// dim3 dimBlock;
		// dim3 dimGrid;
		Wind_Cuda *wind;
	public:
		
		Windmath_Cuda();
		Windmath_Cuda(int nx, int ny, int mbc_1, int mbc_2, int mbc_3, int mbc_4);
		~Windmath_Cuda();

		int get_nx()const { return nx; };
		int get_ny()const { return ny; };
		int get_mbc_1()const { return mbc_1; };
		int get_mbc_2()const { return mbc_2; };
		int get_mbc_3()const { return mbc_3; };
		int get_mbc_4()const { return mbc_4; };
		void set_mbc(int x1, int x2, int x3, int x4) { mbc_1 = x1; mbc_2 = x2; mbc_3 = x3, mbc_4 = x4; };
		void set_nxny(int x1, int x2) { nx = x1; ny = x2; };
		void set_wind(Wind_Cuda &w) { wind = &w; };
		void set_block(int nx, int ny) {
			// dimBlock.x = 8;
			// dimBlock.y = 8;
			//dimBlock(8, 8); 
			// dimGrid.x = (nx + dimBlock.x - 1) / (dimBlock.x);
			// dimGrid.y = (ny + dimBlock.y - 1) / (dimBlock.y);
			//dimGrid((nx + dimBlock.x - 1) / (dimBlock.x), (ny + dimBlock.y - 1) / (dimBlock.y));	 
		};
		// void wind_solver();
		void wind_start(double &cur_time, double &dt);
		void wind_fill();//填洼工具
		void wind_baseflow(int rain_hours, int json_file);//沟道基流
		void diff_xy(double &dt);
		void diff_yx(double &dt);
		void diff_time_step(double &cur_time, double &dt);//,double **max_u, double **max_v);
		void diff_slope();
		//private://私有化
		void bound(double *&domain_data, double *&boundry_data, double *&z_data, int nx, int ny);
		void qxy(double *&h_data, double *&u_data, double *&v_data, double *&qx_data, double *&qy_data);
		void pre_cul(Wind_Cuda &w, int hours);
		void get_max();
		void out_cul(double *&zs_data,double *&h_data,double *&z_data,double *&vel_data,double *&u_data,double *&v_data,int xncols,int ynrows,int mbc_1,int mbc_2,int mbc_3,int mbc_4);
		// void output(Wind_Cuda &w, FilePath &filePath, double &cur_time, float &temp_time_plt, float out_time_plt, float &temp_time_dem, float out_time_dem);
	};
	__global__ void Diff_w(double **w, double **z, double **sc, double **h, double **hcorr1, double **hcorr2, double **hcorr3, double **ht, double **u, double **v, int nx,  int ny);
	__global__ void Diff_u(double **h, double **u, double **qx, double **sx, double **sy, double **manning, double EPS, int nx,  int ny);
	__global__ void Diff_v(double **h, double **v, double **qy, double **sx, double **sy, double **manning, double EPS, int nx,  int ny);
	__global__ void Diff_qx(double **qx, double **h, double **u, int nx,  int ny);
	__global__ void Diff_qy(double **qy, double **h, double **v, int nx,  int ny);
	__global__ void Diff_infiltration(double **R0,double **R_Discount, double **INF, double **INF_Total, double **INF_tp,
		double **INF_Ks, double **INF_U, double **INF_Os, double **INF_Oi,double **Soil_depth, double dt, int nx,  int ny);
	__global__ void Diff_mathtestx(double dt, double Cellsize, double EPS, double **domain, double **sx, double **h, double **ht, double **qx,
		double **R, double **INF, double **tsxf, double **tsxb, double **hcorr1, double **hcorr2, double **hcorr3,double **R_Discount, int nx,  int ny,
		int mb_1,int mbc_2);
	__global__ void Diff_mathtesty(double dt, double Cellsize, double EPS, double **domain, double **sy, double **h, double **ht, double **qy,
		double **R, double **INF, double **tsyf, double **tsyb, double **hcorr1, double **hcorr2, double **hcorr3,double **R_Discount, int nx,  int ny,
		int mbc_3,int mbc_4);
	__global__ void Diff_x_hcorr(double **hcorr, double **hcorr1, double **hcorr2, double **hcorr3, int nx,  int ny,int mbc_1,int mbc_2);
	__global__ void Diff_y_hcorr(double **hcorr, double **hcorr1, double **hcorr2, double **hcorr3, int nx,  int ny,int mbc_3,int mbc_4);
	__global__ void Diff_hx(double **h, double **ht, double **hcorr, int nx,  int ny,int mbc_1,int mbc_2);
	__global__ void Diff_hy(double **h, double **ht, double **hcorr, int nx,  int ny,int mbc_3,int mbc_4);
	__global__ void Diff_flow(double **flow, double **h, double **u, double **v,double **h_max,double **hvel_max, double **vel_max,double Cellsize, int nx,  int ny);
	__global__ void Slope_initial(double **tsxb, double **tsyb, double **tsxf, double **tsyf, double **tsxc, double **tsyc, double **sx, double **sy, int nx,  int ny);
	__global__ void Slope_ts(double Cellsize, double **tsxb, double **tsyb, double **tsxf, double **tsyf, double **tsxc, double **tsyc, double **w, double **sx, double **sy,
				int nx,  int ny);
	__global__ void Slope_where(double **tsxb, double **tsyb, double **tsxf, double **tsyf, double **tsxc, double **tsyc, double **sx, double **sy, int nx,  int ny);
	__global__ void Time_uv(double EPS, double **h, double **u, double **v, double **sx, double **sy, double **manning, int nx,  int ny, int mbc_1, int mbc_2, int mbc_3, int mbc_4,double *dev_maxh);
	__global__ void Time_diff_uv(double **u, double **v, double **max_u, double **max_v, int nx,  int ny,int mb_1,int mbc_2,int mbc_3,int mbc_4);
	__global__ void Time_dt(double **max_u, double **max_v, double **h, double **u, double **v,
				int nx,  int ny,int mbc_1,int mbc_2,int mbc_3,int mbc_4,double *dev_maxu,double *dev_maxv);
	__global__ void Diff_interception(double **R,double **R_P,double **V_Ic,double **V_Smax,double **V_K, double **R_Discount,int nx,int ny,double Cellsize,double dt);
	__global__ void Fill(double **z, double **h, int nx, int ny);
	__device__ double atomicMax(double *address, double value);

}
#endif