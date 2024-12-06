#ifndef WIND_CUDA_H_
#define WIND_CUDA_H_
#include <stdarg.h>
#include <iostream>
#include "dem.cuh"
#include "structfile.cuh"
using namespace STRUCTFILE;
namespace WIND_CUDA
{
	class Wind_Cuda : public DEM::Dem
	{

	private:
		/* data */
		const double G = 9.8;
		const double EPS = 1.E-15;
		double area = 0.;
		// int simulated_time=6;
		int hours;
		int nx, ny; //加上边界虚拟节点后的总大小
		int channel_num = 0;
		int group_num = 24;
		int points_num = 1000;
		int *point_num;
		int *point_x;
		int *point_y;

		double **w;	
		double *w_data;
		double **dev_w;
		double *dev_w_data;

		double **sx;
		double *sx_data;
		double **dev_sx;
		double *dev_sx_data;

		double **sy;
		double *sy_data;
		double **dev_sy;
		double *dev_sy_data;

		double **C;
		double *C_data;
		double **dev_C;
		double *dev_C_data;

		double **qx;
		double *qx_data;
		double **dev_qx;
		double *dev_qx_data;


		double **qy;
		double *qy_data;
		double **dev_qy;
		double *dev_qy_data;

		double **R;
		double *R_data;
		double **dev_R;
		double *dev_R_data;

		double **Soil_depth;
		double *Soil_depth_data;
		double **dev_Soil_depth;
		double *dev_Soil_depth_data;

		double **INF;
		double *INF_data;
		double **dev_INF;
		double *dev_INF_data;

		double **INF_Total;
		double *INF_Total_data;
		double **dev_INF_Total;
		double *dev_INF_Total_data;


		double **INF_tp;
		double *INF_tp_data;
		double **dev_INF_tp;
		double *dev_INF_tp_data;

		double **INF_Ks;
		double *INF_Ks_data;
		double **dev_INF_Ks;
		double *dev_INF_Ks_data;

		double **INF_U;
		double *INF_U_data;
		double **dev_INF_U;
		double *dev_INF_U_data;

		double **INF_Os;
		double *INF_Os_data;
		double **dev_INF_Os;
		double *dev_INF_Os_data;

		double **INF_Oi;
		double *INF_Oi_data;
		double **dev_INF_Oi;
		double *dev_INF_Oi_data;

		double **domain;
		double *domain_data;
		double **dev_domain;
		double *dev_domain_data;

		double **boundry;
		double *boundry_data;
		double **dev_boundry;
		double *dev_boundry_data;

		double **manning;
		double *manning_data;
		double **dev_manning;
		double *dev_manning_data;

		double **tsxb;
		double *tsxb_data;
		double **dev_tsxb;
		double *dev_tsxb_data;

		double **tsyb;
		double *tsyb_data;
		double **dev_tsyb;
		double *dev_tsyb_data;

		double **tsxf;
		double *tsxf_data;
		double **dev_tsxf;
		double *dev_tsxf_data;

		double **tsyf;
		double *tsyf_data;
		double **dev_tsyf;
		double *dev_tsyf_data;

		double **tsxc;
		double *tsxc_data;
		double **dev_tsxc;
		double *dev_tsxc_data;

		double **tsyc;
		double *tsyc_data;
		double **dev_tsyc;
		double *dev_tsyc_data;

		double **max_u;
		double *max_u_data;
		double **dev_max_u;
		double *dev_max_u_data;

		double **max_v;
		double *max_v_data;
		double **dev_max_v;
		double *dev_max_v_data;

		double **sc;
		double *sc_data;
		double **dev_sc;
		double *dev_sc_data;

		double **ht;
		double *ht_data;
		double **dev_ht;
		double *dev_ht_data;

		double **hcorr;
		double *hcorr_data;
		double **dev_hcorr;
		double *dev_hcorr_data;

		double **hcorr1;
		double *hcorr1_data;
		double **dev_hcorr1;
		double *dev_hcorr1_data;

		double **hcorr2;
		double *hcorr2_data;
		double **dev_hcorr2;
		double *dev_hcorr2_data;

		double **hcorr3;
		double *hcorr3_data;
		double **dev_hcorr3;
		double *dev_hcorr3_data;
		//后三项取消了注释
		double **x;
		double *x_data;
		// double **dev_x;
		// double *dev_x_data;

		double **y;
		double *y_data;
		// double **dev_y;
		// double *dev_y_data;

		double **surface;
		double *surface_data;
		double **dev_surface;
		double *dev_surface_data;

		double **z;
		double *z_data;
		double **dev_z;
		double *dev_z_data;

		double **h;
		double *h_data;
		double **dev_h;
		double *dev_h_data;

		double **u;
		double *u_data;
		double **dev_u;
		double *dev_u_data;

		double **v;
		double *v_data;
		double **dev_v;
		double *dev_v_data;

		double **h_max;
		double *h_max_data;
		double **dev_h_max;
		double *dev_h_max_data;

		// double **hu_max;
		// double *hu_max_data;
		// double **dev_hu_max;
		// double *dev_hu_max_data;

		// double **hv_max;
		// double *hv_max_data;
		// double **dev_hv_max;
		// double *dev_hv_max_data;

		double **hvel_max;
		double *hvel_max_data;
		double **dev_hvel_max;
		double *dev_hvel_max_data;

		double **vel_max;
		double *vel_max_data;
		double **dev_vel_max;
		double *dev_vel_max_data;

		double **vel;
		double *vel_data;
		// double **dev_vel;
		// double *dev_vel_data;

		double **zs;
		double *zs_data;
		// double **dev_zs;
		// double *dev_zs_data;

		double **flow;
		double *flow_data;
		double **dev_flow;
		double *dev_flow_data;

		double **R_ALL;//所有降雨

		double *dev_maxh;//h、u、v最大值
		double *dev_maxu;
		double *dev_maxv;

		double **landuse;
		double *landuse_data;

		double **R_P;//累计降雨
		double *R_P_data;
		double **dev_R_P;
		double *dev_R_P_data;

		double **V_Cover;
		double *V_Cover_data;

		double **V_LAI;
		double *V_LAI_data;

		double **V_Smax;
		double *V_Smax_data;
		double **dev_V_Smax;
		double *dev_V_Smax_data;

		double **V_K;
		double *V_K_data;
		double **dev_V_K;
		double *dev_V_K_data;

		double **V_Ic;//累计植被截留
		double *V_Ic_data;
		double **dev_V_Ic;
		double *dev_V_Ic_data;

		double **R_Discount;//R折减
		double *R_Discount_data;
		double **dev_R_Discount;
		double *dev_R_Discount_data;

		double **R0;//R折减
		double *R0_data;
		double **dev_R0;
		double *dev_R0_data;

		double *distance;
		double **R_point_value;
		// double **R_ALL_data2;
		// double *R_ALL_data1;
		void writeDem(double **&x1);
		void replace_R0(string file = "R0.txt");
		void replace_INF_Ks(string file = "INF_Ks.txt");
		void replace_INF_U(string file = "INF_U.txt");
		void replace_INF_Os(string file = "INF_Os.txt");
		void replace_INF_Oi(string file = "INF_Oi.txt");
		void replace_Soil_depth(string file = "INF_Ks.txt");

		void replace_manning(string file = "manning.txt");
		void replace_landuse(string file = "landuse.txt");
	public:
		Wind_Cuda(int x1, int x2);
		Wind_Cuda(FilePath &filepath, int &rain_hours);
		~Wind_Cuda();
		//virtual void outDem(char *fileName,outModeWind from);

		virtual void outPlt();

		// void replace_INF_U(char* file="INF_U.txt");
		// void replace_INF_Os(char* file="INF_Os.txt");
		// void replace_INF_Oi(char* file="INF_Oi.txt");
		// void diff(char dir);

		/**
		 * 复制相关的指针数据
		 * tar->str
		*/
		template <class T>
		void swap_Copy_Data(T **&str, T **&tar);

		template <class T>
		void swap_Copy_Cuda_Data(T *&str, T **&tar);
		//内联函数
		double get_G() const { return G; };
		double get_EPS() const { return EPS; };
		int get_hours() { return hours; };
		int get_nx() { return nx; };
		int get_ny() { return ny; };
		int get_channel_num() { return channel_num; };

		int get_group_num() {return group_num ;};
		int *get_point_num() {return point_num; };
		int *get_point_x() {return point_x; };
		int *get_point_y() {return point_y; };
		double get_area(){return area;};

		double **&get_w() { return w; };
		double *&get_w_data() { return w_data; };
		double **&get_dev_w() { return dev_w; };
		double *&get_dev_w_data() { return dev_w_data; };

		double **&get_sx() { return sx; };
		double *&get_sx_data() { return sx_data; };
		double **&get_dev_sx() { return dev_sx; };
		double *&get_dev_sx_data() { return dev_sx_data; };

		double **&get_sy() { return sy; };
		double *&get_sy_data() { return sy_data; };
		double **&get_dev_sy() { return dev_sy; };
		double *&get_dev_sy_data() { return dev_sy_data; };

		double **&get_C() { return C; };
		double **&get_dev_C() { return dev_C; };

		double **&get_qx() { return qx; };
		double *&get_qx_data() { return qx_data; };
		double **&get_dev_qx() { return dev_qx; };
		double *&get_dev_qx_data() { return dev_qx_data; };

		double **&get_qy() { return qy; };
		double *&get_qy_data() { return qy_data; };
		double **&get_dev_qy() { return dev_qy; };
		double *&get_dev_qy_data() { return dev_qy_data; };

		double **&get_R() { return R; };
		double *&get_R_data() { return R_data; };
		double **&get_dev_R() { return dev_R; };
		double *&get_dev_R_data() { return dev_R_data; };

		double **&get_Soil_depth() { return Soil_depth; };
		double *&get_Soil_depth_data() { return Soil_depth_data; };
		double **&get_dev_Soil_depth() { return dev_Soil_depth; };
		double *&get_dev_Soil_depth_data() { return dev_Soil_depth_data; };

		double **&get_INF() { return INF; };
		double *&get_INF_data() { return INF_data; };
		double **&get_dev_INF() { return dev_INF; };
		double *&get_dev_INF_data() { return dev_INF_data; };

		double **&get_INF_Total() { return INF_Total; };
		double *&get_INF_Total_data() { return INF_Total_data; };
		double **&get_dev_INF_Total() { return dev_INF_Total; };
		double *&get_dev_INF_Total_data() { return dev_INF_Total_data; };

		double **&get_INF_tp() { return INF_tp; };
		double *&get_INF_tp_data() { return INF_tp_data; };
		double **&get_dev_INF_tp() { return dev_INF_tp; };
		double *&get_dev_INF_tp_data() { return dev_INF_tp_data; };

		double **&get_INF_Ks() { return INF_Ks; };
		double **&get_dev_INF_Ks() { return dev_INF_Ks; };

		double **&get_INF_U() { return INF_U; };
		double **&get_dev_INF_U() { return dev_INF_U; };

		double **&get_INF_Os() { return INF_Os; };
		double **&get_dev_INF_Os() { return dev_INF_Os; };

		double **&get_INF_Oi() { return INF_Oi; };
		double **&get_dev_INF_Oi() { return dev_INF_Oi; };

		double **&get_domain() { return domain; };
		double *&get_domain_data() { return domain_data; };
		double **&get_dev_domain() { return dev_domain; };
		double *&get_dev_domain_data() { return dev_domain_data; };

		double **&get_boundry() { return boundry; };
		double **&get_dev_boundry() { return dev_boundry; };

		double **&get_manning() { return manning; };
		double *&get_manning_data() { return manning_data; };
		double **&get_dev_manning() { return dev_manning; };
		double *&get_dev_manning_data() { return dev_manning_data; };

		double **&get_tsxb() { return tsxb; };
		double *&get_tsxb_data() { return tsxb_data; };
		double **&get_dev_tsxb() { return dev_tsxb; };
		double *&get_dev_tsxb_data() { return dev_tsxb_data; };

		double **&get_tsyb() { return tsyb; };
		double *&get_tsyb_data() { return tsyb_data; };
		double **&get_dev_tsyb() { return dev_tsyb; };
		double *&get_dev_tsyb_data() { return dev_tsyb_data; };

		double **&get_tsxf() { return tsxf; };
		double *&get_tsxf_data() { return tsxf_data; };
		double **&get_dev_tsxf() { return dev_tsxf; };
		double *&get_dev_tsxf_data() { return dev_tsxf_data; };

		double **&get_tsyf() { return tsyf; };
		double *&get_tsyf_data() { return tsyf_data; };
		double **&get_dev_tsyf() { return dev_tsyf; };
		double *&get_dev_tsyf_data() { return dev_tsyf_data; };

		double **&get_tsxc() { return tsxc; };
		double *&get_tsxc_data() { return tsxc_data; };
		double **&get_dev_tsxc() { return dev_tsxc; };
		double *&get_dev_tsxc_data() { return dev_tsxc_data; };

		double **&get_tsyc() { return tsyc; };
		double *&get_tsyc_data() { return tsyc_data; };
		double **&get_dev_tsyc() { return dev_tsyc; };
		double *&get_dev_tsyc_data() { return dev_tsyc_data; };

		double **&get_max_u() { return max_u; };
		double *&get_max_u_data() { return max_u_data; };
		double **&get_dev_max_u() { return dev_max_u; };
		double *&get_dev_max_u_data() { return dev_max_u_data; };

		double **&get_max_v() { return max_v; };
		double *&get_max_v_data() { return max_v_data; };
		double **&get_dev_max_v() { return dev_max_v; };
		double *&get_dev_max_v_data() { return dev_max_v_data; };

		double **&get_sc() { return sc; };
		double **&get_dev_sc() { return dev_sc; };

		double **&get_ht() { return ht; };
		double **&get_dev_ht() { return dev_ht; };

		double **&get_hcorr() { return hcorr; };
		double **&get_dev_hcorr() { return dev_hcorr; };
		
		double **&get_hcorr1() { return hcorr1; };
		double **&get_dev_hcorr1() { return dev_hcorr1; };

		double **&get_hcorr2() { return hcorr2; };
		double **&get_dev_hcorr2() { return dev_hcorr2; };

		double **&get_hcorr3() { return hcorr3; };
		double **&get_dev_hcorr3() { return dev_hcorr3; };

		double **&get_x() { return x; };
		double *&get_x_data() { return x_data; };

		double **&get_y() { return y; };
		double *&get_y_data() { return y_data; };

		double **&get_surface() { return surface; };
		double **&get_dev_surface() { return dev_surface; };

		double **&get_z() { return z; };
		double *&get_z_data() { return z_data; };
		double **&get_dev_z() { return dev_z; };
		double *&get_dev_z_data() { return dev_z_data; };

		double **&get_h() { return h; };
		double *&get_h_data() { return h_data; };
		double **&get_dev_h() { return dev_h; };
		double *&get_dev_h_data() { return dev_h_data; };

		double **&get_u() { return u; };
		double *&get_u_data() { return u_data; };
		double **&get_dev_u() { return dev_u; };
		double *&get_dev_u_data() { return dev_u_data; };

		double **&get_v() { return v; };
		double *&get_v_data() { return v_data; };
		double **&get_dev_v() { return dev_v; };
		double *&get_dev_v_data() { return dev_v_data; };

		double **&get_h_max() { return h_max; };
		double *&get_h_max_data() { return h_max_data; };
		double **&get_dev_h_max() { return dev_h_max; };
		double *&get_dev_h_max_data() { return dev_h_max_data; };

		// double **&get_hu_max() { return hu_max; };
		// double *&get_hu_max_data() { return hu_max_data; };
		// double **&get_dev_hu_max() { return dev_hu_max; };
		// double *&get_dev_hu_max_data() { return dev_hu_max_data; };

		// double **&get_hv_max() { return hv_max; };
		// double **&get_dev_hv_max() { return dev_hv_max; };

		double **&get_hvel_max() { return hvel_max; };
		double *&get_hvel_max_data() { return hvel_max_data; };
		double **&get_dev_hvel_max() { return dev_hvel_max; };
		double *&get_dev_hvel_max_data() { return dev_hvel_max_data; };

		double **&get_vel_max() { return vel_max; };
		double *&get_vel_max_data() { return vel_max_data; };
		double **&get_dev_vel_max() { return dev_vel_max; };
		double *&get_dev_vel_max_data() { return dev_vel_max_data; };

		double **&get_vel() { return vel; };
		double *&get_vel_data() { return vel_data; };

		double **&get_zs() { return zs; };
		double *&get_zs_data() { return zs_data; };

		double **&get_flow() { return flow; };
		double *&get_flow_data() { return flow_data; };
		double **&get_dev_flow() { return dev_flow; };
		double *&get_dev_flow_data() { return dev_flow_data; };

		double **&get_R_ALL() { return R_ALL; };

		double *&get_dev_maxh() { return dev_maxh; };
		double *&get_dev_maxu() { return dev_maxu; };
		double *&get_dev_maxv() { return dev_maxv; };

		double **&get_landuse() { return landuse; };
		double *&get_landuse_data() { return landuse_data; };

		double **&get_R_P() { return R_P; };
		double *&get_R_P_data() { return R_P_data; };
		double **&get_dev_R_P() { return dev_R_P; };
		double *&get_dev_R_P_data() { return dev_R_P_data; };

		double **&get_V_Smax() { return V_Smax; };
		double *&get_V_Smax_data() { return V_Smax_data; };
		double **&get_dev_V_Smax() { return dev_V_Smax; };
		double *&get_dev_V_Smax_data() { return dev_V_Smax_data; };

		double **&get_V_K() { return V_K; };
		double *&get_V_K_data() { return V_K_data; };
		double **&get_dev_V_K() { return dev_V_K; };
		double *&get_dev_V_K_data() { return dev_V_K_data; };

		double **&get_V_Ic() { return V_Ic; };
		double *&get_V_Ic_data() { return V_Ic_data; };
		double **&get_dev_V_Ic() { return dev_V_Ic; };
		double *&get_dev_V_Ic_data() { return dev_V_Ic_data; };

		double **&get_R_Discount() { return R_Discount; };
		double *&get_R_Discount_data() { return R_Discount_data; };
		double **&get_dev_R_Discount() { return dev_R_Discount; };
		double *&get_dev_R_Discount_data() { return dev_R_Discount_data; };

		double **&get_R0() { return R0; };
		double *&get_R0_data() { return R0_data; };
		double **&get_dev_R0() { return dev_R0; };
		double *&get_dev_R0_data() { return dev_R0_data; };

		void set_area(double k) {area = k; };
		void set_channel_num(int k) {channel_num = k; };
		void set_group_num(int k) {group_num = k; };
		void set_point_num(int i,int k) { point_num[i] = k; };
		void set_point_x(int i,int k) { point_x[i] = k; };
		void set_point_y(int i,int k) { point_y[i] = k; };
		void set_nxny();
		void set_R(int hour, double *&R_data, double **R_ALL);
		Wind_Cuda &operator=(WIND_CUDA::Wind_Cuda &obj) { return obj; };

		void set_surface(double **x1) { surface = x1; };
		void set_z(double **x1) { swap_Copy_Data(z, x1); };
		void set_h(double **x1) { swap_Copy_Data(h, x1); };
		void set_u(double **x1) { swap_Copy_Data(u, x1); };
		void set_v(double **x1) { swap_Copy_Data(v, x1); };
		//需要重写
		void set_h_max(double **&h_max, double **&h);
		void set_hu_max(double **&hu_max, double **&hv_max, double **&hvel_max,double **&h, double **&u, double **&v); //hu,hv,hvel
		void set_hours(int hour) { hours = hour; };
		void set_xy(double *&x_data, double *&y_data);
		void out_Plt(float cur_time, string filename, string n, int count, ...);
		void out_Point(double *flow_data, double cur_time, string filename,int points, int group);
		void set_R_ALL(int &hours, double area, int R_type);
		void replace_H(string file = "H_initial.txt");
		void replace_U(string file = "U_initial.txt");
		void replace_V(string file = "V_initial.txt");
	};
} // namespace WIND
#endif