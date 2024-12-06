#include <iostream>
#include <fstream>
#include <math.h>
#include <stdarg.h>
#include "wind_cuda.cuh"
#include "structfile.cuh"
#include <cuda_runtime.h>
// #include <cuda_runtime.h>
#ifdef linux  
#include <unistd.h>  
#include <dirent.h>  
#endif  
#ifdef WIN32  
#include <direct.h>  //_mkdir fun
#include <io.h>  //_access fun
#endif 
#define NDEBUG
#define CHECK(res) if(res!=cudaSuccess){exit(-1);}
using namespace WIND_CUDA;
using namespace std;
using namespace DEM;
using namespace STRUCTFILE;
Wind_Cuda::Wind_Cuda(int rows, int cols) : Dem(rows, cols)//type1:
{
#ifndef NDEBUG
	std::cerr << "Wind(int rows, int cols)" << std::endl;
#endif
	set_mbc(0, get_mbc_1(), get_mbc_2(), get_mbc_3(), get_mbc_4(), get_xncols(), get_ynrows());//重新计算mbc_1/2/3/4，并在dem中调用修改
	int nx = get_xncols() + get_mbc_1() + get_mbc_2();
	int ny = get_ynrows() + get_mbc_3() + get_mbc_4();
	cudaError_t res;

	w = (double**)malloc(nx * sizeof(double*));
	w_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_w), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_w_data), nx*ny * sizeof(double)); CHECK(res);

	sx = (double**)malloc(nx * sizeof(double*));
	sx_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_sx), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_sx_data), nx*ny * sizeof(double)); CHECK(res);

	sy = (double**)malloc(nx * sizeof(double*));
	sy_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_sy), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_sy_data), nx*ny * sizeof(double)); CHECK(res);

	C = (double**)malloc(nx * sizeof(double*));
	C_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_C), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_C_data), nx*ny * sizeof(double)); CHECK(res);

	qx = (double**)malloc(nx * sizeof(double*));
	qx_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_qx), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_qx_data), nx*ny * sizeof(double)); CHECK(res);

	qy = (double**)malloc(nx * sizeof(double*));
	qy_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_qy), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_qy_data), nx*ny * sizeof(double)); CHECK(res);

	R = (double**)malloc(nx * sizeof(double*));
	R_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_R), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_R_data), nx*ny * sizeof(double)); CHECK(res);

	Soil_depth = (double**)malloc(nx * sizeof(double*));
	Soil_depth_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_Soil_depth), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_Soil_depth_data), nx*ny * sizeof(double)); CHECK(res);

	INF = (double**)malloc(nx * sizeof(double*));
	INF_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_data), nx*ny * sizeof(double)); CHECK(res);

	INF_Total = (double**)malloc(nx * sizeof(double*));
	INF_Total_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_Total), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_Total_data), nx*ny * sizeof(double)); CHECK(res);

	INF_tp = (double**)malloc(nx * sizeof(double*));
	INF_tp_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_tp), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_tp_data), nx*ny * sizeof(double)); CHECK(res);

	INF_Ks = (double**)malloc(nx * sizeof(double*));
	INF_Ks_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_Ks), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_Ks_data), nx*ny * sizeof(double)); CHECK(res);

	INF_U = (double**)malloc(nx * sizeof(double*));
	INF_U_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_U), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_U_data), nx*ny * sizeof(double)); CHECK(res);

	INF_Os = (double**)malloc(nx * sizeof(double*));
	INF_Os_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_Os), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_Os_data), nx*ny * sizeof(double)); CHECK(res);

	INF_Oi = (double**)malloc(nx * sizeof(double*));
	INF_Oi_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_Oi), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_Oi_data), nx*ny * sizeof(double)); CHECK(res);

	domain = (double**)malloc(nx * sizeof(double*));
	domain_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_domain), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_domain_data), nx*ny * sizeof(double)); CHECK(res);

	boundry = (double**)malloc(nx * sizeof(double*));
	boundry_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_boundry), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_boundry_data), nx*ny * sizeof(double)); CHECK(res);

	manning = (double**)malloc(nx * sizeof(double*));
	manning_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_manning), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_manning_data), nx*ny * sizeof(double)); CHECK(res);

	tsxb = (double**)malloc(nx * sizeof(double*));
	tsxb_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsxb), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsxb_data), nx*ny * sizeof(double)); CHECK(res);

	tsyb = (double**)malloc(nx * sizeof(double*));
	tsyb_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsyb), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsyb_data), nx*ny * sizeof(double)); CHECK(res);

	tsxf = (double**)malloc(nx * sizeof(double*));
	tsxf_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsxf), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsxf_data), nx*ny * sizeof(double)); CHECK(res);

	tsyf = (double**)malloc(nx * sizeof(double*));
	tsyf_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsyf), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsyf_data), nx*ny * sizeof(double)); CHECK(res);

	tsxc = (double**)malloc(nx * sizeof(double*));
	tsxc_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsxc), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsxc_data), nx*ny * sizeof(double)); CHECK(res);

	tsyc = (double**)malloc(nx * sizeof(double*));
	tsyc_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsyc), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsyc_data), nx*ny * sizeof(double)); CHECK(res);

	max_u = (double**)malloc(nx * sizeof(double*));
	max_u_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_max_u), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_max_u_data), nx*ny * sizeof(double)); CHECK(res);

	max_v = (double**)malloc(nx * sizeof(double*));
	max_v_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_max_v), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_max_v_data), nx*ny * sizeof(double)); CHECK(res);

	sc = (double**)malloc(nx * sizeof(double*));
	sc_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_sc), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_sc_data), nx*ny * sizeof(double)); CHECK(res);

	ht = (double**)malloc(nx * sizeof(double*));
	ht_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_ht), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_ht_data), nx*ny * sizeof(double)); CHECK(res);

	hcorr1 = (double**)malloc(nx * sizeof(double*));
	hcorr1_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hcorr1), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hcorr1_data), nx*ny * sizeof(double)); CHECK(res);

	hcorr2 = (double**)malloc(nx * sizeof(double*));
	hcorr2_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hcorr2), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hcorr2_data), nx*ny * sizeof(double)); CHECK(res);

	hcorr3 = (double**)malloc(nx * sizeof(double*));
	hcorr3_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hcorr3), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hcorr3_data), nx*ny * sizeof(double)); CHECK(res);

	x = (double**)malloc(nx * sizeof(double*));
	x_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_x), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_x_data), nx*ny * sizeof(double)); CHECK(res);

	y = (double**)malloc(nx * sizeof(double*));
	y_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_y), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_y_data), nx*ny * sizeof(double)); CHECK(res);

	surface = (double**)malloc(nx * sizeof(double*));
	surface_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_surface), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_surface_data), nx*ny * sizeof(double)); CHECK(res);

	z = (double**)malloc(nx * sizeof(double*));
	z_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_z), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_z_data), nx*ny * sizeof(double)); CHECK(res);

	h = (double**)malloc(nx * sizeof(double*));
	h_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_h), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_h_data), nx*ny * sizeof(double)); CHECK(res);

	u = (double**)malloc(nx * sizeof(double*));
	u_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_u), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_u_data), nx*ny * sizeof(double)); CHECK(res);

	v = (double**)malloc(nx * sizeof(double*));
	v_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_v), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_v_data), nx*ny * sizeof(double)); CHECK(res);

	h_max = (double**)malloc(nx * sizeof(double*));
	h_max_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_h_max), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_h_max_data), nx*ny * sizeof(double)); CHECK(res);

	// hu_max = (double**)malloc(nx * sizeof(double*));
	// hu_max_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_hu_max), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_hu_max_data), nx*ny * sizeof(double)); CHECK(res);

	// hv_max = (double**)malloc(nx * sizeof(double*));
	// hv_max_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_hv_max), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_hv_max_data), nx*ny * sizeof(double)); CHECK(res);

	hvel_max = (double**)malloc(nx * sizeof(double*));
	hvel_max_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hvel_max), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hvel_max_data), nx*ny * sizeof(double)); CHECK(res);

	vel_max = (double**)malloc(nx * sizeof(double*));
	vel_max_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_vel_max), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_vel_max_data), nx*ny * sizeof(double)); CHECK(res);

	zs = (double**)malloc(nx * sizeof(double*));
	zs_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_zs), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_zs_data), nx*ny * sizeof(double)); CHECK(res);

	for (int i = 0; i < nx; i++)
	{
		w[i] = dev_w_data + i * ny;
		sx[i] = dev_sx_data + i * ny;
		sy[i] = dev_sy_data + i * ny;
		C[i] = dev_C_data + i * ny;
		qx[i] = dev_qx_data + i * ny;
		qy[i] = dev_qy_data + i * ny;
		R[i] = dev_R_data + i * ny;
		Soil_depth[i] = dev_Soil_depth_data + i * ny;
		INF[i] = dev_INF_data + i * ny;

		INF_Total[i] = dev_INF_Total_data + i * ny;
		INF_tp[i] = dev_INF_tp_data + i * ny;
		INF_Ks[i] = dev_INF_Ks_data + i * ny;
		INF_U[i] = dev_INF_U_data + i * ny;
		INF_Os[i] = dev_INF_Os_data + i * ny;
		INF_Oi[i] = dev_INF_Oi_data + i * ny;
		domain[i] = dev_domain_data + i * ny;
		boundry[i] = dev_boundry_data + i * ny;
		manning[i] = dev_manning_data + i * ny;
		tsxb[i] = dev_tsxb_data + i * ny;
		tsyb[i] = dev_tsyb_data + i * ny;
		tsxf[i] = dev_tsxf_data + i * ny;
		tsyf[i] = dev_tsyf_data + i * ny;
		tsxc[i] = dev_tsxc_data + i * ny;
		tsyc[i] = dev_tsyc_data + i * ny;
		max_u[i] = dev_max_u_data + i * ny;
		max_v[i] = dev_max_v_data + i * ny;
		sc[i] = dev_sc_data + i * ny;
		ht[i] = dev_ht_data + i * ny;
		hcorr1[i] = dev_hcorr1_data + i * ny;
		hcorr2[i] = dev_hcorr2_data + i * ny;
		hcorr3[i] = dev_hcorr3_data + i * ny;

		// x[i] = dev_x_data + i * ny;
		// y[i] = dev_y_data + i * ny;
		surface[i] = dev_surface_data + i * ny;
		z[i] = dev_z_data + i * ny;
		h[i] = dev_h_data + i * ny;
		u[i] = dev_u_data + i * ny;
		v[i] = dev_v_data + i * ny;
		h_max[i] = dev_h_max_data + i * ny;
	}
	res = cudaMemcpy((void*)(dev_w), (void*)(w), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sx), (void*)(sx), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sy), (void*)(sy), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_C), (void*)(C), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_qx), (void*)(qx), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_qy), (void*)(qy), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R), (void*)(R), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_Soil_depth), (void*)(Soil_depth), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF), (void*)(INF), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Total), (void*)(INF_Total), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_tp), (void*)(INF_tp), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Ks), (void*)(INF_Ks), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_U), (void*)(INF_U), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Os), (void*)(INF_Os), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Oi), (void*)(INF_Oi), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_domain), (void*)(domain), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_boundry), (void*)(boundry), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_manning), (void*)(manning), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxb), (void*)(tsxb), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyb), (void*)(tsyb), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxf), (void*)(tsxf), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyf), (void*)(tsyf), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxc), (void*)(tsxc), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyc), (void*)(tsyc), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_max_u), (void*)(max_u), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_max_v), (void*)(max_v), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sc), (void*)(sc), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_ht), (void*)(ht), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr1), (void*)(hcorr1), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr2), (void*)(hcorr2), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr3), (void*)(hcorr3), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_x), (void*)(x), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_y), (void*)(y), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_surface), (void*)(surface), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_z), (void*)(z), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_h), (void*)(h), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_u), (void*)(u), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_v), (void*)(v), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_h_max), (void*)(h_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_hu_max), (void*)(hu_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_hv_max), (void*)(hv_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hvel_max), (void*)(hvel_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_vel_max), (void*)(vel_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_zs), (void*)(zs), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);

}
Wind_Cuda::Wind_Cuda(FilePath &filepath, int &rain_hours) : Dem(filepath) //type2:变量内存分配及变量读取初始化
{
#ifndef NDEBUG
	std::cerr << "Wind(filePath &filepath)" << std::endl;
	std::cerr << "before set_mbc" << std::endl;
	std::cerr << this->get_xncols() << std::endl;
	std::cerr << this->get_ynrows() << std::endl;
	std::cerr << this->get_mbc_1() << std::endl;
	std::cerr << this->get_mbc_2() << std::endl;
	std::cerr << this->get_mbc_3() << std::endl;
	std::cerr << this->get_mbc_4() << std::endl;
#endif
	set_mbc(0, get_mbc_1(), get_mbc_2(), get_mbc_3(), get_mbc_4(), get_xncols(), get_ynrows());//重新计算mbc_1/2/3/4
	set_hours(rain_hours);//
	set_nxny();
#ifndef NDEBUG
	std::cerr << "xncols: " << nx << std::endl;
	std::cerr << "ynrows: " << ny << std::endl;
#endif
	cudaError_t res;
	point_num = (int*)malloc(group_num * sizeof(int));
	point_x = (int*)malloc(group_num * points_num * sizeof(int));
	point_y = (int*)malloc(group_num * points_num * sizeof(int));

	w = (double**)malloc(nx * sizeof(double*));
	w_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_w), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_w_data), nx*ny * sizeof(double)); CHECK(res);

	sx = (double**)malloc(nx * sizeof(double*));
	sx_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_sx), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_sx_data), nx*ny * sizeof(double)); CHECK(res);

	sy = (double**)malloc(nx * sizeof(double*));
	sy_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_sy), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_sy_data), nx*ny * sizeof(double)); CHECK(res);

	C = (double**)malloc(nx * sizeof(double*));
	C_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_C), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_C_data), nx*ny * sizeof(double)); CHECK(res);

	qx = (double**)malloc(nx * sizeof(double*));
	qx_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_qx), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_qx_data), nx*ny * sizeof(double)); CHECK(res);

	qy = (double**)malloc(nx * sizeof(double*));
	qy_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_qy), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_qy_data), nx*ny * sizeof(double)); CHECK(res);

	R = (double**)malloc(nx * sizeof(double*));
	R_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_R), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_R_data), nx*ny * sizeof(double)); CHECK(res);

	Soil_depth = (double**)malloc(nx * sizeof(double*));
	Soil_depth_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_Soil_depth), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_Soil_depth_data), nx*ny * sizeof(double)); CHECK(res);

	INF = (double**)malloc(nx * sizeof(double*));
	INF_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_data), nx*ny * sizeof(double)); CHECK(res);

	INF_Total = (double**)malloc(nx * sizeof(double*));
	INF_Total_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_Total), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_Total_data), nx*ny * sizeof(double)); CHECK(res);

	INF_tp = (double**)malloc(nx * sizeof(double*));
	INF_tp_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_tp), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_tp_data), nx*ny * sizeof(double)); CHECK(res);

	INF_Ks = (double**)malloc(nx * sizeof(double*));
	INF_Ks_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_Ks), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_Ks_data), nx*ny * sizeof(double)); CHECK(res);

	INF_U = (double**)malloc(nx * sizeof(double*));
	INF_U_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_U), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_U_data), nx*ny * sizeof(double)); CHECK(res);

	INF_Os = (double**)malloc(nx * sizeof(double*));
	INF_Os_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_Os), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_Os_data), nx*ny * sizeof(double)); CHECK(res);

	INF_Oi = (double**)malloc(nx * sizeof(double*));
	INF_Oi_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_INF_Oi), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_INF_Oi_data), nx*ny * sizeof(double)); CHECK(res);

	domain = (double**)malloc(nx * sizeof(double*));
	domain_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_domain), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_domain_data), nx*ny * sizeof(double)); CHECK(res);

	boundry = (double**)malloc(nx * sizeof(double*));
	boundry_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_boundry), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_boundry_data), nx*ny * sizeof(double)); CHECK(res);

	manning = (double**)malloc(nx * sizeof(double*));
	manning_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_manning), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_manning_data), nx*ny * sizeof(double)); CHECK(res);

	tsxb = (double**)malloc(nx * sizeof(double*));
	tsxb_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsxb), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsxb_data), nx*ny * sizeof(double)); CHECK(res);

	tsyb = (double**)malloc(nx * sizeof(double*));
	tsyb_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsyb), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsyb_data), nx*ny * sizeof(double)); CHECK(res);

	tsxf = (double**)malloc(nx * sizeof(double*));
	tsxf_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsxf), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsxf_data), nx*ny * sizeof(double)); CHECK(res);

	tsyf = (double**)malloc(nx * sizeof(double*));
	tsyf_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsyf), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsyf_data), nx*ny * sizeof(double)); CHECK(res);

	tsxc = (double**)malloc(nx * sizeof(double*));
	tsxc_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsxc), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsxc_data), nx*ny * sizeof(double)); CHECK(res);

	tsyc = (double**)malloc(nx * sizeof(double*));
	tsyc_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_tsyc), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_tsyc_data), nx*ny * sizeof(double)); CHECK(res);

	max_u = (double**)malloc(nx * sizeof(double*));
	max_u_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_max_u), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_max_u_data), nx*ny * sizeof(double)); CHECK(res);

	max_v = (double**)malloc(nx * sizeof(double*));
	max_v_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_max_v), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_max_v_data), nx*ny * sizeof(double)); CHECK(res);

	sc = (double**)malloc(nx * sizeof(double*));
	sc_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_sc), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_sc_data), nx*ny * sizeof(double)); CHECK(res);

	ht = (double**)malloc(nx * sizeof(double*));
	ht_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_ht), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_ht_data), nx*ny * sizeof(double)); CHECK(res);

	hcorr = (double**)malloc(nx * sizeof(double*));
	hcorr_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hcorr), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hcorr_data), nx*ny * sizeof(double)); CHECK(res);

	hcorr1 = (double**)malloc(nx * sizeof(double*));
	hcorr1_data= (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hcorr1), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hcorr1_data), nx*ny * sizeof(double)); CHECK(res);

	hcorr2 = (double**)malloc(nx * sizeof(double*));
	hcorr2_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hcorr2), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hcorr2_data), nx*ny * sizeof(double)); CHECK(res);

	hcorr3 = (double**)malloc(nx * sizeof(double*));
	hcorr3_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hcorr3), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hcorr3_data), nx*ny * sizeof(double)); CHECK(res);

	x = (double**)malloc(nx * sizeof(double*));
	x_data= (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_x), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_x_data), nx*ny * sizeof(double)); CHECK(res);

	y = (double**)malloc(nx * sizeof(double*));
	y_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_y), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_y_data), nx*ny * sizeof(double)); CHECK(res);

	surface = (double**)malloc(nx * sizeof(double*));
	surface_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_surface), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_surface_data), nx*ny * sizeof(double)); CHECK(res);

	z = (double**)malloc(nx * sizeof(double*));
	z_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_z), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_z_data), nx*ny * sizeof(double)); CHECK(res);

	h = (double**)malloc(nx * sizeof(double*));
	h_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_h), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_h_data), nx*ny * sizeof(double)); CHECK(res);

	u = (double**)malloc(nx * sizeof(double*));
	u_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_u), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_u_data), nx*ny * sizeof(double)); CHECK(res);

	v = (double**)malloc(nx * sizeof(double*));
	v_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_v), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_v_data), nx*ny * sizeof(double)); CHECK(res);

	h_max = (double**)malloc(nx * sizeof(double*));
	h_max_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_h_max), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_h_max_data), nx*ny * sizeof(double)); CHECK(res);

	// hu_max = (double**)malloc(nx * sizeof(double*));
	// hu_max_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_hu_max), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_hu_max_data), nx*ny * sizeof(double)); CHECK(res);

	// hv_max = (double**)malloc(nx * sizeof(double*));
	// hv_max_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_hv_max), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_hv_max_data), nx*ny * sizeof(double)); CHECK(res);

	hvel_max = (double**)malloc(nx * sizeof(double*));
	hvel_max_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_hvel_max), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_hvel_max_data), nx*ny * sizeof(double)); CHECK(res);

	vel_max = (double**)malloc(nx * sizeof(double*));
	vel_max_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_vel_max), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_vel_max_data), nx*ny * sizeof(double)); CHECK(res);

	vel = (double**)malloc(nx * sizeof(double*));
	vel_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_vel), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_vel_data), nx*ny * sizeof(double)); CHECK(res);

	zs = (double**)malloc(nx * sizeof(double*));
	zs_data = (double*)malloc(nx*ny * sizeof(double));
	// res = cudaMalloc((&dev_zs), nx * sizeof(double*)); CHECK(res);
	// res = cudaMalloc((&dev_zs_data), nx*ny * sizeof(double)); CHECK(res);

	flow = (double**)malloc(nx * sizeof(double*));
	flow_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_flow), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_flow_data), nx*ny * sizeof(double)); CHECK(res);

	res = cudaMalloc((&dev_maxh), sizeof(double)); CHECK(res);
	res = cudaMalloc((&dev_maxu), sizeof(double)); CHECK(res);
	res = cudaMalloc((&dev_maxv), sizeof(double)); CHECK(res);

	landuse = (double**)malloc(nx * sizeof(double*));
	landuse_data= (double*)malloc(nx*ny * sizeof(double));

	V_Cover = (double**)malloc(nx * sizeof(double*));
	V_Cover_data = (double*)malloc(nx*ny * sizeof(double));

	V_LAI = (double**)malloc(nx * sizeof(double*));
	V_LAI_data = (double*)malloc(nx*ny * sizeof(double));

	R_P = (double**)malloc(nx * sizeof(double*));
	R_P_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_R_P), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_R_P_data), nx*ny * sizeof(double)); CHECK(res);

	V_Smax = (double**)malloc(nx * sizeof(double*));
	V_Smax_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_V_Smax), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_V_Smax_data), nx*ny * sizeof(double)); CHECK(res);

	V_K = (double**)malloc(nx * sizeof(double*));
	V_K_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_V_K), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_V_K_data), nx*ny * sizeof(double)); CHECK(res);

	R_Discount = (double**)malloc(nx * sizeof(double*));
	R_Discount_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_R_Discount), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_R_Discount_data), nx*ny * sizeof(double)); CHECK(res);

	R0 = (double**)malloc(nx * sizeof(double*));
	R0_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_R0), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_R0_data), nx*ny * sizeof(double)); CHECK(res);

	V_Ic = (double**)malloc(nx * sizeof(double*));
	V_Ic_data = (double*)malloc(nx*ny * sizeof(double));
	res = cudaMalloc((&dev_V_Ic), nx * sizeof(double*)); CHECK(res);
	res = cudaMalloc((&dev_V_Ic_data), nx*ny * sizeof(double)); CHECK(res);

	R_ALL = new double *[rain_hours];
	for(int i = 0; i < rain_hours; i++)
	{	R_ALL[i] = new double [nx*ny]; }
	for(int i = 0; i < rain_hours; i++)
	{	
		for(int j = 0; j < nx*ny; j++)
		{R_ALL[i][j] = 0.; }
	}

	for (int i = 0; i < nx; i++)
	{
		w[i] = dev_w_data + i * ny;
		sx[i] = dev_sx_data + i * ny;
		sy[i] = dev_sy_data + i * ny;
		C[i] = dev_C_data + i * ny;
		qx[i] = dev_qx_data + i * ny;
		qy[i] = dev_qy_data + i * ny;
		R[i] = dev_R_data + i * ny;
		Soil_depth[i] = dev_Soil_depth_data + i * ny;
		INF[i] = dev_INF_data + i * ny;
		
		INF_Total[i] = dev_INF_Total_data + i * ny;
		INF_tp[i] = dev_INF_tp_data + i * ny;
		INF_Ks[i] = dev_INF_Ks_data + i * ny;
		INF_U[i] = dev_INF_U_data + i * ny;
		INF_Os[i] = dev_INF_Os_data + i * ny;
		INF_Oi[i] = dev_INF_Oi_data + i * ny;
		domain[i] = dev_domain_data + i * ny;
		boundry[i] = dev_boundry_data + i * ny;
		manning[i] = dev_manning_data + i * ny;
		tsxb[i] = dev_tsxb_data + i * ny;
		tsyb[i] = dev_tsyb_data + i * ny;
		tsxf[i] = dev_tsxf_data + i * ny;
		tsyf[i] = dev_tsyf_data + i * ny;
		tsxc[i] = dev_tsxc_data + i * ny;
		tsyc[i] = dev_tsyc_data + i * ny;
		max_u[i] = dev_max_u_data + i * ny;
		max_v[i] = dev_max_v_data + i * ny;
		sc[i] = dev_sc_data + i * ny;
		ht[i] = dev_ht_data + i * ny;
		hcorr[i] = dev_hcorr_data + i * ny;
		hcorr1[i] = dev_hcorr1_data + i * ny;
		hcorr2[i] = dev_hcorr2_data + i * ny;
		hcorr3[i] = dev_hcorr3_data + i * ny;

		// x[i] = dev_x_data + i * ny;
		// y[i] = dev_y_data + i * ny;
		surface[i] = dev_surface_data + i * ny;
		z[i] = dev_z_data + i * ny;
		h[i] = dev_h_data + i * ny;
		u[i] = dev_u_data + i * ny;
		v[i] = dev_v_data + i * ny;
		h_max[i] = dev_h_max_data + i * ny;
		// hu_max[i] = dev_hu_max_data + i * ny;
		// hv_max[i] = dev_hv_max_data + i * ny;
		hvel_max[i] = dev_hvel_max_data + i * ny;
		vel_max[i] = dev_vel_max_data + i * ny;
		// zs[i] = dev_zs_data + i * ny;
		R_P[i] = dev_R_P_data + i * ny;
		flow[i] = dev_flow_data + i * ny;
		V_Smax[i] = dev_V_Smax_data + i * ny;
		V_K[i] = dev_V_K_data + i * ny;
		V_Ic[i] = dev_V_Ic_data + i * ny;
		R_Discount[i] = dev_R_Discount_data + i * ny;
		R0[i] = dev_R0_data + i * ny;
	}
	res = cudaMemcpy((void*)(dev_w), (void*)(w), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sx), (void*)(sx), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sy), (void*)(sy), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_C), (void*)(C), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_qx), (void*)(qx), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_qy), (void*)(qy), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R), (void*)(R), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_Soil_depth), (void*)(Soil_depth), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF), (void*)(INF), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Total), (void*)(INF_Total), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_tp), (void*)(INF_tp), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Ks), (void*)(INF_Ks), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_U), (void*)(INF_U), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Os), (void*)(INF_Os), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Oi), (void*)(INF_Oi), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_domain), (void*)(domain), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_boundry), (void*)(boundry), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_manning), (void*)(manning), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxb), (void*)(tsxb), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyb), (void*)(tsyb), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxf), (void*)(tsxf), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyf), (void*)(tsyf), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxc), (void*)(tsxc), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyc), (void*)(tsyc), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_max_u), (void*)(max_u), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_max_v), (void*)(max_v), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sc), (void*)(sc), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_ht), (void*)(ht), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr), (void*)(hcorr), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr1), (void*)(hcorr1), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr2), (void*)(hcorr2), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr3), (void*)(hcorr3), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_x), (void*)(x), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_y), (void*)(y), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_surface), (void*)(surface), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_z), (void*)(z), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_h), (void*)(h), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_u), (void*)(u), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_v), (void*)(v), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_h_max), (void*)(h_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_hu_max), (void*)(hu_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_hv_max), (void*)(hv_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hvel_max), (void*)(hvel_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_vel_max), (void*)(vel_max), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_zs), (void*)(zs), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_flow), (void*)(flow), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R_P), (void*)(R_P), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_Ic), (void*)(V_Ic), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_Smax), (void*)(V_Smax), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_K), (void*)(V_K), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R_Discount), (void*)(R_Discount), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R0), (void*)(R0), nx * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	double novalue;
	novalue= get_NODATA_VALUE();
	for (int i = 0; i < nx*ny; i++)
	{
		w_data[i] = 0.;
		sx_data[i] = 0.;
		sy_data[i] = 0.;
		C_data[i] = 0.;
		qx_data[i] = 0.;
		qy_data[i] = 0.;
		R_data[i] = 0.;
		Soil_depth_data[i] = 0.2;
		INF_data[i] = 0.;
		INF_Total_data[i] = 0.;
		INF_tp_data[i] = 0.;
		// INF_Ks_data[i] = 1.e-6;
		// INF_U_data[i] = 0.5;
		// INF_Os_data[i] = 0.454;
		// INF_Oi_data[i] = 0.22;
		domain_data[i] = 0.;
		boundry_data[i] = 0.;
		// manning_data[i] = 0.1;
		tsxb_data[i] = 0.;
		tsyb_data[i] = 0.;
		tsxf_data[i] = 0.;
		tsyf_data[i] = 0.;
		tsxc_data[i] = 0.;
		tsyc_data[i] = 0.;
		max_u_data[i] = 0.;
		max_v_data[i] = 0.;
		sc_data[i] = 1.0;
		ht_data[i] = 0.;
		hcorr_data[i] = 0.;
		hcorr1_data[i] = 0.;
		hcorr2_data[i] = 0.;
		hcorr3_data[i] = 0.;
		// x_data[i] = 0.;
		// y_data[i] = 0.;
		surface_data[i] = 0.;
		z_data[i] = novalue;
		h_data[i] = 0.;
		u_data[i] = 0.;
		v_data[i] = 0.;
		h_max_data[i] = 0.;
		// hu_max_data[i] = 0.;
		// hv_max_data[i] = 0.;
		hvel_max_data[i] = 0.;
		vel_max_data[i] = 0.;
		vel_data[i] = 0.;
		zs_data[i] = 0.;
		flow_data[i] = 0.;
		R_P_data[i] = 0.;
		V_Smax_data[i] = 0.;
		V_K_data[i] = 0.;
		V_Ic_data[i] = 0.;
		R_Discount_data[i] = 0.;
		R0_data[i] = 0.;
		landuse_data[i] = 5.;
	}
	for(int i = 0; i < group_num; i++)
	{
		point_num[i] = -99;
		for(int j = 0; j < points_num; j++)
		{
			point_x[ j + i*points_num] = -99;
			point_y[ j + i*points_num] = -99;
		}
	}
	for (int i = 0; i < get_xncols(); i++)
	{
		for (int j = 0; j <= get_ynrows(); j++)
		{
			z_data[j+get_mbc_3() + (i+get_mbc_1()) * (get_ynrows()+get_mbc_3()+get_mbc_4())] = get_Data()[i][j];
		}
	}
	res = cudaMemcpy((void*)(dev_w_data), (void*)(w_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sx_data), (void*)(sx_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sy_data), (void*)(sy_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_C_data), (void*)(C_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_qx_data), (void*)(qx_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_qy_data), (void*)(qy_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R_data), (void*)(R_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_Soil_depth_data), (void*)(Soil_depth_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_data), (void*)(INF_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Total_data), (void*)(INF_Total_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_tp_data), (void*)(INF_tp_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_INF_Ks_data), (void*)(INF_Ks_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_INF_U_data), (void*)(INF_U_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_INF_Os_data), (void*)(INF_Os_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_INF_Oi_data), (void*)(INF_Oi_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_domain_data), (void*)(domain_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_boundry_data), (void*)(boundry_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_manning_data), (void*)(manning_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxb_data), (void*)(tsxb_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyb_data), (void*)(tsyb_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxf_data), (void*)(tsxf_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyf_data), (void*)(tsyf_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsxc_data), (void*)(tsxc_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_tsyc_data), (void*)(tsyc_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_max_u_data), (void*)(max_u_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_max_v_data), (void*)(max_v_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_sc_data), (void*)(sc_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_ht_data), (void*)(ht_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr_data), (void*)(hcorr_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr1_data), (void*)(hcorr1_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr2_data), (void*)(hcorr2_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hcorr3_data), (void*)(hcorr3_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_x_data), (void*)(x_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_y_data), (void*)(y_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_surface_data), (void*)(surface_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_z_data), (void*)(z_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_h_data), (void*)(h_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_u_data), (void*)(u_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_v_data), (void*)(v_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_h_max_data), (void*)(h_max_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_hu_max_data), (void*)(hu_max_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_hv_max_data), (void*)(hv_max_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_hvel_max_data), (void*)(hvel_max_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_vel_max_data), (void*)(vel_max_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_zs_data), (void*)(zs_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_flow_data), (void*)(flow_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R_P_data), (void*)(R_P_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_Smax_data), (void*)(V_Smax_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_K_data), (void*)(V_K_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_Ic_data), (void*)(V_Ic_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R_Discount_data), (void*)(R_Discount_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R0_data), (void*)(R0_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
#ifndef NDEBUG
	std::cerr << "REPLACE OK!" << std::endl;
#endif
	replace_landuse(filepath.inlandusePath.data());
	replace_R0(filepath.inR0Path.data());
	// replace_INF_Ks(filepath.inINF_KsPath.data());
	// replace_INF_U(filepath.inINF_UPath.data());
	// replace_INF_Os(filepath.inINF_OsPath.data());
	// replace_INF_Oi(filepath.inINF_OiPath.data());
	replace_Soil_depth(filepath.inSoil_depthPath.data());
	// replace_H(filepath.inH_initialPath.data());
	// replace_U(filepath.inU_initialPath.data());
	// replace_V(filepath.inV_initialPath.data());
	// replace_manning(filepath.inmanningPath.data());
	int land_type=0 ;
	int type_kind=0 ;
	for (int i = 0; i < nx*ny; i++)
	{
		land_type = floor(landuse_data[i]) ;
		type_kind = floor(landuse_data[i]) ;
		land_type = land_type/10 ;
		type_kind = type_kind%10 ;
		switch(land_type)
		{
			case 0://针对某特殊研究区域
				switch(type_kind)
				{
					case 1://建筑物,不计算植被截留
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;
						V_K_data[i]=0. ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.45;
						manning_data[i] = 0.085;
					// manning_data[i] = 0.085;//system
						break;
					case 2://有作物梯田
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 10. ;
						V_Smax_data[i] = 0.935 + 0.498 * V_LAI_data[i] - 0.00575 * pow(V_LAI_data[i],2) ;//按照农作物
						V_K_data[i]=1 - exp(-V_Cover_data[i] * V_LAI_data[i]) ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.335;
						break;
					case 3://无作物耕地
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;
						V_K_data[i]=0. ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.114;
						break;
					case 4://开阔地
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;
						V_K_data[i]=0. ;
					
						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.12;
						break;
					case 5://稀疏植被
						V_Cover_data[i] = 0.8;
						V_LAI_data[i] = 20. ;
						V_Smax_data[i] = 0.1713 * V_LAI_data[i] ;//按照阔叶林
						V_K_data[i]=1 - exp(-V_Cover_data[i] * V_LAI_data[i]) ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.225;
						break;
					case 6://草地
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 10. ;
						V_Smax_data[i] = 0.59 * pow(V_LAI_data[i] ,0.88);//按照阔叶林
						V_K_data[i]=1 - exp(-V_Cover_data[i] * V_LAI_data[i]) ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.25;
						break;
					case 7://浓密草丛、植被
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 30. ;
						V_Smax_data[i] = 0.2856 * V_LAI_data[i] ;//按照阔叶林
						V_K_data[i]=1 - exp(-V_Cover_data[i] * V_LAI_data[i]) ;

						INF_Ks_data[i] = 1.e-7;
						INF_U_data[i] = 0.2;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.35;
						manning_data[i] = 0.485;
						break;
					case 8://道路
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;//按照阔叶林
						V_K_data[i]=0. ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.45;
						manning_data[i] = 0.085;
						break;
					case 0://其他无法识别
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;//按照草地
						V_K_data[i]=0. ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.25;
						break;
					default:
						cout<<"incorrect landuse type:"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 1://耕地
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 10. ;
						V_Smax_data[i] = 0.935 + 0.498 * V_LAI_data[i] - 0.00575 * pow(V_LAI_data[i],2) ;//按照农作物
						V_K_data[i]=1 - exp(- V_Cover_data[i] * V_LAI_data[i]) ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						// manning_data[i] = 0.1;
						manning_data[i] = 0.1;
						break;
					default:
						cout<<"incorrect landuse 10 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 2://林地
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 30. ;
						V_Smax_data[i] = 0.2856 * V_LAI_data[i] ;//按照阔叶林
						V_K_data[i]=1 - exp(-V_Cover_data[i] * V_LAI_data[i]) ;

						INF_Ks_data[i] = 1.0e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.12;
						break;
					default:
						cout<<"incorrect landuse 20 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 3://草地
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 10. ;
						V_Smax_data[i] = 0.59 * pow(V_LAI_data[i] ,0.88);//按照草地
						V_K_data[i]=1 - exp(-V_Cover_data[i] * V_LAI_data[i]) ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.1;
						break;
					default:
						cout<<"incorrect landuse 30 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 4://灌木地
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] = 20. ;
						V_Smax_data[i] = 0.1713 * V_LAI_data[i] ;//按照阔叶林1/3计算
						V_K_data[i]=1 - exp(-V_Cover_data[i] * V_LAI_data[i]) ;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.2;
						manning_data[i] = 0.09;
						break;
					default:
						cout<<"incorrect landuse 40 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 5://湿地
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0.8 ;
						V_LAI_data[i] =0.;
						V_Smax_data[i] = 0. ;//不计算植被拦截量
						V_K_data[i]=0.;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.45;
						manning_data[i] = 0.02;
						break;
					default:
						cout<<"incorrect landuse 50 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 6://水体
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0. ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;//不计算植被拦截量
						V_K_data[i]=0.;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.45;
						manning_data[i] = 0.015;
						break;
					default:
						cout<<"incorrect landuse 60 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 7://苔原
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0. ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;//不计算植被拦截量
						V_K_data[i]=0.;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.4;
						manning_data[i] = 0.05;
						break;
					default:
						cout<<"incorrect landuse 70 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 8://人造地表
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0. ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;//不计算植被拦截量
						V_K_data[i]=0.;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.45;
						manning_data[i] = 0.015;
						break;
					default:
						cout<<"incorrect landuse 80type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 9://裸地
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0. ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;//不计算植被拦截量
						V_K_data[i]=0.;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.45;
						manning_data[i] = 0.03;
						break;
					default:
						cout<<"incorrect landuse 90 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 10://冰川和永久积雪
				switch(type_kind)
				{
					case 0:
						V_Cover_data[i] = 0. ;
						V_LAI_data[i] = 0. ;
						V_Smax_data[i] = 0. ;//不计算植被拦截量
						V_K_data[i]=0.;

						INF_Ks_data[i] = 1.e-6;
						INF_U_data[i] = 0.1;
						INF_Os_data[i] = 0.45;
						INF_Oi_data[i] = 0.45;
						manning_data[i] = 0.015;
						break;
					default:
						cout<<"incorrect landuse 100 type_kind"<<type_kind<<endl;
						exit(EXIT_FAILURE);
				}
				break;
			default:
				V_Cover_data[i] = 0. ;
				V_LAI_data[i] = 0. ;
				V_Smax_data[i] = 0. ;//按照阔叶林
				V_K_data[i]=0. ;

				INF_Ks_data[i] = 1.e-6;
				INF_U_data[i] = 0.1;
				INF_Os_data[i] = 0.45;
				INF_Oi_data[i] = 0.45;
				manning_data[i] = 0.015;
		}
		INF_Os_data[i] = INF_Os_data[i] - INF_Oi_data[i];
		manning_data[i] = 1./manning_data[i];
		if(INF_Os_data[i]<=0.){INF_Os_data[i] = EPS;}
	}

	// res = cudaMemcpy((void*)(dev_R_data), (void*)(R_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_Soil_depth_data), (void*)(Soil_depth_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Ks_data), (void*)(INF_Ks_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_U_data), (void*)(INF_U_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_INF_Os_data), (void*)(INF_Os_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_R_P_data), (void*)(R_P_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_Smax_data), (void*)(V_Smax_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_K_data), (void*)(V_K_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_V_Ic_data), (void*)(V_Ic_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);

	// res = cudaMemcpy((void*)(dev_h_data), (void*)(h_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_u_data), (void*)(u_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// res = cudaMemcpy((void*)(dev_v_data), (void*)(v_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	res = cudaMemcpy((void*)(dev_manning_data), (void*)(manning_data), nx*ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// outDem("./output/H_initial", h_data, 0.);
	// outDem("./output/U_initial", u_data, 0.);
	// outDem("./output/V_initial", v_data, 0.);
	// outDem("./output/INF_Ks", INF_Ks_data, 0.);
	// outDem("./output/INF_U", INF_U_data, 0.);
	// outDem("./output/INF_Os", INF_Os_data, 0.);
	// outDem("./output/INF_Oi", INF_Oi_data, 0.);
	// outDem("./output/Soil_depth", Soil_depth_data, 0.);
	// outDem("./output/manning", manning_data, 0.);
	// exit(EXIT_FAILURE);
}
Wind_Cuda::~Wind_Cuda()
{
#ifndef NDEBUG
	std::cerr << "~Wind()" << std::endl;
	std::cerr << __func__ << std::endl;
#endif
	delete[] point_num;
	delete[] point_x;
	delete[] point_y;
	// delete[] point_data;
	// cudaFree((void*)dev_point_data);

	delete[] w;
	delete[] w_data;
	cudaFree((void*)dev_w);
	cudaFree((void*)dev_w_data);

	delete[] sx;
	delete[] sx_data;
	cudaFree((void*)dev_sx);
	cudaFree((void*)dev_sx_data);

	delete[] sy;
	delete[] sy_data;
	cudaFree((void*)dev_sy);
	cudaFree((void*)dev_sy_data);

	delete[] C;
	delete[] C_data;
	cudaFree((void*)dev_C);
	cudaFree((void*)dev_C_data);

	delete[] qx;
	delete[] qx_data;
	cudaFree((void*)dev_qx);
	cudaFree((void*)dev_qx_data);

	delete[] qy;
	delete[] qy_data;
	cudaFree((void*)dev_qy);
	cudaFree((void*)dev_qy_data);

	delete[] R;
	delete[] R_data;
	cudaFree((void*)dev_R);
	cudaFree((void*)dev_R_data);

	delete[] Soil_depth;
	delete[] Soil_depth_data;
	cudaFree((void*)dev_Soil_depth);
	cudaFree((void*)dev_Soil_depth_data);

	delete[] INF;
	delete[] INF_data;
	cudaFree((void*)dev_INF);
	cudaFree((void*)dev_INF_data);

	delete[] INF_Total;
	delete[] INF_Total_data;
	cudaFree((void*)dev_INF_Total);
	cudaFree((void*)dev_INF_Total_data);

	delete[] INF_tp;
	delete[] INF_tp_data;
	cudaFree((void*)dev_INF_tp);
	cudaFree((void*)dev_INF_tp_data);

	delete[] INF_Ks;
	delete[] INF_Ks_data;
	cudaFree((void*)dev_INF_Ks);
	cudaFree((void*)dev_INF_Ks_data);

	delete[] INF_U;
	delete[] INF_U_data;
	cudaFree((void*)dev_INF_U);
	cudaFree((void*)dev_INF_U_data);

	delete[] INF_Os;
	delete[] INF_Os_data;
	cudaFree((void*)dev_INF_Os);
	cudaFree((void*)dev_INF_Os_data);

	delete[] INF_Oi;
	delete[] INF_Oi_data;
	cudaFree((void*)dev_INF_Oi);
	cudaFree((void*)dev_INF_Oi_data);

	delete[] domain;
	delete[] domain_data;
	cudaFree((void*)dev_domain);
	cudaFree((void*)dev_domain_data);

	delete[] boundry;
	delete[] boundry_data;
	cudaFree((void*)dev_boundry);
	cudaFree((void*)dev_boundry_data);

	delete[] manning;
	delete[] manning_data;
	cudaFree((void*)dev_manning);
	cudaFree((void*)dev_manning_data);

	delete[] tsxb;
	delete[] tsxb_data;
	cudaFree((void*)dev_tsxb);
	cudaFree((void*)dev_tsxb_data);

	delete[] tsyb;
	delete[] tsyb_data;
	cudaFree((void*)dev_tsyb);
	cudaFree((void*)dev_tsyb_data);

	delete[] tsxf;
	delete[] tsxf_data;
	cudaFree((void*)dev_tsxf);
	cudaFree((void*)dev_tsxf_data);

	delete[] tsyf;
	delete[] tsyf_data;
	cudaFree((void*)dev_tsyf);
	cudaFree((void*)dev_tsyf_data);

	delete[] tsxc;
	delete[] tsxc_data;
	cudaFree((void*)dev_tsxc);
	cudaFree((void*)dev_tsxc_data);

	delete[] tsyc;
	delete[] tsyc_data;
	cudaFree((void*)dev_tsyc);
	cudaFree((void*)dev_tsyc_data);

	delete[] max_u;
	delete[] max_u_data;
	cudaFree((void*)dev_max_u);
	cudaFree((void*)dev_max_u_data);

	delete[] max_v;
	delete[] max_v_data;
	cudaFree((void*)dev_max_v);
	cudaFree((void*)dev_max_v_data);

	delete[] sc;
	delete[] sc_data;
	cudaFree((void*)dev_sc);
	cudaFree((void*)dev_sc_data);

	delete[] ht;
	delete[] ht_data;
	cudaFree((void*)dev_ht);
	cudaFree((void*)dev_ht_data);

	delete[] hcorr;
	delete[] hcorr_data;
	cudaFree((void*)dev_hcorr);
	cudaFree((void*)dev_hcorr_data);

	delete[] hcorr1;
	delete[] hcorr1_data;
	cudaFree((void*)dev_hcorr1);
	cudaFree((void*)dev_hcorr1_data);

	delete[] hcorr2;
	delete[] hcorr2_data;
	cudaFree((void*)dev_hcorr2);
	cudaFree((void*)dev_hcorr2_data);

	delete[] hcorr3;
	delete[] hcorr3_data;
	cudaFree((void*)dev_hcorr3);
	cudaFree((void*)dev_hcorr3_data);

	delete[] x;
	delete[] x_data;
	// cudaFree((void*)dev_x);
	// cudaFree((void*)dev_x_data);

	delete[] y;
	delete[] y_data;
	// cudaFree((void*)dev_y);
	// cudaFree((void*)dev_y_data);

	delete[] surface;
	delete[] surface_data;
	cudaFree((void*)dev_surface);
	cudaFree((void*)dev_surface_data);

	delete[] z;
	delete[] z_data;
	cudaFree((void*)dev_z);
	cudaFree((void*)dev_z_data);

	delete[] h;
	delete[] h_data;
	cudaFree((void*)dev_h);
	cudaFree((void*)dev_h_data);

	delete[] u;
	delete[] u_data;
	cudaFree((void*)dev_u);
	cudaFree((void*)dev_u_data);

	delete[] v;
	delete[] v_data;
	cudaFree((void*)dev_v);
	cudaFree((void*)dev_v_data);

	delete[] h_max;
	delete[] h_max_data;
	cudaFree((void*)dev_h_max);
	cudaFree((void*)dev_h_max_data);

	// delete[] hu_max;
	// delete[] hu_max_data;
	// cudaFree((void*)dev_hu_max);
	// cudaFree((void*)dev_hu_max_data);

	// delete[] hv_max;
	// delete[] hv_max_data;
	// cudaFree((void*)dev_hv_max);
	// cudaFree((void*)dev_hv_max_data);

	delete[] hvel_max;
	delete[] hvel_max_data;
	cudaFree((void*)dev_hvel_max);
	cudaFree((void*)dev_hvel_max_data);

	delete[] vel_max;
	delete[] vel_max_data;
	cudaFree((void*)dev_vel_max);
	cudaFree((void*)dev_vel_max_data);

	delete[] vel;
	delete[] vel_data;

	delete[] zs;
	delete[] zs_data;

	delete[] V_Cover;
	delete[] V_Cover_data;
	
	delete[] V_LAI;
	delete[] V_LAI_data;

	delete[] flow;
	delete[] flow_data;
	cudaFree((void*)dev_flow);
	cudaFree((void*)dev_flow_data);
	// for(int i = 0; i < 9; i++)
	// {delete[] R_ALL[i];}

	delete[] R_P;
	delete[] R_P_data;
	cudaFree((void*)dev_R_P);
	cudaFree((void*)dev_R_P_data);

	delete[] R0;
	delete[] R0_data;
	cudaFree((void*)dev_R0);
	cudaFree((void*)dev_R0_data);

	delete[] V_Smax;
	delete[] V_Smax_data;
	cudaFree((void*)dev_V_Smax);
	cudaFree((void*)dev_V_Smax_data);

	delete[] V_K;
	delete[] V_K_data;
	cudaFree((void*)dev_V_K);
	cudaFree((void*)dev_V_K_data);

	delete[] V_Ic;
	delete[] V_Ic_data;
	cudaFree((void*)dev_V_Ic);
	cudaFree((void*)dev_V_Ic_data);

	delete[] R_Discount;
	delete[] R_Discount_data;
	cudaFree((void*)dev_R_Discount);
	cudaFree((void*)dev_R_Discount_data);

	delete[] R_ALL;

	cudaFree((void*) dev_maxh);
	cudaFree((void*) dev_maxu);
	cudaFree((void*) dev_maxv);
}

void Wind_Cuda::set_R(int hour, double *&R_data, double **R_ALL)//一维到一维
{
	int xncols=get_xncols();
	int ynrows=get_ynrows();
	int mbc_1=get_mbc_1();
	int mbc_2=get_mbc_2();
	int mbc_3=get_mbc_3();
	int mbc_4=get_mbc_4();
	
	for (int i = 0; i < xncols; i++)
	{
		for (int j = 0; j <= ynrows; j++)
		{
			R_data[j+mbc_3 + (i+mbc_1) * ny] = R_ALL[hour-1][j+mbc_3 + (i+mbc_1) * ny];
		}
	}
}

void Wind_Cuda::set_R_ALL(int &hours, double area, int R_type)
{
#ifndef NDEBUG
	cout << "set_R_ALL" << endl;
	cout << __func__ << endl;
#endif
	cudaError_t res;
	int xncols=get_xncols();
	int ynrows=get_ynrows();
	int mbc_1=get_mbc_1();
	int mbc_2=get_mbc_2();
	int mbc_3=get_mbc_3();
	int mbc_4=get_mbc_4();
	double Cellsize = get_Cellsize() * get_Cellsize();
	std::cout<<"mbc_1"<<mbc_1<<endl;
	std::cout<<"mbc_2"<<mbc_2<<endl;
	std::cout<<"mbc_3"<<mbc_3<<endl;
	std::cout<<"mbc_4"<<mbc_4<<endl;
	double no_value = get_NODATA_VALUE();

	if( R_type == 1)
	{
		hours = 1;
		FILE *Rain_file;
		string filename = "./input/R/R1.txt";
		const char *fname = filename.data();

		if((Rain_file = fopen(fname, "rb")) == NULL)//R如果读取失败，停止执行
		{std::cout << "Failed to open file:" << filename << std::endl;exit(EXIT_FAILURE);}

		for(int i=1;i > 0;i++)//依次打开降雨文件，记录文件总数
		{
			hours=i;
			filename = "./input/R/R" + std::to_string(hours) + ".txt";
			fname = filename.data();
			if((Rain_file = fopen(fname, "rb")) == NULL){i=-1;hours -= 1;}
		}
		cout<<"Files:"<<hours<<endl;

		for(int i=1; i <= hours; i++)//依次读取降雨文件内容
		{
			filename = "./input/R/R" + to_string(i) + ".txt";
			Dem dem(filename); //实例化一个DEM
			this->swap_Copy_Cuda_Data(R_ALL[i-1], dem.get_Data());
		} 
	}
	else if(R_type ==2)
	{
		FILE *R_file;
		string R_name = "./input/R.txt";
		const char *R_fname = R_name.data();

		cout << "in set_R_ALL:"<<endl;
		if((R_file = fopen(R_fname, "rb")) == NULL)//查询分布式降雨文件R.txt，没有停止执行
		{std::cout << "Failed to open file:" << R_name << std::endl;exit(EXIT_FAILURE);}
		else{std::cout<<"get R.txt"<<endl;}

		string R_point = "point";//分布式降雨
		const char *str = R_point.data();//分布式降雨
		int point_num = 9;//分布式降雨
		int time_num = 218;//分布式降雨
		int title_x[10];//分布式降雨
		int title_y[10];//分布式降雨
		R_point_value = new double *[point_num];//分布式降雨
		for( int i = 0; i < point_num ;i++)//
		{
			R_point_value[i] = new double [time_num];
			title_x[i] = 0;
			title_y[i] = 0;
		}
		for(int i = 0; i < point_num; i++)//分布式降雨
		{	
			for(int j = 0; j < time_num; j++)
			{R_point_value[i][j] = 0.; }
		}
		fscanf(R_file ,"%d", &point_num);//分布式降雨
		fscanf(R_file ,"%d", &time_num);//分布式降雨
		hours = time_num ;
		std::cout<<"point_num"<<point_num<<"time_num"<<time_num<<endl;//分布式降雨
			for( int i = 0; i < point_num ;i++)//获取降雨监测点位置
			{
				fscanf(R_file ,"%s", &str);
				printf("%s: ", &str);
				fscanf(R_file ,"%d" , &title_x[i]);
				printf("%d\n", title_x[i]);
				fscanf(R_file ,"%d" , &title_y[i]);
				printf("%d\n", title_y[i]);
			}
			for(int j = 0; j < time_num ;j++)//读取降雨监测点时间序列降雨数据
			{
				for(int i = 0; i < point_num ;i++)
				{
					fscanf(R_file, "%lf", &R_point_value[i][j]);
					std::cout<<"R_point_value["<<i<<"]["<<j<<"]"<<R_point_value[i][j]<<endl;
				}
			}
		double dis_per = 0.;
		double R_total = 0.;
		distance = new double [point_num] ;

		for(int Time = 0; Time < time_num ;Time++)//遍历时间序列
		{
			for (int j = ynrows + mbc_3 - 1;j >= mbc_3;  j--)//遍历所有网格点
			{
				for (int i = mbc_1; i < xncols + mbc_1; i++)
				{
					dis_per = 0.;
					R_total = 0.;
					for(int p = 0; p < point_num; p++)//遍历监测点,如果位置位于监测点，数据特殊处理
					{
						distance[p] = sqrt( pow(abs(i - title_x[p]- mbc_1),2) + pow(abs(j - title_y[p] - mbc_3 ),2) );//计算点的距离
						dis_per = dis_per + 1. / pow(distance[p], 2);//计算总权
					}
					for(int p=0; p < point_num; p++)//根据反距离计算权重
					{
						distance[p] = 1. / pow(distance[p] ,2)/ dis_per ;
					}
					for(int p=0; p < point_num; p++)//根据反距离插值给每个点赋值
					{
							R_total = R_total + R_point_value[p][Time] * distance[p] ;
					}
					if(domain_data[j + i * (ynrows + mbc_3 + mbc_4)] < 0.5)
					{
						R_ALL[Time][j + i * (ynrows + mbc_3 + mbc_4)] = R_total;
					}
					else
					{
						R_ALL[Time][j + i * (ynrows + mbc_3 + mbc_4)] = no_value;
					}
					for(int p = 0; p < point_num; p++)//遍历监测点,如果位置位于监测点，数据特殊处理
					{
						distance[p] = sqrt( pow(abs(i - title_x[p]- mbc_1),2) + pow(abs(j - title_y[p] - mbc_3 ),2) );//计算点的距离
						if(distance[p] < EPS)//锁定监测点特殊赋值
						{
							R_ALL[Time][j + i * (ynrows + mbc_3 + mbc_4)] = R_point_value[p][Time];
						}
					}
				}
			}
		}
	}
	else if(R_type ==3)
	{
		FILE *R_file;
		string R_name = "./input/R.txt";
		const char *R_fname = R_name.data();
		int time_num = 72;//分布式降雨

		cout << "in set_R_ALL:"<<endl;
		if((R_file = fopen(R_fname, "rb")) == NULL)//查询分布式降雨文件R.txt，没有停止执行
		{std::cout << "Failed to open file:" << R_name << std::endl;exit(EXIT_FAILURE);}
		else{std::cout<<"get R.txt"<<endl;}
		
		char str[100];
		fscanf(R_file ,"%d", &time_num);//降雨时序数
		hours = time_num ;//设置计算总时间

		//double title[time_num];
		double *title;
		title = (double*)malloc(time_num * sizeof(double));


		for (int i = 0; i <time_num; i++)//读取降雨时序
		{
			fscanf(R_file, "%s", &str);
			fscanf(R_file, "%lf", &title[i]);
			cout<<"R"<<i<<":"<<title[i]<<endl;
		}

			for (int k = 0; k <= time_num;  k++)//遍历所有网格点。赋值降雨
			{
				for (int j = ynrows + mbc_3 - 1;j >= mbc_3;  j--)//遍历所有网格点
				{
					for (int i = mbc_1; i < xncols + mbc_1; i++)
					{
						R_ALL[k][j + i * (ynrows + mbc_3 + mbc_4)] = title[k];
					}
				}
			}
	}
	for(int i=0; i<nx*ny; i++)//设置地形区域为计算范围
	{
		if(z_data[i]>9000 || z_data[i] < 0)
		{domain_data[i] = 1.; }
		else
		{area = area + 1;}
		if(R0_data[i] ==1.)
		{channel_num = channel_num + 1;}
	}
	area = area * Cellsize/1000000.;
	cout<<"channel_num:"<<channel_num<<endl;
	cout<<"area:"<<area<<" km2"<<endl;
	set_area(area);
	res = cudaMemcpy((void*)(dev_domain_data), (void*)(domain_data), nx *ny * sizeof(double*), cudaMemcpyHostToDevice); CHECK(res);
	// delete[] distance;
	// delete[] R_point_value;
}

void Wind_Cuda::replace_landuse(string filepath)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *file = fopen(filepath.c_str(), "r"))
	{
		fclose(file);
	}
	else
	{
#ifndef NDEBUG
		cout << "R File read failed!!" << endl;
#endif
		return;
	}
	cout << "in  replace Landuse.txt" << endl;
	Dem dem(filepath); //实例化一个DEM
	//R=dem.get_Data_COPY();
	this->swap_Copy_Cuda_Data(landuse_data, dem.get_Data());
	cout << "Get Landuse.txt" << endl;
}

void Wind_Cuda::replace_R0(string filepath)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *file = fopen(filepath.c_str(), "r"))
	{
		fclose(file);
	}
	else
	{
		cout << "channel file read failed!!" << endl;
		return;
	}
	Dem dem(filepath); //实例化一个DEM
	//R=dem.get_Data_COPY();
	this->swap_Copy_Cuda_Data(R0_data, dem.get_Data());
}
void Wind_Cuda::replace_Soil_depth(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
// #ifndef NDEBUG
		cout << "Soil_depth File read failed!!" << endl;
// #endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	//v=dem.get_Data_COPY();
	this->swap_Copy_Cuda_Data(Soil_depth_data, dem.get_Data());
	cout << "Get Soil_depth.txt" << endl;
}
void Wind_Cuda::replace_INF_Ks(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
#ifndef NDEBUG
		cout << "INF_Ks File read failed!!" << endl;
#endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	//v=dem.get_Data_COPY();
	this->swap_Copy_Cuda_Data(INF_Ks_data, dem.get_Data());
}
void Wind_Cuda::replace_INF_U(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
#ifndef NDEBUG
		cout << "INF_U File read failed!!" << endl;
#endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	//v=dem.get_Data_COPY();
	this->swap_Copy_Cuda_Data(INF_U_data, dem.get_Data());
}
void Wind_Cuda::replace_INF_Os(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
#ifndef NDEBUG
		cout << "INF_Os File read failed!!" << endl;
#endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	this->swap_Copy_Cuda_Data(INF_Os_data, dem.get_Data());
}
void Wind_Cuda::replace_INF_Oi(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
#ifndef NDEBUG
		cout << "INF_Oi File read failed!!" << endl;
#endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	v=dem.get_Data_COPY();
	this->swap_Copy_Cuda_Data(INF_Oi_data, dem.get_Data());
}
void Wind_Cuda::replace_H(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
#ifndef NDEBUG
		cout << "h_iniital File read failed!!" << endl;
#endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	//h=dem.get_Data_COPY();
	swap_Copy_Cuda_Data(h_data, dem.get_Data());
	cout << "Get h_iniital" << endl;
}
void Wind_Cuda::replace_U(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
#ifndef NDEBUG
		cout << "u_iniital file read failed!!" << endl;
#endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	//u=dem.get_Data_COPY();
	swap_Copy_Cuda_Data(u_data, dem.get_Data());	
	cout << "Get u_iniital" << endl;
}
void Wind_Cuda::replace_V(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
#ifndef NDEBUG
		cout << "v_iniital File read failed!!" << endl;
#endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	//v=dem.get_Data_COPY();
	swap_Copy_Cuda_Data(v_data, dem.get_Data());
	cout << "Get v_iniital" << endl;
}
void Wind_Cuda::replace_manning(string file)
{
#ifndef NDEBUG
	cout << "Instatiation" << endl;
	cout << __func__ << endl;
#endif
	if (FILE *filepath = fopen(file.c_str(), "r"))
	{
		fclose(filepath);
	}
	else
	{
#ifndef NDEBUG
		cout << "manning File read failed!!" << endl;
#endif
		return;
	}
	Dem dem(file); //实例化一个DEM
	//v=dem.get_Data_COPY();
	swap_Copy_Cuda_Data(manning_data, dem.get_Data());
	cout << "Get v_iniital" << endl;
}

void Wind_Cuda::outPlt()
{
}
/**
* 复制相关的指针数据
* tar->str
*/
template <class T>
void Wind_Cuda::swap_Copy_Data(T **&str, T **&tar)
{
	try
	{
		for (int i = get_mbc_1(); i < get_xncols() + get_mbc_1(); i++)
		{
			for (int j = get_mbc_3(); j < get_ynrows() + get_mbc_3(); j++)
			{
				str[i][j] = tar[i - get_mbc_1()][j - get_mbc_3()];
			}
		}
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}
}
template <class T>
void Wind_Cuda::swap_Copy_Cuda_Data(T *&str, T **&tar)//二维到一维
{
	int xncols=get_xncols();
	int ynrows=get_ynrows();
	int mbc_1=get_mbc_1();
	int mbc_2=get_mbc_2();
	int mbc_3=get_mbc_3();
	int mbc_4=get_mbc_4();
	try
	{
		cudaError_t res;
		for (int i = 0; i < xncols; i++)
		{
			for (int j = 0; j < ynrows; j++)
			{
				str[j+mbc_3 + (i+mbc_1) * ny] = tar[i][j];
			}
		}
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}
}
void Wind_Cuda::set_nxny()
{
	nx = get_xncols() + get_mbc_1() + get_mbc_2();
	ny = get_ynrows() + get_mbc_3() + get_mbc_4();
}
void Wind_Cuda::set_xy(double *&x_data, double *&y_data)
{
	int ynrows=get_ynrows();
	int xncols=get_xncols();
	int mbc_1=get_mbc_1();
	int mbc_2=get_mbc_2();
	int mbc_3=get_mbc_3();
	int mbc_4=get_mbc_4();
	double xll=get_xllcorner();
	double yll=get_yllcorner();
	double Cellsize=get_Cellsize();
	for (int i = mbc_1; i < xncols + mbc_1; i++)
	{
		for (int j = mbc_3; j < ynrows + mbc_3; j++)
		{
			x_data[j + i * (ynrows + mbc_3 + mbc_4)] = xll + (i - mbc_1)* Cellsize;
			y_data[j + i * (ynrows + mbc_3 + mbc_4)] = yll + (j - mbc_3)* Cellsize;
		}
	}
}
void Wind_Cuda::out_Plt(float cur_time, string filename, string ss, int count, ...)
{
#ifndef NDEBUG
	cout << "OUT PLT FILE" << endl;
	cout << __func__ << endl;
#endif
	// int fileOK=1;
	int ynrows=get_ynrows();
	int xncols=get_xncols();
	int mbc_1=get_mbc_1();
	int mbc_2=get_mbc_2();
	int mbc_3=get_mbc_3();
	int mbc_4=get_mbc_4();
	const char *fname = filename.data();
	FILE *fp; //打开文件
	if ((fp= fopen(fname, "a"))== NULL)
	{
		// std::cout << fileOK << std::endl;
		std::cout << "Failed to open file:" << fname << std::endl;
		exit(EXIT_FAILURE);
	}
	// std::cout <<"VARIABLES= "<< ss.data() << std::endl;
	fprintf(fp, "  VARIABLES =");
	fprintf(fp, ss.data());
	fprintf(fp, "\n");
	fprintf(fp, "ZONE T=\"");
	fprintf(fp, "%14.6f\"", cur_time);
	fprintf(fp, " i=");
	fprintf(fp, "%6d", get_xncols());
	fprintf(fp, ",j=");
	fprintf(fp, "%6d", get_ynrows());
	fprintf(fp, "\n");
	// int c = 0;
	for (int j = ynrows + mbc_3 - 1; j >= mbc_3; j--)
	{
		for (int i = mbc_1; i < xncols +mbc_1; i++)
		{
			va_list args;
			va_start(args, count);
			for (int k = 0; k < count; k++)
			{
				// std::cout << (va_arg(args, double **))[i][j];
				// if((va_arg(args, double **))[i][j] == get_NODATA_VALUE())
				//     fprintf(fp, "%22.8lf", get_NODATA_VALUE());
				// else
				fprintf(fp, "%14.6lf", (va_arg(args, double *))[j + i * (ynrows + mbc_3 + mbc_4)]);
				// fprintf(fp, "  ");
				// c++;
			}
			fprintf(fp, "\n");
			va_end(args);
			// std::cout << "    " << c << std::endl;
			// std::cout << b.get_h()[i][j] << std::endl;
		}

		// std::cout << *b.get_h()[i] << std::endl;
	}
	fclose(fp);
}
void Wind_Cuda::out_Point(double *flow_data, double cur_time, string filename, int points, int group)
{
	double Q = 0.;
	// stringstream ss;
	// string group_string;
	// ss << group;
	// ss >> group_string;
	filename = filename + std::to_string(group) + ".txt";
	const char *fname = filename.data();
	FILE* fp;
	if ((fp= fopen(fname, "a"))== NULL)
	{std::cout << "Failed to open  file:" << fname << std::endl; exit(1);}
	// std::cout<<"in out_point function"<<endl;
	// std::cout<<"filename"<<filename<<endl;
	// std::cout<<"group"<<group<<endl;
	// std::cout<<"points"<<points<<endl;
	for(int i = 0 ; i < points ; i++)
	{	
		// std::cout<<"points:"<<points<<endl;
			// std::cout<<"point_x:"<<point_x[i + (group-1) * points_num]<<endl;
			// std::cout<<"point_y:"<<point_y[i + (group-1) * points_num]<<endl;
		if(point_x[i + (group-1) * points_num] > 0 & point_y[i + (group-1) * points_num] > 0)
		{
			// std::cout<<flow_data[ (point_y[i + (group-1) * points_num] - 1 + get_mbc_1())*ny + (point_x[i + (group-1) * points_num] - 1 + get_mbc_3())  ]<<endl;
			Q = Q + flow_data[ (point_y[i + (group-1) * points_num] - 1 + get_mbc_1())*ny + (point_x[i + (group-1) * points_num] - 1 + get_mbc_3())  ];
		}
	}
	fprintf(fp, "%-20.8lf  ", cur_time);
	fprintf(fp, "%-20.8lf\n", Q);
	fclose(fp);
}
