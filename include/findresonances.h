
#pragma once
#include <armadillo>
#include <map>
#include <cstring>



using namespace arma;
using namespace std;

const std::complex<double> IU(0, 1);
const double PI = 3.1415926;


struct FsDetaild
{
	vec frequency;  
	vec decay_constant;
	vec Q_factor;
	vec amplitude;
	vec phase;
	vec error_estimate;
};


struct FdmMatrix
{
	cx_mat U0;
	cx_mat U1;
	cx_mat U2;
};

FdmMatrix fdmMatrixGChenGuo(
	const cx_vec & c,
	const int M,
	const cx_vec & z );

//# 输入
//# c : 时域信号          复数
//# t   :时间序列         
//# freq_window : 频率窗  {100 ， 200}
//# Jnum : 初始基函数数量    
//# rc : 消除特征值的百分比，小于该值去掉特征分量
void findresonances(
	const cx_vec& complex_signal,
	const vec& time_series,
	const vec& freqWin,
	uword j_basis_count=0,
	double removal_criteria=1e-6,
	uword max_iterations=100,
	double error_threshold = 1e-10);