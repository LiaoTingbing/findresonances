

#include "../include/findresonances.h"

FdmMatrix fdmMatrixGChenGuo(const cx_vec& complex_signal,
	const int M,
	const cx_vec& z)
{
	uword jNum = z.n_elem;

	cx_mat U0(jNum, jNum, fill::zeros);
	cx_mat U1(jNum, jNum, fill::zeros);
	cx_mat U2(jNum, jNum, fill::zeros);


	vec nv1 = linspace(0, M, M + 1);
	vec nvd = linspace(0, M * 2, M * 2 + 1);

	cx_vec c10 = complex_signal.rows(0, M);				//M个
	cx_vec c30 = complex_signal.rows(M + 1, 2 * M + 1);			//M个
	cx_vec c11 = complex_signal.rows(1, M + 1);					//M个
	cx_vec c31 = complex_signal.rows(M + 2, 2 * M + 2);			//M个
	cx_vec c12 = complex_signal.rows(2, M + 2);				//M个
	cx_vec c32 = complex_signal.rows(M + 3, 2 * M + 3);		//M个

	cx_vec G10(jNum, fill::zeros);
	cx_vec G30(jNum, fill::zeros);
	cx_vec G11(jNum, fill::zeros);
	cx_vec G31(jNum, fill::zeros);
	cx_vec G12(jNum, fill::zeros);
	cx_vec G32(jNum, fill::zeros);

	for (size_t i = 0; i < jNum; i++)
	{
		cx_vec ulk = pow(
			cx_vec(nv1.n_elem, fill::value(z(i))),
			-nv1);

		G10(i) = dot(ulk, c10);
		G30(i) = dot(ulk, c30);
		G11(i) = dot(ulk, c11);
		G31(i) = dot(ulk, c31);
		G12(i) = dot(ulk, c12);
		G32(i) = dot(ulk, c32);

	}
	for (size_t i = 0; i < jNum; i++)
	{
		cx_double u1 = z(i);
		cx_double u1K = pow(u1, -M);

		cx_double G10i = G10(i);
		cx_double G11i = G11(i);
		cx_double G12i = G12(i);
		cx_double G30i = G30(i);
		cx_double G31i = G31(i);
		cx_double G32i = G32(i);

		for (size_t j = i + 1; j < jNum; j++)
		{
			cx_double u2 = z(j);
			cx_double u2K = pow(u2, -M);
			cx_double coff = 1.0 / (u1 - u2);

			U0(i, j) = coff * (u1 * G10(j) - u2 * G10i - u1K * G30(j) + u2K * G30i);
			U1(i, j) = coff * (u1 * G11(j) - u2 * G11i - u1K * G31(j) + u2K * G31i);
			U2(i, j) = coff * (u1 * G12(j) - u2 * G12i - u1K * G32(j) + u2K * G32i);

			U0(j, i) = coff * (u1 * G10(j) - u2 * G10i - u1K * G30(j) + u2K * G30i);
			U1(j, i) = coff * (u1 * G11(j) - u2 * G11i - u1K * G31(j) + u2K * G31i);
			U2(j, i) = coff * (u1 * G12(j) - u2 * G12i - u1K * G32(j) + u2K * G32i);
		}
	}

	// 2M+1长度
	vec Kd = M - abs(M - nvd) + 1;

	cx_vec	cd0 = complex_signal.rows(0, 2 * M) % Kd;
	cx_vec	cd1 = complex_signal.rows(0 + 1, 2 * M + 1) % Kd;
	cx_vec	cd2 = complex_signal.rows(0 + 2, 2 * M + 2) % Kd;

	// 计算对角

	for (size_t i = 0; i < jNum; i++)
	{
		cx_vec ud = pow(
			cx_vec(nvd.n_elem, fill::value(z(i))), -nvd);

		U0(i, i) = dot(ud, cd0);
		U1(i, i) = dot(ud, cd1);
		U2(i, i) = dot(ud, cd2);

	}
	return { U0 , U1 , U2 };
}


void findresonances(
	const cx_vec& input_signal,
	const vec& time_series,
	const vec& freqWin,
	uword j_basis_count,
	double removal_criteria,
	uword max_iterations,
	double error_threshold)
{
	//c.print();
	double time_step = time_series(1) - time_series(0);
	uword signal_length = input_signal.n_elem;
	uword signal_half_length = (signal_length - 4) / 2;
	if (j_basis_count == 0)
	{
		j_basis_count = ceil(signal_length * time_step * 2.0 * PI * (freqWin(1) - freqWin(0)) / 4.0 / PI);
	}
	cx_vec z = exp(-IU * time_step * 2.0 * PI * linspace(freqWin(0), freqWin(1), j_basis_count));

	cx_vec zNxet(size(z), fill::ones);
	cx_mat v_singular_matrix_reduced, P;
	mat singular_values_matrix_reduced;
	cx_vec u;	//特征值
	FdmMatrix fdm_matrix;
	for (size_t i = 0; i < max_iterations; i++)
	{
		if (i > 0)
		{
			z = zNxet;
		}
		fdm_matrix = fdmMatrixGChenGuo(input_signal, signal_half_length, z);

		cx_mat u_singular_matrix;
		vec singular_values;
		cx_mat v_singular_matrix;
		svd(u_singular_matrix,
			singular_values,
			v_singular_matrix, fdm_matrix.U0);			// 奇异值分解
		//ds.print();

		//singular_values.print();
		//	最大主成分
		uvec significant_indices = find(abs(singular_values) > ((max(abs(singular_values)) * removal_criteria)));

		// 去除低特征值分量
		cx_mat u_singular_matrix_reduced = u_singular_matrix.cols(significant_indices);		//共轭正交矩阵
		vec singular_values_reduced = singular_values.rows(significant_indices);
		singular_values_matrix_reduced = diagmat(pow(singular_values_reduced, -0.5));	//对角矩阵
		v_singular_matrix_reduced = v_singular_matrix.cols(significant_indices);				//共轭正交矩阵


		// U0 和U1是对称矩阵
		cx_mat sMatrix = singular_values_matrix_reduced * u_singular_matrix_reduced.t() * fdm_matrix.U1 * v_singular_matrix_reduced * singular_values_matrix_reduced;

		//P;	//特征向量
		eig_gen(u, P, sMatrix);

		zNxet = u / abs(u);

		if (z.n_elem == zNxet.n_elem and norm(z - zNxet) < error_threshold)
		{
			break;
		}
	}

	cx_mat B = v_singular_matrix_reduced * singular_values_matrix_reduced * P;

	vec frequency = real(IU * log(u) / time_step / 2 / PI);
	vec decay_constant = -imag(IU * log(u) / time_step);
	vec Q_factor = 2 * PI * abs(frequency) / 2 / decay_constant;
	//# z是输入，u是输出特征值
	//	# 误差计算和归一化
	vec	R(size(u));
	for (size_t i = 0; i < u.n_elem; i++)
	{
		cx_mat normNum = sqrt(B.col(i).st() * fdm_matrix.U0 * B.col(i));
		B.col(i) = B.col(i) / normNum(0);
		cx_double u0 = u(i);
		cx_mat uGuest = sqrt(B.col(i).st() * fdm_matrix.U2 * B.col(i));
		if (real(u0 / uGuest(0)) < 0.0)
		{
			uGuest(0) = -uGuest(0);
		}
		R(i) = abs(log(uGuest(0) / u0)) / abs(log(u0));
	}
	// 计算幅值
	cx_vec cz(z.n_elem);
	vec nv1 = linspace(0, signal_half_length - 1, signal_half_length);
	cx_vec c10 = input_signal.rows(0, signal_half_length - 1);
	for (size_t i = 0; i < z.n_elem; i++)
	{
		cz(i) = dot(c10, pow(cx_vec(nv1.n_elem, fill::value(z(i))), -nv1));
	}

	cx_vec complexAmplitude(u.n_elem);
	for (size_t i = 0; i < u.n_elem; i++)
	{
		complexAmplitude(i) = pow(dot(B.col(i), cz), 2);
	}
	vec amplitude = abs(complexAmplitude);
	vec error_estimate = R;
	vec phase = arg(complexAmplitude);



	//frequency.st().print("frequency");
	//decay_constant.st().print("decay_constant");
	//Q_factor.st().print("Q_factor");
	//amplitude.st().print("amplitude");
	//phase.st().print("phase");
	//error_estimate.st().print("error_estimate");

	//cout << "\tfrequency\t" << "decay_constant\t" << "Q_factor\t"
	//	<< "amplitude\t" << "phase\t" << "error_estimate\n";
	//frequency.st().print();
	//decay_constant.st().print();
	//Q_factor.st().print();
	//amplitude.st().print();
	//phase.st().print();
	//error_estimate.st().print();

	mat fs_detailed(frequency.n_elem , 6 , fill::zeros);
	fs_detailed.col(0) = frequency;
	fs_detailed.col(1) = decay_constant;
	fs_detailed.col(2) = Q_factor;
	fs_detailed.col(3) = amplitude;
	fs_detailed.col(4) = phase;
	fs_detailed.col(5) = error_estimate;
	fs_detailed.print();


	map<string, vec> fs;
	fs["frequency"] = frequency;
	fs["decay_constant"] = decay_constant;
	fs["Q_factor"] = Q_factor;
	fs["amplitude"] = amplitude;
	fs["phase"] = phase;
	fs["error_estimate"] = error_estimate;
 

	 


}




